"""
Joint Gaussian decomposition of positive CCP4 omit-map density.

Requires numpy, scipy, pandas, matplotlib, gemmi.
Fits 1..K shared-width Gaussian components + constant background in a local
Cartesian spherical ROI, ranks by BIC, and checks center identifiability.
The CSV is MODEL-derived, not a validated atomic assignment.
"""

import argparse
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gemmi
from scipy.ndimage import maximum_filter
from scipy.optimize import least_squares
from scipy.optimize import linear_sum_assignment


def read_ccp4(path):
    mp=gemmi.read_ccp4_map(str(path), setup=False)
    arr=np.array(mp.grid, dtype=np.float32, copy=True)
    axes=np.array(mp.axis_positions(),dtype=int)
    starts=np.array([mp.header_i32(i) for i in (5,6,7)],dtype=float)
    sampling=np.array([mp.header_i32(i) for i in (8,9,10)],dtype=float)
    mat=np.array(mp.grid.unit_cell.orth.mat,dtype=float)
    if np.any(sampling<=0): raise ValueError('Invalid CCP4 sampling')
    return arr, axes, starts, sampling, mat


def to_cart(indices, axes, starts, sampling, orth):
    f=(np.asarray(indices,dtype=float)[...,axes]+starts[axes])/sampling
    return f @ orth.T


def extract_roi(rho, axes, starts, sampling, orth, center, radius):
    center=np.asarray(center,float)
    base=to_cart(np.zeros(3),axes,starts,sampling,orth)
    steps=np.stack([to_cart(np.eye(3)[i],axes,starts,sampling,orth)-base for i in range(3)])
    # Bounding grid indices via inverse lattice step matrix (works with oblique cells).
    stepmat=steps.T
    est=np.linalg.solve(stepmat, center-base)
    extent=np.ceil(radius*np.linalg.norm(np.linalg.inv(stepmat),axis=1)).astype(int)+2
    lo=np.maximum(0,np.floor(est).astype(int)-extent)
    hi=np.minimum(rho.shape,np.ceil(est).astype(int)+extent+1)
    if np.any(hi<=lo): raise ValueError('ROI center lies outside the sampled map')
    idx=np.indices(tuple(hi-lo)).reshape(3,-1).T+lo
    xyz=to_cart(idx,axes,starts,sampling,orth)
    d=np.linalg.norm(xyz-center,axis=1)
    val=rho[tuple(idx.T)].astype(float)
    keep=(d<=radius)&np.isfinite(val)
    if keep.sum()<30: raise ValueError('Too few finite density voxels in ROI')
    return xyz[keep],val[keep],float(np.min(np.linalg.norm(steps,axis=1)))


def seeds(xyz, y, max_count, separation, min_amplitude=0):
    # Greedy high-density initialization, with spatial nonmaximum suppression.
    order=np.argsort(y)[::-1]
    chosen=[]
    for j in order:
        if y[j]<=min_amplitude: break
        if all(np.linalg.norm(xyz[j]-xyz[k])>=separation for k in chosen): chosen.append(j)
        if len(chosen)>=max_count: break
    return chosen


def prediction(p, xyz, k):
    b=p[-1]
    out=np.full(len(xyz),b,dtype=float)
    for j in range(k):
        amp,*rest=p[5*j:5*j+5]
        ctr=np.array(rest[:3]); width=rest[3]
        out+=amp*np.exp(-0.5*np.sum(((xyz-ctr)/width)**2,axis=1))
    return out


def one_fit(xyz,y,k,starts,width_bounds=(0.3,2.0),n_restarts=3,rng=None):
    rng=np.random.default_rng(42) if rng is None else rng
    bg=float(np.quantile(y,0.2))
    span=np.ptp(xyz,axis=0)
    lower=[];upper=[]
    ampmax=max(float(np.ptp(y))*8,1.)
    for j in range(k):
        lower.extend([0,*list(np.min(xyz,axis=0)),width_bounds[0]])
        upper.extend([ampmax,*list(np.max(xyz,axis=0)),width_bounds[1]])
    lower.append(float(np.min(y)-np.ptp(y)*2-1))
    upper.append(float(np.max(y)+np.ptp(y)*2+1))
    best=None
    for restart in range(n_restarts):
        p0=[]
        for j in range(k):
            ctr=xyz[starts[j]].copy()
            if restart: ctr+=rng.normal(scale=0.12,size=3)
            ctr=np.clip(ctr,np.min(xyz,axis=0)+1e-5,np.max(xyz,axis=0)-1e-5)
            p0.extend([np.clip(y[starts[j]]-bg,0.01,ampmax-1e-5),*ctr,0.8])
        p0.append(bg)
        fit=least_squares(lambda p:prediction(p,xyz,k)-y,np.clip(p0,np.array(lower)+1e-8,np.array(upper)-1e-8),
                          bounds=(lower,upper),max_nfev=450,method='trf',x_scale='jac')
        rss=float(np.sum(fit.fun**2))
        if best is None or rss<best[0]: best=(rss,fit)
    return best[1]


def analyze_arrays(xyz,y,max_components=4,min_sep=0.6,width_bounds=(0.3,2.0),restarts=3):
    xyz=np.asarray(xyz,float);y=np.asarray(y,float)
    n=len(y); chosen=seeds(xyz,y,max_components,min_sep,float(np.quantile(y,0.2)))
    if not chosen: raise ValueError('No positive candidate peaks in selected ROI')
    rng=np.random.default_rng(2026); fits=[]; rows=[]
    for k in range(1,min(max_components,len(chosen))+1):
        fit=one_fit(xyz,y,k,chosen,width_bounds,restarts,rng)
        rss=float(np.dot(fit.fun,fit.fun)); pnum=5*k+1
        # Common Gaussian noise variance estimated by RSS/n; relative BIC only.
        bic=n*np.log(max(rss/n,1e-14))+pnum*np.log(n)
        aic=n*np.log(max(rss/n,1e-14))+2*pnum
        r2=1-rss/np.sum((y-y.mean())**2) if np.ptp(y)>0 else np.nan
        rows.append(dict(K=k,N=n,RSS=rss,RMSE=np.sqrt(rss/n),R2=r2,BIC=bic,AIC=aic,success=fit.success))
        fits.append(fit)
    comparison=pd.DataFrame(rows); ibest=int(comparison.BIC.argmin());fit=fits[ibest];k=int(comparison.iloc[ibest].K)
    pars=fit.x; components=[]
    # Approximate local curvature errors; not valid for systematic map errors.
    dof=max(n-len(pars),1)
    try:
        cov=np.linalg.pinv(fit.jac.T@fit.jac)*(np.dot(fit.fun,fit.fun)/dof)
        errors=np.sqrt(np.maximum(np.diag(cov),0))
    except np.linalg.LinAlgError: errors=np.full(len(pars),np.nan)
    for j in range(k):
        amp,x,yc,z,w=pars[j*5:j*5+5]
        ea,ex,ey,ez,ew=errors[j*5:j*5+5]
        components.append(dict(component=j+1,amplitude=amp,x_A=x,y_A=yc,z_A=z,
                               sigma_A=w,FWHM_A=2.354820045*w,center_se_x_A=ex,
                               center_se_y_A=ey,center_se_z_A=ez,
                               integrated_gaussian=amp*(2*np.pi)**1.5*w**3))
    comps=pd.DataFrame(components)
    pairs=[]
    for a,b in itertools.combinations(range(k),2):
        ca=comps.iloc[a];cb=comps.iloc[b]
        d=np.linalg.norm(np.array([ca.x_A,ca.y_A,ca.z_A])-np.array([cb.x_A,cb.y_A,cb.z_A]))
        pairs.append(dict(i=a+1,j=b+1,distance_A=d,d_over_sum_sigma=d/(ca.sigma_A+cb.sigma_A),
                          d_over_mean_fwhm=d/(0.5*(ca.FWHM_A+cb.FWHM_A))))
    return comparison,comps,pd.DataFrame(pairs),fit


def plot_fit(xyz,y,comp,fit,path):
    k=len(comp); pred=prediction(fit.x,xyz,k)
    center=np.mean(comp[['x_A','y_A','z_A']].to_numpy(),axis=0)
    axis=np.argmax(np.ptp(xyz,axis=0));sel=np.argsort(xyz[:,axis]); coord=xyz[sel,axis]
    fig,axes=plt.subplots(1,2,figsize=(11,4.2))
    axes[0].scatter(y,pred,s=4,alpha=.35); lo=min(y.min(),pred.min()); hi=max(y.max(),pred.max())
    axes[0].plot([lo,hi],[lo,hi],'k--',lw=1);axes[0].set(xlabel='Observed density',ylabel='Model density',title='Voxel-wise fit')
    axes[1].scatter(xyz[:,axis],y-pred,s=4,alpha=.35)
    axes[1].axhline(0,color='k',ls='--',lw=1)
    axes[1].set(xlabel=['X','Y','Z'][axis]+' (Å)',ylabel='Observed − model',title='Residuals')
    fig.tight_layout();fig.savefig(path,dpi=170);plt.close(fig)


def synthetic_benchmark(output,distances=(0.7,1.0,1.4,1.8,2.2,2.8),width=0.8,noise=0.08):
    rng=np.random.default_rng(5)
    u=np.arange(-3.2,3.21,0.35)
    xyz=np.stack(np.meshgrid(u,u,u,indexing='ij'),-1).reshape(-1,3)
    records=[]
    for d in distances:
        true=np.array([[-d/2,0,0],[d/2,0,0]])
        r2=np.sum((xyz[:,None,:]-true[None,:,:])**2,axis=2)
        y=np.exp(-r2/(2*width*width)).sum(axis=1)+rng.normal(0,noise,size=len(xyz))
        compare,comp,pairs,fit=analyze_arrays(xyz,y,2,0.45,(0.4,1.4),2)
        best=int(compare.loc[compare.BIC.idxmin(),'K'])
        sep=np.nan
        if best==2:
            sep=float(pairs.iloc[0].distance_A)
        records.append(dict(true_separation_A=d,selected_K=best,recovered_separation_A=sep,
                            delta_BIC_1_minus_2=float(compare.iloc[0].BIC-compare.iloc[1].BIC),
                            recovered_two=(best==2),noise_sd=noise,truth_sigma_A=width))
    data=pd.DataFrame(records);data.to_csv(output,index=False)
    return data


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('mapfile',nargs='?',help='CCP4 map (.map / .ccp4)')
    ap.add_argument('--center',nargs=3,type=float,metavar=('X','Y','Z'),help='ROI center in Cartesian Å; default strongest positive voxel')
    ap.add_argument('--radius',type=float,default=4.5,help='ROI sphere radius (Å)')
    ap.add_argument('--max-components',type=int,default=4)
    ap.add_argument('--min-separation',type=float,default=0.6,help='Candidate seed separation (Å), NOT a fitted distance constraint')
    ap.add_argument('--min-width',type=float,default=0.35)
    ap.add_argument('--max-width',type=float,default=1.8)
    ap.add_argument('--restarts',type=int,default=3)
    ap.add_argument('--output-prefix',default='multi_mn')
    ap.add_argument('--synthetic',action='store_true',help='Run synthetic two-center resolution benchmark')
    args=ap.parse_args()
    if args.synthetic:
        df=synthetic_benchmark(args.output_prefix+'_synthetic.csv')
        print(df.to_string(index=False));return
    if not args.mapfile: ap.error('mapfile is required unless --synthetic')
    rho,axes,starts,sampling,orth=read_ccp4(args.mapfile)
    if args.center is None:
        ix=np.unravel_index(np.nanargmax(rho),rho.shape)
        center=to_cart(np.array(ix),axes,starts,sampling,orth)
    else: center=np.array(args.center)
    print('ROI center (Å):',center,'map shape:',rho.shape)
    xyz,y,spacing=extract_roi(rho,axes,starts,sampling,orth,center,args.radius)
    print('Finite ROI voxels:',len(y),'grid spacing approx (Å):',spacing)
    comparison,comp,pairs,fit=analyze_arrays(xyz,y,args.max_components,args.min_separation,
                                               (args.min_width,args.max_width),args.restarts)
    prefix=args.output_prefix
    comparison.to_csv(prefix+'_models.csv',index=False);comp.to_csv(prefix+'_centers.csv',index=False)
    pairs.to_csv(prefix+'_pairs.csv',index=False)
    plot_fit(xyz,y,comp,fit,prefix+'_diagnostics.png')
    print('\nModel comparison:\n',comparison.to_string(index=False))
    print('\nBIC-selected centers (Cartesian Å):\n',comp.to_string(index=False))
    print('\nPair separations:\n',pairs.to_string(index=False) if len(pairs) else '(one component)')
    print('\nNote: K, BIC, and covariance SE do not alone prove physical resolvability.')

if __name__=='__main__': main()

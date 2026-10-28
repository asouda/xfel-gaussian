"""
Sequential fitting/subtraction of positive difference-map peaks.

Dependencies: pip install gemmi numpy scipy pandas matplotlib
Usage:
 python sequential_ccp4_gaussian.py omit_OEC_15A.map --n-peaks 6 --output-prefix OEC_seq

Read CCP4 without altering its extent, orientation, or origin. Fit a local
spherical 3D Gaussian + constant baseline at the highest remaining peak;
subtract only the Gaussian (NOT the local baseline) from the entire sampled
map; repeat. Model-derived centers do not establish elemental identities.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
from scipy.optimize import least_squares
from scipy.ndimage import maximum_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FWHM_FACTOR = 2 * np.sqrt(2 * np.log(2))

class DensityMap:
    def __init__(self, path):
        self.ccp4 = gemmi.read_ccp4_map(str(path), setup=False)
        self.original = np.array(self.ccp4.grid, dtype=np.float64, copy=True)
        self.axis = np.asarray(self.ccp4.axis_positions(), dtype=int)
        self.start = np.asarray([self.ccp4.header_i32(i) for i in (5, 6, 7)], dtype=float)
        self.sampling = np.asarray([self.ccp4.header_i32(i) for i in (8, 9, 10)], dtype=float)
        if np.any(self.sampling <= 0):
            raise ValueError('Invalid CCP4 sampling intervals')
        orth = np.array(self.ccp4.grid.unit_cell.orth.mat, dtype=float)
        # The 3 array-index directions are not necessarily x,y,z crystallographic directions.
        self.origin = self.to_cart(np.zeros(3), orth)
        self.step = np.stack([self.to_cart(np.eye(3)[i], orth) - self.origin
                              for i in range(3)], axis=1)
        self.invstep = np.linalg.inv(self.step)
        self.voxel_A3 = abs(np.linalg.det(self.step))
        self.residual = self.original.copy()

    def to_cart(self, idx, orth=None):
        if orth is None:
            orth = np.array(self.ccp4.grid.unit_cell.orth.mat, dtype=float)
        fractional = (np.asarray(idx, dtype=float)[..., self.axis] + self.start[self.axis]) / self.sampling
        return fractional @ orth.T

    def to_index(self, xyz):
        return (np.asarray(xyz) - self.origin) @ self.invstep.T

    def local(self, center, radius, values=None):
        """Return Cartesian positions, densities and indices within a real-space sphere."""
        values = self.residual if values is None else values
        c = self.to_index(center)
        reach = radius * np.linalg.norm(self.invstep, axis=1) + 2
        lo = np.maximum(0, np.floor(c - reach).astype(int))
        hi = np.minimum(values.shape, np.ceil(c + reach).astype(int) + 1)
        if np.any(lo >= hi):
            raise ValueError('Center outside map')
        idx = np.indices(tuple(hi-lo)).reshape(3, -1).T + lo
        xyz = self.to_cart(idx)
        y = values[tuple(idx.T)]
        keep = (np.linalg.norm(xyz - center, axis=1) <= radius) & np.isfinite(y)
        return xyz[keep], y[keep], idx[keep]

    def subtract(self, amp, center, sigma, radius_factor=5.0):
        # Update a compact region; beyond 5 sigma, component is negligible.
        xyz, _, idx = self.local(center, radius_factor*sigma)
        component = gaussian(xyz, amp, center, sigma)
        self.residual[tuple(idx.T)] -= component

    def save_map(self, path):
        # Preserve the CCP4 extent, starts, axes, and crystallographic frame.
        a = np.array(self.ccp4.grid, copy=False)
        a[:] = np.where(np.isfinite(self.residual), self.residual, 0).astype(np.float32)
        self.ccp4.write_ccp4_map(str(path))


def gaussian(xyz, amp, center, sigma):
    return amp * np.exp(-0.5 * np.sum(((xyz - center)/sigma)**2, axis=-1))


def robust_noise(y):
    z = y[np.isfinite(y)]
    if len(z) < 20: raise ValueError('Not enough valid map voxels')
    median = np.median(z)
    mad = np.median(np.abs(z - median))
    return max(1.4826*mad, 1e-7)


def choose_peak(rho, original_mask, exclusion, local_max_size=3):
    # Exclude edges where local Gaussian parameters may be poorly constrained.
    valid = np.isfinite(rho) & original_mask & ~exclusion
    candidates = valid & (rho == maximum_filter(np.where(valid, rho, -np.inf),
                                                  size=local_max_size, mode='constant', cval=-np.inf))
    if not np.any(candidates): return None
    ind = np.argmax(np.where(candidates, rho, -np.inf))
    idx = np.unravel_index(ind, rho.shape)
    if not np.isfinite(rho[idx]): return None
    return np.array(idx), float(rho[idx])


def fit_peak(dmap, seed, radius, sigma_bounds, shift_limit, max_nfev=300):
    xyz, y, _ = dmap.local(seed, radius)
    if len(y) < 35: raise RuntimeError('Insufficient local valid voxels')
    bg = float(np.quantile(y, .2))
    amp = max(float(np.max(y)-bg), 1e-3)
    initial = np.r_[amp, seed, np.clip(radius/3, *sigma_bounds), bg]
    lo = np.r_[0., seed-shift_limit, sigma_bounds[0], np.min(y)-np.ptp(y)-1.]
    hi = np.r_[max(10*amp, 1.), seed+shift_limit, sigma_bounds[1], np.max(y)+np.ptp(y)+1.]
    initial = np.clip(initial, lo+1e-7, hi-1e-7)
    def pred(p): return gaussian(xyz, p[0], p[1:4], p[4])+p[5]
    res = least_squares(lambda p: pred(p)-y, initial, bounds=(lo, hi),
                        method='trf', loss='linear', x_scale='jac', max_nfev=max_nfev)
    p = res.x
    rss = float(np.dot(res.fun, res.fun))
    tot = float(np.dot(y-y.mean(), y-y.mean()))
    r2 = 1-rss/tot if tot > 0 else float('nan')
    return dict(amp=p[0], xyz=p[1:4], sigma=p[4], bg=p[5], r2=r2,
                rmse=np.sqrt(rss/len(y)), n=len(y), success=res.success)


def run(args):
    dm = DensityMap(args.map)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    valid = np.isfinite(dm.original)
    reference_noise = robust_noise(dm.original)
    threshold = args.threshold_sigma * reference_noise
    excluded = np.zeros(dm.original.shape, bool)
    history = []
    if args.save_steps: dm.save_map(f'{prefix}_step00.map')
    print(f'Grid: {dm.original.shape}; robust map sigma={reference_noise:.5g}; '
          f'positive peak threshold={threshold:.5g}')
    for k in range(1, args.n_peaks+1):
        peak = choose_peak(dm.residual, valid, excluded)
        if peak is None: break
        idx, height = peak
        if height < threshold:
            print(f'Stopping: strongest remaining eligible peak {height:.4f} below threshold.')
            break
        seed = dm.to_cart(idx)
        try:
            f = fit_peak(dm, seed, args.fit_radius, (args.sigma_min, args.sigma_max),
                         args.max_shift)
        except (RuntimeError, ValueError) as exc:
            print(f'Skipping seed {idx.tolist()}: {exc}')
            f = None
        if f is not None and f['success'] and f['amp'] > 0:
            center = f['xyz']
            # Gaussian subtraction only; retaining local background prevents map-wide drift.
            dm.subtract(f['amp'], center, f['sigma'])
            row = dict(iteration=k, seed_x_A=seed[0],seed_y_A=seed[1],seed_z_A=seed[2],
                       x_A=center[0], y_A=center[1], z_A=center[2],
                       peak_before=height, amplitude=f['amp'], sigma_A=f['sigma'],
                       FWHM_A=FWHM_FACTOR*f['sigma'], background=f['bg'],
                       integrated_gaussian=f['amp']*(2*np.pi)**1.5*f['sigma']**3,
                       R2=f['r2'], RMSE=f['rmse'], n_voxels=f['n'])
            history.append(row)
            print(f'#{k}: center=({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}) Å '
                  f'amp={f["amp"]:.3f}, width={f["sigma"]:.3f} Å, R2={f["r2"]:.3f}')
            if args.save_steps: dm.save_map(f'{prefix}_step{k:02d}.map')
        # Avoid repeatedly selecting same site, including failed/unstable fits.
        marker = f['xyz'] if f is not None else seed
        _, _, around = dm.local(marker, args.min_separation, dm.original)
        excluded[tuple(around.T)] = True
    df = pd.DataFrame(history)
    df.to_csv(f'{prefix}_centers.csv', index=False)
    dm.save_map(f'{prefix}_residual.map')
    if len(df):
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
        axs[0].plot(df.iteration, df.amplitude, 'o-')
        axs[0].set(xlabel='Subtraction iteration',ylabel='Fitted Gaussian amplitude',title='Sequential peak amplitudes')
        axs[1].plot(df.iteration, df.R2, 'o-')
        axs[1].set(xlabel='Subtraction iteration', ylabel='Local fit R²',title='Fit quality')
        axs[1].set_ylim(min(-.1, float(df.R2.min())-.05), 1.05)
        fig.tight_layout(); fig.savefig(f'{prefix}_diagnostics.png',dpi=180);plt.close(fig)
    print(f'Saved {len(df)} components; results: {prefix}_centers.csv and {prefix}_residual.map')
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('map',help='CCP4 or MAP file (full or cropped)')
    p.add_argument('--output-prefix',default='sequential_peaks')
    p.add_argument('--n-peaks',type=int,default=6)
    p.add_argument('--threshold-sigma',type=float,default=2.0,
                   help='Positive amplitude threshold in robust original-map sigma')
    p.add_argument('--fit-radius',type=float,default=2.2,help='Local spherical fit radius (Å)')
    p.add_argument('--sigma-min',type=float,default=.35)
    p.add_argument('--sigma-max',type=float,default=1.8)
    p.add_argument('--max-shift',type=float,default=1.2,
                   help='Max fitted displacement from seed in each Cartesian axis (Å)')
    p.add_argument('--min-separation',type=float,default=.7,
                   help='Seed exclusion around accepted fitted centers (Å)')
    p.add_argument('--save-steps',action='store_true',help='Export CCP4 residual after each subtraction')
    args=p.parse_args()
    if not (args.fit_radius>0 and args.sigma_min>0 and args.sigma_max>args.sigma_min):
        p.error('Invalid fitting radius or width bounds')
    run(args)

if __name__=='__main__': main()

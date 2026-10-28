"""
Joint Gaussian decomposition of positive CCP4 omit-map density.

Requires numpy, scipy, pandas, matplotlib, gemmi, and pyyaml.
Fits 1..K Gaussian components plus a constant background in a local Cartesian
spherical ROI, ranks the models by BIC, and checks center identifiability.
Each component is either isotropic or a fully rotated ellipsoid (full covariance).
Centers are either refined or held at the detected local maxima.

The CSV is MODEL-derived, not a validated atomic assignment. The PDB file
labels the first four centers MN1–MN4 and later centers O1, O2, and so on
only so the coordinates can be opened beside a structure.

Usage:
    python -m xfel_gaussian.multi_peaks_gaussian config/peaks.yaml

The YAML file may be a flat mapping or a combined file with a ``multi_peaks``
section and a shared ``map`` entry. See config/peaks.yaml.
"""

import argparse
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gemmi
import yaml
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular
from scipy.spatial import cKDTree

FWHM_FACTOR = 2 * np.sqrt(2 * np.log(2))

DEFAULTS = dict(
    center=None,
    radius=4.5,
    max_components=4,
    min_separation=0.6,
    min_width=0.35,
    max_width=1.8,
    restarts=3,
    output_prefix='multi_mn',
    synthetic=False,
    fix_centers=False,
    full_covariance=False,
)


def load_config(path):
    """Read a YAML mapping and remember its path so relative files resolve beside it."""
    path = Path(path)
    with path.open() as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f'{path} must contain a YAML mapping')
    data['_config_path'] = str(path.resolve())
    return data


def section_config(cfg, name):
    """Return this program's mapping from a flat file or a combined multi-section file."""
    other = 'sequential_peaks' if name == 'multi_peaks' else 'multi_peaks'
    base = {'_config_path': cfg.get('_config_path')}
    if name in cfg and isinstance(cfg[name], dict):
        base.update(cfg[name])
        base.setdefault('map', cfg.get('map'))
        return base
    if other in cfg and isinstance(cfg.get(other), dict):
        raise KeyError(f'Configuration has {other!r} but not {name!r}')
    flat = {key: value for key, value in cfg.items() if key != name}
    flat['_config_path'] = cfg.get('_config_path')
    return flat


def _config_dir(cfg):
    raw = cfg.get('_config_path')
    return Path(raw).resolve().parent if raw else Path.cwd()


def _resolve_path(cfg, value):
    path = Path(value)
    if not path.is_absolute():
        path = _config_dir(cfg) / path
    return path.resolve()


def _as_bool(value, key):
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in ('true', 'yes', '1', 'false', 'no', '0'):
        return value.strip().lower() in ('true', 'yes', '1')
    raise ValueError(f'{key} must be a boolean, got {value!r}')


def resolve_config(cfg):
    """Apply defaults to the multi-peak section and check the keys."""
    raw = section_config(cfg, 'multi_peaks')
    unknown = set(raw) - set(DEFAULTS) - {'map', '_config_path'}
    if unknown:
        raise ValueError(f'Unknown multi_peaks settings: {", ".join(sorted(unknown))}')
    resolved = dict(DEFAULTS)
    for key, value in raw.items():
        if value is None and key in DEFAULTS:
            continue
        resolved[key] = value
    resolved['_config_path'] = raw.get('_config_path')
    resolved['map'] = raw.get('map')
    resolved['fix_centers'] = _as_bool(resolved['fix_centers'], 'fix_centers')
    resolved['full_covariance'] = _as_bool(resolved['full_covariance'], 'full_covariance')
    resolved['synthetic'] = _as_bool(resolved['synthetic'], 'synthetic')
    if resolved['min_width'] <= 0 or resolved['max_width'] <= resolved['min_width']:
        raise ValueError('min_width and max_width must satisfy 0 < min_width < max_width')
    if resolved['radius'] <= 0 or resolved['max_components'] < 1 or resolved['restarts'] < 1:
        raise ValueError('radius, max_components, and restarts must be positive')
    if not resolved['synthetic']:
        if not resolved.get('map'):
            raise ValueError('map is required unless synthetic is true')
        resolved['map'] = _resolve_path(resolved, resolved['map'])
        if not resolved['map'].is_file():
            raise FileNotFoundError(resolved['map'])
    resolved['output_prefix'] = str(_resolve_path(resolved, resolved['output_prefix']))
    if resolved['center'] is not None:
        center = np.asarray(resolved['center'], dtype=float)
        if center.shape != (3,):
            raise ValueError('center must be three Cartesian coordinates in Å')
        resolved['center'] = center
    return resolved


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


def estimate_spacing(xyz):
    if len(xyz) < 2:
        return 0.5
    distances, _ = cKDTree(xyz).query(xyz, k=2)
    nearest = distances[:, 1]
    nearest = nearest[np.isfinite(nearest) & (nearest > 0)]
    if len(nearest) == 0:
        return 0.5
    return float(np.median(nearest))


def local_maximum_mask(xyz, y, spacing):
    """True where a sample is at least as high as every voxel within one grid cell."""
    cutoff = 1.8 * float(spacing)
    groups = cKDTree(xyz).query_ball_point(xyz, cutoff)
    mask = np.ones(len(y), dtype=bool)
    for i, group in enumerate(groups):
        if y[i] < np.max(y[group]):
            mask[i] = False
    return mask


def cholesky_matrix(values):
    """Lower-triangular factor from [L00, L10, L20, L11, L21, L22]."""
    return np.array([
        [values[0], 0.0, 0.0],
        [values[1], values[3], 0.0],
        [values[2], values[4], values[5]],
    ], dtype=float)


def mahalanobis2(xyz, center, cholesky):
    delta = np.asarray(xyz, dtype=float) - np.asarray(center, dtype=float)
    local = solve_triangular(cholesky, np.swapaxes(delta, -1, 0), lower=True, check_finite=False)
    return np.sum(local * local, axis=0)


def principal_frame(cholesky):
    """Descending principal widths, a right-handed eigenframe, and Σ = LL^T."""
    covariance = cholesky @ cholesky.T
    values, vectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    order = np.argsort(values)[::-1]
    values = np.clip(values[order], 0.0, None)
    vectors = vectors[:, order]
    for i in range(3):
        dominant = int(np.argmax(np.abs(vectors[:, i])))
        if vectors[dominant, i] < 0:
            vectors[:, i] *= -1
    if np.linalg.det(vectors) < 0:
        vectors[:, 2] *= -1
    return np.sqrt(values), vectors, covariance


def component_stride(fix_centers, full_covariance):
    """Free parameters per component: amplitude, optional center, then shape."""
    return 1 + (0 if fix_centers else 3) + (6 if full_covariance else 1)


def prediction(parameters, xyz, k, centers, fix_centers, full_covariance):
    background = parameters[-1]
    out = np.full(len(xyz), background, dtype=float)
    stride = component_stride(fix_centers, full_covariance)
    for j in range(k):
        block = parameters[stride * j:stride * (j + 1)]
        amplitude = block[0]
        cursor = 1
        if fix_centers:
            center = centers[j]
        else:
            center = block[cursor:cursor + 3]
            cursor += 3
        if full_covariance:
            quadratic = mahalanobis2(xyz, center, cholesky_matrix(block[cursor:cursor + 6]))
        else:
            quadratic = np.sum(((xyz - center) / block[cursor]) ** 2, axis=1)
        out += amplitude * np.exp(-0.5 * quadratic)
    return out


def _shape_record(cholesky):
    widths, axes, covariance = principal_frame(cholesky)
    equivalent = float(np.prod(widths) ** (1.0 / 3.0))
    record = dict(
        sigma_A=equivalent,
        sigma_1_A=float(widths[0]), sigma_2_A=float(widths[1]), sigma_3_A=float(widths[2]),
        FWHM_A=float(FWHM_FACTOR * equivalent),
        fwhm_1_A=float(FWHM_FACTOR * widths[0]),
        fwhm_2_A=float(FWHM_FACTOR * widths[1]),
        fwhm_3_A=float(FWHM_FACTOR * widths[2]),
        cov_xx=float(covariance[0, 0]), cov_yy=float(covariance[1, 1]), cov_zz=float(covariance[2, 2]),
        cov_xy=float(covariance[0, 1]), cov_xz=float(covariance[0, 2]), cov_yz=float(covariance[1, 2]),
        L00=float(cholesky[0, 0]), L10=float(cholesky[1, 0]), L20=float(cholesky[2, 0]),
        L11=float(cholesky[1, 1]), L21=float(cholesky[2, 1]), L22=float(cholesky[2, 2]),
        condition_number=float(np.linalg.cond(covariance)),
        integrated_scale=float(np.prod(np.diag(cholesky))),
    )
    for axis, label in enumerate(('v1', 'v2', 'v3')):
        record[f'{label}_x'] = float(axes[0, axis])
        record[f'{label}_y'] = float(axes[1, axis])
        record[f'{label}_z'] = float(axes[2, axis])
    return record


def one_fit(xyz, y, k, starts, width_bounds=(0.3, 2.0), n_restarts=3, rng=None,
            fix_centers=False, full_covariance=False):
    rng = np.random.default_rng(42) if rng is None else rng
    bg = float(np.quantile(y, 0.2))
    ampmax = max(float(np.ptp(y)) * 8, 1.)
    seed_centers = np.asarray(xyz[list(starts[:k])], dtype=float)
    smin, smax = width_bounds
    lower, upper = [], []
    for _ in range(k):
        lower.append(0.0)
        upper.append(ampmax)
        if not fix_centers:
            lower.extend(np.min(xyz, axis=0))
            upper.extend(np.max(xyz, axis=0))
        if full_covariance:
            lower.extend([smin, -smax, -smax, smin, -smax, smin])
            upper.extend([smax, smax, smax, smax, smax, smax])
        else:
            lower.append(smin)
            upper.append(smax)
    lower.append(float(np.min(y) - np.ptp(y) * 2 - 1))
    upper.append(float(np.max(y) + np.ptp(y) * 2 + 1))
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    best = None
    stride = component_stride(fix_centers, full_covariance)
    nfev = max(450, 60 * (stride * k + 1))
    for restart in range(n_restarts):
        p0 = []
        for j in range(k):
            amplitude = np.clip(y[starts[j]] - bg, 0.01, ampmax - 1e-5)
            p0.append(amplitude)
            if not fix_centers:
                center = seed_centers[j].copy()
                if restart:
                    center = center + rng.normal(scale=0.12, size=3)
                center = np.clip(center, np.min(xyz, axis=0) + 1e-5, np.max(xyz, axis=0) - 1e-5)
                p0.extend(center)
            width0 = float(np.clip(0.8 + (rng.normal(scale=0.05) if restart else 0.0), smin, smax))
            if full_covariance:
                p0.extend([width0, 0.0, 0.0, width0, 0.0, width0])
            else:
                p0.append(width0)
        p0.append(bg)
        p0 = np.clip(np.asarray(p0, dtype=float), lower + 1e-8, upper - 1e-8)
        fit = least_squares(
            lambda p: prediction(p, xyz, k, seed_centers, fix_centers, full_covariance) - y,
            p0, bounds=(lower, upper), max_nfev=nfev, method='trf', x_scale='jac',
        )
        rss = float(np.sum(fit.fun ** 2))
        if best is None or rss < best[0]:
            best = (rss, fit)
    return best[1]


def _describe_components(parameters, errors, k, seed_centers, fix_centers, full_covariance):
    stride = component_stride(fix_centers, full_covariance)
    background = float(parameters[-1])
    components = []
    for j in range(k):
        block = parameters[stride * j:stride * (j + 1)]
        err = errors[stride * j:stride * (j + 1)]
        cursor = 1
        if fix_centers:
            center = seed_centers[j]
            center_se = np.full(3, np.nan)
        else:
            center = block[cursor:cursor + 3]
            center_se = err[cursor:cursor + 3]
            cursor += 3
        if full_covariance:
            factor = cholesky_matrix(block[cursor:cursor + 6])
        else:
            factor = np.diag([block[cursor], block[cursor], block[cursor]])
        shape = _shape_record(factor)
        shape.update(dict(
            component=j + 1,
            amplitude=float(block[0]),
            x_A=float(center[0]), y_A=float(center[1]), z_A=float(center[2]),
            center_se_x_A=float(center_se[0]), center_se_y_A=float(center_se[1]), center_se_z_A=float(center_se[2]),
            center_fixed=bool(fix_centers),
            full_covariance=bool(full_covariance),
            background=background,
            integrated_gaussian=float(block[0] * (2 * np.pi) ** 1.5 * shape.pop('integrated_scale')),
        ))
        components.append(shape)
    return components


def analyze_arrays(xyz, y, max_components=4, min_sep=0.6, width_bounds=(0.3, 2.0), restarts=3,
                   fix_centers=False, full_covariance=False, spacing=None):
    xyz = np.asarray(xyz, float)
    y = np.asarray(y, float)
    n = len(y)
    if spacing is None:
        spacing = estimate_spacing(xyz)
    if fix_centers:
        pool = np.flatnonzero(local_maximum_mask(xyz, y, spacing))
        if len(pool) == 0:
            raise ValueError('No local maxima in the selected ROI')
        chosen_local = seeds(xyz[pool], y[pool], max_components, min_sep, float(np.quantile(y, 0.2)))
        chosen = [int(pool[j]) for j in chosen_local]
    else:
        chosen = seeds(xyz, y, max_components, min_sep, float(np.quantile(y, 0.2)))
    if not chosen:
        raise ValueError('No positive candidate peaks in selected ROI')
    rng = np.random.default_rng(2026)
    fits, rows = [], []
    stride = component_stride(fix_centers, full_covariance)
    for k in range(1, min(max_components, len(chosen)) + 1):
        fit = one_fit(xyz, y, k, chosen, width_bounds, restarts, rng, fix_centers, full_covariance)
        rss = float(np.dot(fit.fun, fit.fun))
        pnum = stride * k + 1
        # Common Gaussian noise variance estimated by RSS/n; relative BIC only.
        bic = n * np.log(max(rss / n, 1e-14)) + pnum * np.log(n)
        aic = n * np.log(max(rss / n, 1e-14)) + 2 * pnum
        r2 = 1 - rss / np.sum((y - y.mean()) ** 2) if np.ptp(y) > 0 else np.nan
        rows.append(dict(
            K=k, N=n, n_parameters=pnum, RSS=rss, RMSE=np.sqrt(rss / n), R2=r2, BIC=bic, AIC=aic,
            success=bool(fit.success), n_at_bound=int(np.count_nonzero(fit.active_mask)),
            fix_centers=bool(fix_centers), full_covariance=bool(full_covariance),
        ))
        fits.append(fit)
    comparison = pd.DataFrame(rows)
    ibest = int(comparison.BIC.argmin())
    fit = fits[ibest]
    k = int(comparison.iloc[ibest].K)
    seed_centers = np.asarray(xyz[chosen[:k]], dtype=float)
    dof = max(n - len(fit.x), 1)
    try:
        curvature = np.linalg.pinv(fit.jac.T @ fit.jac) * (np.dot(fit.fun, fit.fun) / dof)
        errors = np.sqrt(np.maximum(np.diag(curvature), 0))
    except np.linalg.LinAlgError:
        errors = np.full(len(fit.x), np.nan)
    components = _describe_components(fit.x, errors, k, seed_centers, fix_centers, full_covariance)
    comps = pd.DataFrame(components)
    pairs = []
    for a, b in itertools.combinations(range(k), 2):
        ca, cb = comps.iloc[a], comps.iloc[b]
        distance = np.linalg.norm(np.array([ca.x_A, ca.y_A, ca.z_A]) - np.array([cb.x_A, cb.y_A, cb.z_A]))
        pairs.append(dict(
            i=a + 1, j=b + 1, distance_A=distance,
            d_over_sum_sigma=distance / (ca.sigma_A + cb.sigma_A),
            d_over_mean_fwhm=distance / (0.5 * (ca.FWHM_A + cb.FWHM_A)),
        ))
    return comparison, comps, pd.DataFrame(pairs), fit


def plot_fit(xyz, y, comp, fit, path):
    k = len(comp)
    fix_centers = bool(comp['center_fixed'].iloc[0])
    full_covariance = bool(comp['full_covariance'].iloc[0])
    centers = comp[['x_A', 'y_A', 'z_A']].to_numpy() if fix_centers else None
    pred = prediction(fit.x, xyz, k, centers, fix_centers, full_covariance)
    axis = int(np.argmax(np.ptp(xyz, axis=0)))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].scatter(y, pred, s=4, alpha=.35)
    lo = min(y.min(), pred.min())
    hi = max(y.max(), pred.max())
    axes[0].plot([lo, hi], [lo, hi], 'k--', lw=1)
    axes[0].set(xlabel='Observed density', ylabel='Model density', title='Voxel-wise fit')
    axes[1].scatter(xyz[:, axis], y - pred, s=4, alpha=.35)
    axes[1].axhline(0, color='k', ls='--', lw=1)
    axes[1].set(xlabel=['X', 'Y', 'Z'][axis] + ' (Å)', ylabel='Observed − model', title='Residuals')
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def synthetic_benchmark(output, distances=(0.7, 1.0, 1.4, 1.8, 2.2, 2.8), width=0.8, noise=0.08,
                        fix_centers=False, full_covariance=False):
    rng = np.random.default_rng(5)
    u = np.arange(-3.2, 3.21, 0.35)
    xyz = np.stack(np.meshgrid(u, u, u, indexing='ij'), -1).reshape(-1, 3)
    records = []
    for d in distances:
        true = np.array([[-d / 2, 0, 0], [d / 2, 0, 0]])
        r2 = np.sum((xyz[:, None, :] - true[None, :, :]) ** 2, axis=2)
        y = np.exp(-r2 / (2 * width * width)).sum(axis=1) + rng.normal(0, noise, size=len(xyz))
        compare, comp, pairs, fit = analyze_arrays(
            xyz, y, 2, 0.45, (0.4, 1.4), 2,
            fix_centers=fix_centers, full_covariance=full_covariance, spacing=0.35,
        )
        best = int(compare.loc[compare.BIC.idxmin(), 'K'])
        sep = float(pairs.iloc[0].distance_A) if best == 2 else np.nan
        records.append(dict(
            true_separation_A=d, selected_K=best, recovered_separation_A=sep,
            delta_BIC_1_minus_2=float(compare.iloc[0].BIC - compare.iloc[1].BIC),
            recovered_two=(best == 2), noise_sd=noise, truth_sigma_A=width,
            fix_centers=fix_centers, full_covariance=full_covariance,
        ))
    data = pd.DataFrame(records)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(output, index=False)
    return data


def write_centers_pdb(path, coordinates, serial_start=43042):
    """Write fitted Cartesian centers as fixed-column HETATM records.

    Records are 80 columns. The first four centers are named MN1–MN4 and
    later centers O1, O2, .... Trailing spaces are part of the record.
    """
    lines = []
    for i, xyz in enumerate(coordinates):
        if i < 4:
            name, element = f"MN{i + 1}", "MN"
        else:
            name, element = f"O{i - 3}", "O"
        x, y, z = (float(value) for value in xyz)
        lines.append(
            f"HETATM{serial_start + i:5d} {name:<4s}BOEY A{601:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{27.27:6.2f}          {element:<4s}\n"
        )
    Path(path).write_text("".join(lines))


def run(cfg, verbose=True):
    """Fit the multi-peak model described by a YAML mapping. Returns result tables."""
    cfg = resolve_config(cfg)
    prefix = Path(cfg['output_prefix'])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if cfg['synthetic']:
        frame = synthetic_benchmark(
            str(prefix) + '_synthetic.csv',
            fix_centers=cfg['fix_centers'], full_covariance=cfg['full_covariance'],
        )
        if verbose:
            print(frame.to_string(index=False))
        return {'synthetic': frame, 'synthetic_csv': str(prefix) + '_synthetic.csv'}

    rho, axes, starts, sampling, orth = read_ccp4(cfg['map'])
    if cfg['center'] is None:
        ix = np.unravel_index(np.nanargmax(rho), rho.shape)
        center = to_cart(np.array(ix), axes, starts, sampling, orth)
    else:
        center = np.asarray(cfg['center'], dtype=float)
    if verbose:
        mode = 'full covariance' if cfg['full_covariance'] else 'isotropic'
        where = 'fixed at detected maxima' if cfg['fix_centers'] else 'refined'
        print(f'ROI center (Å): {center}  map shape: {rho.shape}')
        print(f'Model: {mode}; centers {where}')
    xyz, y, spacing = extract_roi(rho, axes, starts, sampling, orth, center, cfg['radius'])
    if verbose:
        print(f'Finite ROI voxels: {len(y)}  grid spacing approx (Å): {spacing:.4g}')
    comparison, comps, pairs, fit = analyze_arrays(
        xyz, y, cfg['max_components'], cfg['min_separation'],
        (cfg['min_width'], cfg['max_width']), cfg['restarts'],
        fix_centers=cfg['fix_centers'], full_covariance=cfg['full_covariance'], spacing=spacing,
    )
    comparison.to_csv(str(prefix) + '_models.csv', index=False)
    comps.to_csv(str(prefix) + '_centers.csv', index=False)
    centers_pdb = str(prefix) + '_centers.pdb'
    write_centers_pdb(centers_pdb, comps[['x_A', 'y_A', 'z_A']].to_numpy())
    pairs.to_csv(str(prefix) + '_pairs.csv', index=False)
    diagnostics = str(prefix) + '_diagnostics.png'
    plot_fit(xyz, y, comps, fit, diagnostics)
    if verbose:
        print('\nModel comparison:\n', comparison.to_string(index=False))
        print('\nBIC-selected centers (Cartesian Å):\n', comps.to_string(index=False))
        print('\nPair separations:\n', pairs.to_string(index=False) if len(pairs) else '(one component)')
        print('\nNote: K, BIC, and covariance SE do not alone prove physical resolvability.')
        print(f'Centers PDB: {centers_pdb}')
    return {
        'comparison': comparison,
        'components': comps,
        'pairs': pairs,
        'diagnostics': diagnostics,
        'models_csv': str(prefix) + '_models.csv',
        'centers_csv': str(prefix) + '_centers.csv',
        'centers_pdb': centers_pdb,
        'pairs_csv': str(prefix) + '_pairs.csv',
        'roi_center': np.asarray(center, dtype=float),
        'background': float(fit.x[-1]),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('config', help='YAML configuration file')
    args = parser.parse_args(argv)
    return run(load_config(args.config))

if __name__ == '__main__':
    main()

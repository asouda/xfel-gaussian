"""
Sequential fitting and subtraction of positive difference-map peaks.

Dependencies: gemmi, numpy, scipy, pandas, matplotlib, pyyaml.
Usage:
    python -m xfel_gaussian.sequential_peaks_gaussian config/peaks.yaml

Read CCP4 without altering its extent, orientation, or origin. At the highest
remaining local maximum, fit a local 3D Gaussian plus a constant baseline.
The Gaussian is isotropic or a fully rotated ellipsoid. Its center is either
refined or held at that detected maximum. Subtract only the Gaussian (not the
local baseline) from the sampled map, then repeat.

The YAML file may be a flat mapping or a combined file with a
``sequential_peaks`` section and a shared ``map`` entry. See config/peaks.yaml.
Model-derived centers do not establish elemental identities. The PDB file
labels the first four centers MN1–MN4 and later centers O1, O2, and so on
only so the coordinates can be opened beside a structure.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
import yaml
from scipy.optimize import least_squares
from scipy.ndimage import maximum_filter
from scipy.linalg import solve_triangular
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

FWHM_FACTOR = 2 * np.sqrt(2 * np.log(2))

DEFAULTS = dict(
    output_prefix='sequential_peaks',
    n_peaks=6,
    threshold_sigma=2.0,
    fit_radius=2.2,
    sigma_min=0.35,
    sigma_max=1.8,
    max_shift=1.2,
    min_separation=0.7,
    save_steps=False,
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
    if name in cfg and isinstance(cfg[name], dict):
        out = {'_config_path': cfg.get('_config_path')}
        out.update(cfg[name])
        out.setdefault('map', cfg.get('map'))
        return out
    if other in cfg and isinstance(cfg.get(other), dict):
        raise KeyError(f'Configuration has {other!r} but not {name!r}')
    flat = dict(cfg)
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
    """Apply defaults to the sequential-peak section and check the keys."""
    raw = section_config(cfg, 'sequential_peaks')
    unknown = set(raw) - set(DEFAULTS) - {'map', '_config_path'}
    if unknown:
        raise ValueError(f'Unknown sequential_peaks settings: {", ".join(sorted(unknown))}')
    resolved = dict(DEFAULTS)
    for key, value in raw.items():
        if value is None and key in DEFAULTS:
            continue
        resolved[key] = value
    resolved['_config_path'] = raw.get('_config_path')
    if not raw.get('map'):
        raise ValueError('map is required')
    resolved['map'] = _resolve_path(resolved, raw['map'])
    if not Path(resolved['map']).is_file():
        raise FileNotFoundError(resolved['map'])
    resolved['output_prefix'] = str(_resolve_path(resolved, resolved['output_prefix']))
    for key in ('save_steps', 'fix_centers', 'full_covariance'):
        resolved[key] = _as_bool(resolved[key], key)
    if not (resolved['fit_radius'] > 0 and resolved['sigma_min'] > 0 and resolved['sigma_max'] > resolved['sigma_min']):
        raise ValueError('Invalid fitting radius or width bounds')
    if resolved['n_peaks'] < 1:
        raise ValueError('n_peaks must be positive')
    return resolved


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

    def subtract(self, amp, center, sigma, radius_factor=5.0, cholesky=None):
        # Update a compact region. sigma is the longest principal width.
        xyz, _, idx = self.local(center, radius_factor * float(sigma))
        if cholesky is None:
            component = gaussian(xyz, amp, center, sigma)
        else:
            component = amp * np.exp(-0.5 * mahalanobis2(xyz, center, cholesky))
        self.residual[tuple(idx.T)] -= component

    def save_map(self, path):
        # Preserve the CCP4 extent, starts, axes, and crystallographic frame.
        a = np.array(self.ccp4.grid, copy=False)
        a[:] = np.where(np.isfinite(self.residual), self.residual, 0).astype(np.float32)
        self.ccp4.write_ccp4_map(str(path))


def gaussian(xyz, amp, center, sigma):
    return amp * np.exp(-0.5 * np.sum(((xyz - center)/sigma)**2, axis=-1))


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


def shape_record(cholesky):
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
        subtraction_sigma=float(widths[0]),
    )
    for axis, label in enumerate(('v1', 'v2', 'v3')):
        record[f'{label}_x'] = float(axes[0, axis])
        record[f'{label}_y'] = float(axes[1, axis])
        record[f'{label}_z'] = float(axes[2, axis])
    return record


def model_density(xyz, parameters, seed, fix_centers, full_covariance):
    """Evaluate one Gaussian. Returns density, center, and Cholesky factor or None."""
    amplitude = parameters[0]
    cursor = 1
    if fix_centers:
        center = np.asarray(seed, dtype=float)
    else:
        center = parameters[cursor:cursor + 3]
        cursor += 3
    if full_covariance:
        factor = cholesky_matrix(parameters[cursor:cursor + 6])
        cursor += 6
        signal = amplitude * np.exp(-0.5 * mahalanobis2(xyz, center, factor))
    else:
        width = parameters[cursor]
        cursor += 1
        factor = None
        signal = gaussian(xyz, amplitude, center, width)
    return signal + parameters[cursor], center, factor


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


def fit_peak(dmap, seed, radius, sigma_bounds, shift_limit, max_nfev=300,
             fix_centers=False, full_covariance=False):
    xyz, y, _ = dmap.local(seed, radius)
    if len(y) < 35:
        raise RuntimeError('Insufficient local valid voxels')
    seed = np.asarray(seed, dtype=float)
    bg = float(np.quantile(y, .2))
    amp = max(float(np.max(y) - bg), 1e-3)
    smin, smax = sigma_bounds
    s0 = float(np.clip(radius / 3, smin, smax))
    initial = [amp]
    lower = [0.0]
    upper = [max(10 * amp, 1.0)]
    if not fix_centers:
        initial.extend(seed)
        lower.extend(seed - shift_limit)
        upper.extend(seed + shift_limit)
    if full_covariance:
        initial.extend([s0, 0.0, 0.0, s0, 0.0, s0])
        lower.extend([smin, -smax, -smax, smin, -smax, smin])
        upper.extend([smax, smax, smax, smax, smax, smax])
    else:
        initial.append(s0)
        lower.append(smin)
        upper.append(smax)
    initial.append(bg)
    lower.append(float(np.min(y) - np.ptp(y) - 1.))
    upper.append(float(np.max(y) + np.ptp(y) + 1.))
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    initial = np.clip(np.asarray(initial, dtype=float), lower + 1e-8, upper - 1e-8)
    nfev = max(max_nfev, 900) if full_covariance else max_nfev

    def predicted(parameters):
        values, _, _ = model_density(xyz, parameters, seed, fix_centers, full_covariance)
        return values

    res = least_squares(lambda p: predicted(p) - y, initial, bounds=(lower, upper),
                        method='trf', loss='linear', x_scale='jac', max_nfev=nfev)
    values, center, factor = model_density(xyz, res.x, seed, fix_centers, full_covariance)
    if factor is None:
        width = res.x[1 if fix_centers else 4]
        factor = np.diag([width, width, width])
    shape = shape_record(factor)
    rss = float(np.dot(res.fun, res.fun))
    tot = float(np.dot(y - y.mean(), y - y.mean()))
    r2 = 1 - rss / tot if tot > 0 else float('nan')
    shape.update(dict(
        amp=float(res.x[0]),
        xyz=np.asarray(center, dtype=float),
        sigma=shape.pop('subtraction_sigma'),
        bg=float(res.x[-1]),
        r2=float(r2),
        rmse=float(np.sqrt(rss / len(y))),
        n=int(len(y)),
        success=bool(res.success),
        n_at_bound=int(np.count_nonzero(res.active_mask)),
        cholesky=factor,
        integrated_gaussian=float(res.x[0] * (2 * np.pi) ** 1.5 * shape.pop('integrated_scale')),
    ))
    return shape


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
    cfg = resolve_config(cfg)
    dm = DensityMap(cfg['map'])
    prefix = Path(cfg['output_prefix'])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    valid = np.isfinite(dm.original)
    reference_noise = robust_noise(dm.original)
    threshold = cfg['threshold_sigma'] * reference_noise
    excluded = np.zeros(dm.original.shape, bool)
    history = []
    if cfg['save_steps']:
        dm.save_map(f'{prefix}_step00.map')
    if verbose:
        mode = 'full covariance' if cfg['full_covariance'] else 'isotropic'
        where = 'fixed at detected maxima' if cfg['fix_centers'] else 'refined'
        print(f'Grid: {dm.original.shape}; robust map sigma={reference_noise:.5g}; '
              f'positive peak threshold={threshold:.5g}')
        print(f'Model: {mode}; centers {where}')
    for k in range(1, cfg['n_peaks'] + 1):
        peak = choose_peak(dm.residual, valid, excluded)
        if peak is None:
            break
        idx, height = peak
        if height < threshold:
            if verbose:
                print(f'Stopping: strongest remaining eligible peak {height:.4f} below threshold.')
            break
        seed = dm.to_cart(idx)
        try:
            fitted = fit_peak(
                dm, seed, cfg['fit_radius'], (cfg['sigma_min'], cfg['sigma_max']), cfg['max_shift'],
                fix_centers=cfg['fix_centers'], full_covariance=cfg['full_covariance'],
            )
        except (RuntimeError, ValueError) as exc:
            if verbose:
                print(f'Skipping seed {idx.tolist()}: {exc}')
            fitted = None
        if fitted is not None and fitted['success'] and fitted['amp'] > 0:
            center = fitted['xyz']
            # Gaussian subtraction only; retaining local background prevents map-wide drift.
            dm.subtract(fitted['amp'], center, fitted['sigma'], cholesky=fitted['cholesky'])
            row = dict(
                iteration=k,
                seed_x_A=float(seed[0]), seed_y_A=float(seed[1]), seed_z_A=float(seed[2]),
                x_A=float(center[0]), y_A=float(center[1]), z_A=float(center[2]),
                peak_before=height, amplitude=fitted['amp'], background=fitted['bg'],
                center_fixed=bool(cfg['fix_centers']), full_covariance=bool(cfg['full_covariance']),
                R2=fitted['r2'], RMSE=fitted['rmse'], n_voxels=fitted['n'],
                n_at_bound=fitted['n_at_bound'], success=fitted['success'],
                integrated_gaussian=fitted['integrated_gaussian'],
            )
            for key, value in fitted.items():
                if key in ('amp', 'xyz', 'sigma', 'bg', 'r2', 'rmse', 'n', 'success', 'n_at_bound',
                           'cholesky', 'integrated_gaussian'):
                    continue
                row[key] = value
            history.append(row)
            if verbose:
                print(f'#{k}: center=({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}) Å '
                      f'amp={fitted["amp"]:.3f}, principal σ='
                      f'({fitted["sigma_1_A"]:.3f}, {fitted["sigma_2_A"]:.3f}, {fitted["sigma_3_A"]:.3f}) Å, '
                      f'R2={fitted["r2"]:.3f}')
            if cfg['save_steps']:
                dm.save_map(f'{prefix}_step{k:02d}.map')
        marker = fitted['xyz'] if fitted is not None else seed
        _, _, around = dm.local(marker, cfg['min_separation'], dm.original)
        excluded[tuple(around.T)] = True
    df = pd.DataFrame(history)
    centers_csv = f'{prefix}_centers.csv'
    centers_pdb = f'{prefix}_centers.pdb'
    residual_map = f'{prefix}_residual.map'
    diagnostics = f'{prefix}_diagnostics.png'
    df.to_csv(centers_csv, index=False)
    coordinates = df[['x_A', 'y_A', 'z_A']].to_numpy() if len(df) else []
    write_centers_pdb(centers_pdb, coordinates)
    dm.save_map(residual_map)
    if len(df):
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
        axs[0].plot(df.iteration, df.amplitude, 'o-')
        axs[0].set(xlabel='Subtraction iteration', ylabel='Fitted Gaussian amplitude', title='Sequential peak amplitudes')
        axs[1].plot(df.iteration, df.R2, 'o-')
        axs[1].set(xlabel='Subtraction iteration', ylabel='Local fit R²', title='Fit quality')
        axs[1].set_ylim(min(-.1, float(df.R2.min()) - .05), 1.05)
        fig.tight_layout()
        fig.savefig(diagnostics, dpi=180)
        plt.close(fig)
    if verbose:
        print(f'Saved {len(df)} components; results: {centers_csv}, {centers_pdb}, and {residual_map}')
    df.attrs['centers_csv'] = centers_csv
    df.attrs['centers_pdb'] = centers_pdb
    df.attrs['residual_map'] = residual_map
    df.attrs['diagnostics'] = diagnostics
    return df


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('config', help='YAML configuration file')
    args = parser.parse_args(argv)
    return run(load_config(args.config))

if __name__ == '__main__':
    main()

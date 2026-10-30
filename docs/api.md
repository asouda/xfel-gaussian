# API

The public entry points are the two modules. Names that start with an underscore are internal. Both modules accept a filesystem path through `load_config` and a mapping through `run`.

```python
from xfel_gaussian.multi_peaks_gaussian import load_config, run as run_multi
from xfel_gaussian.sequential_peaks_gaussian import run as run_sequential

cfg = load_config("config/peaks.yaml")
multi = run_multi(cfg)          # dict of tables and paths
history = run_sequential(cfg)   # DataFrame, paths on history.attrs
```

Command-line wrappers take the same YAML path:

```bash
python -m xfel_gaussian.multi_peaks_gaussian config/peaks.yaml
python -m xfel_gaussian.sequential_peaks_gaussian config/peaks.yaml
```

`verbose` defaults to true and prints the ROI, the model comparison, and the centers. Pass `verbose=False` to keep only the return value. The console script name `xfel-gaussian` is declared in `pyproject.toml` as `xfel_gaussian:main`. The package initializer does not define `main`, so the working interfaces are the two modules above.

## Configuration

`load_config(path)` reads a YAML mapping and stores the resolved path under `_config_path`. `run` uses that path to resolve `map` and `output_prefix`.

`resolve_config(cfg)` copies the defaults, reads either the module's section or a flat mapping, checks the keys, and returns a plain dictionary. Relative paths in the returned dictionary are absolute. Call it directly when a later step needs the resolved map path without fitting.

`section_config(cfg, name)` returns the `multi_peaks` or `sequential_peaks` mapping, including the shared `map` key. A file that contains only the other section raises `KeyError`.

## Joint module

`xfel_gaussian.multi_peaks_gaussian`

### `run(cfg, verbose=True)`

Fits the joint model and writes the [output files](outputs.md). Returns the dictionary documented there. `cfg` is a mapping from `load_config` or an equivalent dictionary.

### Geometry and design matrix

`read_ccp4(path)` returns `arr, axes, starts, sampling, orth`. `axes` is gemmi's `axis_positions()` tuple. `starts` is header words 5–7. `sampling` is header words 8–10. `orth` is the orthogonalization matrix.

`to_cart(indices, axes, starts, sampling, orth)` maps stored indices to Cartesian Å, including a batch of indices of shape `(..., 3)`.

`extract_roi(rho, axes, starts, sampling, orth, center, radius)` returns `xyz`, `y`, and the shortest grid-step length for the finite voxels inside the sphere.

`estimate_spacing(xyz)` is the median nearest-neighbor distance, used when `analyze_arrays` is called without a spacing. The CLI passes the shortest grid step instead.

`local_maximum_mask(xyz, y, spacing)` is true where a sample is at least as high as every sample within `1.8 * spacing`.

`seeds(xyz, y, max_count, separation, min_amplitude=0)` is the greedy high-to-low selection with a minimum separation.

### Model

`component_stride(fix_centers, full_covariance)` is the number of free parameters per component, excluding the shared background.

`cholesky_matrix(values)` builds \(L\) from `[L00, L10, L20, L11, L21, L22]`.

`mahalanobis2(xyz, center, cholesky)` returns \(\lVert L^{-1}(\mathbf r-\mathbf r_0)\rVert^2\) for each row of `xyz`.

`prediction(parameters, xyz, k, centers, fix_centers, full_covariance)` evaluates the joint density, including the background stored as the last parameter. `centers` is used only when centers are fixed.

`principal_frame(cholesky)` returns descending principal widths, the right-handed eigenframe, and \(\Sigma\).

`one_fit(xyz, y, k, starts, width_bounds=(0.3, 2.0), n_restarts=3, rng=None, fix_centers=False, full_covariance=False)` runs the bounded least squares for a single \(K\). `starts` is a sequence of seed indices into `xyz`. The return value is the `scipy.optimize.OptimizeResult` with the lowest residual sum of squares.

`analyze_arrays(xyz, y, max_components=4, min_sep=0.6, width_bounds=(0.3, 2.0), restarts=3, fix_centers=False, full_covariance=False, spacing=None)` fits \(K = 1, \ldots\) and returns `comparison`, `components`, `pairs`, and the selected `OptimizeResult`.

### Figures, benchmark, and PDB

`plot_fit(xyz, y, comp, fit, path)` writes the two-panel diagnostic figure. `comp` is the component DataFrame of the selected model.

`residual_map_path(map_path)` is the input path with `_residual` inserted before the suffix.

`gaussian_residual(rho, axes, starts, sampling, orth, components, radius_factor=5.0)` returns the map minus each component's Gaussian. The shared background is not removed. Each component is subtracted inside `radius_factor` times `sigma_1_A`.

`write_map_like(source, values, destination)` writes `values` into a copy of the CCP4 header and grid at `source`. Non-finite samples become zero.

`synthetic_benchmark(output, distances=(0.7, 1.0, 1.4, 1.8, 2.2, 2.8), width=0.8, noise=0.08, fix_centers=False, full_covariance=False)` writes and returns the two-peak recovery table.

`write_centers_pdb(path, coordinates, serial_start=43042)` writes the HETATM records described in [Outputs](outputs.md).

`main(argv=None)` parses a single config-path argument and calls `run`.

## Sequential module

`xfel_gaussian.sequential_peaks_gaussian`

### `DensityMap(path)`

Holds the original grid, a residual copy, and the Cartesian mapping.

- `to_cart(idx)` and `to_index(xyz)` convert stored indices and Cartesian Å.
- `local(center, radius, values=None)` returns Cartesian coordinates, densities, and stored indices inside the sphere. `values` defaults to the residual.
- `subtract(amp, center, sigma, radius_factor=5.0, cholesky=None)` removes a Gaussian from the residual inside `radius_factor * sigma`. A Cholesky factor selects the ellipsoidal form.
- `save_map(path)` writes the residual into the original CCP4 header and grid.

`voxel_A3` is computed on the grid step and is not used by the fit.

### One peak

`robust_noise(y)` returns \(1.4826\) times the median absolute deviation of the finite values, and at least \(10^{-7}\).

`choose_peak(rho, original_mask, exclusion, neighborhood_size=3)` returns `(index, height)` of the strongest eligible local maximum, or `None`. `neighborhood_size` is the odd side length of the cubic window.

`spline_maximum_index(rho, idx, neighborhood_size, spline_order)` returns the continuous stored index of the maximum of a spline interpolant of the neighborhood around `idx`. Order 3 is cubic. Orders below 2, a window that is not entirely finite, and a spline that cannot be built all return `idx`.

`gaussian(xyz, amp, center, sigma)` is the isotropic component without a background.

`model_density(xyz, parameters, seed, fix_centers, full_covariance)` returns the predicted density, the center, and the Cholesky factor or `None` for an isotropic parameter vector.

`fit_peak(dmap, seed, radius, sigma_bounds, shift_limit, max_nfev=300, fix_centers=False, full_covariance=False)` fits one sphere and returns a dictionary of shape descriptors plus `amp`, `xyz`, `sigma` (longest principal width), `bg`, `r2`, `rmse`, `n`, `success`, `n_at_bound`, `cholesky`, and `integrated_gaussian`.

`shape_record(cholesky)` is the covariance, width, and axis dictionary used to build that result. It includes `subtraction_sigma`, the longest principal width, which `fit_peak` then renames to `sigma`.

`cholesky_matrix`, `mahalanobis2`, `principal_frame`, and `write_centers_pdb` match the joint module.

### `run(cfg, verbose=True)`

Runs the subtraction loop and returns the centers DataFrame. Output paths are on `df.attrs` under `centers_csv`, `centers_pdb`, `residual_map`, `prefix_residual_map`, and `diagnostics`. `residual_map` is the file beside the input.

`main(argv=None)` parses a single config-path argument and calls `run`.

## Defaults

Library defaults, before `config/peaks.yaml` overrides them:

| Joint key | Default | Sequential key | Default |
| --- | --- | --- | --- |
| `center` | `None` | `n_peaks` | `6` |
| `radius` | `4.5` | `threshold_sigma` | `2.0` |
| `max_components` | `4` | `fit_radius` | `2.2` |
| `min_separation` | `0.6` | `sigma_min` | `0.35` |
| `min_width` | `0.35` | `sigma_max` | `1.8` |
| `max_width` | `1.8` | `max_shift` | `1.2` |
| `restarts` | `3` | `min_separation` | `0.7` |
| `output_prefix` | `multi_mn` | `output_prefix` | `sequential_peaks` |
| `synthetic` | `False` | `save_steps` | `False` |
| `fix_centers` | `False` | `fix_centers` | `False` |
| `full_covariance` | `False` | `full_covariance` | `False` |
| | | `neighborhood_size` | `3` |
| | | `spline_seed` | `False` |
| | | `spline_order` | `3` |

The example `config/peaks.yaml` turns fixed centers and the full covariance on for both workflows, raises the joint component limit to 9, and turns sequential spline seeding on (`spline_seed: true`, order 3, neighborhood side 3). The file, not this table, is the configuration of a run.

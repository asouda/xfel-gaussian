# Configuration

This folder holds the YAML file read by both fitting workflows. The worked example is `peaks.yaml`. One file can configure both programs: a shared `map` entry, a `multi_peaks` section, and a `sequential_peaks` section. Either program also accepts a flat file that contains only its own keys.

Paths are resolved relative to the YAML file, not the shell's working directory. In `peaks.yaml`, `../data/...` and `../results/...` therefore point at the `data` and `results` directories next to `config`. Absolute paths are left unchanged. The parent directory of `output_prefix` is created if needed.

`map` is a CCP4 `.map` or `.ccp4` file. It is required for the sequential workflow, and for the joint workflow unless `synthetic` is true. The file must exist when the program starts. The grid is read as stored: the sample is not expanded to a full unit cell and is not treated as periodic.

Unknown keys are an error. Boolean values may be real YAML booleans or the strings `true`/`yes`/`1` and `false`/`no`/`0`. A key set to `null` keeps that key's default when the key is one of the defaults listed below.

## `multi_peaks`

Joint fit of K = 1, ..., Kmax positive Gaussians plus one constant background inside a Cartesian sphere. The BIC minimum selects K.

| Key | Default | Meaning |
| --- | --- | --- |
| `center` | `null` | ROI center in Cartesian Å, three numbers. `null` uses the strongest positive voxel in the map. |
| `radius` | `4.5` | Sphere radius in Å. Must be positive. |
| `max_components` | `4` | Largest K to fit. Must be at least 1. BIC may select a smaller model. |
| `min_separation` | `0.6` | Minimum distance in Å between seed points. Refined centers are not held apart by this value. |
| `min_width` | `0.35` | Lower bound on an isotropic width, or on a Cholesky diagonal, in Å. |
| `max_width` | `1.8` | Upper bound on an isotropic width, or on the absolute value of a Cholesky entry, in Å. Must exceed `min_width`. |
| `restarts` | `3` | Bounded least-squares starts. The first start is unperturbed. Later starts jitter a free center and the initial width. Must be at least 1. |
| `output_prefix` | `multi_mn` | Path prefix for the CSV, PDB, and diagnostic figure. |
| `synthetic` | `false` | If true, ignore the map and run the built-in two-peak recovery. `fix_centers`, `full_covariance`, and `output_prefix` still apply. The benchmark uses its own grid, width bounds `(0.4, 1.4)`, and two restarts. |
| `fix_centers` | `false` | Hold each center on its detected seed. Seeds are then local maxima. If false, centers are free inside the ROI bounding box. |
| `full_covariance` | `false` | Fit a rotated ellipsoid through six Cholesky entries. If false, each component has one isotropic width. |

`peaks.yaml` overrides several of these defaults (`max_components: 9`, `restarts: 2`, `fix_centers: true`, `full_covariance: true`, and an output prefix under `results/`).

## `sequential_peaks`

Fit one positive Gaussian and a local baseline at the strongest remaining local maximum, subtract only the Gaussian, and repeat.

| Key | Default | Meaning |
| --- | --- | --- |
| `output_prefix` | `sequential_peaks` | Path prefix for the CSV, PDB, residual map, and diagnostic figure. |
| `n_peaks` | `6` | Maximum number of subtraction attempts. The loop also stops when no eligible peak remains above threshold. |
| `threshold_sigma` | `2.0` | Stop when the strongest remaining peak is below this multiple of the robust scale of the original map. |
| `fit_radius` | `2.2` | Radius in Å of the sphere used for each one-peak fit. Must be positive. |
| `sigma_min` | `0.35` | Lower width bound in Å, applied as `min_width` is in the joint fit. |
| `sigma_max` | `1.8` | Upper width bound in Å. Must exceed `sigma_min`. |
| `max_shift` | `1.2` | If centers are refined, each center stays within this distance in Å of its seed. Ignored when `fix_centers` is true. |
| `min_separation` | `0.7` | Radius in Å of the exclusion sphere placed around a seed after each attempt, so the next peak is a different site. |
| `save_steps` | `false` | Write a CCP4 map of the residual before subtraction and after each accepted peak. |
| `fix_centers` | `false` | Hold the center at the detected voxel. If false, the center may move by up to `max_shift`. |
| `full_covariance` | `false` | Same shape choice as in the joint fit. |

## Which section is read

`multi_peaks_gaussian` reads `multi_peaks` and the shared `map`. `sequential_peaks_gaussian` reads `sequential_peaks` and the shared `map`. A file that contains one of those sections and not the other is an error for the missing workflow. A file that contains neither section is treated as a flat configuration for whichever program is launched.

The commands, the output files, and the meaning of BIC, Cholesky widths, and fixed centers are documented in the MkDocs site (`uv run mkdocs serve` from the repository root).

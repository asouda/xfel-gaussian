# Outputs

Paths are `{output_prefix}` plus a suffix. The prefix is resolved next to the YAML file and may contain directories. The example configuration writes `results/multi_peaks_*` and `results/sequential_peaks_*` from the repository root. `results/` is gitignored.

Centers are Cartesian ångströms in the map's orthogonal frame. Density, amplitude, and background are in the map's stored units. Widths are in Å, covariance entries are in Å\(^2\), and `integrated_gaussian` is in map-unit·Å\(^3\).

!!! warning "Labels are not atoms"
    The center tables are model-derived. The PDB atom names `MN1`–`MN4` and then `O1`, `O2`, ... exist so the coordinates can be opened beside a structure. They are not an assignment of manganese or oxygen.

## Residual maps

After the selected peaks have been fitted, both workflows write the input density minus those Gaussians. The joint shared background and the sequential local baselines stay in the map. Each Gaussian is removed inside a sphere of radius five times its longest principal width. Voxels outside every such sphere are unchanged. Non-finite samples are written as zero. The CCP4 grid, starts, axis order, and cell are copied from the input.

The same array is written twice:

- Beside the input. `directory/name.map` becomes `directory/name_residual.map`. A `.ccp4` input keeps that suffix. Both workflows use this path, so the file matches the workflow that finished last.
- `{prefix}_residual.map`, next to that workflow's other products. The other workflow does not replace this copy.

An empty sequential search still writes both files, and they match the input. `synthetic: true` does not read a map and does not write a residual.

## Joint fit

A map run writes the tables and figure below, plus the residual maps above. `synthetic: true` writes only the synthetic CSV and skips the rest.

### `{prefix}_models.csv`

One row per component count \(K = 1, \ldots, K_{\max}\). The selected model is the row with the smallest `BIC`. An exact tie keeps the smaller \(K\).

| Column | Meaning |
| --- | --- |
| `K` | Number of Gaussians. |
| `N` | Number of finite voxels in the sphere. |
| `n_parameters` | Free parameters, including the shared background. |
| `RSS` | Residual sum of squares. |
| `RMSE` | \(\sqrt{\mathrm{RSS}/N}\). |
| `R2` | Coefficient of determination inside the sphere. NaN if the densities are constant. |
| `BIC` | \(N\ln(\mathrm{RSS}/N) + p\ln N\), up to a constant shared by every row. \(\mathrm{RSS}/N\) is floored at \(10^{-14}\). |
| `AIC` | The same residual term with penalty \(2p\). Reported, not used for selection. |
| `success` | Solver success flag for the start that won. A false flag does not drop the row. |
| `n_at_bound` | Parameters sitting on a bound. |
| `fix_centers` | Whether centers were held fixed. |
| `full_covariance` | Whether the shape was a full covariance. |

### `{prefix}_centers.csv`

One row per component of the BIC-selected model, in seed order. The shared `background` is repeated on every row. Isotropic fits fill the principal-axis columns with three equal widths, a diagonal covariance, and condition number 1.

| Column | Meaning |
| --- | --- |
| `component` | 1-based index in seed order. |
| `amplitude` | \(A\), map units, non-negative. |
| `x_A`, `y_A`, `z_A` | Center in Å. |
| `center_se_x_A`, `center_se_y_A`, `center_se_z_A` | Approximate standard errors in Å from the Jacobian. NaN when the center is fixed. |
| `center_fixed`, `full_covariance` | The switches used for this fit. |
| `background` | Shared baseline, map units. |
| `sigma_A` | Equivalent width \((s_1 s_2 s_3)^{1/3}\). |
| `sigma_1_A`, `sigma_2_A`, `sigma_3_A` | Principal widths in Å, longest first. |
| `FWHM_A`, `fwhm_1_A`, `fwhm_2_A`, `fwhm_3_A` | \(2\sqrt{2\ln 2}\) times the equivalent width and the three principal widths. |
| `cov_xx`, `cov_yy`, `cov_zz`, `cov_xy`, `cov_xz`, `cov_yz` | Unique entries of \(\Sigma\). |
| `L00`, `L10`, `L20`, `L11`, `L21`, `L22` | Cholesky entries that were optimized, or \(\sigma\) on the diagonal when the component is isotropic. |
| `v1_x`, `v1_y`, `v1_z`, and `v2_*`, `v3_*` | Principal-axis direction cosines, right-handed, longest axis first. |
| `condition_number` | \(\lambda_{\max}/\lambda_{\min}\) of \(\Sigma\). |
| `integrated_gaussian` | \(A(2\pi)^{3/2}\ell_{00}\ell_{11}\ell_{22}\). |

### `{prefix}_pairs.csv`

One row per unordered pair of selected centers. The file has no pair rows when BIC selects a single component.

| Column | Meaning |
| --- | --- |
| `i`, `j` | 1-based component indices. |
| `distance_A` | Euclidean separation in Å. |
| `d_over_sum_sigma` | Distance divided by the sum of the two equivalent widths. |
| `d_over_mean_fwhm` | Distance divided by the average of the two equivalent FWHMs. |

### `{prefix}_diagnostics.png`

Left: observed voxel density against the joint model. Right: residual against the Cartesian axis along which the ROI points span the largest range. The horizontal axis of the residual panel is that coordinate, not a principal axis.

### `{prefix}_synthetic.csv`

Written only when `synthetic: true`.

| Column | Meaning |
| --- | --- |
| `true_separation_A` | Separation of the two simulated centers. |
| `selected_K` | 1 or 2, the BIC choice. |
| `recovered_separation_A` | Distance between fitted centers when \(K = 2\). NaN when \(K = 1\). |
| `delta_BIC_1_minus_2` | \(\mathrm{BIC}(K=1) - \mathrm{BIC}(K=2)\). Positive when the two-component model wins. |
| `recovered_two` | True when BIC selects \(K = 2\). |
| `noise_sd` | Standard deviation of the added noise. |
| `truth_sigma_A` | True isotropic width, \(0.8\) Å. |
| `fix_centers`, `full_covariance` | Switches passed into the benchmark. |

## Sequential fit

### `{prefix}_centers.csv`

One row per accepted peak, in subtraction order. A skipped or non-positive fit does not appear, and `iteration` keeps the attempt number, so the column can skip. An empty search writes a CSV with no data rows and still writes the other files.

| Column | Meaning |
| --- | --- |
| `iteration` | Search attempt, counting from 1. Skipped sites are not rows, so the numbers can jump. |
| `seed_x_A`, `seed_y_A`, `seed_z_A` | Cartesian position of the seed. The detected voxel when `spline_seed` is false, and the spline maximum when `spline_seed` is true. |
| `x_A`, `y_A`, `z_A` | Fitted center. Equal to the seed when `fix_centers` is true. |
| `peak_before` | Residual density at the detected voxel before this fit. The threshold uses this value. |
| `amplitude` | Fitted \(A\). |
| `background` | Baseline of this peak only. It is not subtracted. |
| `center_fixed`, `full_covariance` | The switches used. |
| `spline_seed`, `spline_order`, `neighborhood_size` | How the seed was placed. `spline_order` is used when `spline_seed` is true. Orders below 2 leave the seed on the detected voxel. |
| `R2`, `RMSE` | Fit quality inside that peak's sphere. |
| `n_voxels` | Finite voxels in the sphere. |
| `n_at_bound` | Parameters sitting on a bound. |
| `success` | Solver success flag. Accepted rows are successes with positive amplitude. |
| `integrated_gaussian` | Same integral as in the joint table. |
| shape columns | `sigma_A`, `sigma_1_A` … `sigma_3_A`, `FWHM_A`, `fwhm_1_A` … `fwhm_3_A`, the six `cov_*` entries, the six `L*` entries, `condition_number`, and `v1_*`, `v2_*`, `v3_*`, with the same meanings as above. |

The longest principal width is the radius scale of the subtraction. It is not given a separate column; it is `sigma_1_A`.

### `{prefix}_centers.pdb`

See below. The sequential and joint writers use the same record layout.

### `{prefix}_stepNN.map`

Written only when `save_steps` is true. `step00` is the residual before any subtraction. `step01` and later follow each accepted peak. Failed attempts do not produce a step.

### `{prefix}_diagnostics.png`

Written when at least one peak was accepted. Left: fitted amplitude against iteration. Right: local \(R^2\) against iteration.

## PDB records

Both workflows write one `HETATM` record per center, 80 columns, serial numbers starting at 43042. The file ends with an `END` record.

- Columns 13–16 are the atom name: `MN1`, `MN2`, `MN3`, `MN4` for the first four centers, then `O1`, `O2`, and so on.
- The four characters `BOEY` are written starting at column 17, so the alternate-location column is `B` and the residue name is `OEY`. Chain is `A` and the residue sequence is 601.
- Occupancy is 1.00 and the temperature factor is 27.27. Both are constants.
- Columns 31–54 are the Cartesian center in Å, in the standard orthogonal \(x\), \(y\), \(z\) fields.

Open the file to see the centers. Do not read the residue name or the atom name as chemistry.

## Notebook products

The standalone notebooks write their own CSV files in the working directory, typically `notebooks/` if that is where the kernel was started.

- `6w1v_fitted_gaussian_peaks.csv` from the axis-aligned tutorial: sign, grid index, Cartesian center, signed amplitude, three Cartesian widths, FWHM, background, signed integral, \(R^2\), RMSE, and the convergence flag.
- `6w1v_fitted_full_covariance_peaks.csv` from the full-covariance tutorial: the detected extremum, the fixed center, signed amplitude, background, principal widths and FWHMs, the six covariance entries, the Cholesky factor, direction cosines, condition number, signed integral, and fit diagnostics.

Those filenames are set inside the notebooks and do not use `output_prefix`. The module notebook does not invent a third table. It displays the files the two library `run` functions write.

## Return values

The joint `run` returns a dictionary. A map run contains `comparison`, `components`, and `pairs` as DataFrames, plus `diagnostics`, `models_csv`, `centers_csv`, `centers_pdb`, `pairs_csv`, `residual_map`, `prefix_residual_map`, `roi_center`, and `background`. `residual_map` is the file beside the input. `prefix_residual_map` is `{prefix}_residual.map`. A synthetic run contains `synthetic` and `synthetic_csv` only.

The sequential `run` returns the centers DataFrame. Paths are stored on `df.attrs` under `centers_csv`, `centers_pdb`, `residual_map`, `prefix_residual_map`, and `diagnostics`. `residual_map` is the file beside the input and `prefix_residual_map` is `{prefix}_residual.map`. The diagnostics path is set even when the figure was not written because no peak was accepted.

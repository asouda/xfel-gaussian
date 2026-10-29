# Notebooks

Tutorial notebooks for the Gaussian peak models. They are meant to be read and run from this folder. The two standalone notebooks reimplement a single-peak fit so the model, the CCP4 geometry, and the plots are visible in one file. They do not import `xfel_gaussian`. The third notebook calls the installed package on `config/peaks.yaml`.

Maps are not in the repository. For the standalone notebooks, set `MAP_FILE` to a `.map` or `.ccp4` path, or place a map in `data/` at the repository root (the notebooks look for `../data/*.map` and then `../data/*.ccp4`). The module notebook uses the `map` path in the YAML file. A local crop is a better first run than a full unit-cell grid.

Start Jupyter with the project environment from the repository root:

```bash
uv sync
uv run --with jupyterlab jupyter lab notebooks
```

`ipykernel` is already a project dependency, so an editor that uses the project virtual environment can open the notebooks without JupyterLab.

Difference density is signed, spatially correlated, and affected by the map coefficients, sharpening, truncation, and the omitted model. A fitted width is a descriptor of the map feature. It is not a calibrated atomic displacement parameter. The Gaussian integral is not an electron count. The strongest peak is not a named atom until its Cartesian position is compared with a structural model in the same frame.

## `xfel_gaussian_fitting_tutorial.ipynb`

Axis-aligned anisotropic Gaussian for one extremum at a time, using 6W1V omit-OEC difference density as the worked example.

```text
Δρ(r) = b + A exp(-1/2 * ((x-x0)^2/sx^2 + (y-y0)^2/sy^2 + (z-z0)^2/sz^2))
```

The eight parameters are a signed amplitude, a center, three positive widths along Cartesian x, y, and z, and a constant baseline. Positive and negative extrema are both detected. A candidate must be a grid extremum among its 26 neighbors, pass a threshold in units of the map-wide standard deviation, and lie at least `MIN_DISTANCE` Å from an already selected extremum. The center is refined inside the fit sphere. The notebook writes `6w1v_fitted_gaussian_peaks.csv` and shows a stored-array section through each fitted peak.

This is the special case of a diagonal covariance whose axes stay on x, y, and z. A tilted feature needs the full-covariance notebook or the library's `full_covariance` switch.

## `xfel_gaussian_full_covariance_tutorial.ipynb`

The same detection workflow with a fully rotated ellipsoid. The center is held at the detected extremum. The quadratic form is a full covariance Σ = L Lᵀ, with L the Cholesky factor that is actually optimized. Principal widths and axis directions are reported after the fit by diagonalizing Σ.

The notebook first recovers a synthetic ellipsoid whose principal axes are not aligned with x, y, and z, then fits the map. Sections show the 1σ and 2σ elliptical cuts of the ellipsoid (r − r0)ᵀ Σ⁻¹ (r − r0) = n². Results are written to `6w1v_fitted_full_covariance_peaks.csv`.

## `xfel_gaussian_module_calculation.ipynb`

Runs both library workflows from one configuration. It loads `../config/peaks.yaml` (or `config/peaks.yaml` if that relative path is the one that exists), calls `multi.run` and `sequential.run`, and displays the BIC table, the component table, and the diagnostic figures those functions write.

For each selected component it also draws the section layout of the full-covariance tutorial: stored-array sections of the observed map, the model, and the residual, then Cartesian xy, xz, and yz planes with the 1σ and 2σ ellipses. The joint-model panels use every component plus the shared background. Each sequential panel uses that peak's Gaussian plus the local background fitted with it, drawn on the original map.

Library fits constrain the amplitude to be non-negative. The signed single-peak models live in the other two notebooks.

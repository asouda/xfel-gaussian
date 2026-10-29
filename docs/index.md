# XFEL Gaussian

`xfel-gaussian` fits three-dimensional Gaussians to positive peaks in a CCP4 map. The intended maps are omit and difference electron-density maps, where an unmodeled atom appears as a localized positive feature. The same code reads any CCP4 grid whose peaks of interest are positive.

Two workflows share `config/peaks.yaml`:

- The **joint fit** places one to \(K_{\max}\) positive Gaussians, plus a single constant background, inside a Cartesian sphere. It fits every component count and keeps the model with the lowest Bayesian information criterion (BIC).
- The **sequential fit** finds the strongest remaining local maximum, fits one Gaussian and a local baseline there, subtracts only the Gaussian, and repeats.

Each Gaussian is either isotropic or a fully rotated ellipsoid. The ellipsoid is parameterized by a Cholesky factor so the covariance stays positive definite. Centers are either free parameters or held on the detected local maxima.

```bash
uv run python -m xfel_gaussian.multi_peaks_gaussian config/peaks.yaml
uv run python -m xfel_gaussian.sequential_peaks_gaussian config/peaks.yaml
```

!!! warning "What a fit is"
    A center, a width, and a Gaussian integral describe the map feature inside the fit window. They are not a chemical assignment, an atomic displacement parameter, or an electron count. Atom names in the PDB file exist so the coordinates can be opened beside a structure.

## Where to read next

| Chapter | What it covers |
| --- | --- |
| [Installation](installation.md) | Python 3.14, `uv sync`, and building this site |
| [Configuration](configuration.md) | The YAML file in `config/` |
| [Workflows](workflows.md) | The joint and sequential algorithms, and how the notebooks differ |
| [Theory](theory.md) | The Gaussian models, CCP4 geometry, BIC, and the limits of the fit |
| [Outputs](outputs.md) | CSV columns, the PDB labels, residual maps, and figures |
| [Notebooks](notebooks.md) | The three tutorials in `notebooks/` |
| [API](api.md) | Functions a caller uses from Python |

The package author recorded in the project metadata is Alexander Soudackov. The source is released under the GNU General Public License v3.

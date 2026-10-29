# xfel-gaussian

Gaussian models for positive peaks in CCP4 omit and difference electron-density maps. Two workflows share one YAML file:

- **Joint fit** (`xfel_gaussian.multi_peaks_gaussian`): several positive Gaussians and one background in a spherical region of interest. The number of components is chosen by the Bayesian information criterion (BIC).
- **Sequential fit** (`xfel_gaussian.sequential_peaks_gaussian`): fit the strongest remaining local maximum, subtract that Gaussian from the map, and repeat.

Each component is either an isotropic Gaussian or a fully rotated ellipsoid. Centers are either refined or held at the detected local maxima.

The fitted centers describe map features. Names written into the PDB output are viewing labels, and the Gaussian integral is not a calibrated electron count. The [theory chapter](docs/theory.md) states what the parameters mean.

## Installation

The package requires Python 3.14. From a clone of this repository:

```bash
uv sync
```

`uv sync` creates `.venv` and installs the package plus its dependencies (`gemmi`, NumPy, SciPy, pandas, Matplotlib, PyYAML).

Maps are not included. `data/` is ignored by git. Put a CCP4 `.map` or `.ccp4` file on disk and point `map` in [`config/peaks.yaml`](config/peaks.yaml) at it. Paths in that file are relative to the file itself, so the example entry `../data/9z35_omitOEC_withCa_extended_10A.map` means `data/9z35_omitOEC_withCa_extended_10A.map` from the repository root.

## Usage

From the repository root:

```bash
uv run python -m xfel_gaussian.multi_peaks_gaussian config/peaks.yaml
uv run python -m xfel_gaussian.sequential_peaks_gaussian config/peaks.yaml
```

The same calls from Python:

```python
from xfel_gaussian.multi_peaks_gaussian import load_config, run as run_multi
from xfel_gaussian.sequential_peaks_gaussian import run as run_sequential

cfg = load_config("config/peaks.yaml")
multi = run_multi(cfg)
sequential = run_sequential(cfg)
```

`multi` is a dictionary of tables and output paths. `sequential` is a pandas `DataFrame` of accepted peaks. With the example configuration, files are written under `results/`. That directory is gitignored.

Set `synthetic: true` under `multi_peaks` to skip the map and run a two-peak recovery check. Details of every key are in [`config/README.md`](config/README.md).

## Documentation

The documentation site includes the theoretical background, both workflows, the output columns, and the Python API.

```bash
uv sync --group docs
uv run mkdocs serve
```

The site is then at <http://127.0.0.1:8000>. `uv run mkdocs build` writes a static site to `site/`.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/xfel_gaussian/` | Joint and sequential fitting modules |
| [`config/`](config/README.md) | Shared YAML configuration |
| [`notebooks/`](notebooks/README.md) | Tutorials for the library and for the related single-peak models |
| `docs/` | MkDocs sources |

## License

GNU General Public License v3. See [`LICENSE`](LICENSE).

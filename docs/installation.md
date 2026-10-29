# Installation

## Requirements

The package requires Python 3.14 or newer. Runtime dependencies are [gemmi](https://gemmi.readthedocs.io/) (CCP4 I/O and unit-cell orthogonalization), NumPy, SciPy, pandas, Matplotlib, and PyYAML. The repository pins them in `uv.lock` and declares them in `pyproject.toml`.

Install [uv](https://docs.astral.sh/uv/) if it is not already on the machine. The `.python-version` file in the repository root asks uv for Python 3.14.

## Install the package

Clone the repository and synchronize the environment from the repository root:

```bash
git clone git@github.com:asouda/xfel-gaussian.git
cd xfel-gaussian
uv sync
```

`uv sync` creates `.venv`, installs the locked dependencies, and installs `xfel-gaussian` itself from `src/`. Later commands use that environment through `uv run`, or through `.venv` after it is activated.

Maps are not part of the clone. The `data/` directory is listed in `.gitignore`. Copy a CCP4 `.map` or `.ccp4` file onto disk and set the `map` key in `config/peaks.yaml`. Paths inside the YAML file are relative to that file. The example value `../data/9z35_omitOEC_withCa_extended_10A.map` is `data/9z35_omitOEC_withCa_extended_10A.map` when the shell is at the repository root.

## Run a fit

```bash
uv run python -m xfel_gaussian.multi_peaks_gaussian config/peaks.yaml
uv run python -m xfel_gaussian.sequential_peaks_gaussian config/peaks.yaml
```

Each command reads its own section of the YAML file. The joint fit writes comparison and center tables under the `multi_peaks.output_prefix` path. The sequential fit writes centers, a residual CCP4 map, and a diagnostic figure under `sequential_peaks.output_prefix`. With the example file those prefixes land in `results/`, which is gitignored.

From Python the same entry points are `run` after `load_config`:

```python
from xfel_gaussian.multi_peaks_gaussian import load_config, run as run_multi
from xfel_gaussian.sequential_peaks_gaussian import run as run_sequential

cfg = load_config("config/peaks.yaml")
multi = run_multi(cfg)
sequential = run_sequential(cfg)
```

`load_config` records the YAML path. Relative map and output paths are then resolved next to the file. A dictionary passed to `run` without that bookkeeping resolves relative paths from the process working directory.

Set `synthetic: true` in the `multi_peaks` section to run the built-in two-peak recovery instead of reading a map. The sequential workflow always requires a map file.

The key reference is the [configuration](configuration.md) chapter. The shape of each result table is the [outputs](outputs.md) chapter.

## Notebooks

`ipykernel` is installed with the package. Launch a notebook UI against the same environment with:

```bash
uv run --with jupyterlab jupyter lab notebooks
```

The two standalone tutorials need a map in `data/` or a `MAP_FILE` path. `xfel_gaussian_module_calculation.ipynb` reads `config/peaks.yaml`. See [Notebooks](notebooks.md).

## Build this documentation

The documentation dependencies are a uv dependency group, so a plain `uv sync` does not install them.

```bash
uv sync --group docs
uv run mkdocs serve
```

The preview is served at <http://127.0.0.1:8000>. A static build is:

```bash
uv run mkdocs build
```

MkDocs writes HTML to `site/`. That directory is gitignored. `mkdocs.yml` sets `strict: true`, so a broken internal link fails the build. The configuration and notebook chapters are included from `config/README.md` and `notebooks/README.md`, which keeps the site and those folder descriptions on the same text.

## License

The package is distributed under the GNU General Public License v3. The full text is `LICENSE` in the repository root.

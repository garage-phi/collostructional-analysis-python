# Collostructional Analysis for Python

[![DOI](https://zenodo.org/badge/1138207212.svg)](https://doi.org/10.5281/zenodo.18599761)

A Python implementation of Stefan Th. Gries’s Collostructional Analysis methods, maintaining numerical consistency with the original R scripts.

## Project Structure

* **`core/`**: Main analysis logic and algorithm implementation.
    * Refer to the `README.md` within this directory for detailed usage instructions.
* **`simulator/`**: Interactive tools for visualizing and comparing association measures.
    * Refer to the `README.md` within this directory for detailed usage instructions.
* **`validation/`**: Procedures for numerical verification against Gries (2024) *Coll.analysis v4.1*. 
    * Refer to the comments at the beginning of `validator_for_core.py` for execution instructions.
* **`tests`** : Contains test scripts for pytest.
    * Required test CSVs (shared with `validator_for_core.py`) must be placed in `assets/` to run these tests.
* **`assets/`**: data directory for testing.
    * Required CSVs are not bundled due to copyright restrictions.

---

## Installation & Setup

This project supports both local environments and Google Colab.

### Using `uv`

```bash
uv sync
uv sync --extra notebook  # Includes simulator support

```

- If you only need the core analysis: use the default installation.
- If you want interactive visualizations (simulator): install with the `notebook` extra.


### Using `pip`

```bash
pip install .
pip install ".[notebook]" # Includes simulator support

```

### Google Colab Entry

This project is optimized for Google Colab. You can upload the `core/collostructional_analysis.py` script directly to your session or paste the code into a cell to get started immediately.


---

## Features by Module

### 1. Core Analysis (`core/`)

* **Comprehensive Analysis**: Supports Simple, Distinctive, and Co-varying collexeme analyses.
* **High Compatibility**: Includes custom Fisher-Yates Exact test methods to match the results of the original R scripts.
* **Signed Metrics**: Toggle `signed_metrics=True` to represent Repulsion as negative values for LLR and FYE.

### 2. Interactive Simulation (`simulator/`)

* **Metric Behavioral Visualization**: Explore how metrics like PMI and LOR react to changes in frequency and ratio.
* **Theoretical Case Studies**: Demonstrations of "Ranking Flips" between PMI and LOR, and scale variance/invariance between different measures.

### 3. Verification (`validation/`) and Testing (`tests/`)

* **Reproducibility Checks**: Automated comparisons to ensure results stay within predefined tolerance of reference R outputs.
    * Numerical comparisons are performed with predefined relative and absolute tolerances, to account for floating-point and implementation-level differences.
* **Environment-Specific Testing**: Provides both a standalone validator in `validation/` (compatible with Google Colab) and scripts in `tests/` designed for pytest.

---

## Implementation & Theoretical Notes

For in-depth information regarding the technical and theoretical aspects of this project, please refer to the following:

* **Algorithm & Implementation Details**: See `core/README.md` for notes on the implementation of FYE, LOR, and other metrics.
* **Metric Behavior & Case Studies**: See `simulator/README.md` for detailed discussions on PMI-LOR behavior and ranking inconsistencies.

---

## References


- Gries, Stefan Th. 2019. 15 years of collostructions. <i>International Journal of Corpus Linguistics</i>, 24(3), 385–412.
- Gries, Stefan Th. 2022. What do (most of) our dispersion measures measure (most)? Dispersion? <i>Journal of Second Language Studies</i>, 5(2), 171–205.
- Gries, Stefan Th. 2023. Overhauling Collostructional Analysis: Towards More Descriptive Simplicity and More Explanatory Adequacy. <i>Cognitive Semantics</i> 9(3): 351–386.
- Gries, Stefan Th. 2024. Coll.analysis 4.1. A script for R to compute perform collostructional analyses. <https://www.stgries.info/teaching/groningen/index.html>.
- Stefanowitsch, Anatol and Stefan Th. Gries. 2003. Collostructions: investigating the interaction between words and constructions. <i>International Journal of Corpus Linguistics</i> 8(2): 209–243.

## Acknowledgments

I would like to acknowledge Anatol Stefanowitsch, Stefan Th. Gries, and their collaborators for their foundational and pioneering work in collostructural analysis since 2003.

Special thanks are also due to Stefan Th. Gries for the continuous development of the original R scripts, including the latest 2024 update (v4.1), which served as the foundation for this Python implementation.

---


## Metadata

- **Maintainer:** [yz-rrr](https://github.com/yz-rrr)
- **ORCID:** [0009-0009-5953-3964](https://orcid.org/0009-0009-5953-3964)
- **License:** [MIT License](LICENSE).
- **Citation:** See [CITATION.cff](CITATION.cff) or the "Cite this repository" button on the sidebar.



---

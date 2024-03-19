
[![codecov](https://codecov.io/gh/cfe-lab/CFEIntact/branch/master/graph/badge.svg?token=OCYKUD7QET)](https://codecov.io/gh/cfe-lab/CFEIntact)
[![types - Mypy](https://img.shields.io/badge/types-Mypy-blue.svg)](https://github.com/python/mypy)
[![flake8 checked](https://img.shields.io/badge/flake8-checked-blueviolet.svg)](https://github.com/PyCQA/flake8)
[![License - MIT](https://img.shields.io/badge/license-MIT-9400d3.svg)](https://spdx.org/licenses/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/cfe-lab/CFEIntact/pulls)

# CFEIntact &middot; ![Python version](https://img.shields.io/badge/python-3.6+-blue.svg)

CFEIntact is an automated proviral intactness checker for HIV-1 consensus sequences, developed in Python. It enables researchers and clinicians to analyze HIV-1 sequences for intactness, aiding in the understanding of viral replication and persistence.

- **Automated Analysis:** Simplifies the process of checking HIV-1 consensus sequences for intactness with automated tools.
- **Extensible:** Designed to be extendable with additional checks and features as the field of HIV research evolves.
- **Open Source:** Encourages contributions and improvements from the community to make the tool more accurate and robust.

## Installation

CFEIntact is built with Python 3 and requires several dependencies, including MAFFT and BLAST+, for sequence alignment and analysis.

### Prerequisites

- `Python3` for running the main program.
- `pip` for managing Python packages.
- `MAFFT` for sequence alignment.
- `BLAST+` for sequence analysis when certain checks are enabled.

Clone the repository and install CFEIntact:

```shell
git clone https://github.com/cfe-lab/CFEIntact
cd CFEIntact
python3 -m pip install .
```

## Running CFEIntact

To analyze a set of FASTA sequences for HIV-1 proviral intactness:

```shell
cfeintact check --subtype HXB2 sequences.fasta
```

Note: currently it is recommended to use HXB2 even if your sequences are not subtype B.
This is because many well known reference subtypes contain defects in them.

## Documentation

Refer to the project's [website](https://cfe-lab.github.io/CFEIntact) for comprehensive documentation, including setup, usage examples, and development guidelines.

## Project Background

CFEIntact originated as a fork of [HIVIntact](https://github.com/ramics/HIVIntact), initially developed by Imogen Wright et al. This project owes its foundational concept to the original work and appreciates the efforts put in by the original team in advancing the field of HIV research.

However, over time, CFEIntact has undergone significant modifications and enhancements to better meet our evolving needs. These changes have led the project to diverge in ways that the original HIVIntact paper no longer accurately describes the forked version. CFEIntact has introduced new features, optimizations, and methodologies that are distinct from its origin, catering to a broader spectrum of analyses.

We encourage users and contributors to view CFEIntact as an extension and evolution of the ideas first introduced in HIVIntact, adjusted and expanded upon to serve contemporary research purposes.

## Contributing

We welcome contributions to CFEIntact! Whether you're fixing a bug, adding a new feature, or improving the documentation, your help makes CFEIntact better for everyone.

### License

CFEIntact is [MIT licensed](./LICENSE).

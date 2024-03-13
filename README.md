
[![codecov](https://codecov.io/gh/cfe-lab/CFEIntact/branch/master/graph/badge.svg?token=OCYKUD7QET)](https://codecov.io/gh/cfe-lab/CFEIntact)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/cfe-lab/CFEIntact/pulls)

# CFEIntact &middot; ![Python version](https://img.shields.io/badge/python-3.6+-blue.svg)

CFEIntact is an automated proviral intactness checker for HIV-1 consensus sequences, developed in Python. It enables researchers and clinicians to analyze HIV-1 sequences for intactness, aiding in the understanding of viral replication and persistence.

* **Automated Analysis:** Simplifies the process of checking HIV-1 consensus sequences for intactness with automated tools.
* **Extensible:** Designed to be extendable with additional checks and features as the field of HIV research evolves.
* **Open Source:** Encourages contributions and improvements from the community to make the tool more accurate and robust.

## Installation

CFEIntact is built with Python 3 and requires several dependencies, including MAFFT and BLAST+, for sequence alignment and analysis.

### Prerequisites

- Python 3 installed on your system.
- `pip` for managing Python packages.
- MAFFT for sequence alignment.
- BLAST+ for sequence analysis when certain checks are enabled.

### Install CFEIntact

Clone the repository and install CFEIntact using `virtualenv`:

```shell
git clone https://github.com/cfe-lab/CFEIntact
cd CFEIntact
python3 -m pip install .
```

## Running CFEIntact

To analyze a set of FASTA sequences for HIV-1 proviral intactness:

```shell
proviral intact --subtype B sequences.fasta
```

## Documentation

Refer to the project's [website](https://cfe-lab.github.io/CFEIntact) for comprehensive documentation, including setup, usage examples, and development guidelines.

## Contributing

We welcome contributions to CFEIntact! Whether you're fixing a bug, adding a new feature, or improving the documentation, your help makes CFEIntact better for everyone.

### License

CFEIntact is [MIT licensed](./LICENSE).

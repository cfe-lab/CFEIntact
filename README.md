
[![codecov](https://codecov.io/gh/cfe-lab/CFEIntact/branch/master/graph/badge.svg?token=OCYKUD7QET)](https://codecov.io/gh/cfe-lab/CFEIntact)
[![types - Mypy](https://img.shields.io/badge/types-Mypy-blue.svg)](https://github.com/python/mypy)
[![flake8 checked](https://img.shields.io/badge/flake8-checked-blueviolet.svg)](https://github.com/PyCQA/flake8)
[![License - MIT](https://img.shields.io/badge/license-MIT-9400d3.svg)](https://spdx.org/licenses/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/cfe-lab/CFEIntact/pulls)

# CFEIntact &middot; ![Python version](https://img.shields.io/badge/python-3.6+-blue.svg)

CFEIntact is [BCCfE](https://github.com/cfe-lab)'s version of [HIVIntact](https://github.com/ramics/HIVIntact) -- an automated proviral intactness checker for HIV-1 consensus sequences, developed in Python. It enables researchers and clinicians to analyze HIV-1 sequences for intactness, aiding in the understanding of viral replication and persistence.

- **Automated Analysis:** Simplifies the process of checking HIV-1 consensus sequences for intactness with automated tools.
- **Extensible:** Designed to be extendable with additional checks and features as the field of HIV research evolves.
- **Free software:** Encourages contributions and improvements from the community to make the tool more accurate and robust.

## Installation

The easiest way to install and run CFEIntact is through Docker. The pre-built Docker image is available on Docker Hub.

Ensure Docker is installed on your system. You can download and install Docker from [the official Docker website](https://www.docker.com/get-started).

Once Docker is installed, run:

```shell
docker run --rm cfelab/cfeintact version
```

This command should print something like `1.23.2`.

## Running CFEIntact

Assuming you have a file called `sequences.fasta` in the current directory.
Then to analyze a set of FASTA sequences that are included in that file, run:

```shell
docker run --rm -v .:/w cfelab/cfeintact check sequences.fasta
```

This command mounts the current directory to `/w` in the container, making your local files accessible. Be sure to adjust the command for your specific analysis (e.g., `sequences.fasta` is your input file).

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

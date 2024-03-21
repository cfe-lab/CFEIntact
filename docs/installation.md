
# CFEIntact Installation Guide

Welcome to the CFEIntact installation guide! This page aims to provide detailed instructions for installing CFEIntact on your system, guaranteeing smooth execution and operation of the program. CFEIntact, an automated proviral intactness checker, is primarily developed in Python and requires certain prerequisites for its proper function.

## Prerequisites

Before installing CFEIntact, ensure that you have the following software preinstalled on your system:

1. **Python3**: CFEIntact is developed in Python and requires Python3 to run. You can download it from the [official Python website](https://www.python.org/downloads/). After installation, verify the Python version by executing `python3 --version` in your system's command line.

2. **pip**: pip, a recursive acronym for "pip installs packages", is a Python package installer. It usually comes pre-installed with Python. You can verify its installation by running the `python3 -m pip --version` command. 

3. **mafft**: Multiple Alignment using Fast Fourier Transform (mafft) underpins high-speed and accurate alignment of multiple sequences, a critical operation in CFEIntact. To install mafft, follow this [link](https://mafft.cbrc.jp/alignment/software/source.html). 

4. **blastn**: Part of the BLAST+ suite from NCBI, blastn cross-references nucleotide sequences to locate regions of similarity. This sequence analysis and comparison tool is crucial for CFEIntact. Follow [this guide](https://www.ncbi.nlm.nih.gov/books/NBK569861/) for detailed installation instructions.

These prerequisites are essential for CFEIntact to operate as expected. For successful program execution, please ensure all these tools are installed and working correctly.

## Installing CFEIntact

Now that the prerequisites are set, let's proceed with installing CFEIntact. As a Python-based program, CFEIntact can be installed with pip through the following command:

```shell
python3 -m pip install 'git+https://github.com/cfe-lab/CFEIntact'
```

## Platform-specific Instructions

The installation process for CFEIntact will be largely the same across different platforms thanks to Python's cross-platform support. However, ensure you have the necessary permissions to install software on your system or use `sudo` (for Unix-based systems) if required. For Windows users, you may need to run the command prompt as an Administrator for the installation process. Mac users, depending on your system's security and privacy settings, may need to allow the system to install programs from identified developers or use the `sudo` command.

We hope this guide assists you in installing CFEIntact effectively. If you're experiencing unexpected issues or need further assistance, don't hesitate to reach out via our [GitHub issues](https://github.com/cfe-lab/CFEIntact/issues). Your feedback helps us continually improve and extend the usability of CFEIntact.

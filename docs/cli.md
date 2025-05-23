---
title: Interacting with CFEIntact
---

CFEIntact is executed through Docker, and it is a command-line only program.

Before running CFEIntact, ensure that you have [installed it](installation.md).

# Basic Usage

To check consensus sequences for intactness, use the `check` command followed by the required and optional arguments. Here is the basic syntax:

```bash
docker run -v .:/w cfelab/cfeintact check [OPTIONS] INPUT_FILE
```

Argument `INPUT_FILE` denotes the path to the sequence file in FASTA format that you wish to analyze. Make sure the file is located in the directory you're mounting, i.e., your current working directory.

# Options

The `check` command provides several options that allow you to customize the analysis based on your requirements:

- `--subtype`: Specifies the HIV subtype for the analysis. It is a required parameter.
- `--check-packaging-signal` / `--ignore-packaging-signal`: Enables or disables the check for the Packaging Signal Region (PSI). Enabled by default.
- `--check-rre` / `--ignore-rre`: Enables or disables the check for the Rev Response Element (RRE). Enabled by default.
- `--check-major-splice-donor-site` / `--ignore-major-splice-donor-site`: Enables or disables the check for the Major Splice Donor (MSD) site. Enabled by default.
- `--check-hypermut` / `--ignore-hypermut`: Enables or disables the check for hypermutations. Enabled by default.
- `--check-long-deletion` / `--ignore-long-deletion`: Enables or disables the check for large deletions. Enabled by default.
- `--check-nonhiv` / `--ignore-nonhiv`: Enables or disables the check for non-HIV sequences. Enabled by default.
- `--check-scramble` / `--ignore-scramble`: Enables or disables the check for scrambled sequences. Enabled by default.
- `--check-internal-inversion` / `--ignore-internal-inversion`: Enables or disables the check for internal inversions. Enabled by default.
- `--check-unknown-nucleotides` / `--ignore-unknown-nucleotides`: Enables or disables the check for unknown nucleotides in the sequences. Enabled by default.
- `--check-small-orfs` / `--ignore-small-orfs`: Enables or disables the analysis of small Open Reading Frames (ORFs). Enabled by default.
- `--check-distance` / `--ignore-distance`: Enables or disables the distance-based analysis. This includes the `InsertionInOrf` error, as well as `SequenceDivergence` one. Enabled by default.
- `--output-csv` / `--output-json`: Chooses between CSV and JSON output format. CSV is selected by default.
- `--output`: Specifies the directory where output files will be stored. Defaults to the current working directory (which is `/w` in the case of the docker version).

# Examples

Run an analysis on `sequences.fasta` for subtype 'B', checking all aspects except for genetic distance:

```bash
docker run -v .:/w cfelab/cfeintact check --subtype B --ignore-distance sequences.fasta
```

Check sequences in `input.fasta` for subtype 'A' with default settings and save outputs in JSON format:

```bash
docker run -v .:/w cfelab/cfeintact check --subtype A --output-json input.fasta
```

# Understanding the Output

CFEIntact generates several output files based on the analysis, including intact sequences, nonintact sequences, subtypes, ORFs information, holistic info, and error details. The format of these files can be JSON or CSV, as specified by the user. These files provide a comprehensive report on the analysis, including detected defects, ORF analysis, subtype information, and more.

Consult the generated output files in the specified `--output` folder to review the analysis results and take necessary actions based on the detected genetic anomalies.

Refer to [the Inputs and Outputs page](io.md) to learn more.

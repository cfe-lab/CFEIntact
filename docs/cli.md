
# Running CFEIntact

CFEIntact is a comprehensive tool designed for analyzing HIV sequences to check for intactness and detect various genetic anomalies such as hypermutations, large deletions, non-HIV fragments, scrambling, and inversions. It also assesses key genetic regions, including the Packaging Signal Region (PSI), the Rev Response Element (RRE), and the Major Splice Donor (MSD) site. This document will guide you on how to run the CFEIntact tool and understand the options available for a customized analysis.

## Prerequisites

Before running CFEIntact, ensure that you have [installed it](installation.md).

## Running the Tool

CFEIntact is executed through a command-line interface (CLI).

### Basic Usage

To check consensus sequences for intactness, use the `check` command followed by the required and optional arguments. Here is the basic syntax:

```bash
cfeintact check [OPTIONS] INPUT_FILE
```

`INPUT_FILE` denotes the path to the sequence file in FASTA format that you wish to analyze.

### Options

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
- `--output-csv` / `--output-json`: Chooses between CSV and JSON output format. CSV is selected by default.
- `--working-folder`: Specifies the directory where output files will be stored. Defaults to the current working directory.

### Examples

Run an analysis on `sequences.fasta` for subtype 'A', checking all aspects except small ORFs, and output the results in JSON format:

```bash
cfeintact check --subtype A --ignore-small-orfs --output-json sequences.fasta
```

Check sequences in `input.fasta` for subtype 'B' with default settings and save outputs in CSV format in a specific directory:

```bash
cfeintact check --subtype B --output-csv --working-folder /path/to/output_folder input.fasta
```

## Understanding the Output

CFEIntact generates several output files based on the analysis, including intact sequences, nonintact sequences, subtypes, ORFs information, holistic info, and error details. The format of these files can be JSON or CSV, as specified by the user. These files provide a comprehensive report on the analysis, including detected defects, ORF analysis, subtype information, and more.

Consult the generated output files in the specified working folder to review the analysis results and take necessary actions based on the detected genetic anomalies.

Refer to [the Inputs and Outputs page](io.md) to learn more.

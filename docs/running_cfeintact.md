---
title: Running CFEIntact
---

Once your data are prepared and organized into a single FASTA file, it's time to run CFEIntact and analyze your HIV-1 sequences. This guide walks you through the essential steps to launch the analysis using Docker.

---

## Quick Start Command

The simplest way to run CFEIntact is via Docker. Open your terminal in the directory where your FASTA file is located (for example, "sequences.fasta") and run the following command:

```shell
docker run --rm -v .:/w cfelab/cfeintact check --subtype B sequences.fasta
```

### Breaking It Down

- **docker run --rm**: Runs the container and automatically removes it after execution.
- **-v .:/w**: Mounts your current directory into the container at `/w`, allowing the tool to access your data.
- **cfelab/cfeintact**: The Docker image for CFEIntact.
- **check**: Tells CFEIntact to analyze your sequences.
- **--subtype B**: Specifies the HIV subtype ("B" in this case). Adjust this option based on your knowledge of the sequence subtype.
- **sequences.fasta**: The name of your input FASTA file containing the viral sequences.

---

## Additional Options

CFEIntact offers a range of command-line options to enable or disable specific checks, such as examining the packaging signal, RRE, hypermutation, and more. For a quick run, the default options are typically sufficient. However, for more detailed analyses, you can explore the following:

- **--check-packaging-signal / --ignore-packaging-signal**
- **--check-rre / --ignore-rre**
- **--check-major-splice-donor-site / --ignore-major-splice-donor-site**
- **--check-hypermut / --ignore-hypermut**
- **--check-long-deletion / --ignore-long-deletion**
- **--check-nonhiv / --ignore-nonhiv**
- **--check-scramble / --ignore-scramble**
- **--check-internal-inversion / --ignore-internal-inversion**
- **--check-unknown-nucleotides / --ignore-unknown-nucleotides**
- **--check-small-orfs / --ignore-small-orfs**

You may combine these options according to the level of detail needed for your research. For the quick start, the example command covers the essentials.

---

## What Happens Next?

When you run the command:

- **Processing:**
  CFEIntact reads your FASTA file, aligns sequences against appropriate reference subtypes, and evaluates each sequence for defects such as deletions, frame shifts, and hypermutations.

- **Output Files:**
  The tool generates several outputs in your mounted directory including:
  - **regions.csv or regions.json**: Details the detected open reading frames.
  - **holistic.csv or holistic.json**: Provides overall analysis information for each sequence.
  - **defects.csv or defects.json**: Lists any issues detected in the sequences.
  - **subtypes.fasta**: Contains reference subtype sequences used in the alignment.
  - Optionally, a BLAST output file might be produced for further investigation.

---

Next: [output interpretation](io.html).

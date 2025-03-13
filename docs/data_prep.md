---
title: Data Preparation
---

Before running CFEIntact, you need to prepare your sequence file so that the analysis runs flawlessly. The core requirement is a properly formatted FASTA file containing the viral sequences you wish to analyze. Once you have this file ready, you can quickly move on to processing your data through CFEIntact. Below are the essential steps and background information to help you get started.

---

## Get Your Data Ready

- **Prepare a Single FASTA File:**
  Ensure that you have a FASTA file containing all the HIV--1 consensus sequences. For example:

  ```FASTA
  >Sample_1
  ATGCGTACGATCGATCGT
  >Sample_2
  ATGCTAGCTAGCTAGCTATACGATCGAT
  ```

- **Clean and Quality--Check the Data:**
  If your sequences come from high--throughput platforms, consider running a quality control step to trim low--quality regions or remove ambiguous characters. This helps ensure that the analysis focuses on high--confidence sequence regions.

- **Data Organization:**
  If you have multiple FASTA files, combine them into a single file. A consolidated file minimizes complications when mounting data into the Docker container with CFEIntact.

Once your FASTA file is ready and validated, you're set to run CFEIntact.

---

## Example Data

For practical understanding and experimentation, download [this FASTA file](https://raw.githubusercontent.com/cfe-lab/CFEIntact/refs/heads/master/tests/data-small.fasta) for a sample input that work for CFEIntact. Using these samples, you can become familiar with the structure and usage of the input files.

---

Next: [running CFEIntact](running_cfeintact.html).

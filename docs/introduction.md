---
title: Introduction
---

This guide is designed to help people without programming expertise analyze proviral sequences using CFEIntact. By leveraging a Docker container, we simplify the technical setup, allowing them to focus on the biological insights from their data. This page provides a high level overview, before the dive into practical matters of using the program.

---

# Background

HIV/AIDS research sits at the intersection of modern genomics and clinical medicine, where understanding the minutiae of viral structure can lead to breakthroughs in treatment and prevention. CFEIntact operates in this vibrant domain by analyzing the genetic integrity of HIV--1 proviral sequences. In essence, the tool examines whether the viral genomes embedded within infected cells remain "intact" or have incurred defects that would render them nonfunctional.

### The biology behind proviral intactness

When HIV infects a host cell, it integrates its genomic material into the host's DNA as a "provirus." Although antiretroviral therapies can control active viral replication, the reservoir of these integrated proviruses often persists in a latent state. However, not all integrated sequences are created equal. Many are defective; they can contain large deletions, frameshift mutations, internal stop codons, or other structural anomalies that prevent the production of infectious virus. Conversely, a subset of these proviruses remains intact --- capable of reactivation and potentially driving viral rebound if therapy is interrupted.

![image](https://s7d1.scene7.com/is/image/CENODS/09705-scicon5-hiv?&wid=400)

Understanding the composition of proviral reservoir is crucial. It informs both the evaluation of long--term therapy success and strategic efforts aimed at eradicating the latent reservoir, a key barrier to achieving a cure.

### The bioinformatics challenge

The complexity of the HIV genome --- with its overlapping reading frames, high mutability, and error--prone replication --- demands sophisticated bioinformatics tools. Traditional sequence alignment techniques are enhanced by specialized algorithms that can detect subtle genetic deviations. CFEIntact leverages those algorithms to examine Open Reading Frames, evaluate insertion and deletion events, and detect hypermutation signatures induced by host defense enzymes.

In performing this analysis, CFEIntact converts messy biological data into actionable insights. By providing comprehensive reports across multiple regions of the HIV genome --- from structural genes like gag, pol, and env to smaller regulatory ORFs --- researchers gain a deeper perspective on viral diversity and the functional viability of these sequences.

---

# Why CFEIntact?

Here are several reasons why CFEIntact is a particularly well--suited tool for the problems highlighted above:

- CFEIntact detects a diverse range of defects including large deletions, out--of--frame insertions/deletions, internal stop codons, as well as mutations in key functional domains like the Packaging Signal, Rev Response Element, and Major Splice Donor site.
- It examines almost every critical region of the HIV genome --- from the major structural genes (gag, pol, env) to the smaller regulatory open reading frames (vif, vpr, tat, vpu, rev, nef).
- CFEIntact calculates useful metrics like genetic distance and indel impact.
- With a rich set of command--line options, users can easily toggle individual checks. For example, one can enable or disable tests for hypermutation, packaging signal integrity, revision of the donor site, as well as checks on sequence divergence.
- This modularity means that CFEIntact can be tailored to the specific needs of different research projects or experimental setups.
- Comprehensive reports are generated, including:
  - A holistic summary that covers overall sequence integrity, alignment coverage, and subtype inference.
  - Detailed breakdowns of identified defects.
  - Annotations of the detected open reading frames and specific regions affected by insertions, deletions, or other irregularities.
- Additionally, a BLAST file and a subtype FASTA file are produced, allowing for further inspection and cross-referencing against known HIV reference sequences.
- The tool outputs its findings in both CSV and JSON formats, making integration with downstream pipelines and further statistical analysis straightforward.
- With a dedicated Dockerfile and accompanying development container configuration, users can quickly deploy the tool in a consistent, reproducible environment.
- The availability of both `pip`--based installation and Docker deployment means that users can easily integrate CFEIntact into laboratory workflows without the hassle of complex dependency issues.
- The codebase is thoroughly tested using GitHub Actions with multiple test suites, ensuring the tool remains robust as it evolves.
- Linting and type checking via `flake8` and `mypy` help maintain code quality, making it a reliable platform for cutting--edge HIV research.

CFEIntact exemplifies how modern bioinformatics tools can distill complex genetic data into clear, meaningful narratives about viral integrity. It is a product of interdisciplinary design --- drawing on molecular biology, computer science, and statistical analysis --- to meet one of the most pressing challenges in HIV research. The tool stands as a testament to the progress made when innovative technology is harnessed in service of understanding life at its most intricate and impactful level.

---

# Setup and Environment

CFEIntact is a command-line program --- it is designed to be interacted with from within a terminal.
If you are not experienced with that, see the following guides:
- For Windows: [a beginner's guide to the Windows Command Prompt](https://www.makeuseof.com/tag/a-beginners-guide-to-the-windows-command-line).
- For POSIX (GNU, BSD, MacOS, ...): [command line for beginners](https://www.freecodecamp.org/news/command-line-for-beginners).

To make CFEIntact accessible and easy to use, we have encapsulated it within a Docker container. This ensures all dependencies and configurations are pre--installed, minimizing compatibility issues across different systems. Docker is a platform to develop, ship, and run applications in isolated environments.

By following the rest of this guide, you will be up and running with proviral sequence analysis quickly, even with minimal technical background.

---

Next: [installation](quick_install.html).

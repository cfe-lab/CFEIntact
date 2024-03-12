
# CFEIntact workflow

## Preparation

Before analyzing a sequence, CFEIntact does some initial preprocessing that is used later in the analysis.

### Alignment to Reference

For each input sequence, CFEIntact uses the `mafft` software to align it to its subtype sequence.
This operation is repeated with the reverse complement (RC) of the input sequence to determine if the fit is better.
If the RC provides a better alignment, CFEIntact uses it instead.
This ensures that the direction in which the original sequence is read does not affect the analysis.

The alignment is global, and it never fails.


### ORF Detection

(TODO: describe the logic)

Detection never fails, but it can output ORFs that have a length of 0.

The outputs from this procedure are used to produce the `orfs.json` file.


### BLAST Analysis

Optionally, CFEIntact calls the NCBI's `blastn` program to obtain alignment data that is region-based,
as opposed to the global alignment provided by `mafft`.

## Analysis steps

CFEIntact performs multiple analysis steps that are independent of each other.
Most of these analyses are optional and controlled by passing a command line option to the main program.

Table below lists all of the independent steps:

| Name                      | Enabled by Default? | Command Line Option                       |
| --------------------------| --------------------| ------------------------------------------|
| PSI check                 | Yes                 | `--exclude-packaging-signal`              |
| RRE check                 | Yes                 | `--exclude-rre`                           |
| MSD check                 | Yes                 | `--ignore-major-splice-donor-site`        |
| Hypermutation check       | No                  | `--run-hypermut`                          |
| Large deletion check      | No                  | `--check-long-deletion`                   |
| NonHIV check              | No                  | `--check-nonhiv`                          |
| Scramble check            | No                  | `--check-scramble`                        |
| Inversion check           | No                  | `--check-internal-inversion`              |
| Large ORFs analysis       | Yes                 |                                           |
| Small ORFs analysis       | No                  | `--include-small-orfs`                    |


Each step works on a single input sequence and has a set of potential errors it can detect for that sequence.

We describe the logic of each step below.

### PSI check

Determines presence and possible intactness of HIV Packaging Signal Region.

Based on the alignment, CFEIntact locates the PSI region in the input sequence and checks its length.
For lengths smaller than the [tolerable limit](cutoffs.md), an error with code `PackagingSignalDeletion` is reported.

### RRE check

Determines presence and possible intactness of HIV Rev Response Element.

The analysis performed is the same as the PSI check, but the error code is `RevResponseElementDeletion`.

### MSD check

Determines whether the Major Splice Donor site is mutated.

Based on the alignment, CFEIntact locates the region that is expected to contain the MSD subsequence.
If the found subsequence is anything but `G` followed by `T`,
an error with `MajorSpliceDonorSiteMutated` error code is reported.

### Hypermutation check

APOBEC3G/F hypermutation scan and test based on Rose and Korber, Bioinformatics (2000).
Briefly, scans reference for APOBEC possible signatures and non-signatures and performs
fisher test based on ratio of G->A in ref -> query at these signatures.

If there is enough evidence that the sequence is hypermutated (p-value <0.05)
this step outputs the `APOBECHypermutationDetected` error code.

### Large deletion check

If the input sequence is shorter than 8000 nucleotide bases,
an error with code `LongDeletion` is produced.

### NonHIV check

(TODO: describe the logic)

This analysis step is copied from [HIVSeqinR software](https://github.com/guineverelee/HIVSeqinR).

It outputs the `NonHIV` error code.

### Scramble check

(TODO: describe the logic)

This analysis step is copied from [HIVSeqinR software](https://github.com/guineverelee/HIVSeqinR).

It outputs the `Scramble` error code.

### Inversion check

(TODO: describe the logic)

This analysis step is copied from [HIVSeqinR software](https://github.com/guineverelee/HIVSeqinR).

It outputs the `InternalInversion` error code.

### Large ORFs analysis

Large ORFs are `gag`, `pol` and `env`.

Using data from the ORF detection procedure,
CFEIntact goes through each of the listed ORFs and checks two things:

- their lengths
- and possible out-of-frame indels.

The length check is based on comparing it to predefined length limits.
Go to [cuttofs page](cutoffs.md) now to learn the limits that CFEIntact adhears to.
If the length is too long, the error code is `InsertionInOrf`.
If the length is too short, then CFEIntact also checks if there is an internal stop codon in the analyzed ORF.
Depending on that, the code is either `DeletionInOrf` or `InternalStopInOrf`.
Notably, internal stop codons that do not make the resulting protein too short are ignored.

When out-of-frame indels detection is run,
the assumption is that they are common but not all of them render the respective ORF dysfunctional.
So we try to estimate the impact that detected indels have.
This is done with the following algorithm:

1. Redo the alignment, not global this time, but only for the current ORF.
2. Iterate through each position in the alignment, keeping the index of "current frame":

    - if insertion is encountered, change the current frame by +1,
    - if deletion is encountered, change the current frame by -1,

    counting all nucleotides that are not in the initial frame as we go.

3. Compare the amount of out of frame nucleotides to a predefined limit.

If the impact is large enough, meaning that too many nucleotides got frame shifted,
then output an error with code `FrameshiftInOrf`.

(TODO: show example)

### Small ORFs analysis

Small ORFs are `vif`, `vpr`, `tat`, `vpu`, `rev` and `nef`.
The analysis performed for them is the same as the Large ORFs analysis,
as well as the possible error codes that this step outputs.

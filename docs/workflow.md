
# CFEIntact workflow

## Preparation

Before analyzing a sequence, CFEIntact does some initial preprocessing that is used later in the analysis.

### BLAST Analysis

CFEIntact calls the NCBI's `blastn` program to obtain alignment data that is region-based, as opposed to global alignment.

The subtype of the sequence is determined at this point as well - BLAST tries to align the sequence to every reference subtype sequence specified via the `--subtype` CLI option.

This step is optional, depending on the command line arguments.

### Alignment to Reference

CFEIntact uses BioPython's global pairwise aligner to align the sequence to its subtype sequence.
This operation is repeated with the reverse complement (RC) of the input sequence to determine if the fit is better.
If the RC provides a better alignment, CFEIntact uses it instead.
This ensures that the direction in which the original sequence is read does not affect the analysis.

The alignment is global, and it never fails.

### ORF Detection

The logic for Open Reading Frame (ORF) detection operates under the principle of identifying gene segments within the HIV genome that have the potential to code for proteins. The steps generally involve:

1. **Mapping ORFs** to known HIV genes (e.g., gag, pol, env) based on their positions within the sequence alignment to the reference.
2. **Scanning the sequence** for start codons that may indicate the beginning of an ORF.
3. **Tracing the sequence** from the start codon to the nearest stop codon (e.g., "TAA", "TAG", "TGA") without encountering other stop codons in between, as this would indicate a potential ORF.

Detection never fails, but it can output ORFs that have a length of 0.

The outputs from this procedure are used to produce the `orfs.json` file.

## Analysis steps

CFEIntact performs multiple analysis steps that are independent of each other.
Most of these analyses are optional and controlled by passing a command line option to the main program.

Table below lists all of the independent steps:

| Name                      | Command Line Option                |
|---------------------------|------------------------------------|
| PSI check                 | `--ignore-packaging-signal`        |
| RRE check                 | `--ignore-rre`                     |
| MSD check                 | `--ignore-major-splice-donor-site` |
| Hypermutation check       | `--ignore-hypermut`                |
| Large deletion check      | `--ignore-long-deletion`           |
| NonHIV check              | `--ignore-nonhiv`                  |
| Scramble check            | `--ignore-scramble`                |
| Inversion check           | `--ignore-internal-inversion`      |
| Large ORFs analysis       |                                    |
| Small ORFs analysis       | `--ignore-small-orfs`              |

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
this step outputs the `APOBECHypermutation` error code.

### Large deletion check

If the input sequence is shorter than 8000 nucleotide bases,
an error with code `LongDeletion` is produced.

### NonHIV check

To determine whether the sequence might be non-HIV or significantly divergent from known HIV sequences, the following logic is applied:

1. **Calculate the coverage** of the input sequence by the known HIV sequences using the alignment data from the BLAST analysis. Coverage is measured as the percentage of the input sequence that aligns with known HIV sequences.
2. **Threshold determination**: If the coverage falls below a certain threshold (e.g., 80%), it suggests that a significant portion of the sequence does not align with known HIV sequences, potentially indicating a non-HIV origin or considerable divergence.
3. **Report NonHIV**: If the coverage is below the threshold, a `NonHIV` error is reported, suggesting that the sequence may not be HIV or may be a highly divergent strain.

This analysis step is copied from [HIVSeqinR software](https://github.com/guineverelee/HIVSeqinR).

### Scramble check

For the scramble check, the logic is aimed at detecting if segments of the HIV genome appear in an unusual order, suggesting potential recombination events or errors in sequence assembly:

1. **Analyze BLAST alignment**: Review the alignment start and end positions of each segment from the BLAST analysis.
2. **Expect sequential alignment**: In a non-scrambled, intact HIV genome, segments should align sequentially without overlap or significant gaps.
3. **Detect discrepancies**: If segments are found to align out of expected order or with unexpected overlaps or gaps, it suggests the sequence may be scrambled.
4. **Report Scramble**: If evidence of scrambling is detected, a `Scramble` error is reported, indicating potential recombination or assembly issues.

This analysis step is copied from [HIVSeqinR software](https://github.com/guineverelee/HIVSeqinR).

### Inversion check

The inversion check is designed to identify sequences that may have undergone an inversion, a genetic event where a segment of the genome is reversed end to end:

1. **Orientation analysis**: Using BLAST alignment data, assess the orientation of each aligned segment compared to the reference sequence. Each segment has an associated strand orientation - forward or reverse.
2. **Detect mixed orientations**: In an intact HIV genome, all segments should align in the same orientation. The presence of segments that align in both forward and reverse orientations within the same sequence suggests possible inversion.
3. **Report Internal Inversion**: If evidence of internal inversion is detected, an `InternalInversion` error is reported, indicating potential genetic rearrangement within the sequence.

This analysis step is copied from [HIVSeqinR software](https://github.com/guineverelee/HIVSeqinR).

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

3. Check if the indel impact is nonzero (meaning out-of-frame indels were detected).
4. Check if the amino acid sequence distance exceeds the distance threshold for that ORF.

A `FrameshiftInOrf` error is reported if and only if:
- The indel impact is nonzero (out-of-frame indels are present), **AND**
- The sequence distance crosses the distance threshold.

This combined approach ensures that frameshifts are only reported when they occur alongside significant sequence divergence, reducing false positives from minor frame variations in otherwise intact sequences.

### Small ORFs analysis

Small ORFs are `vif`, `vpr`, `tat`, `vpu`, `rev` and `nef`.
The analysis performed for them is the same as the Large ORFs analysis,
as well as the possible error codes that this step outputs.

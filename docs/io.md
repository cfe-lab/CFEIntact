
# Inputs and Outputs

The only input is a path to a `.fasta` file that you want to analyze.
This file should include the genomic sequences of interest.

Once the analysis is complete, CFEIntact will generate four output files:

| Filename        | Description                                                                                                                      |
|-----------------|----------------------------------------------------------------------------------------------------------------------------------|
| errors.csv      | This file contains associations between sequences and their identified defects.                                                  |
| orfs.csv        | This file contains associations between sequences and their identified Open Reading Frames (ORFs).                               |
| holistic.csv    | This file contains comprehensive analysis details related to the sequences as a whole, rather than being specific to each ORF.   |
| intact.fasta    | This file contains a list of sequences with no fatal defects identified. These sequences are considered putative intact genomes. |
| nonintact.fasta | This file contains a list of sequences with identified defects. These sequences are not considered putative intact genomes.      |
| blast.csv       | This file contains output from BLASTN software. It is produced conditionally.                                                    |

If you pass the `--output-json` option to CFEIntact, the output format will be `.json` instead of `.csv`.

Here is an example of the contents of the `errors.csv` file:

| sequence_name | error                       | message                                                                                                             | orf |
|---------------|-----------------------------|---------------------------------------------------------------------------------------------------------------------|-----|
| KX505501.1    | DeletionInOrf               | ORF pol at 2084-5096 can have maximum deletions 30, got 2892.                                                       | pol |
| KX505501.1    | RevResponseElementDeletion  | Query Sequence exceeds maximum deletion tolerance in RRE. Contains 35 deletions with max tolerance of 20 deletions. |     |
| MN691959      | InternalStopInOrf           | Smaller ORF vpu at 6060-6309 contains an internal stop codon at 6100.                                               | vpu |
| MK114856.1    | APOBECHypermutationDetected | Query sequence shows evidence of APOBEC3F/G-mediated hypermutation (p = 3.639064030015132e-65).                     |     |
| MK116110.1    | PackagingSignalDeletion     | Query Sequence exceeds maximum deletion tolerance in PSI. Contains 93 deletions with max tolerance of 10 deletions. |     |

Here is an example of the contents of the `holistic.csv` file:

| seqid      | intact | qlen | hypermutation_probablility | inferred_subtype                     | blast_matched_qlen | blast_sseq_coverage | blast_qseq_coverage | blast_sseq_orfs_coverage | orfs_start | orfs_end | blast_n_conseqs |
|------------|--------|------|----------------------------|--------------------------------------|--------------------|---------------------|---------------------|--------------------------|------------|----------|-----------------|
| KX505501.1 | False  | 1997 | 0.71                       | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455 | 1997               | 0.25                | 1.22                | 0.18                     | 789        | 8793     | 4               |
| MN691959   | False  | 9493 | 0.20                       | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455 | 9493               | 1.08                | 1.11                | 1                        | 789        | 8793     | 3               |
| MN692074   | False  | 4178 | 0.36                       | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455 | 4178               | 0.50                | 1.17                | 0.41                     | 789        | 8793     | 4               |
| MN692145   | True   | 9689 | 0.17                       | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455 | 9689               | 1.13                | 1.13                | 1                        | 789        | 8793     | 3               |

### Field descriptions

- `seqid`: The identifier or name for the sequence (same as in the input FASTA file)
- `intact`: Whether the sequence is considered to be intact (True) or not (False)
- `qlen`: Length of the _query_ sequence
- `hypermutation_probablility`: The probability that the sequence shows evidence of hypermutation
- `inferred_subtype`: The suspected subtype of the sequence based on analysis
- `blast_matched_qlen`: The length of the matched subsequence from the BLASTN software
- `blast_sseq_coverage`: The coverage of the subject sequence from BLASTN output
- `blast_qseq_coverage`: The coverage of the query sequence from BLASTN output
- `blast_sseq_orfs_coverage`: The percentage coverage of ORFs in the subject sequence from BLASTN
- `orfs_start`: The starting point of the ORFs
- `orfs_end`: The ending point of the ORFs
- `blast_n_conseqs`: The number of consecutive sequences from BLASTN results

Here is an example of the contents of the `orfs.csv` file:

| seqid      | name | start | end  | orientation | distance            | protein | aminoacids | nucleotides | subtype_start | subtype_end | subtype_aminoacids | subtype_nucleotides |
|------------|------|-------|------|-------------|---------------------|---------|------------|-------------|---------------|-------------|--------------------|---------------------|
| KX505501.1 | env  | 0     | 1823 | forward     | 0.7623480451210163  | MGAR... | GLSG*...   | GGTCT...    | 6223          | 8793        | MRVKE...           | ATGAG...            |
| KX505501.1 | vif  | 0     | 1823 | forward     | 0.7647696476964769  | MGAR... | GLSG*...   | GGTCT...    | 5040          | 5618        | MENR...            | ATGG...             |
| MN691959   | env  | 6070  | 8655 | forward     | 0.1405525502318391  | MRVK... | MRVK...    | ATGAG...    | 6223          | 8793        | MRVKE...           | ATGAG...            |
| MN691959   | vif  | 4890  | 5468 | forward     | 0.09157509157509158 | MENR... | MENR...    | ATGG...     | 5040          | 5618        | MENR...            | ATGG...             |

### Field descriptions

- `seqid`: The identifier for the sequence from which the ORF (Open Reading Frame) is derived.
- `name`: The name of the gene associated with the ORF.
- `start`: The start position of the ORF within the query sequence.
- `end`: The end position of the ORF within the query sequence.
- `orientation`: Indicates if the ORF is in the forward or reverse orientation with respect to the original sequence.
- `distance`: A numerical value representing some form of calculated metric or distance related to the ORF. It is likely related to the real genetic distance.
- `protein`: A representation of the sequence encoded by the ORF.  Due to space constraints, a snippet or indicator (like "GLSG*...") is provided in this example.
- `aminoacids`: The full or partial amino acid sequence encoded by the ORF.
- `nucleotides`: The nucleotide sequence of the ORF.
- `subtype_start`: The start position of the ORF within the identified subtype sequence.
- `subtype_end`: The end position of the ORF within the identified subtype sequence.
- `subtype_aminoacids`: The amino acid sequence relevant to the identified subtype.
- `subtype_nucleotides`: The nucleotide sequence relevant to the identified subtype.

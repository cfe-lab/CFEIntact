
# Inputs and Outputs

CFEIntact generates several output files based on the analysis, including intact sequences, nonintact sequences, subtypes, ORFs information, holistic info, and error details. The format of these files can be JSON or CSV, as specified by the user. These files provide a comprehensive report on the analysis, including detected defects, ORF analysis, subtype information, and more.

The only input is a path to a `.fasta` file that you want to analyze.
This file should include the genomic sequences of interest.

Once the analysis is complete, CFEIntact will generate these output files:

- `defects.csv`
- `regions.csv`
- `holistic.csv`
- `subtypes.fasta`
- `blast.csv`

If you pass the `--output-json` option to CFEIntact, the output format will be `.json` instead of `.csv`.

## `defects.csv`

This file contains associations between sequences and their identified defects.
Here is an example of the contents of the `defects.csv` file:

| qseqid     | code                       | message                                                                                                             | region    |
|------------|----------------------------|---------------------------------------------------------------------------------------------------------------------|-----------|
| KX505501.1 | InternalStopInOrf          | ORF 'pol' at 1629-1927 contains an internal stop codon at 1746.                                                     | pol       |
| KX505501.1 | RevResponseElementDeletion | Query Sequence exceeds maximum deletion tolerance in RRE. Contains 35 deletions with max tolerance of 20 deletions. |           |
| MN691959   | DeletionInOrf              | ORF 'tat_exon2' exeeds maximum deletion tolerance. Contains 45 deletions with max tolerance of 0 deletions.         | tat_exon2 |
| MK114856.1 | APOBECHypermutation        | Query sequence shows evidence of APOBEC3F/G-mediated hypermutation (p = 3.639064030015132e-65).                     |           |
| MK116110.1 | PackagingSignalDeletion    | Query Sequence exceeds maximum deletion tolerance in PSI. Contains 93 deletions with max tolerance of 10 deletions. |           |

Here, and below `qseqid` stands for "Query Sequence Id", which is the same sequence name as in the input `.fasta` file.

## `regions.csv`

This file contains associations between sequences and their identified Open Reading Frames (ORFs).
Here is an example of the contents of the `regions.csv` file:

| qseqid     | name | start | end  | orientation | distance            | protein | aminoacids | nucleotides | subtype_start | subtype_end | subtype_aminoacids | subtype_nucleotides |
|------------|------|-------|------|-------------|---------------------|---------|------------|-------------|---------------|-------------|--------------------|---------------------|
| KX505501.1 | env  | 0     | 1823 | forward     | 0.7623480451210163  | MGAR... | GLSG*...   | GGTCT...    | 6223          | 8793        | MRVKE...           | ATGAG...            |
| KX505501.1 | vif  | 0     | 1823 | forward     | 0.7647696476964769  | MGAR... | GLSG*...   | GGTCT...    | 5040          | 5618        | MENR...            | ATGG...             |
| MN691959   | env  | 6070  | 8655 | forward     | 0.1405525502318391  | MRVK... | MRVK...    | ATGAG...    | 6223          | 8793        | MRVKE...           | ATGAG...            |
| MN691959   | vif  | 4890  | 5468 | forward     | 0.09157509157509158 | MENR... | MENR...    | ATGG...     | 5040          | 5618        | MENR...            | ATGG...             |

### Field descriptions

- `qseqid`: The identifier for the sequence from which the ORF (Open Reading Frame) is derived.
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

## `holistic.csv`

This file contains comprehensive analysis details related to the sequences as a whole, rather than being specific to each ORF.
Here is an example of the contents of the `holistic.csv` file:

| qseqid     | intact | qlen | inferred_subtype                         | hypermutation_probablility | is_reverse_complement | blast_matched_qlen | blast_sseq_coverage | blast_qseq_coverage | blast_sseq_orfs_coverage | orfs_start | orfs_end | blast_n_conseqs |
|------------|--------|------|------------------------------------------|----------------------------|-----------------------|--------------------|---------------------|---------------------|--------------------------|------------|----------|-----------------|
| KX505501.1 | False  | 1997 | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.CfE | 0.71                       | False                 | 1997               | 0.25                | 1.22                | 0.18                     | 789        | 8793     | 4               |
| MN691959   | False  | 9493 | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.CfE | 0.20                       | False                 | 9493               | 1.08                | 1.11                | 1                        | 789        | 8793     | 3               |
| MN692074   | False  | 4178 | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.CfE | 0.36                       | False                 | 4178               | 0.50                | 1.17                | 0.41                     | 789        | 8793     | 4               |
| MN692145   | True   | 9689 | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.CfE | 0.17                       | False                 | 9689               | 1.13                | 1.13                | 1                        | 789        | 8793     | 3               |

### Field descriptions

- `qseqid`: The identifier or name for the sequence (same as in the input `.fasta` file)
- `intact`: Whether the query sequence is considered to be intact (True) or not (False)
- `qlen`: Length of the _query_ sequence
- `inferred_subtype`: The suspected subtype of the sequence based on analysis. Same as `sseqid` in `blast.csv`.
- `hypermutation_probablility`: The probability that the sequence shows evidence of hypermutation
- `is_reverse_complement`: Whether the query sequence has been reverse complemented to better fit the reference sequence
- `blast_matched_qlen`: The length of the matched subsequence from the BLASTN software
- `blast_sseq_coverage`: The coverage of the subject sequence from BLASTN output
- `blast_qseq_coverage`: The coverage of the query sequence from BLASTN output
- `blast_sseq_orfs_coverage`: The percentage coverage of ORFs in the subject sequence from BLASTN
- `orfs_start`: The starting point of the ORFs
- `orfs_end`: The ending point of the ORFs
- `blast_n_conseqs`: The number of consecutive sequences from BLASTN results

## `subtypes.fasta`

This file contains a list of reference sequences.
The names are exactly those `inferred_subtype`s found in `holistic.csv`.
Here is an example of the contents of the `subtypes.fasta` file:

```fasta
>Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.CfE
TGGAAGGGCTAATTC...GACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA
```

## `blast.csv`

This file contains output from BLASTN software.
It is produced conditionally.

Here is an example of the contents of the `blast.csv` file:

| qseqid      | sseqid                                      | sgi | qlen | slen | length | qstart | qend | sstart | send | evalue  | bitscore | pident  | nident | sstrand |
|-------------|---------------------------------------------|-----|------|------|--------|--------|------|--------|------|---------|----------|---------|--------|---------|
| KX505501.1  | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.CfE    | 0   | 1997 | 9718 | 1751   | 1      | 1746 | 455    | 2202 | 0.0     | 2186     | 93.946  | 1645   | plus    |
| KX505501.1  | Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.CfE    | 0   | 1997 | 9718 | 251    | 1747   | 1997 | 301    | 550  | 1.84e-85| 312      | 93.625  | 235    | plus    |

### Field descriptions

- `qseqid` - The identifier or name for the sequence (same as in the input FASTA file)
- `sseqid` - Means subject Seq-id. In our case, it's the subtype's sequence id and its the same as `inferred_subtype` of `holistic.csv`
- `sgi` - Means subject GI
- `qlen` - Length of the _query_ sequence
- `slen` - Means subject sequence length
- `length` - Means alignment length
- `qstart` - Means start of alignment in query
- `qend` - Means end of alignment in query
- `sstart` - Means start of alignment in subject
- `send` - Means end of alignment in subject
- `evalue` - Means expect value
- `bitscore` - Means bit score
- `pident` - Means percentage of identical matches
- `nident` - Means number of identical matches
- `sstrand` - Means subject strand

Consult [BLASTn documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/) for more details on these values.

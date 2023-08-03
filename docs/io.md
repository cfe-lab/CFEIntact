
# Inputs and Outputs

The only input is a path to a `.fasta` file that you want to analyze.
This file should include the genomic sequences of interest.

Once the analysis is complete, HIVIntact will generate four output files:

| Filename         | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| errors.json      | This file contains associations between sequences and their identified defects. |
| orfs.json        | This file contains associations between sequences and their identified Open Reading Frames (ORFs). |
| intact.fasta     | This file contains a list of sequences with no fatal defects identified. These sequences are considered putative intact genomes. |
| nonintact.fasta  | This file contains a list of sequences with identified defects. These sequences are not considered putative intact genomes. |

If you pass the `--output-csv` option to HIVIntact, the output format will be `.csv` instead of `.json`.

Here is an example of the contents of the `errors.csv` file:

| sequence_name | error                          | message                                                                                             |
|---------------|--------------------------------|-----------------------------------------------------------------------------------------------------|
| KX505501.1    | DeletionInOrf                  | ORF pol at 2084-5096 can have maximum deletions 30, got 2892                                    |
| KX505501.1    | RevResponseElementDeletion     | Query Sequence exceeds maximum deletion tolerance in RRE. Contains 35 deletions with max tolerance of 20 deletions. |
| MN691959      | InternalStopInOrf              | Smaller ORF vpu at 6060-6309 contains an internal stop codon                                     |
| MK114856.1    | APOBECHypermutationDetected    | Query sequence shows evidence of APOBEC3F/G-mediated hypermutation (p = 3.639064030015132e-65).  |
| MK116110.1    | PackagingSignalDeletion        | Query Sequence exceeds maximum deletion tolerance in PSI. Contains 93 deletions with max tolerance of 10 deletions. |


## Inputs and outputs

HIVIntact expects a filepath to a `.fasta` file to be provided.
This file should contain the sequences that we want to analyze.

Then after the run is finished, 4 files will be produced:

   [Filename]       [Description]
   errors.json      File containing associations of sequences to their list of identified defects.
   orfs.json        File containing associations of sequences to their list of identified ORFs.
   intact.fasta     File containing a list of files with no fatal defects identified.
   nonintact.fasta  File containing a list of files with no fatal defects.

If `--output-csv` option is passed to HIVIntact,
then the output format is `.csv` instead of `.json`.

Below is an example `errors.csv` contents:

      sequence_name,error,message
      KX505501.1,DeletionInOrf,"ORF pol at 2084-5096 can have maximum deletions 30, got 2892"
      KX505501.1,RevResponseElementDeletion,Query Sequence exceeds maximum deletion tolerance in RRE. Contains 35 deletions with max tolerance of 20 deletions.
      MN691959,InternalStopInOrf,Smaller ORF vpu at 6060-6309 contains an internal stop codon
      MK114856.1,APOBECHypermutationDetected,Query sequence shows evidence of APOBEC3F/G-mediated hypermutation (p = 3.639064030015132e-65).
      MK116110.1,PackagingSignalDeletion,Query Sequence exceeds maximum deletion tolerance in PSI. Contains 93 deletions with max tolerance of 10 deletions.

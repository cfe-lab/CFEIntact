# Constants for intactness
# In format start, end, maximum_allowed_deletion
# HXB2 coordinates, must be converted for other subtypes

# All coordinates one-based here.
# Note that because of the defective insertion in VPR of HXB2,
# the coordinates after that insertion are shifted by -1 (automatically).

from typing import List, Tuple

# (name, start, end, max_deletions, max_insertions, max_distance)
ORFDefinition = Tuple[str, int, int, int, int, float]
ORFsDefinition = List[ORFDefinition]

DEFAULT_FORWARD_ORFs: ORFsDefinition \
    = [("gag", 790, 2292, 30, 45, 0.36329365079365084),
       ("pol", 2085, 5096, 87, 51, 0.30019762845849796),
       ("env", 6225, 8795, 150, 123, 0.8311602209944755)]
DEFAULT_REVERSE_ORFS: ORFsDefinition = []
DEFAULT_SMALL_FORWARD_ORFS: ORFsDefinition \
    = [("vif", 5041, 5619, 15, 12, 0.59375),
       ("vpr", 5559, 5850, 16, 8, 0.5625),
       ("tat_exon1", 5831, 6045, 3, 9, 0.8949367088607598),
       ("rev_exon1", 5970, 6045, 3, 9, 1.1806451612903224),
       ("vpu", 6062, 6310, 15, 24, 1.117777777777778),
       ("tat_exon2", 8377, 8469, 24, 15, 1.4210526315789473),
       ("rev_exon2", 8378, 8653, 27, 39, 0.7577319587628868),
       ("nef", 8797, 9417, 27, 48, 0.8023584905660375)]
DEFAULT_SMALL_REVERSE_ORFS: ORFsDefinition = []

DEFAULT_ORF_LENGTH = 1000
DEFAULT_SMALL_ORF_LENGTH = 1
DEFAULT_ERROR_BAR = 1

# Major Splice Donor Site Location -- HXB2 coordinates only
DEFAULT_MSD_SITE_LOCUS = 743
DEFAULT_MSD_SEQUENCE = "GT"

# Packaging Signal Location -- HXB2 coordinates only
DEFAULT_PSI_LOCUS = (680, 809)
PSI_ERROR_TOLERANCE = 10  # 10 Nucleotide deletion tolerance

# Rev Response Element Location -- HXB2 coordinates only
DEFAULT_RRE_LOCUS = (7755, 8020)
RRE_ERROR_TOLERANCE = 20

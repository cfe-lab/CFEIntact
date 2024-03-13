# Constants for intactness
# In format start, end, maximum_allowed_deletion
# HXB2 coordinates, must be converted for other subtypes

# All coordinates one-based here.
# Note that because of the defective insertion in VPR of HXB2,
# the coordinates after that insertion are shifted by -1 (automatically).

from typing import List, Tuple

ORFDefinition = Tuple[str, int, int, int]
ORFsDefinition = List[ORFDefinition]

# (start, end, error_bar, breaks_intactness)
DEFAULT_FORWARD_ORFs: ORFsDefinition \
    = [("gag", 790, 2292, 30),
       ("pol", 2085, 5096, 30),
       ("env", 6225, 8795, 100)]
DEFAULT_REVERSE_ORFS: ORFsDefinition = []
DEFAULT_SMALL_FORWARD_ORFS: ORFsDefinition \
    = [("vif", 5041, 5619, 30),
       ("vpr", 5559, 5850, 30),
       ("tat_exon1", 5831, 6045, 30),
       ("rev_exon1", 5970, 6045, 30),
       ("vpu", 6062, 6310, 30),
       ("tat_exon2", 8377, 8469, 30),
       ("rev_exon2", 8378, 8653, 30),
       ("nef", 8797, 9417, 30)]
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

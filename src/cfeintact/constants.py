# Constants for intactness
# In format start, end, maximum_allowed_deletion
# HXB2 coordinates, must be converted for other subtypes

# All coordinates one-based here.
# Note that because of the defective insertion in VPR of HXB2,
# the coordinates after that insertion are shifted by -1 (automatically).

from typing import List, Tuple
from fractions import Fraction

# (name, start, end, max_deletions, max_insertions, max_distance)
ORFDefinition = Tuple[str, int, int, int, int, Fraction]
ORFsDefinition = List[ORFDefinition]

DEFAULT_FORWARD_ORFs: ORFsDefinition \
    = [("gag", 790, 2292, 42, 207, Fraction("0.6135702022840934")),
       ("pol", 2085, 5096, 93, 24, Fraction("0.7523414583707471")),
       ("env", 6225, 8795, 54, 459, Fraction("0.6862261858364505"))]
DEFAULT_REVERSE_ORFS: ORFsDefinition = []
DEFAULT_SMALL_FORWARD_ORFS: ORFsDefinition \
    = [("vif", 5041, 5619, 12, 9, Fraction("0.5596656062735806")),
       ("vpr", 5559, 5850, 6, 6, Fraction("0.4875248737180463")),
       ("tat_exon1", 5831, 6045, 0, 0, Fraction("0.5802194034967433")),
       ("rev_exon1", 5970, 6045, 0, 0, Fraction("0.5982053838484547")),
       ("vpu", 6062, 6310, 6, 24, Fraction("0.6351840675920339")),
       ("tat_exon2", 8377, 8469, 0, 15, Fraction("0.6561922365988909")),
       ("rev_exon2", 8378, 8653, 7, 9, Fraction("0.5499850344208321")),
       ("nef", 8797, 9417, 48, 54, Fraction("0.6839129079656255"))]
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

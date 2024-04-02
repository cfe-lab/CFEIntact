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
    = [("gag", 790, 2292, 39, 201, Fraction("0.5981496276061171")),
       ("pol", 2085, 5096, 91, 21, Fraction("0.7442043222003922")),
       ("env", 6225, 8795, 54, 363, Fraction("0.6303960531725368"))]
DEFAULT_REVERSE_ORFS: ORFsDefinition = []
DEFAULT_SMALL_FORWARD_ORFS: ORFsDefinition \
    = [("vif", 5041, 5619, 3, 3, Fraction("0.6257074581813602")),
       ("vpr", 5559, 5850, 6, 3, Fraction("0.4742978272390036")),
       ("tat_exon1", 5831, 6045, 0, 0, Fraction("0.5616889405384466")),
       ("rev_exon1", 5970, 6045, 0, 0, Fraction("0.6713780918727915")),
       ("vpu", 6062, 6310, 6, 24, Fraction("0.5982053838484547")),
       ("tat_exon2", 8377, 8469, 0, 3, Fraction("0.5396509491733008")),
       ("rev_exon2", 8378, 8653, 6, 3, Fraction("0.6868760338334845")),
       ("nef", 8797, 9417, 48, 51, Fraction("0.6840688395346594"))]
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

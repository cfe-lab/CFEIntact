
# Predefined limits

Many decisions in HIVIntact logic rely on predefined "limits".
This page lists all of them,
except for the standard constants such as the 5% limit for p-value.

## ORF lengths

| ORF name    | deletion tolerance | insertion tolerance |
|-------------|--------------------|---------------------|
| gag         | -30                | +90                 |
| pol         | -30                | +90                 |
| env         | -100               | +300                |
| vif         | -30                | +90                 |
| vpr         | -30                | +90                 |
| tat_exon1   | -30                | +90                 |
| rev_exon1   | -30                | +90                 |
| vpu         | -30                | +90                 |
| tat_exon2   | -30                | +90                 |
| rev_exon2   | -30                | +90                 |
| nef         | -30                | +90                 |

## Other

The PSI error tolerance is 10.

The RRE error tolerance is 20.

Minimum tolerable sequence length is 8000.

Maximum tolerable number of nucleotides impacted by frame shifts is 10% of the ORF's length.


# Predefined limits

Many decisions in CFEIntact logic rely on predefined "limits".
This page lists all of them,
except for the standard constants such as the 5% limit for p-value.

## ORF lengths

| ORF name  | deletion tolerance | insertion tolerance |
|-----------|--------------------|---------------------|
| gag       | -12                | +33                 |
| pol       | -21                | +24                 |
| env       | -141               | +72                 |
| vif       | -6                 | +3                  |
| vpr       | -15                | +3                  |
| tat_exon1 | -1                 | +3                  |
| rev_exon1 | -0                 | +3                  |
| vpu       | -9                 | +21                 |
| tat_exon2 | -0                 | +3                  |
| rev_exon2 | -4                 | +3                  |
| nef       | -42                | +45                 |

The values above are based on the following emperical data:

| ORF name    | Median | Minimum | Maximum | Mean   | Mode | Stdev | Variance | Count |
|-------------|--------|---------|---------|--------|------|-------|----------|-------|
| gag         | 1503   | 1491    | 1536    | 1507.9 | 1503 | 8.47  | 71.75    | 1265  |
| pol         | 3012   | 2991    | 3036    | 3013.92| 3012 | 6.42  | 41.22    | 628   |
| env         | 2571   | 2430    | 2643    | 2571.65| 2559 | 31.41 | 986.43   | 1468  |
| vif         | 579    | 573     | 582     | 579.04 | 579  | 0.89  | 0.79     | 6661  |
| vpr         | 291    | 276     | 294     | 290.54 | 291  | 1.87  | 3.48     | 1061  |
| tat_exon1   | 215    | 214     | 218     | 215.01 | 215  | 0.29  | 0.08     | 1019  |
| rev_exon1   | 76     | 76      | 79      | 76.03  | 76   | 0.27  | 0.07     | 1074  |
| vpu         | 246    | 237     | 267     | 246.49 | 246  | 4.5   | 20.25    | 1016  |
| tat_exon2   | 91     | 91      | 94      | 91.04  | 91   | 0.34  | 0.12     | 2146  |
| rev_exon2   | 275    | 271     | 278     | 275.01 | 275  | 0.46  | 0.21     | 2100  |
| nef         | 621    | 579     | 666     | 625.5  | 621  | 13.12 | 172.04   | 1042  |

<div style='width: 100%; text-align: right;'><tiny>(source: study by British Columbia Centre for Excellence in HIV/AIDS)</tiny></div>


The table above presents an analysis of Open Reading Frame (ORF) lengths for different genes. The analysis summarizes the lengths by providing statistical measures such as the median, minimum, maximum, mean, mode, standard deviation, variance, and the count of the samples.


Statistical measures help in understanding the distribution of the data. For instance, the mean gives an idea of the average ORF length, while the median provides a measure of the center of the distribution. The mode indicates the most frequently observed length, and the variance or standard deviation shows how spread apart the lengths are from the mean.


Each gene has different limits for deletion and insertion tolerance, as shown in the first table. Insights into the length distribution further our understanding of how these tolerances affect gene functionality. For example, if a gene has a high variance, it suggests that the gene can tolerate more variations in its length, hence a higher tolerance limit.

## Other

The PSI error tolerance is 10.

The RRE error tolerance is 20.

Minimum tolerable sequence length is 8000.

Maximum tolerable number of nucleotides impacted by frame shifts is 10% of the ORF's length.

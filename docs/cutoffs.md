
# Predefined limits

Many decisions in CFEIntact logic rely on predefined "limits".
This page lists all of them,
except for the standard constants such as the 5% limit for p-value.

## ORF lengths

| ORF name  | deletion tolerance | insertion tolerance |
|-----------|--------------------|---------------------|
| gag       | -39                | +201                |
| pol       | -91                | +21                 |
| env       | -54                | +363                |
| vif       | -3                 | +3                  |
| vpr       | -6                 | +3                  |
| tat_exon1 | -0                 | +0                  |
| rev_exon1 | -0                 | +0                  |
| vpu       | -6                 | +24                 |
| tat_exon2 | -0                 | +3                  |
| rev_exon2 | -6                 | +3                  |
| nef       | -48                | +51                 |

The values above are based on the following emperical data:

| ORF name  | Median | Minimum | Maximum |    Mean | Mode | Stdev | Count |
|-----------+--------+---------+---------+---------+------+-------+-------|
| gag       |   1508 |    1469 |    1709 | 1509.93 | 1502 | 18.79 |  1042 |
| pol       |   3011 |    2920 |    3032 | 3009.14 | 3011 | 22.94 |  1048 |
| env       |   2576 |    2522 |    2939 | 2581.93 | 2549 | 45.63 |  1042 |
| vif       |    578 |     575 |     581 |  577.98 |  578 |  0.99 |  1046 |
| vpr       |    290 |     284 |     293 |  290.01 |  290 |  1.28 |  1060 |
| tat_exon1 |    214 |     214 |     214 |     214 |  214 |   0.0 |  1061 |
| rev_exon1 |     75 |      75 |      75 |      75 |   75 |   0.0 |  1064 |
| vpu       |    245 |     239 |     269 |   247.2 |  245 |  7.38 |  1055 |
| tat_exon2 |     92 |      92 |      95 |   92.16 |   92 |  0.68 |  1057 |
| rev_exon2 |    275 |     269 |     278 |  275.07 |  275 |  0.53 |  1053 |
| nef       |    629 |     581 |     680 |  632.23 |  620 | 19.49 |  1050 |

<div style='width: 100%; text-align: right;'><tiny>(source: study by British Columbia Centre for Excellence in HIV/AIDS)</tiny></div>


The table above presents an analysis of Open Reading Frame (ORF) lengths for different genes. The analysis summarizes the lengths by providing statistical measures such as the median, minimum, maximum, mean, mode, standard deviation, variance, and the count of the samples.


Statistical measures help in understanding the distribution of the data. For instance, the mean gives an idea of the average ORF length, while the median provides a measure of the center of the distribution. The mode indicates the most frequently observed length, and the variance or standard deviation shows how spread apart the lengths are from the mean.


Each gene has different limits for deletion and insertion tolerance, as shown in the first table. Insights into the length distribution further our understanding of how these tolerances affect gene functionality. For example, if a gene has a high variance, it suggests that the gene can tolerate more variations in its length, hence a higher tolerance limit.

## Other

The PSI error tolerance is 10.

The RRE error tolerance is 20.

Minimum tolerable sequence length is 8000.

Maximum tolerable number of nucleotides impacted by frame shifts is 10% of the ORF's length.


# Predefined limits

Many decisions in CFEIntact logic rely on predefined "limits".
This page lists all of them,
except for the standard constants such as the 5% limit for p-value.

## ORF lengths

| ORF name  | Deletion tolerance | Insertion tolerance |
|-----------|--------------------|---------------------|
| gag       | -30                | +45                 |
| pol       | -87                | +51                 |
| env       | -150               | +123                |
| vif       | -15                | +12                 |
| vpr       | -16                | +8                  |
| tat_exon1 | -3                 | +9                  |
| rev_exon1 | -3                 | +9                  |
| vpu       | -15                | +24                 |
| tat_exon2 | -24                | +15                 |
| rev_exon2 | -27                | +39                 |
| nef       | -27                | +48                 |

The values above are based on the following emperical data:

| ORF name  | Median | Minimum | Maximum | Mean    | Stdev | Count |
|-----------|--------|---------|---------|---------|-------|-------|
| gag       | 1508   | 1466    | 1715    | 1509.82 | 22.38 | 1405  |
| pol       | 3011   | 2918    | 3035    | 3008.63 | 22.12 | 1410  |
| env       | 2576   | 2522    | 3035    | 2585.8  | 53.29 | 1408  |
| vif       | 578    | 566     | 587     | 578.15  | 1.97  | 1412  |
| vpr       | 290    | 284     | 296     | 289.78  | 1.53  | 1414  |
| tat_exon1 | 214    | 214     | 214     | 214     | 0     | 1428  |
| rev_exon1 | 75     | 75      | 75      | 75      | 0     | 1434  |
| vpu       | 245    | 239     | 269     | 246.95  | 6.71  | 1410  |
| tat_exon2 | 92     | 92      | 107     | 92.76   | 2.92  | 1420  |
| rev_exon2 | 275    | 268     | 284     | 276.32  | 3.34  | 1427  |
| nef       | 629    | 581     | 683     | 631.04  | 18.2  | 1408  |

  <div style='width: 100%; text-align: right;'><tiny>(source: study by British Columbia Centre for Excellence in HIV/AIDS)</tiny></div>

The table above presents an analysis of Open Reading Frame (ORF) lengths for different genes. The analysis summarizes the lengths by providing statistical measures such as the median, minimum, maximum, mean, mode, standard deviation, variance, and the count of the samples.

Then deletion and insertion tolerance is calculated as `|Median - Minimum|` and `|Median - Maximum|`.

Statistical measures help in understanding the distribution of the data. For instance, the mean gives an idea of the average ORF length, while the median provides a measure of the center of the distribution. The mode indicates the most frequently observed length, and the variance or standard deviation shows how spread apart the lengths are from the mean.


Each gene has different limits for deletion and insertion tolerance, as shown in the first table. Insights into the length distribution further our understanding of how these tolerances affect gene functionality. For example, if a gene has a high variance, it suggests that the gene can tolerate more variations in its length, hence a higher tolerance limit.

## Distance

The distance tolerance values for each ORF name as presented below are
crucial for determining the functional viability of variations within
the gene sequences. These tolerances demonstrate the degree to which
each ORF can deviate from its typical genetic make-up before the
changes are considered to significantly impact the gene's function or
structural integrity.

| ORF name  | Distance tolerance     |
|-----------|-----------------------|
| gag       | 0.36329365079365084    |
| pol       | 0.30019762845849796    |
| env       | 0.8311602209944755     |
| vif       | 0.59375                |
| vpr       | 0.5625                 |
| tat_exon1 | 0.8949367088607598     |
| rev_exon1 | 1.1806451612903224     |
| vpu       | 1.117777777777778      |
| tat_exon2 | 1.4210526315789473     |
| rev_exon2 | 0.7577319587628868     |
| nef       | 0.8023584905660375     |

The following table complements the distance tolerance data by
providing statistical measures on observed distances. These include
the mean, median, mode, minimum, maximum, standard deviation (Stdev),
and count metrics for each ORF. Such statistics support the determined
distance tolerances by reflecting the variability and the general
distribution of distances within the pool of observed instances.

| ORF name  | Median | Minimum | Maximum | Mean | Stdev | Count |
|-----------|--------|---------|---------|------|-------|-------|
| gag       | 0.23   | 0.15    | 0.98    | 0.26 | 0.1   | 1413  |
| pol       | 0.17   | 0.11    | 1.88    | 0.24 | 0.28  | 1408  |
| env       | 0.6    | 0.5     | 1.36    | 0.64 | 0.15  | 1411  |
| vif       | 0.36   | 0.19    | 0.79    | 0.37 | 0.08  | 1420  |
| vpr       | 0.31   | 0.16    | 0.59    | 0.31 | 0.09  | 1421  |
| tat_exon1 | 0.47   | 0.21    | 0.86    | 0.47 | 0.11  | 1419  |
| rev_exon1 | 0.46   | 0.23    | 0.92    | 0.49 | 0.14  | 1424  |
| vpu       | 0.78   | 0.55    | 1.08    | 0.77 | 0.1   | 1408  |
| tat_exon2 | 0.58   | 0.19    | 1.18    | 0.61 | 0.17  | 1418  |
| rev_exon2 | 0.48   | 0.16    | 0.76    | 0.48 | 0.12  | 1422  |
| nef       | 0.56   | 0.39    | 1.34    | 0.57 | 0.12  | 1409  |

   <div style='width: 100%; text-align: right;'><tiny>( source: https://github.com/cfe-lab/hivstats )</tiny></div>

The distance measures are an essential component of the genetic
analysis, providing insights into the physical distances within the
genetic material that can lead to functional or structural changes. By
understanding these tolerances and the distribution of natural
variances, researchers can better predict and analyse the impact of
genetic mutations or variations within the specified ORFs.

## Other

The PSI error tolerance is 10.

The RRE error tolerance is 20.

Minimum tolerable sequence length is 8000.

Maximum tolerable number of nucleotides impacted by frame shifts is 10% of the ORF's length.

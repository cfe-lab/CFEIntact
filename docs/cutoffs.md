
# Predefined limits

Many decisions in CFEIntact logic rely on predefined "limits".
This page lists all of them,
except for the standard constants such as the 5% limit for p-value.

## ORF lengths

| ORF name  | Deletion tolerance | Insertion tolerance |
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

| ORF name  | Mean    | Minimum | Maximum | Median | Mode | Stdev | Count |
|-----------|---------|---------|---------|--------|------|-------|-------|
| gag       | 1509.93 | 1469    | 1709    | 1508   | 1502 | 18.79 | 1042  |
| pol       | 3009.14 | 2920    | 3032    | 3011   | 3011 | 22.94 | 1048  |
| vpu       | 247.2   | 239     | 269     | 245    | 245  | 7.38  | 1055  |
| vif       | 577.98  | 575     | 581     | 578    | 578  | 0.99  | 1046  |
| vpr       | 290.01  | 284     | 293     | 290    | 290  | 1.28  | 1060  |
| tat_exon1 | 214     | 214     | 214     | 214    | 214  | 0.0   | 1061  |
| tat_exon2 | 92.16   | 92      | 95      | 92     | 92   | 0.68  | 1057  |
| rev_exon1 | 75      | 75      | 75      | 75     | 75   | 0.0   | 1064  |
| rev_exon2 | 275.07  | 269     | 278     | 275    | 275  | 0.53  | 1053  |
| env       | 2581.93 | 2522    | 2939    | 2576   | 2549 | 45.63 | 1042  |
| nef       | 632.23  | 581     | 680     | 629    | 620  | 19.49 | 1050  |

  <div style='width: 100%; text-align: right;'><tiny>(source: study by British Columbia Centre for Excellence in HIV/AIDS)</tiny></div>


The table above presents an analysis of Open Reading Frame (ORF) lengths for different genes. The analysis summarizes the lengths by providing statistical measures such as the median, minimum, maximum, mean, mode, standard deviation, variance, and the count of the samples.

Then deletion and insertion tolerance is calculated as `|Median - Minimum|` and `|Median - Maximum|`.
Note that `InsertionInOrf` error is only reported when the distance measure is also crossed.

Statistical measures help in understanding the distribution of the data. For instance, the mean gives an idea of the average ORF length, while the median provides a measure of the center of the distribution. The mode indicates the most frequently observed length, and the variance or standard deviation shows how spread apart the lengths are from the mean.


Each gene has different limits for deletion and insertion tolerance, as shown in the first table. Insights into the length distribution further our understanding of how these tolerances affect gene functionality. For example, if a gene has a high variance, it suggests that the gene can tolerate more variations in its length, hence a higher tolerance limit.

## Distance

The distance tolerance values for each ORF name as presented below are
crucial for determining the functional viability of variations within
the gene sequences. These tolerances demonstrate the degree to which
each ORF can deviate from its typical genetic make-up before the
changes are considered to significantly impact the gene's function or
structural integrity.

| ORF name  | Distance tolerance |
|-----------|--------------------|
| gag       | 0.5981496276061171 |
| pol       | 0.7442043222003922 |
| env       | 0.6303960531725368 |
| vif       | 0.6257074581813602 |
| vpr       | 0.4742978272390036 |
| tat_exon1 | 0.5616889405384466 |
| rev_exon1 | 0.6713780918727915 |
| vpu       | 0.5982053838484547 |
| tat_exon2 | 0.5396509491733008 |
| rev_exon2 | 0.6868760338334845 |
| nef       | 0.6840688395346594 |

The following table complements the distance tolerance data by
providing statistical measures on observed distances. These include
the mean, median, mode, minimum, maximum, standard deviation (Stdev),
and count metrics for each ORF. Such statistics support the determined
distance tolerances by reflecting the variability and the general
distribution of distances within the pool of observed instances.

| ORF name  | Mean | Median | Mode | Stdev | Minimum | Maximum | Count |
|-----------|------|--------|------|-------|---------|---------|-------|
| gag       | 0.28 | 0.27   | 0.26 | 0.06  | 0.20    | 0.60    | 1046  |
| pol       | 0.24 | 0.21   | 0.19 | 0.11  | 0.15    | 0.74    | 1046  |
| vpu       | 0.55 | 0.55   | 0.51 | 0.03  | 0.47    | 0.63    | 1042  |
| vif       | 0.36 | 0.36   | 0.37 | 0.05  | 0.23    | 0.63    | 1051  |
| vpr       | 0.33 | 0.34   | 0.38 | 0.06  | 0.23    | 0.47    | 1051  |
| tat_exon1 | 0.42 | 0.43   | 0.48 | 0.06  | 0.25    | 0.56    | 1058  |
| tat_exon2 | 0.47 | 0.48   | 0.44 | 0.07  | 0.24    | 0.67    | 1048  |
| rev_exon1 | 0.43 | 0.43   | 0.43 | 0.07  | 0.27    | 0.60    | 1053  |
| rev_exon2 | 0.42 | 0.42   | 0.41 | 0.07  | 0.21    | 0.54    | 1054  |
| env       | 0.5  | 0.49   | 0.45 | 0.05  | 0.44    | 0.69    | 1043  |
| nef       | 0.47 | 0.47   | 0.45 | 0.04  | 0.40    | 0.68    | 1043  |

   <div style='width: 100%; text-align: right;'><tiny>(source: study by British Columbia Centre for Excellence in HIV/AIDS)</tiny></div>

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

Codebook created on 2023-03-08 at 2023-03-08 09:07:53
================

A codebook contains documentation and metadata describing the contents,
structure, and layout of a data file.

## Dataset description

The data contains 6329 cases and 14 variables.

## Codebook

| name        | type    |    n | missing | unique |   mean | median |   mode |     sd | min |    max |  range | skew | skew_2se |  kurt | kurt_2se |
|:------------|:--------|-----:|--------:|-------:|-------:|-------:|-------:|-------:|----:|-------:|-------:|-----:|---------:|------:|---------:|
| ppnr        | numeric | 6329 |     0.0 |     95 | 274.81 | 281.00 | 281.00 | 196.15 | 2.0 | 603.00 | 601.00 | 0.02 |     0.27 | -1.49 |   -12.11 |
| triggerid   | numeric | 6329 |     0.0 |     76 |  34.84 |  35.00 |  35.00 |  19.28 | 2.0 |  77.00 |  75.00 | 0.01 |     0.17 | -1.18 |    -9.62 |
| day         | integer | 6329 |     0.0 |     31 |  15.12 |  15.00 |  15.00 |   8.97 | 1.0 |  31.00 |  30.00 | 0.08 |     1.35 | -1.32 |   -10.68 |
| ANGRY       | numeric | 5721 |     0.1 |    101 |  14.33 |   8.00 |   8.00 |  17.44 | 1.0 | 100.00 |  99.00 | 2.15 |    33.26 |  4.89 |    37.74 |
| SAD         | numeric | 5719 |     0.1 |    100 |  18.01 |  11.00 |  11.00 |  20.35 | 1.0 | 100.00 |  99.00 | 1.71 |    26.41 |  2.49 |    19.21 |
| ANX         | numeric | 5720 |     0.1 |     98 |  13.46 |   8.00 |   8.00 |  16.63 | 1.0 | 100.00 |  99.00 | 2.35 |    36.27 |  5.99 |    46.24 |
| DEPRE       | numeric | 5720 |     0.1 |    101 |  17.03 |  10.00 |  10.00 |  20.47 | 1.0 | 100.00 |  99.00 | 1.89 |    29.24 |  3.26 |    25.18 |
| RUMI.       | numeric | 5719 |     0.1 |    101 |  27.29 |  18.00 |  18.00 |  26.07 | 1.0 | 100.00 |  99.00 | 1.01 |    15.62 | -0.02 |    -0.18 |
| DIST.       | numeric | 5719 |     0.1 |    101 |  29.41 |  20.00 |  20.00 |  26.43 | 1.0 | 100.00 |  99.00 | 0.81 |    12.51 | -0.53 |    -4.11 |
| RFLCT.      | numeric | 5719 |     0.1 |    101 |  23.20 |  16.00 |  16.00 |  22.11 | 1.0 | 100.00 |  99.00 | 1.15 |    17.80 |  0.56 |     4.29 |
| REAPP.      | numeric | 5719 |     0.1 |     99 |  18.27 |  12.00 |  12.00 |  19.07 | 1.0 | 100.00 |  99.00 | 1.51 |    23.39 |  1.84 |    14.23 |
| SUPP.       | numeric | 5721 |     0.1 |    101 |  23.89 |  15.00 |  15.00 |  24.36 | 1.0 | 100.00 |  99.00 | 1.24 |    19.13 |  0.55 |     4.28 |
| SOCSH.      | numeric | 5719 |     0.1 |    101 |  21.49 |  12.00 |  12.00 |  23.90 | 1.0 | 100.00 |  99.00 | 1.32 |    20.37 |  0.71 |     5.52 |
| person_CESD | numeric | 6329 |     0.0 |     25 |   0.95 |   0.95 |   0.95 |   0.28 | 0.5 |   2.05 |   1.55 | 1.11 |    18.01 |  1.97 |    16.02 |

### Legend

-   **Name**: Variable name
-   **type**: Data type of the variable
-   **missing**: Proportion of missing values for this variable
-   **unique**: Number of unique values
-   **mean**: Mean value
-   **median**: Median value
-   **mode**: Most common value (for categorical variables, this shows
    the frequency of the most common category)
-   **mode_value**: For categorical variables, the value of the most
    common category
-   **sd**: Standard deviation (measure of dispersion for numerical
    variables
-   **v**: Agresti’s V (measure of dispersion for categorical variables)
-   **min**: Minimum value
-   **max**: Maximum value
-   **range**: Range between minimum and maximum value
-   **skew**: Skewness of the variable
-   **skew_2se**: Skewness of the variable divided by 2\*SE of the
    skewness. If this is greater than abs(1), skewness is significant
-   **kurt**: Kurtosis (peakedness) of the variable
-   **kurt_2se**: Kurtosis of the variable divided by 2\*SE of the
    kurtosis. If this is greater than abs(1), kurtosis is significant.

This codebook was generated using the [Workflow for Open Reproducible
Code in Science (WORCS)](https://osf.io/zcvbs/)

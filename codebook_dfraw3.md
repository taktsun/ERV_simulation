Codebook created on 2023-03-08 at 2023-03-08 09:07:58
================

A codebook contains documentation and metadata describing the contents,
structure, and layout of a data file.

## Dataset description

The data contains 14098 cases and 14 variables.

## Codebook

| name        | type    |     n | missing | unique |   mean | median |  mode |    sd |   min |   max | range | skew | skew_2se |  kurt | kurt_2se |
|:------------|:--------|------:|--------:|-------:|-------:|-------:|------:|------:|------:|------:|------:|-----:|---------:|------:|---------:|
| ppnr        | numeric | 14098 |    0.00 |    200 | 256.19 |  256.0 | 256.0 | 95.57 | 101.0 | 801.0 | 700.0 | 0.82 |    19.76 |  3.54 |    42.91 |
| triggerid   | numeric | 14098 |    0.00 |     80 |  38.15 |   38.0 |  38.0 | 20.50 |   1.0 |  80.0 |  79.0 | 0.01 |     0.29 | -1.16 |   -14.07 |
| day         | integer | 14098 |    0.00 |     31 |  14.04 |   13.0 |  13.0 |  8.74 |   1.0 |  31.0 |  30.0 | 0.29 |     6.99 | -1.11 |   -13.44 |
| kwaad       | numeric | 12293 |    0.13 |    102 |  11.87 |    7.0 |   7.0 | 16.30 |   0.0 | 100.0 | 100.0 | 2.62 |    59.38 |  7.96 |    90.08 |
| droevig     | numeric | 12293 |    0.13 |    102 |  13.11 |    8.0 |   8.0 | 17.13 |   0.0 | 100.0 | 100.0 | 2.30 |    52.04 |  5.78 |    65.43 |
| angstig     | numeric | 12293 |    0.13 |    100 |  10.27 |    7.0 |   7.0 | 13.76 |   0.0 | 100.0 | 100.0 | 2.81 |    63.49 |  9.96 |   112.76 |
| depressief  | numeric | 12293 |    0.13 |    101 |  12.62 |    8.0 |   8.0 | 16.72 |   0.0 | 100.0 | 100.0 | 2.32 |    52.57 |  5.93 |    67.14 |
| Rumi_past   | numeric | 12293 |    0.13 |    101 |  18.35 |   10.0 |  10.0 | 22.28 |   0.0 | 100.0 | 100.0 | 1.64 |    37.18 |  1.98 |    22.41 |
| Rumi_future | numeric | 12293 |    0.13 |    102 |  31.12 |   21.0 |  21.0 | 28.58 |   0.0 | 100.0 | 100.0 | 0.65 |    14.82 | -0.81 |    -9.22 |
| Dist        | numeric | 12293 |    0.13 |    102 |  24.36 |   12.0 |  12.0 | 27.05 |   0.0 | 100.0 | 100.0 | 1.19 |    26.94 |  0.35 |     3.95 |
| Reap        | numeric | 12293 |    0.13 |    100 |  15.09 |    9.0 |   9.0 | 17.86 |   0.0 | 100.0 | 100.0 | 1.75 |    39.58 |  2.93 |    33.14 |
| Supp        | numeric | 12293 |    0.13 |    102 |  20.14 |   10.0 |  10.0 | 24.60 |   0.0 | 100.0 | 100.0 | 1.60 |    36.15 |  1.72 |    19.43 |
| Sosu        | numeric | 12293 |    0.13 |    102 |  18.02 |   10.0 |  10.0 | 21.47 |   0.0 | 100.0 | 100.0 | 1.58 |    35.79 |  1.86 |    21.03 |
| person_CESD | numeric | 14098 |    0.00 |     23 |   0.92 |    0.9 |   0.9 |  0.23 |   0.5 |   1.7 |   1.2 | 0.77 |    18.56 |  0.38 |     4.56 |

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

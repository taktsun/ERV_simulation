Codebook created on 2023-03-08 at 2023-03-08 08:54:04
================

A codebook contains documentation and metadata describing the contents,
structure, and layout of a data file.

## Dataset description

The data contains 5040 cases and 13 variables.

## Codebook

| name        | type    |    n | missing | unique |  mean | median |  mode |    sd |  min |   max | range |  skew | skew_2se |  kurt | kurt_2se |
|:------------|:--------|-----:|--------:|-------:|------:|-------:|------:|------:|-----:|------:|------:|------:|---------:|------:|---------:|
| ppnr        | numeric | 5040 |    0.00 |     70 | 35.50 |  35.50 | 35.50 | 20.21 | 1.00 | 70.00 |  69.0 |  0.00 |     0.00 | -1.20 |    -8.71 |
| triggerid   | numeric | 5040 |    0.00 |     72 | 35.50 |  35.50 | 35.50 | 20.78 | 0.00 | 71.00 |  71.0 |  0.00 |     0.00 | -1.20 |    -8.71 |
| day         | integer | 5040 |    0.00 |     31 | 17.80 |  21.00 | 21.00 |  9.61 | 1.00 | 31.00 |  30.0 | -0.37 |    -5.32 | -1.29 |    -9.33 |
| M_nerv      | numeric | 3824 |    0.24 |      8 |  1.64 |   1.00 |  1.00 |  1.57 | 0.00 |  6.00 |   6.0 |  0.68 |     8.53 | -0.48 |    -3.00 |
| M_nieder    | numeric | 3826 |    0.24 |      8 |  1.17 |   1.00 |  1.00 |  1.44 | 0.00 |  6.00 |   6.0 |  1.22 |    15.45 |  0.80 |     5.06 |
| M_bekuem    | numeric | 3825 |    0.24 |      8 |  1.38 |   1.00 |  1.00 |  1.53 | 0.00 |  6.00 |   6.0 |  0.96 |    12.16 |  0.07 |     0.45 |
| M_rum2      | numeric | 3812 |    0.24 |      8 |  1.98 |   2.00 |  2.00 |  1.75 | 0.00 |  6.00 |   6.0 |  0.43 |     5.46 | -0.94 |    -5.91 |
| M_rum1      | numeric | 3817 |    0.24 |      8 |  1.62 |   1.00 |  1.00 |  1.57 | 0.00 |  6.00 |   6.0 |  0.74 |     9.29 | -0.36 |    -2.26 |
| M_distr2    | numeric | 3814 |    0.24 |      8 |  1.45 |   1.00 |  1.00 |  1.74 | 0.00 |  6.00 |   6.0 |  1.02 |    12.88 | -0.07 |    -0.46 |
| M_distr1    | numeric | 3812 |    0.24 |      8 |  1.19 |   1.00 |  1.00 |  1.54 | 0.00 |  6.00 |   6.0 |  1.26 |    15.92 |  0.69 |     4.38 |
| M_refl2     | numeric | 3812 |    0.24 |      8 |  2.66 |   3.00 |  3.00 |  1.55 | 0.00 |  6.00 |   6.0 | -0.01 |    -0.10 | -0.59 |    -3.74 |
| M_refl1     | numeric | 3813 |    0.24 |      8 |  2.29 |   2.00 |  2.00 |  1.54 | 0.00 |  6.00 |   6.0 |  0.22 |     2.77 | -0.56 |    -3.56 |
| person_CESD | numeric | 5040 |    0.00 |     28 |  1.40 |   1.35 |  1.35 |  0.40 | 0.65 |  2.85 |   2.2 |  0.70 |    10.18 |  0.88 |     6.35 |

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

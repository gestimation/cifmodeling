# Data from a prostate cancer trial in Byer & Green (1980)

Anonymized data from a randomized clinical trial of prostate cancer
published in Byer & Green (1980).

## Usage

``` r
data(prostate)
```

## Format

A data frame with 502 observations and 18 variables, including:

- dtime:

  Follow-up time in months.

- status:

  Event status ("alive", "dead - prostatic ca", "dead - other ca",
  "dead - heart or vascular", "dead - cerebrovascular").

- rx:

  Treatment assignment to diethylstilbestrol (DES) or a placebo.

- age:

  Age at baseline (years).

- wt:

  Weight in pounds.

- pf:

  Performance status.

- hx:

  History of cardiovascular disease.

- sbp:

  Systolic blood pressure.

- dbp:

  Diastolic blood pressure.

- ekg:

  Electrocardiogram category.

- hg:

  Hemoglobin level.

- sz:

  Size of the primary tumor.

- sg:

  Stage/grade of disease.

- ap:

  Serum acid phosphatase.

- bm:

  Bone metastases indicator.

- stage:

  Clinical stage.

- sdate:

  Start date.

- patno:

  Patient number.

## Source

Byer, D. P. & Green, S. B. (1980), 'Prognostic variables for survival in
a randomized comparison of treatments for prostatic cancer', Bulletin du
Cancer 67, 477-488

## Details

The dataset records follow-up for cause of death together with treatment
assignment and baseline characteristics. It is used in the package
documentation to illustrate stratified cumulative incidence analyses.

## Examples

``` r
data(prostate)
head(prostate)
#>   patno stage              rx dtime                 status age  wt
#> 1     1     3 0.2 mg estrogen    72                  alive  75  76
#> 2     2     3 0.2 mg estrogen     1        dead - other ca  54 116
#> 3     3     3 5.0 mg estrogen    40 dead - cerebrovascular  69 102
#> 4     4     3 0.2 mg estrogen    20 dead - cerebrovascular  75  94
#> 5     5     3         placebo    65                  alive  67  99
#> 6     6     3 0.2 mg estrogen    24    dead - prostatic ca  71  98
#>                     pf hx sbp dbp                           ekg       hg sz sg
#> 1      normal activity  0  15   9                  heart strain 13.79883  2  8
#> 2      normal activity  0  13   7 heart block or conduction def 14.59961 42 NA
#> 3      normal activity  1  14   8                  heart strain 13.39844  3  9
#> 4 in bed < 50% daytime  1  14   7                        benign 17.59766  4  8
#> 5      normal activity  0  17  10                        normal 13.39844 34  8
#> 6      normal activity  0  19  10                        normal 15.09961 10 11
#>          ap bm      sdate
#> 1 0.2999878  0 1977-08-10
#> 2 0.6999512  0 1977-09-21
#> 3 0.2999878  0 1978-01-12
#> 4 0.8999023  0 1978-03-19
#> 5 0.5000000  0 1978-03-22
#> 6 0.5999756  0 1978-06-14
```

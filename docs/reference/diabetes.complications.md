# Data from a cohort study of patients with type 2 diabetes

Anonymized data from a cohort study of patients with type 2 diabetes
followed for ocular and macro-vascular complications.

## Usage

``` r
data(diabetes.complications)
```

## Format

A data frame with 978 observations and 19 variables:

- t:

  Follow-up time in years.

- epsilon:

  Event type indicator (0 = censored, 1 = diabetic retinopathy, 2 =
  macro-vascular complication).

- fruit:

  Fruit intake (g/day).

- fruitq:

  Quartile of fruit intake.

- fruitq1:

  Binary indicator for low fruit intake.

- strata:

  Stratum used for inverse probability of censoring weights.

- age:

  Age at baseline (years).

- sex:

  Sex coded as 0 = woman, 1 = man.

- bmi:

  Body mass index at baseline.

- hba1c:

  Hemoglobin A1c (%).

- diabetes_duration:

  Duration of diabetes (years).

- drug_oha:

  Indicator for oral hypoglycemic agent use.

- drug_insulin:

  Indicator for insulin use.

- sbp:

  Systolic blood pressure (mmHg).

- ldl:

  Low-density lipoprotein cholesterol (mg/dL).

- hdl:

  High-density lipoprotein cholesterol (mg/dL).

- tg:

  Triglycerides (mg/dL).

- current_smoker:

  Indicator for current smoking status.

- alcohol_drinker:

  Indicator for current alcohol drinking.

- ltpa:

  Leisure-time physical activity (METs).

## Source

Anonymized data supplied with the package for documentation and
demonstration purposes.

## Details

The variables include follow-up time, cause-specific event indicators,
exposure indicators for fruit intake, censoring strata, and a set of
covariates used in the package vignettes.

## Examples

``` r
data(diabetes.complications)
str(diabetes.complications)
#> 'data.frame':    978 obs. of  20 variables:
#>  $ t                : num  8.62 8.51 7.79 8.91 8.94 ...
#>  $ epsilon          : int  0 0 0 0 0 1 0 1 0 1 ...
#>  $ strata           : int  1 4 3 1 2 3 1 3 3 1 ...
#>  $ fruit            : num  75 26.8 64.3 5.35 211.05 ...
#>  $ fruitq1          : int  0 1 0 1 0 1 0 0 0 0 ...
#>  $ age              : int  45 68 63 49 55 61 56 64 67 58 ...
#>  $ sex              : int  0 0 0 0 0 0 0 1 1 1 ...
#>  $ bmi              : num  21.5 18.3 23.9 22.9 18.7 23.4 20.1 28.6 25.6 22.8 ...
#>  $ hba1c            : num  6.97 8.02 6.89 7.24 8.28 ...
#>  $ diabetes_duration: num  4.2 2.9 14.3 4.2 16.3 8.9 13 6.1 20.3 13.3 ...
#>  $ drug_oha         : int  0 0 1 1 1 1 0 1 0 0 ...
#>  $ drug_insulin     : int  0 1 0 0 0 0 1 0 0 1 ...
#>  $ sbp              : int  124 128 164 126 136 146 118 136 136 142 ...
#>  $ ldl              : num  187.3 87.6 74.6 95.7 50.5 ...
#>  $ hdl              : num  58.1 57.2 35 34.7 55.7 36.7 57.7 52.2 60.5 55.4 ...
#>  $ tg               : num  123 71 252 83 139 181 71 57 91 115 ...
#>  $ current_smoker   : int  0 1 1 1 1 0 1 0 0 0 ...
#>  $ alcohol_drinker  : int  0 0 1 0 0 0 0 0 0 0 ...
#>  $ ltpa             : num  52.5 11.03 4.38 9.38 12.38 ...
#>  $ fruitq           : Factor w/ 4 levels "Q1","Q2","Q3",..: 2 1 2 1 4 1 4 2 2 2 ...
```

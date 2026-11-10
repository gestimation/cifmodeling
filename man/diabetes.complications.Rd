\name{diabetes.complications}
\alias{diabetes.complications}
\docType{data}
\title{Data from a cohort study of patients with type 2 diabetes}

\format{
A data frame with 978 observations and 19 variables:
\describe{
  \item{t}{Follow-up time in years.}
  \item{epsilon}{Event type indicator (0 = censored, 1 = diabetic retinopathy, 2 = macro-vascular complication).}
  \item{fruit}{Fruit intake (g/day).}
  \item{fruitq}{Quartile of fruit intake.}
  \item{fruitq1}{Binary indicator for low fruit intake.}
  \item{strata}{Stratum used for inverse probability of censoring weights.}
  \item{age}{Age at baseline (years).}
  \item{sex}{Sex coded as 0 = woman, 1 = man.}
  \item{bmi}{Body mass index at baseline.}
  \item{hba1c}{Hemoglobin A1c (\%).}
  \item{diabetes_duration}{Duration of diabetes (years).}
  \item{drug_oha}{Indicator for oral hypoglycemic agent use.}
  \item{drug_insulin}{Indicator for insulin use.}
  \item{sbp}{Systolic blood pressure (mmHg).}
  \item{ldl}{Low-density lipoprotein cholesterol (mg/dL).}
  \item{hdl}{High-density lipoprotein cholesterol (mg/dL).}
  \item{tg}{Triglycerides (mg/dL).}
  \item{current_smoker}{Indicator for current smoking status.}
  \item{alcohol_drinker}{Indicator for current alcohol drinking.}
  \item{ltpa}{Leisure-time physical activity (METs).}
}}
\source{
Anonymized example data from the Japan Diabetes Complications Study.
The original cohort analysis of fruit intake and incident diabetic retinopathy
was reported by Tanaka S, Yoshimura Y, Kawasaki R, Kamada C, Tanaka S,
Horikawa C, Ohashi Y, Araki A, Ito H, Akanuma Y, Yamada N, Yamashita H,
Sone H, for the Japan Diabetes Complications Study Group.
Fruit intake and incident diabetic retinopathy with type 2 diabetes.
\emph{Epidemiology}. 2013, 24(2), 204--211.
\doi{10.1097/EDE.0b013e318281725e}
}

\usage{data(diabetes.complications)}

\description{
Anonymized data from a cohort study of patients with type 2 diabetes followed for ocular and macro-vascular complications.
}

\details{
The variables include follow-up time, cause-specific event indicators, exposure indicators for fruit intake, censoring strata, and a set of covariates used in the package vignettes.
}

\examples{
data(diabetes.complications)
str(diabetes.complications)
}

\keyword{datasets}

\name{prostate}
\alias{prostate}
\docType{data}
\title{Prostate cancer trial data from Byer & Green (1980)}
\format{A data frame with 502 observations and 16 variables, including:
\describe{
  \item{dtime}{Follow-up time in days.}
  \item{status}{Event status ("alive", "dead - prostatic ca", "dead - other ca", "dead - heart or vascular", "dead - cerebrovascular").}
  \item{rx}{Treatment assignment to diethylstilbestrol (DES) or a placebo.}
  \item{age}{Age at baseline (years).}
  \item{wt}{Weight in pounds.}
  \item{pf}{Performance status.}
  \item{hx}{History of cardiovascular disease.}
  \item{sbp}{Systolic blood pressure.}
  \item{dbp}{Diastolic blood pressure.}
  \item{ekg}{Electrocardiogram category.}
  \item{hg}{Hemoglobin level.}
  \item{sz}{Size of the primary tumor.}
  \item{sg}{Stage/grade of disease.}
  \item{ap}{Serum acid phosphatase.}
  \item{bm}{Bone metastases indicator.}
  \item{stage}{Clinical stage.}
  \item{sdate}{Start date.}
  \item{patno}{Patient number.}
}}
\source{Byer, D. P. & Green, S. B. (1980), 'Prognostic variables for survival in a randomized comparison of treatments for prostatic cancer', Bulletin du Cancer 67, 477-488}
\usage{data(prostate)}
\description{Anonymized data from a randomized clinical trial of prostate cancer published in Byer & Green (1980). }
\details{The dataset records follow-up for cause of death together with treatment assignment and baseline characteristics. It is used in the package documentation to illustrate stratified cumulative incidence analyses.}
\examples{
data(prostate)
head(prostate)
}
\keyword{datasets}

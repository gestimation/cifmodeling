# Package index

## Nonparametric estimation and visualization

- [`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
  : Calculate the Kaplan-Meier estimator and the Aalen-Johansen
  estimator
- [`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md)
  : Generate a survival/CIF curve with marks that represent censoring,
  competing risks and intercurrent events
- [`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md)
  : Arrange multiple survival and CIF plots in a panel display
- [`cifflowchart()`](https://gestimation.github.io/cifmodeling/reference/cifflowchart.md)
  : Create a flowchart of exclusions, treatment groups, and outcome
  status

## Direct polytomous regression

- [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
  : Fit coherent regression models of CIFs using polytomous log odds
  products
- [`coef(`*`<polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`vcov(`*`<polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`nobs(`*`<polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`summary(`*`<polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`print(`*`<summary.polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`effect_label.polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`tidy(`*`<polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`glance(`*`<polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  [`augment(`*`<polyreg>`*`)`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  : Methods for polyreg objects

## Data

- [`diabetes.complications`](https://gestimation.github.io/cifmodeling/reference/diabetes.complications.md)
  : Data from a cohort study of patients with type 2 diabetes
- [`prostate`](https://gestimation.github.io/cifmodeling/reference/prostate.md)
  : Data from a prostate cancer trial in Byer & Green (1980)

## Helpers

- [`extract_time_to_event()`](https://gestimation.github.io/cifmodeling/reference/extract_time_to_event.md)
  : Extract per-stratum event times from a formula and data
- [`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
  : Create a survival or competing-risks response

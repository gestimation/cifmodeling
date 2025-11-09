# cifmodeling 0.7.0

* Finalize API: argument names switched to dot.case; outcome.type values standardized to lower-hyphen.

# cifmodeling 0.6.0

* Initial CRAN submission.

# cifmodeling 0.4.0

* Bug fixes
  * `order.strata` is now applied via discrete scale limits so legend and fill ordering remain stable, including when `panel.per.variable = TRUE`.
  * Numeric stratification variables are coerced to factors when they have fewer than nine unique values and are split at the median (Below/Above median) when they have nine or more unique values.
  * `cifplot(panel.per.event = TRUE)` keeps both event-specific `ggplot` objects in `attr(x, "plots")` while returning the combined panel.
  * Default y-axis labels consistently use "Survival" for survival outcomes and "Cumulative incidence" phrasing for competing-risk outcomes.

# cifmodeling 0.2.0

* Maintenance updates.

# cifmodeling 0.1.0

* Initial CRAN submission.

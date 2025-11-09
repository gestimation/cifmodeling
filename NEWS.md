<<<<<<< HEAD
# cifmodeling 0.6.0
=======
# cifmodeling 0.7.0

# cifmodeling 0.4.0

* Bug fixes
  * `order.strata` is now applied via discrete scale limits so legend and fill
    ordering remain stable, including when `printEachVar = TRUE`.
  * Numeric stratification variables are coerced to factors when they have fewer
    than nine unique values and are split at the median (Below/Above median)
    when they have nine or more unique values.
  * `cifplot(printEachEvent = TRUE)` keeps both event-specific `ggplot`
    objects in `attr(x, "plots")` while returning the combined panel.
  * Default y-axis labels consistently use "Survival" for survival outcomes and
    "Cumulative incidence" phrasing for competing-risk outcomes.

# cifmodeling 0.2.0

3# cifmodeling 0.1.0
>>>>>>> 8cb48b52ceae88e76eaa96eb4e5a44b4956cdada

* Initial CRAN submission.

#' These arguments are shared by [cifplot()] and [cifpanel()].
#'
#' @name cif-visual-arguments
#' @keywords internal
#'
#' @param add.risktable Logical; if `TRUE`, adds a numbers-at-risk table under the plot.
#'   Default `TRUE`. **Note:** when a panel mode is active, tables are suppressed.
#' @param add.estimate.table Logical; if `TRUE`, adds a table of estimates and CIs.
#'   Default `FALSE`. **Note:** when a panel mode is active, tables are suppressed.
#' @param symbol.risk.table Character specifying the symbol used in the risk table to denote
#'   strata: `"square"`, `"circle"`, or `"triangle"` (default `"square"`).
#' @param font.size.risk.table Numeric font size for texts in risk / estimate tables (default `3`).
#' @param label.strata Character vector or named character vector specifying labels for strata.
#'   Names (if present) must match the (re-ordered) underlying strata levels.
#'   **Note:** when any of the panel modes is active
#'   (`panel.per.variable = TRUE`, `panel.per.event = TRUE`, `panel.censoring = TRUE`,
#'   or `panel.mode = "auto"` and it actually dispatches to a panel),
#'   strata labels are suppressed to avoid duplicated legends across sub-plots.
#' @param level.strata Optional character vector giving the full set of expected strata levels.
#'   When provided, both `order.strata` and `label.strata` are validated against it
#'   before application.
#' @param order.strata Optional character vector specifying the display order of strata
#'   in the legend/number-at-risk table. Specify the levels of strata. Levels not listed are dropped.
#' @param legend.position Character specifying the legend position:
#'   `"top"`, `"right"`, `"bottom"`, `"left"`, or `"none"` (default `"top"`).
#' @param type.y Character string specifying the y-scale. For survival/CIF curves,
#'  `"surv"` implies survival probabilities and `"risk"` implies CIF
#'  (1-survival in simple survival settings). Specify `"cumhaz"` to plot cumulative hazard
#'  or `"cloglog"` to generate a complementary log-log plot.
#' If `NULL`, a default is chosen from `outcome.type` or the survfit object.
#' @param label.x Character x-axis label (default `"Time"`).
#' @param label.y Character y-axis label (default is chosen automatically from `outcome.type`
#'   and `type.y`, e.g. "Survival", "Cumulative incidence" or "Cumulative hazard").
#' @param limits.x Numeric length-2 vector specifying x-axis limits. If `NULL`, it is
#'   set from the fitted object (typically `c(0, max(time))`).
#' @param limits.y Numeric length-2 vector specifying y-axis limits. If `NULL`, it is
#'   set to `c(0, 1)` for probability-type outcomes.
#' @param breaks.x Numeric vector of x-axis breaks (default `NULL`).
#' @param breaks.y Numeric vector of y-axis breaks (default `NULL`).
#' @param use.coord.cartesian Logical; if `TRUE`, uses `ggplot2::coord_cartesian()` for zooming
#'   instead of changing the scale limits (default `FALSE`).
#' @param add.conf Logical; if `TRUE`, adds a CI ribbon
#'   (via `ggsurvfit::add_confidence_interval()`). Default `TRUE`.
#' @param add.censor.mark Logical; if `TRUE`, draws censoring marks on each curve
#'   (via `ggsurvfit::add_censor_mark()`). Default `TRUE`.
#' @param shape.censor.mark Integer point shape used for censoring marks (default `3`).
#' @param size.censor.mark Numeric point size used for censoring marks (default `2`).
#' @param add.competing.risk.mark Logical; if `TRUE`, draws time marks for the competing event
#'   (event 2). If no times are supplied via `competing.risk.time`, the function tries to
#'   extract them automatically from the data. Default `FALSE`.
#' @param competing.risk.time A **named list** of numeric vectors. Each name must correspond to a
#'   strata label, and its numeric vector gives the times at which the competing event occurred
#'   in that stratum. Typically left as `list()` and filled internally.
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default `16`).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default `2`).
#' @param add.intercurrent.event.mark Logical; if `TRUE`, overlays user-specified intercurrent-event
#'   times per stratum. Default `FALSE`.
#' @param intercurrent.event.time A **named list** of numeric vectors for intercurrent events
#'   (names must match strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks
#'   (default `1`).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks
#'   (default `2`).
#' @param add.quantile Logical; if `TRUE`, adds a quantile reference line (via
#'   `ggsurvfit::add_quantile()`). Default `FALSE`.
#' @param level.quantile Numeric quantile level to be shown (default `0.5` for the median).
#' @param rows.columns.panel Optional integer vector `c(nrow, ncol)` controlling
#'   the layout of the panel returned by the panel modes. If `NULL`, an automatic
#'   layout is determined from the number of subplots.
#' @param style Character choosing the base plot style: `"classic"`, `"bold"`,
#' `"framed"`, `"grid"`, `"gray"` or `"ggsurvfit"` (default `"classic"`).
#'   Abbreviations such as `"C"`, `"B"`, `"F"`, or `"G"` are also accepted.
#' @param palette Optional character vector specifying the color palette to use across strata.
#' @param linewidth Optional numeric specifying the line width of curve (default `0.8`).
#' @param linetype Optional logical using different line types of curve (default `FALSE`).
#' @param font.family Character specifying the font family: `"sans"`,  `"serif"`, or
#' `"mono"` (default `"sans"`).
#' @param font.size Integer specifying the base font size (default `12`).
#' @param print.panel Logical. When `TRUE`, panel displays created internally are
#'   printed automatically in interactive sessions; otherwise they are returned
#'   invisibly for further modification (default `FALSE`).
#' @param filename.ggsave Character; if non-`NULL`, save the plot to this file.
#' @param width.ggsave Numeric width passed to `ggplot2::ggsave()` (default `6`).
#' @param height.ggsave Numeric height passed to `ggplot2::ggsave()` (default `6`).
#' @param dpi.ggsave Numeric DPI passed to `ggplot2::ggsave()` (default `300`).
#'
NULL

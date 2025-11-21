#' @title Generate a survival/CIF curve with marks that represent
#' censoring, competing risks and intermediate events
#'
#' @description
#' This function generates a survival or CIF curve from a unified formula–data interface
#' or from an existing survfit object. When a formula is supplied, the LHS is typically
#' `Event()` or `survivai::Surv()`, and the RHS specifies an optional
#' stratification variable. In addition to the curves themselves,
#' [cifplot()] can add numbers-at-risk tables, tables of point estimates and
#' CIs, censoring marks, competing-risk marks, and
#' intercurrent-event marks.
#'
#' For more complex multi-panel displays, [cifplot()] can internally call
#' [cifpanel()] via several “panel modes” (per event, per variable, or
#' censoring-focused). The function returns an object whose `plot`
#' component is a regular ggplot object that can be further modified (compatible with `+` and `%+%`).
#'
#' @inheritParams cif-stat-arguments
#' @inheritParams cif-visual-arguments
#'
#' @param formula_or_fit Either a model formula or a survfit object. When a formula is
#'   supplied, the LHS must be `Event(time, status)` or `Surv(time, status)`.
#'   The RHS specifies an optional stratification variable.
#' @param code.events Optional numeric length-3 vector `c(event1, event2, censoring)`.
#'   When supplied, it overrides `code.event1`, `code.event2`, and `code.censoring`
#'   (primarily used when [cifpanel()] is called or when `panel.per.event = TRUE`).
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
#'   in the legend / risk table. Specify the levels of strata. Levels not listed are dropped.
#' @param legend.position Character specifying the legend position:
#'   `"top"`, `"right"`, `"bottom"`, `"left"`, or `"none"` (default `"top"`).
#' @param panel.per.event Logical. **Explicit panel mode.** If `TRUE` and
#'   `outcome.type == "competing-risk"`, [cifplot()] internally calls [cifpanel()]
#'   to display two event-specific CIFs side-by-side (event 1 and event 2) using
#'   reversed `code.events`. Ignored for non-competing-risk outcomes.
#' @param panel.censoring Logical. **Explicit panel mode.** If `TRUE` and
#'   `outcome.type == "survival"`, [cifplot()] internally calls [cifpanel()]
#'   to display KM curves for `(event, censor)` and `(censor, event)` so that
#'   censoring patterns can be inspected.
#' @param panel.per.variable Logical. **Explicit panel mode.** If `TRUE` and the RHS
#'   of the formula has multiple covariates (e.g. `~ a + b + c`), the function produces
#'   a panel where each variable in the RHS is used once as the stratification factor.
#' @param panel.mode Character specifying **Automatic panel mode.** If `"auto"` and none of
#'   `panel.per.variable`, `panel.per.event`, `panel.censoring` has been set to `TRUE`,
#'   the function chooses a suitable panel mode automatically:
#'   (i) if the formula RHS has 2+ variables, it behaves like `panel.per.variable = TRUE`;
#'   (ii) otherwise, if `outcome.type == "competing-risk"`, it behaves like
#'   `panel.per.event = TRUE`; (iii) otherwise, if `outcome.type == "survival"`, it
#'   behaves like `panel.censoring = TRUE`. If a panel mode is explicitly specified,
#'   `panel.mode` is ignored.
#' @param survfit.info,axis.info,visual.info,panel.info,style.info,print.info,ggsave.info
#'   Internal lists used for programmatic control. Not intended for direct user input.
#'
#' @details
#'
#' ### Typical use cases
#' -   Draw one survival/CIF curve set by exposure groups (e.g., treatment vs control).
#' -   Call `cifpanel()` with a simplified code to create a panel displaying plots of multiple stratified survival/CIF curves or CIF curves for each event type.
#' -   Add CIs and censor/competing-risk/intercurrent-event marks.
#' -   Add a risk table to display the number at risk or the estimated survival probabilities or CIFs and CIs at each point in time.
#'
#' ### Key arguments shared with cifcurve()
#' -   **Outcome type and estimator**
#'       -   `outcome.type = "survival"`: Kaplan-Meier estimator
#'       -   `outcome.type = "competing-risk"`: Aalen-Johansen estimator
#'
#' -   **Confidence intervals**
#'       -   `conf.int` sets the two-sided level (default 0.95)
#'       -   `conf.type` chooses the transformation (`"arcsine-square root"`, `"plain"`, `"log"`, `"log-log"`, `"logit"`, or `"none"`)
#'       -   `error` chooses the estimator for SE (`"greenwood"`, `"tsiatis"` or `"if"` for survival curves and `"delta"`, `"aalen"` or `"if"` for CIFs)
#'
#' ### Key arguments for cifplot()
#' -   **Data visualization**
#'       -   `add.conf` adds CIs on the ggplot2-based plot
#'       -   `add.competing.risk.mark` and `add.intercurrent.event.mark` adds symbols to describe competing risks or intercurrent events in addition to conventional censoring marks with `add.censor.mark`
#'       -   `add.risktable` adds numbers at risk
#'       -   `add.estimate.table` adds time-by-time estimates and CIs
#'       -   `add.quantile` adds a reference line at a chosen quantile level
#'
#' -   **Plot customization**
#'       -   `type.y` chooses y-axis (`"surv"` for survival and `"risk"` for 1-survival/CIF)
#'       -   `limits.x`, `limits.y`, `breaks.x`, `breaks.y`: numeric vectors for axis control
#'       -   `style` specifies the appearance of plot (`"classic"`, `"bold"`, `"framed"`, `"grid"`, `"gray"` or `"ggsurvfit"`)
#'       -   `palette` specifies color of each curve (e.g. `palette=c("blue1", "cyan3", "navy", "deepskyblue3"))`)
#'
#' -   **Panel display**
#'       -   `panel.per.variable` produces multiple survival/CIF curves per stratification variable specified in the formula
#'       -   `panel.per.event` produces CIF curves for each event type
#'       -   `panel.censoring` produces the Kaplan–Meier curves for (event, censor) and (censor, event) so that censoring patterns can be inspected
#'       -   `panel.mode` uses automatic panel mode
#'
#' When \code{panel.per.event = TRUE}, two panels are created with
#' \code{code.events = list(c(e1, e2, c), c(e2, e1, c))}, where
#' \code{code.events = c(e1, e2, c)} is the input coding for event1, event2, and censoring.
#' Common legend is collected by default (\code{legend.collect = TRUE}).
#'
#' Numeric stratification variables are normalized automatically. Columns with
#' fewer than nine distinct numeric values are coerced to factors; columns with
#' nine or more distinct numeric values are split at the median into
#' \dQuote{Below median} and \dQuote{Above median} strata.
#'
#' ### Advanced control not required for typical use
#'
#' The arguments below fine-tune internal estimation and figure appearance.
#' **Most users do not need to change these defaults.**
#'
#' #### Graphical layers
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `add.conf` | Add confidence interval ribbon. | `TRUE` |
#' | `add.risktable` | Add numbers-at-risk table below the plot. | `TRUE` |
#' | `add.estimate.table` | Add estimates and confidence intervals table. | `FALSE` |
#' | `symbol.risk.table` |  Symbol for strata in risk / estimate tables | `"square"` |
#' | `font.size.risk.table` |  Font size for texts in risk / estimate tables | `3` |
#' | `add.censor.mark` | Add censoring marks. | `TRUE` |
#' | `add.competing.risk.mark` | Add marks for event2 of "competing-risk" outcome. | `FALSE` |
#' | `add.intercurrent.event.mark` | Add intercurrent event marks at user-specified times. | `FALSE` |
#' | `add.quantile` | Add quantile reference lines. | `FALSE` |
#' | `level.quantile` | Quantile level for `add.quantile`. | `0.5` |
#'
#' #### Time for marks
#'
#' | Argument | Description|
#' |---|---|
#' | `competing.risk.time` | **Named list** of numeric vectors that contains times of competing risks. Names must match strata labels. Typically created internally|
#' | `intercurrent.event.time` | **Named list** of numeric vectors that contains times of intercurrent events. Names must match strata labels. Typically created by `extract_time_to_event()`. |
#'
#' #### Appearance of marks
#'
#' | Argument | Applies to | Default |
#' |---|---|---|
#' | `shape.censor.mark` | Censoring marks | `3` (cross) |
#' | `size.censor.mark` | Censoring marks | `2` |
#' | `shape.competing.risk.mark` | Competing-risk marks | `16` (filled circle) |
#' | `size.competing.risk.mark` | Competing-risk marks | `2` |
#' | `shape.intercurrent.event.mark` | Intercurrent marks | `1` (circle) |
#' | `size.intercurrent.event.mark` | Intercurrent marks | `2` |
#'
#' #### Panel display
#'
#' | Argument | Description |
#' |---|---|
#' | `panel.per.variable` | One panel per stratification variable |
#' | `panel.per.event` | For `"competing-risk"`, show CIFs of event 1 and event 2 |
#' | `panel.censoring` | For survival, show (event, censor) vs (censor, event) |
#' | `panel.mode` with 2+ stratification variables | Behave like `panel.per.variable` |
#' | `panel.mode` with outcome.type = "competing-risk" | Behave like `panel.per.event` |
#' | `panel.mode` with outcome.type = "survival" | Behave like `panel.censoring` |
#'
#' #### Axes and legend
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `limits.x`, `limits.y` | Axis limits (`c(min, max)`) | Auto |
#' | `breaks.x`, `breaks.y` | Tick breaks for x and y axes | Auto |
#' | `use.coord.cartesian` | For zooming use `coord_cartesian()` | `FALSE` |
#' | `legend.position` |`"top"`, `"right"`, `"bottom"`, `"left"`, `"none"` | `"top"` |
#'
#' #### Export
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `filename.ggsave` | If non-`NULL`, save the plot using `ggsave()` | `NULL` |
#' | `width.ggsave` | Size passed to `ggsave()` | `6`|
#' | `height.ggsave` | Size passed to `ggsave()` | `6`|
#' | `dpi.ggsave` | DPI passed to `ggsave()` | `300` |
#'
#' **Notes**
#' - For CIF displays, set `type.y = "risk"`. For survival scale, use `type.y = NULL` or `= "surv"`.
#' - Event coding can be controlled via `code.event1`, `code.event2`, `code.censoring`.
#'   For ADaM-style data, use `code.event1 = 0`, `code.censoring = 1`.
#' - Per-stratum time lists should have names identical to plotted strata labels.
#'
#' @return
#' A `"cifplot"` object (a list) with at least the following elements:
#'
#' - `plot`: a ggplot object containing the main plot
#' - `patchwork`: reserved for compatibility with panel displays
#'   (typically `NULL` for single-panel plots)
#' - `survfit.info`, `axis.info`, `visual.info`, `panel.info`,
#'   `style.info`, `inset.info`, `print.info`, `ggsave.info`:
#'   internal lists storing the fitted curves and display settings
#' - `version`: a character string giving the cifmodeling version used
#' - `call`: the original function call
#'
#' The object is returned invisibly. When a panel mode is active and
#' `print.panel = TRUE`, the panel is also printed in interactive sessions.
#'
#' @examples
#' data(diabetes.complications)
#' cifplot(Event(t,epsilon) ~ fruitq,
#'         data = diabetes.complications,
#'         outcome.type="competing-risk",
#'         add.risktable = FALSE,
#'         label.y='CIF of diabetic retinopathy',
#'         label.x='Years from registration')
#'
#' @importFrom ggsurvfit ggsurvfit add_confidence_interval add_risktable add_risktable_strata_symbol add_censor_mark add_quantile
#' @importFrom ggplot2 theme_classic theme_bw element_text  element_rect element_blank labs lims geom_point aes
#' ggsave guides scale_color_discrete scale_fill_discrete scale_color_manual
#' scale_fill_manual scale_linetype_manual scale_linetype_discrete scale_shape_discrete scale_shape_manual
#' @importFrom grDevices gray
#' @importFrom patchwork wrap_plots
#'
#' @name cifplot
#' @keywords internal
#' @section Lifecycle:
#' \lifecycle{stable}
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifcurve()] for KM/AJ estimators; [cifpanel()] for display of multiple CIFs; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
cifplot <- function(
    formula_or_fit,
    data                          = NULL,
    weights                       = NULL,
    subset.condition              = NULL,
    na.action                     = na.omit,
    outcome.type                  = c("competing-risk", "survival"),
    code.event1                   = 1,
    code.event2                   = 2,
    code.censoring                = 0,
    code.events                   = NULL,
    error                         = NULL,
    conf.type                     = "arcsine-square root",
    conf.int                      = 0.95,
    type.y                        = NULL,
    label.x                       = "Time",
    label.y                       = NULL,
    label.strata                  = NULL,
    level.strata                  = NULL,
    order.strata                  = NULL,
    limits.x                      = NULL,
    limits.y                      = NULL,
    breaks.x                      = NULL,
    breaks.y                      = NULL,
    use.coord.cartesian           = FALSE,
    add.conf                      = TRUE,
    add.risktable                 = TRUE,
    add.estimate.table            = FALSE,
    symbol.risk.table             = "square",
    font.size.risk.table          = 3,
    add.censor.mark               = TRUE,
    shape.censor.mark             = 3,
    size.censor.mark              = 2,
    add.competing.risk.mark       = FALSE,
    competing.risk.time           = list(),
    shape.competing.risk.mark     = 16,
    size.competing.risk.mark      = 2,
    add.intercurrent.event.mark   = FALSE,
    intercurrent.event.time       = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark  = 2,
    add.quantile                  = FALSE,
    level.quantile                = 0.5,
    panel.per.event               = FALSE,
    panel.censoring               = FALSE,
    panel.per.variable            = FALSE,
    panel.mode                    = "auto",
    rows.columns.panel            = NULL,
    style                         = "classic",
    palette                       = NULL,
    linewidth                     = 0.8,
    linetype                      = FALSE,
    font.family                   = "sans",
    font.size                     = 12,
    legend.position               = "top",
    print.panel                   = FALSE,
    filename.ggsave               = NULL,
    width.ggsave                  = 6,
    height.ggsave                 = 6,
    dpi.ggsave                    = 300,
    survfit.info                  = NULL,
    axis.info                     = NULL,
    visual.info                   = NULL,
    panel.info                    = NULL,
    style.info                    = NULL,
    inset.info                    = NULL,
    print.info                    = NULL,
    ggsave.info                   = NULL,
    ...
) {
  dots <- list(...)
  dots$print.panel <- dots$print.panel %||% print.panel

  print_panel <- isTRUE(panel.per.variable) || isTRUE(panel.per.event) || isTRUE(panel.censoring)
  if (print_panel) {
    font.size <- font.size/2
    if (!is.null(label.strata)) {
      .warn("panel_disables_labelstrata")
    }
    if (isTRUE(add.risktable) || isTRUE(add.estimate.table)) {
      .warn("panel_disables_tables")
    }
    label.strata     <- NULL
    legend.position  <- "none"
    add.risktable     <- FALSE
    add.estimate.table <- FALSE
  }

  infos <- cifplot_build_info(
    error                         = error,
    conf.type                     = conf.type,
    conf.int                      = conf.int,

    type.y                        = type.y,
    label.x                       = label.x,
    label.y                       = label.y,
    level.strata                  = level.strata,
    label.strata                  = label.strata,
    order.strata                  = order.strata,
    limits.x                      = limits.x,
    limits.y                      = limits.y,
    breaks.x                      = breaks.x,
    breaks.y                      = breaks.y,
    use.coord.cartesian           = use.coord.cartesian,

    add.conf         = add.conf,
    add.risktable                  = add.risktable,
    symbol.risk.table              = symbol.risk.table,
    add.estimate.table              = add.estimate.table,
    font.size.risk.table           = font.size.risk.table,
    add.censor.mark                 = add.censor.mark,
    shape.censor.mark             = shape.censor.mark,
    size.censor.mark              = size.censor.mark,
    add.competing.risk.mark          = add.competing.risk.mark,
    competing.risk.time           = competing.risk.time,
    shape.competing.risk.mark     = shape.competing.risk.mark,
    size.competing.risk.mark      = size.competing.risk.mark,
    add.intercurrent.event.mark      = add.intercurrent.event.mark,
    intercurrent.event.time       = intercurrent.event.time,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark,
    size.intercurrent.event.mark  = size.intercurrent.event.mark,
    add.quantile               = add.quantile,
    level.quantile                      = level.quantile,

    panel.per.event                = panel.per.event,
    panel.censoring                = panel.censoring,
    panel.per.variable                  = panel.per.variable,
    rows.columns.panel            = rows.columns.panel,

    style           = style,
    palette         = palette,
    linewidth       = linewidth,
    linetype        = linetype,
    font.family     = font.family,
    font.size       = font.size,
    legend.position = legend.position,

    filename.ggsave = filename.ggsave,
    width.ggsave    = width.ggsave,
    height.ggsave   = height.ggsave,
    dpi.ggsave      = dpi.ggsave,

    survfit.info    = survfit.info,
    axis.info       = axis.info,
    visual.info     = visual.info,
    panel.info      = panel.info,
    style.info      = style.info,
    ggsave.info     = ggsave.info
  )

  survfit.info <- infos$survfit.info
  axis.info    <- infos$axis.info
  visual.info  <- infos$visual.info
  panel.info   <- infos$panel.info
  style.info   <- infos$style.info
  ggsave.info  <- infos$ggsave.info
  inset.info   <- inset.info %||% list()
  print.info   <- print.info %||% list()

  style.info$font.family <- style.info$font.family %||% "sans"
  style.info$font.size   <- style.info$font.size   %||% 12

  level_input <- axis.info$level.strata
  order_input <- axis.info$order.strata

  norm <- normalize_strata_info(
    level.strata = axis.info$level.strata,
    label.strata = axis.info$label.strata,
    order.strata = axis.info$order.strata
  )
  axis.info$level.strata <- norm$level
  axis.info$order.strata <- norm$order_data
  axis.info$label.strata <- norm$label_map

  if (!is.null(axis.info$label.strata)) {
    stopifnot(!is.null(names(axis.info$label.strata)))
  }
  if (!is.null(order_input) && !is.null(level_input)) {
    if (!all(as.character(order_input) %in% as.character(level_input))) {
      warning("order.strata has unknown levels; ignoring order/label application.")
      axis.info$order.strata <- NULL
      axis.info$label.strata <- NULL
    }
  }

  .assert(!(isTRUE(panel.info$panel.per.variable) && isTRUE(panel.info$panel.per.event) && isTRUE(panel.info$panel.censoring)),
          "incompatible_flags",
          which = "panel.per.variable, panel.per.event and panel.censoring")

  if (!inherits(formula_or_fit, "survfit")) {
    outcome.type <- util_check_outcome_type(outcome.type, formula=formula_or_fit, data = data)
  }

  if (!is.null(code.events)) {
    ce <- plot_check_code_events(code.events)
    .assert(length(ce) == 3L, "code_events_len_vec")
    .assert(all(ce == as.integer(ce)), "code_events_integer")
    .assert(ce[1L] != ce[2L], "code_events_distinct")
    code.event1    <- ce[1L]
    code.event2    <- ce[2L]
    code.censoring <- ce[3L]
  }

  if (isTRUE(add.competing.risk.mark) && length(competing.risk.time) == 0) {
    visual.info$competing.risk.time <- extract_time_to_event(
      formula_or_fit, data=data, subset.condition=subset.condition, na.action=na.action, which.event="event2",
      code.event1=code.event1, code.event2=code.event2, code.censoring=code.censoring)
  }

  panel_mode <- plot_decide_panel_mode(
    formula_or_fit = formula_or_fit,
    data           = data,
    outcome.type   = outcome.type,
    panel.info     = panel.info,
    panel.mode     = panel.mode
  )

  panel.info$panel.per.variable     <- identical(panel_mode, "each_var")
  panel.info$panel.per.event   <- identical(panel_mode, "each_event")
  panel.info$panel.censoring   <- identical(panel_mode, "censoring")

  if (isTRUE(panel.info$panel.per.variable)) {
    style.info$legend.position <- "none"
    .assert(!inherits(formula_or_fit, "survfit"), "need_formula_for_panel.per.variable")
    .assert(inherits(formula_or_fit, "formula"),  "formula_must_be")
    .assert(!is.null(data),                       "need_data")

    out_pv <- do.call(
      plot_panel.per.variable,
      c(list(
        formula          = formula_or_fit,
        data             = data,
        weights          = weights,
        subset.condition = subset.condition,
        na.action        = na.action,
        outcome.type     = outcome.type,
        code.event1      = code.event1,
        code.event2      = code.event2,
        code.censoring   = code.censoring,
        survfit.info     = survfit.info,
        axis.info        = axis.info,
        visual.info      = visual.info,
        panel.info       = panel.info,
        style.info       = style.info,
        ggsave.info      = ggsave.info
      ), dots)
    )
    if (!is.null(out_pv)) return(out_pv)
  }

  if (isTRUE(panel.info$panel.per.event)) {
    style.info$legend.position <- "none"
    if (inherits(formula_or_fit, "survfit")) {
      warning("panel.per.event=TRUE requires a formula interface; falling back to single-plot.")
    } else {
      out_pe <- plot_panel.per.event(
        formula            = formula_or_fit,
        data               = data,
        axis.info          = axis.info,
        visual.info        = visual.info,
        panel.info         = panel.info,
        style.info         = style.info,
        ggsave.info        = ggsave.info,
        survfit.info       = survfit.info,
        dots               = dots,
        code.event1        = code.event1,
        code.event2        = code.event2,
        code.censoring     = code.censoring,
        outcome.type       = outcome.type,
        rows.columns.panel = rows.columns.panel,
        subset.condition   = subset.condition,
        na.action          = na.action
      )
      if (!is.null(out_pe)) return(out_pe)
    }
  }

  if (isTRUE(panel.info$panel.censoring)) {
    style.info$legend.position <- "none"
    if (inherits(formula_or_fit, "survfit")) {
      warning("panel.censoring=TRUE requires a formula interface; falling back to single-plot.")
    } else {
      out_pe <- plot_panel.censoring(
        formula            = formula_or_fit,
        data               = data,
        axis.info          = axis.info,
        visual.info        = visual.info,
        panel.info         = panel.info,
        style.info         = style.info,
        ggsave.info        = ggsave.info,
        survfit.info       = survfit.info,
        dots               = dots,
        code.event1        = code.event1,
        code.event2        = code.event2,
        code.censoring     = code.censoring,
        outcome.type       = outcome.type,
        rows.columns.panel = rows.columns.panel,
        subset.condition   = subset.condition,
        na.action          = na.action
      )
      if (!is.null(out_pe)) return(out_pe)
    }
  }

  dots_clean <- plot_make_dots_clean(dots)
  args_single <- c(
    list(
      formula_or_fit   = formula_or_fit,
      data             = data,
      weights          = weights,
      subset.condition = subset.condition,
      na.action        = na.action,
      outcome.type     = outcome.type,
      code.event1      = code.event1,
      code.event2      = code.event2,
      code.censoring   = code.censoring,
      code.events      = NULL,
      survfit.info     = survfit.info,
      axis.info        = axis.info,
      visual.info      = visual.info,
      panel.info       = panel.info,
      style.info       = style.info,
      ggsave.info      = ggsave.info
    ),
    dots_clean
  )

  out_plot <- do.call(cifplot_single, args_single)
  out_plot <- apply_strata_to_plots(
    list(out_plot),
    order_data = axis.info$order.strata,
    label_map  = axis.info$label.strata
  )[[1]]

  survfit.info$formula_or_fit <- survfit.info$formula_or_fit %||% formula_or_fit
  survfit.info$outcome.type   <- survfit.info$outcome.type   %||% outcome.type
  survfit.info$code.event1    <- survfit.info$code.event1    %||% code.event1
  survfit.info$code.event2    <- survfit.info$code.event2    %||% code.event2
  survfit.info$code.censoring <- survfit.info$code.censoring %||% code.censoring
  survfit.info$code.events    <- survfit.info$code.events    %||% code.events
  survfit.info$data.name      <- survfit.info$data.name      %||% deparse(substitute(data))

  print.info <- print.info %||% list()
  print.info$print.panel <- isTRUE(print.panel)
  print.info$engine      <- print.info$engine %||% "cifplot_single"

  ret <- list(
    plot         = out_plot,
    patchwork    = NULL,
    survfit.info = survfit.info,
    axis.info    = axis.info,
    visual.info  = visual.info,
    panel.info   = panel.info,
    style.info   = style.info,
    inset.info   = inset.info %||% list(),
    print.info   = print.info,
    ggsave.info  = ggsave.info,
    version      = utils::packageVersion("cifmodeling"),
    call         = match.call()
  )
  class(ret) <- c("cifplot", class(ret))

  if (interactive() && isTRUE(ret$print.info$print.panel)) {
    print(ret$plot)
    invisible(ret)
  } else {
    ret
  }
}

plot_panel.per.variable <- function(
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = NULL,
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    survfit.info = NULL,
    axis.info = NULL,
    visual.info = NULL,
    panel.info = NULL,
    style.info = NULL,
    ggsave.info = NULL,
    ...
) {
  dots <- list(...)

  if (!inherits(formula, "formula")) {
    warning("panel.per.variable=TRUE requires a formula; falling back to single-plot.")
    return(NULL)
  }
  if (is.null(data)) {
    warning("panel.per.variable=TRUE requires `data`; falling back to single-plot.")
    return(NULL)
  }

  tt   <- terms(formula)
  vars <- attr(tt, "term.labels")

  if (length(vars) <= 1L) {
    return(NULL)
  }

  if (is.null(outcome.type)) {
    outcome.type <- util_check_outcome_type(
      formula      = formula,
      data         = data,
      na.action    = na.action,
      auto_message = FALSE
    )
  }

  lhs_chr  <- deparse(formula[[2L]])
  formulas <- lapply(vars, function(v) as.formula(paste0(lhs_chr, " ~ ", v)))
  K        <- length(formulas)

  if (identical(outcome.type, "competing-risk")) {
    code.events <- rep(list(c(code.event1, code.event2, code.censoring)), K)
  } else {
    code.events <- rep(list(c(code.event1, code.censoring)), K)
  }

  panel.info <- panel.info %||% list()
  if (is.null(panel.info$rows.columns.panel)) {
    ncol <- ceiling(sqrt(K))
    nrow <- ceiling(K / ncol)
    panel.info$rows.columns.panel <- c(nrow, ncol)
  }
  panel.info$panel.per.variable   <- FALSE
  panel.info$panel.per.event <- FALSE
  panel.info$panel.censoring <- FALSE
  if (is.null(panel.info$title.plot)) {
    panel.info$title.plot <- vars
  }

  if (is.null(dots$legend.collect)) dots$legend.collect <- TRUE
  if (is.null(dots$print.panel))    dots$print.panel    <- FALSE

  res <- do.call(
    cifpanel,
    c(list(
      formulas         = formulas,
      data             = data,
      weights          = weights,
      subset.condition = subset.condition,
      na.action        = na.action,
      outcome.type     = outcome.type,
      code.events      = code.events,
      survfit.info     = survfit.info,
      axis.info        = axis.info,
      visual.info      = visual.info,
      panel.info       = panel.info,
      style.info       = style.info,
      ggsave.info      = ggsave.info
    ), dots)
  )

  if (is.list(res) && !is.null(res$patchwork)) {
    return(res)
  }
  res
}

plot_panel.per.event <- function(
    formula,
    data,
    axis.info,
    visual.info,
    panel.info,
    style.info,
    ggsave.info,
    survfit.info,
    dots,
    code.event1,
    code.event2,
    code.censoring,
    outcome.type,
    rows.columns.panel,
    subset.condition = NULL,
    na.action = na.omit
) {
  if (is.null(outcome.type) || outcome.type != "competing-risk") {
    warning("panel.per.event=TRUE is only for COMPETING-RISK; falling back to single-plot.")
    return(NULL)
  }

  ce_panel <- plot_check_code_events(c(code.event1, code.event2, code.censoring))

  if (!is.null(dots$title.plot)) dots$title.plot <- NULL

#  ylabs_vec <- axis.info$label.y
#  if (is.null(ylabs_vec) && !is.null(dots$label.y)) ylabs_vec <- dots$label.y
#  if (is.null(ylabs_vec)) {
#    ylabs_vec <- plot_default_event_y_labels()
#  } else {
#    if (length(ylabs_vec) == 1L) ylabs_vec <- rep(ylabs_vec, 2L)
#    if (length(ylabs_vec)  > 2L) ylabs_vec <- ylabs_vec[1:2]
#  }
  if (!is.null(dots$label.y)) dots$label.y <- NULL
  ylabs_vec <- c("Cumulative incidence of interest", "Cumulative incidence of competing risk")

  axis.info.panel <- panel_modify_list(axis.info, list(
    label.y      = ylabs_vec,
    label.strata = axis.info$label.strata,
    order.strata = axis.info$order.strata,
    level.strata = axis.info$level.strata
  ))

  visual.info.panel <- panel_modify_list(visual.info, list(
  ))

  panel.info.panel <- panel_modify_list(panel.info, list(
    rows.columns.panel = if (is.null(rows.columns.panel)) c(1L, 2L) else rows.columns.panel
  ))

  style_cur     <- style.info$style
  palette_cur   <- style.info$palette
  linewidth_cur <- style.info$linewidth
  linetype_cur  <- style.info$linetype
  ff_cur        <- style.info$font.family
  fs_cur        <- style.info$font.size
  lg_cur        <- style.info$legend.position

  ggsave.info.panel <- panel_modify_list(ggsave.info, list(
  ))

  panel_args  <- list(
    formula           = formula,
    data              = data,
    subset.condition  = subset.condition,
    na.action         = na.action,
    outcome.type      = "competing-risk",
    code.events       = list(ce_panel, c(ce_panel[2L], ce_panel[1L], ce_panel[3L])),
    axis.info         = axis.info.panel,
    visual.info       = visual.info.panel,
    panel.info        = panel.info.panel,
    style.info        = list(
      style           = style_cur,
      palette         = palette_cur,
      linewidth       = linewidth_cur,
      linetype        = linetype_cur,
      font.family     = ff_cur,
      font.size       = fs_cur,
      legend.position = lg_cur
    ),
    ggsave.info       = ggsave.info.panel,
    survfit.info      = survfit.info
  )

  if (is.null(dots$rows.columns.panel)) dots$rows.columns.panel <- panel.info.panel$rows.columns.panel
  if (is.null(dots$legend.collect))     dots$legend.collect     <- TRUE
  if (is.null(dots$print.panel))        dots$print.panel        <- FALSE
  dots$visual.info <- NULL

  panel_out <- do.call(cifpanel, c(panel_args, dots))

  if (is.list(panel_out) && !is.null(panel_out$patchwork)) {
    return(panel_out)
  }
  panel_out
}

plot_panel.censoring <- function(
    formula,
    data,
    axis.info,
    visual.info,
    panel.info,
    style.info,
    ggsave.info,
    survfit.info,
    dots,
    code.event1,
    code.event2,
    code.censoring,
    outcome.type,
    rows.columns.panel,
    subset.condition = NULL,
    na.action = na.omit
) {
  if (is.null(outcome.type) || outcome.type != "survival") {
    warning("panel.censoring=TRUE is only for SURVIVAL outcome; falling back to single-plot.")
    return(NULL)
  }

  if (!is.null(dots$title.plot)) dots$title.plot <- NULL

  ylabs_vec <- c("Survival for event of interest", "Survival with censoring as event")

  axis.info.panel <- panel_modify_list(axis.info, list(
    label.y      = ylabs_vec,
    label.strata = axis.info$label.strata,
    order.strata = axis.info$order.strata,
    level.strata = axis.info$level.strata
  ))

  panel.info.panel <- panel_modify_list(panel.info, list(
    rows.columns.panel = if (is.null(rows.columns.panel)) c(1L, 2L) else rows.columns.panel
  ))

  panel_args  <- list(
    formula           = formula,
    data              = data,
    subset.condition  = subset.condition,
    na.action         = na.action,
    outcome.type      = "survival",
    code.events       = list(
      c(code.event1,    code.censoring),
      c(code.censoring, code.event1)
    ),
    axis.info         = axis.info.panel,
    visual.info       = visual.info,
    panel.info        = panel.info.panel,
    style.info        = style.info,
    ggsave.info       = ggsave.info,
    survfit.info      = survfit.info
  )

  if (is.null(dots$rows.columns.panel)) dots$rows.columns.panel <- panel.info.panel$rows.columns.panel
  if (is.null(dots$legend.collect))     dots$legend.collect     <- TRUE
  if (is.null(dots$print.panel))        dots$print.panel        <- FALSE
  dots$visual.info <- NULL

  panel_out <- do.call(cifpanel, c(panel_args, dots))

  if (is.list(panel_out) && !is.null(panel_out$patchwork)) {
    return(panel_out)
  }
  panel_out
}


cifplot_single <- function(
    formula_or_fit,
    data             = NULL,
    weights          = NULL,
    subset.condition = NULL,
    na.action        = na.omit,
    outcome.type     = c("competing-risk", "survival"),
    code.event1      = 1,
    code.event2      = 2,
    code.censoring   = 0,
    code.events      = NULL,
    survfit.info     = NULL,
    axis.info        = NULL,
    visual.info      = NULL,
    panel.info       = NULL,
    style.info       = NULL,
    ggsave.info      = NULL,
    ...
) {
  n_panel_flags <- sum(isTRUE(panel.info$panel.per.variable), isTRUE(panel.info$panel.per.event), isTRUE(panel.info$panel.censoring))
  .assert(n_panel_flags <= 1,
          "incompatible_flags",
          which = "panel.per.variable, panel.per.event, panel.censoring")

  dots         <- list(...)

  survfit.info <- survfit.info %||% list()
  axis.info    <- axis.info    %||% list()
  visual.info  <- visual.info  %||% list()
  panel.info   <- panel.info   %||% list()
  style.info   <- style.info   %||% list()
  ggsave.info  <- ggsave.info  %||% list()

  if (!is.list(survfit.info)) survfit.info <- list(value = survfit.info)
  if (!is.list(axis.info))    axis.info    <- list(value = axis.info)
  if (!is.list(visual.info))  visual.info  <- list(value = visual.info)
  if (!is.list(panel.info))   panel.info   <- list(value = panel.info)
  if (!is.list(style.info))   style.info   <- list(style = style.info)
  if (!is.list(ggsave.info))  ggsave.info  <- list(value = ggsave.info)

  if (!is.null(dots$style)) {
    style.info$style <- dots$style
    dots$style <- NULL
  }
  if (!is.null(dots$palette)) {
    style.info$palette <- dots$palette
    dots$palette <- NULL
  }
  if (!is.null(dots$linewidth)) {
    style.info$linewidth <- dots$linewidth
    dots$linewidth <- NULL
  }
  if (!is.null(dots$linetype)) {
    style.info$linetype <- dots$linetype
    dots$linetype <- NULL
  }
  if (!is.null(dots$font.family)) {
    style.info$font.family <- dots$font.family
    dots$font.family <- NULL
  }
  if (!is.null(dots$font.size)) {
    style.info$font.size <- dots$font.size
    dots$font.size <- NULL
  }
  if (!is.null(dots$legend.position)) {
    style.info$legend.position <- dots$legend.position
    dots$legend.position <- NULL
  }

  survfit.info <- panel_modify_list(list(
    error     = NULL,
    conf.type = "arcsine-square root",
    conf.int  = 0.95
  ), survfit.info)

  axis.info <- panel_modify_list(list(
    type.y            = NULL,
    label.x           = "Time",
    label.y           = NULL,
    label.strata      = NULL,
    level.strata      = NULL,
    order.strata      = NULL,
    limits.x          = NULL,
    limits.y          = NULL,
    breaks.x          = NULL,
    breaks.y          = NULL,
    use.coord.cartesian = FALSE
  ), axis.info)

  visual.info <- panel_modify_list(list(
    add.conf         = TRUE,
    add.risktable                  = FALSE,
    add.estimate.table              = FALSE,
    symbol.risk.table              = "square",
    font.size.risk.table           = 3,
    add.censor.mark                 = TRUE,
    shape.censor.mark             = 3,
    size.censor.mark              = 2,
    add.competing.risk.mark          = FALSE,
    competing.risk.time           = list(),
    shape.competing.risk.mark     = 16,
    size.competing.risk.mark      = 2,
    add.intercurrent.event.mark      = FALSE,
    intercurrent.event.time       = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark  = 2,
    add.quantile               = FALSE,
    level.quantile                      = 0.5
  ), visual.info)

  panel.info <- panel_modify_list(list(
    panel.per.event     = FALSE,
    panel.censoring     = FALSE,
    panel.per.variable       = FALSE,
    rows.columns.panel = NULL
  ), panel.info)

  style.info <- panel_modify_list(list(
    style           = "classic",
    palette         = NULL,
    linewidth       = 0.8,
    linetype        = FALSE,
    font.family     = "sans",
    font.size       = 12,
    legend.position = "top"
  ), style.info)

  ggsave.info <- panel_modify_list(list(
    filename.ggsave = NULL,
    width.ggsave    = 6,
    height.ggsave   = 6,
    dpi.ggsave      = 300,
    units           = "in"
  ), ggsave.info)

  error     <- survfit.info$error
  conf.type <- survfit.info$conf.type
  conf.int  <- survfit.info$conf.int

  type.y              <- axis.info$type.y
  label.x             <- axis.info$label.x
  label.y             <- axis.info$label.y
  label.strata        <- axis.info$label.strata
  level.strata        <- axis.info$level.strata
  order.strata        <- axis.info$order.strata
  limits.x            <- axis.info$limits.x
  limits.y            <- axis.info$limits.y
  breaks.x            <- axis.info$breaks.x
  breaks.y            <- axis.info$breaks.y
  use.coord.cartesian <- isTRUE(axis.info$use.coord.cartesian)

  add.conf        <- visual.info$add.conf
  add.risktable                 <- visual.info$add.risktable
  add.estimate.table             <- visual.info$add.estimate.table
  symbol.risk.table             <- visual.info$symbol.risk.table
  font.size.risk.table          <- visual.info$font.size.risk.table
  add.censor.mark                <- visual.info$add.censor.mark
  shape.censor.mark            <- visual.info$shape.censor.mark
  size.censor.mark             <- visual.info$size.censor.mark
  add.competing.risk.mark         <- visual.info$add.competing.risk.mark
  competing.risk.time          <- visual.info$competing.risk.time
  shape.competing.risk.mark    <- visual.info$shape.competing.risk.mark
  size.competing.risk.mark     <- visual.info$size.competing.risk.mark
  add.intercurrent.event.mark     <- visual.info$add.intercurrent.event.mark
  intercurrent.event.time      <- visual.info$intercurrent.event.time
  shape.intercurrent.event.mark<- visual.info$shape.intercurrent.event.mark
  size.intercurrent.event.mark <- visual.info$size.intercurrent.event.mark
  add.quantile              <- visual.info$add.quantile
  level.quantile                     <- visual.info$level.quantile

  panel.per.event     <- isTRUE(panel.info$panel.per.event)
  panel.censoring     <- isTRUE(panel.info$panel.censoring)
  panel.per.variable       <- isTRUE(panel.info$panel.per.variable)
  rows.columns.panel <- panel.info$rows.columns.panel

  style           <- style.info$style
  palette         <- style.info$palette
  linewidth       <- style.info$linewidth
  linetype        <- style.info$linetype
  font.family     <- style.info$font.family
  font.size       <- style.info$font.size
  legend.position <- style.info$legend.position

  filename.ggsave <- ggsave.info$filename.ggsave
  width.ggsave    <- ggsave.info$width.ggsave
  height.ggsave   <- ggsave.info$height.ggsave
  dpi.ggsave      <- ggsave.info$dpi.ggsave
  ggsave.units    <- ggsave.info$units %||% "in"

  level_input     <- axis.info$level.strata
  order_input     <- axis.info$order.strata

  norm <- normalize_strata_info(
    level.strata = axis.info$level.strata,
    label.strata = axis.info$label.strata,
    order.strata = axis.info$order.strata
  )

  axis.info$level.strata <- norm$level
  axis.info$order.strata <- norm$order_data
  axis.info$label.strata <- norm$label_map

  if (!is.null(axis.info$label.strata)) {
    stopifnot(!is.null(names(axis.info$label.strata)))
  }
  if (!is.null(order_input) && !is.null(level_input)) {
    if (!all(as.character(order_input) %in% as.character(level_input))) {
      warning("order.strata has unknown levels; ignoring order/label application.")
      axis.info$order.strata <- NULL
      axis.info$label.strata <- NULL
    }
  }
  level.strata <- axis.info$level.strata
  label.strata <- axis.info$label.strata
  order.strata <- axis.info$order.strata

  if (is.null(outcome.type)) {
    outcome.type <- util_check_outcome_type(formula = if (inherits(formula_or_fit,"survfit")) NULL
                                            else formula_or_fit, data = data, na.action = na.action, auto_message = FALSE)
  } else {
    outcome.type <- match.arg(outcome.type, c("competing-risk","survival"))
  }

  if (!inherits(formula_or_fit, "survfit")) {
    if (is.null(data)) stop("When `formula` is a formula, `data` must be provided.")
    norm_inputs <- plot_normalize_formula_data(formula_or_fit, data)
    data_working <- norm_inputs$data
    formula_or_fit <- cifcurve(formula_or_fit, data = data_working, weights = weights, subset.condition = subset.condition, na.action = na.action,
                               outcome.type = outcome.type, code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
                               error = error, conf.type = conf.type, conf.int = conf.int)
    formula_or_fit <- plot_survfit_short_strata_names(formula_or_fit)
  }

  p <- call_ggsurvfit(
    survfit_object = formula_or_fit,
    out_read_surv  = NULL,
    survfit.info   = survfit.info,
    axis.info      = axis.info,
    visual.info    = visual.info,
    panel.info     = panel.info,
    style.info     = style.info,
    ggsave.info    = ggsave.info
  )
  if (!is.null(filename.ggsave)) ggplot2::ggsave(filename.ggsave, plot = p, width = width.ggsave, height = height.ggsave, dpi = dpi.ggsave, units = ggsave.units)
  return(p)
}

#' Plot survival or cumulative incidence curves with ggsurvfit
#'
#' @param survfit_object A \code{survfit} object.
#' @param out_read_surv (optional) List returned by your \code{util_read_surv()} to auto-set x limits.
#' @param conf.type Character transformation for CI on the probability scale (default \code{"arcsine-square root"}).

#' @param add.conf Logical add \code{add_confidence_interval()} to plot. It calls geom_ribbon() (default \code{TRUE}).
#' @param add.risktable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param add.estimate.table Logical add \code{add_risktable(risktable_stats="estimate (conf.low, conf.high)")} to plot (default \code{FALSE}).
#' @param add.censor.mark Logical add \code{add_censor_mark()} to plot. It calls geom_point() (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param add.competing.risk.mark Logical add time marks to describe event2 specified by Event(), usually the competing events. It calls geom_point() (default \code{TRUE}).
#' @param competing.risk.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param add.intercurrent.event.mark Logical overlay user-specified time marks per strata calls geom_point() (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param add.quantile Logical add \code{add_quantile()} to plot. It calls geom_segment() (default \code{TRUE}).
#' @param level.quantile Numeric quantile level passed to \code{add_quantile()} (default \code{0.5}).

#' @param type.y \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param label.x Character x-axis labels (default \code{"Time"}).
#' @param label.y Character y-axis labels (default internally set to \code{"Survival"} or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata.

#' @param limits.x Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,max(out_read_surv$t))}.
#' @param limits.y Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,1)}.
#' @param breaks.x Numeric vectors for axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vectors for axis breaks (default \code{NULL}).
#' @param use.coord.cartesian Logical specify use of coord_cartesian() (default \code{FALSE}).

#' @param style Character plot theme controls (default \code{"classic"}).
#' @param font.family Character plot theme controls (default \code{"sans"}).
#' @param font.size Integer plot theme controls (default \code{14}).
#' @param legend.position Character specify position of legend: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).

#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
call_ggsurvfit <- function(
    survfit_object,
    out_read_surv = NULL,
    survfit.info  = NULL,
    axis.info     = NULL,
    visual.info   = NULL,
    panel.info    = NULL,
    style.info    = NULL,
    ggsave.info   = NULL
){
  survfit.info        <- survfit.info %||% list()
  axis.info           <- axis.info    %||% list()
  visual.info         <- visual.info  %||% list()
  panel.info          <- panel.info   %||% list()
  style.info          <- style.info   %||% list()
  ggsave.info         <- ggsave.info  %||% list()

  error               <- survfit.info$error
  conf.type           <- survfit.info$conf.type
  conf.int            <- survfit.info$conf.int

  type.y              <- axis.info$type.y
  label.x.user        <- axis.info$label.x
  label.y.user        <- axis.info$label.y
  label.strata        <- axis.info$label.strata
  level.strata        <- axis.info$level.strata
  order.strata        <- axis.info$order.strata
  limits.x            <- axis.info$limits.x
  limits.y            <- axis.info$limits.y
  breaks.x            <- axis.info$breaks.x
  breaks.y            <- axis.info$breaks.y
  use.coord.cartesian <- isTRUE(axis.info$use.coord.cartesian)

  add.conf                      <- visual.info$add.conf
  add.risktable                 <- visual.info$add.risktable
  add.estimate.table            <- visual.info$add.estimate.table
  symbol.risk.table             <- visual.info$symbol.risk.table
  font.size.risk.table          <- visual.info$font.size.risk.table
  add.censor.mark               <- visual.info$add.censor.mark
  shape.censor.mark             <- visual.info$shape.censor.mark
  size.censor.mark              <- visual.info$size.censor.mark
  add.competing.risk.mark       <- visual.info$add.competing.risk.mark
  competing.risk.time           <- visual.info$competing.risk.time
  shape.competing.risk.mark     <- visual.info$shape.competing.risk.mark
  size.competing.risk.mark      <- visual.info$size.competing.risk.mark
  add.intercurrent.event.mark   <- visual.info$add.intercurrent.event.mark
  intercurrent.event.time       <- visual.info$intercurrent.event.time
  shape.intercurrent.event.mark <- visual.info$shape.intercurrent.event.mark
  size.intercurrent.event.mark  <- visual.info$size.intercurrent.event.mark
  add.quantile                  <- visual.info$add.quantile
  level.quantile                <- visual.info$level.quantile

  panel.per.event               <- isTRUE(panel.info$panel.per.event)
  panel.censoring               <- isTRUE(panel.info$panel.censoring)
  panel.per.variable            <- isTRUE(panel.info$panel.per.variable)

  style              <- style.info$style
  palette            <- style.info$palette
  linewidth          <- style.info$linewidth
  linetype           <- style.info$linetype
  font.family        <- style.info$font.family
  font.size          <- style.info$font.size
  legend.position    <- style.info$legend.position

  filename.ggsave    <- ggsave.info$filename.ggsave
  width.ggsave       <- ggsave.info$width.ggsave
  height.ggsave      <- ggsave.info$height.ggsave
  dpi.ggsave         <- ggsave.info$dpi.ggsave
  ggsave.units       <- ggsave.info$units %||% "in"

  if (!identical(legend.position, "none")) {
    label.strata.map <- plot_make_label.strata.map(
      survfit_object = survfit_object,
      label.strata   = label.strata,
      level.strata   = level.strata
    )
    res <- plot_reconcile_order_and_labels(
      survfit_object   = survfit_object,
      label.strata.map = label.strata.map,
      level.strata     = level.strata,
      order.strata     = order.strata
    )
    limits_arg          <- res$limits_arg
    label.strata.map    <- res$label.strata.map
    strata_levels_final <- res$strata_levels_final
    strata_labels_final <- res$strata_labels_final
    n_strata_effective  <- length(limits_arg)
  } else {
    limits_arg          <- NULL
    label.strata.map    <- NULL
    strata_levels_final <- NULL
    strata_labels_final <- NULL
    n_strata_effective  <- NULL
  }

  out_cg <- check_ggsurvfit(
    survfit_object  = survfit_object,
    survfit.info    = survfit.info,
    axis.info       = axis.info,
    visual.info     = visual.info,
    style.info      = style.info,
    out_read_surv   = out_read_surv
  )

  p <- out_cg$out_survfit_object +
    ggplot2::labs(
      x = label.x.user %||% "Time",
      y = label.y.user %||% out_cg$label.y
    )

  if (isTRUE(add.conf)) {
    p <- p + add_confidence_interval()
  }
  if (isTRUE(add.censor.mark)) {
    p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)
  }
  if (isTRUE(add.quantile)) {
    p <- p + add_quantile(y_value = level.quantile)
  }

  apply_add_risktable_strata_symbol <- function (p, symbol.risk.table) {
    if (symbol.risk.table=="square") {
      p <- p + add_risktable_strata_symbol(symbol = "\U25A0", size = 14)
    } else if (symbol.risk.table=="circle") {
      p <- p + add_risktable_strata_symbol(symbol = "\U25CF", size = 14)
    } else if (symbol.risk.table=="triangle") {
      p <- p + add_risktable_strata_symbol(symbol = "\U25B2", size = 14)
    }
    p
  }

  if (isTRUE(add.estimate.table) && isTRUE(add.risktable)) {
    p <- p + add_risktable(
      risktable_stats = c("n.risk", "{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"),
      stats_label     = c("No. at risk", "Estimate (95% CI)"),
      theme           = plot_theme_risktable_font(font.family = font.family, plot.title.size = font.size),
      risktable_group = "risktable_stats",
      size            = font.size.risk.table
    )
    p <- apply_add_risktable_strata_symbol(p, symbol.risk.table)
  } else if (isTRUE(add.estimate.table)) {
    p <- p + add_risktable(
      risktable_stats = c("{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"),
      stats_label     = c("Estimate (95% CI)"),
      theme           = plot_theme_risktable_font(font.family = font.family, plot.title.size = font.size),
      size            = font.size.risk.table
    )
    p <- apply_add_risktable_strata_symbol(p, symbol.risk.table)
  } else if (isTRUE(add.risktable)) {
    p <- p + add_risktable(
      risktable_stats = c("n.risk"),
      stats_label     = c("No. at risk"),
      theme           = plot_theme_risktable_font(font.family = font.family, plot.title.size = font.size),
      size            = font.size.risk.table
    )
    p <- apply_add_risktable_strata_symbol(p, symbol.risk.table)
  }

  if (isTRUE(add.competing.risk.mark) && length(competing.risk.time)>0) {
    p <- plot_draw_marks(
      p, survfit_object,
      competing.risk.time, out_cg$type.y,
      shape = shape.competing.risk.mark,
      size  = size.competing.risk.mark
    )
  }
  if (isTRUE(add.intercurrent.event.mark) && length(intercurrent.event.time)>0) {
    p <- plot_draw_marks(
      p, survfit_object,
      intercurrent.event.time, out_cg$type.y,
      shape = shape.intercurrent.event.mark,
      size  = size.intercurrent.event.mark
    )
  }

  x_max <- plot_make_x_max(survfit_object)
  if (isTRUE(use.coord.cartesian)) {
    if (!is.null(breaks.x)) p <- p + ggplot2::scale_x_continuous(breaks = breaks.x)
    if (!is.null(breaks.y)) p <- p + ggplot2::scale_y_continuous(breaks = breaks.y)
    if (!is.null(limits.x) || !is.null(limits.y)) {
      p <- p + ggplot2::coord_cartesian(xlim = limits.x, ylim = limits.y, expand = FALSE)
    } else if (!is.null(limits.x)) {
      p <- p + ggplot2::coord_cartesian(xlim = limits.x, ylim = limits.y, expand = FALSE)
    } else if (!is.null(limits.y)) {
      p <- p + ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = limits.y, expand = FALSE)
    } else {
      p <- p + ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = c(0, 1), expand = FALSE)
    }
  } else {
    if (!is.null(breaks.y)) {
      p <- p + ggplot2::scale_y_continuous(breaks = breaks.y)
    } else if (!is.null(limits.y)) {
      p <- p + ggplot2::lims(y = limits.y)
    } else {
      p <- p + ggplot2::lims(y = c(0, 1))
    }
    if (!is.null(breaks.x)) {
      p <- p + ggplot2::scale_x_continuous(breaks = breaks.x)
    } else if (!is.null(limits.x)) {
      p <- p + ggplot2::lims(x = limits.x)
    } else {
      p <- p + ggplot2::lims(x = c(0, x_max))
    }
  }

  if (!identical(style, "ggsurvfit")) {
    p <- plot_apply_style(
      p,
      style               = style,
      font.family         = font.family,
      font.size           = font.size,
      legend.position     = legend.position,
      n_strata            = n_strata_effective,
      palette_colors      = palette,
      strata_levels_final = strata_levels_final,
      strata_labels_final = strata_labels_final
    )
  }

  p <- plot_apply_all_scales(
    p,
    style               = style,
    palette             = palette,
    n_strata            = n_strata_effective,
    strata_levels_final = strata_levels_final,
    strata_labels_final = strata_labels_final,
    limits_arg          = limits_arg
  )

  p <- p + ggplot2::guides(fill = "none", alpha = "none", shape = "none") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(fill = NA)))

  p <- plot_fix_palette_vector_arg(p)
  return(p)
}


check_ggsurvfit <- function(
    survfit_object,
    survfit.info  = NULL,
    axis.info     = NULL,
    visual.info   = NULL,
    style.info    = NULL,
    out_read_surv = NULL
){
  survfit.info <- survfit.info %||% list()
  axis.info    <- axis.info    %||% list()
  visual.info  <- visual.info  %||% list()
  style.info   <- style.info   %||% list()

  conf.type           <- survfit.info$conf.type
  type.y              <- axis.info$type.y
  label.y             <- axis.info$label.y
  limits.x            <- axis.info$limits.x
  limits.y            <- axis.info$limits.y
  breaks.x            <- axis.info$breaks.x
  breaks.y            <- axis.info$breaks.y
  use.coord.cartesian <- isTRUE(axis.info$use.coord.cartesian)

  add.conf                      <- visual.info$add.conf
  add.censor.mark               <- visual.info$add.censor.mark
  add.competing.risk.mark       <- visual.info$add.competing.risk.mark
  add.intercurrent.event.mark   <- visual.info$add.intercurrent.event.mark
  shape.censor.mark             <- visual.info$shape.censor.mark
  shape.competing.risk.mark     <- visual.info$shape.competing.risk.mark
  shape.intercurrent.event.mark <- visual.info$shape.intercurrent.event.mark

  style     <- style.info$style
  palette   <- style.info$palette
  linewidth <- style.info$linewidth
  linetype  <- style.info$linetype

  if (isTRUE(add.censor.mark) && isTRUE(add.intercurrent.event.mark) &&
      identical(shape.censor.mark, shape.intercurrent.event.mark)) {
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.intercurrent.event.mark")
  }
  if (isTRUE(add.censor.mark) && isTRUE(add.competing.risk.mark) &&
      identical(shape.censor.mark, shape.competing.risk.mark)) {
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.competing.risk.mark")
  }
  if (isTRUE(add.competing.risk.mark) && isTRUE(add.intercurrent.event.mark) &&
      identical(shape.competing.risk.mark, shape.intercurrent.event.mark)) {
    .warn("shape_identical", a = "shape.competing.risk.mark", b = "shape.intercurrent.event.mark")
  }

  is_len2_num <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  is_nondec   <- function(x) all(diff(x) >= 0, na.rm = TRUE)

  if (!is.null(limits.x)) {
    if (!(is.numeric(limits.x) && length(limits.x) == 2L && all(is.finite(limits.x)))) {
      .warn("limits_len2", arg = "limits.x")
    } else if (!all(diff(limits.x) > 0)) {
      .warn("limits_increasing", arg = "limits.x")
    }
    tmax <- suppressWarnings(max(survfit_object$time, na.rm = TRUE))
    if (is.finite(tmax)) {
      if (tmax < limits.x[1] || tmax > limits.x[2]) {
        .warn("limits_x_outside", tmax = signif(tmax, 6), arg = "limits.x",
              a = signif(limits.x[1], 6), b = signif(limits.x[2], 6))
      }
    }
  } else if (!is.null(out_read_surv) && !is.null(out_read_surv$t)) {
    tmax <- suppressWarnings(max(out_read_surv$t, na.rm = TRUE))
    if (!is.finite(tmax) || tmax <= 0) .warn("ors_tmax_bad")
  }

  if (!is.null(limits.y)) {
    if (!(is.numeric(limits.y) && length(limits.y) == 2L && all(is.finite(limits.y)))) {
      .warn("limits_len2", arg = "limits.y")
    } else if (!all(diff(limits.y) > 0)) {
      .warn("limits_increasing", arg = "limits.y")
    }

    surv  <- survfit_object$surv
    upper <- survfit_object$upper
    lower <- survfit_object$lower

    if (identical(type.y, "risk")) {
      surv  <- 1 - surv
      if (!is.null(upper)) upper <- 1 - upper
      if (!is.null(lower)) lower <- 1 - lower
    }

    if (any(surv < limits.y[1] | surv > limits.y[2], na.rm = TRUE))
      .warn("est_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])

    if (isTRUE(add.conf)) {
      if (!is.null(upper) && any(upper < limits.y[1] | upper > limits.y[2], na.rm = TRUE))
        .warn("upper_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])
      if (!is.null(lower) && any(lower < limits.y[1] | lower > limits.y[2], na.rm = TRUE))
        .warn("lower_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])
    }
  }

  check_breaks <- function(bk, nm, lims) {
    if (is.null(bk) || is.function(bk)) return(invisible())
    if (!is.numeric(bk)) {
      warning(sprintf("`%s` should be numeric (or a function).", nm), call. = FALSE); return(invisible())
    }
    if (!is_nondec(bk)) {
      warning(sprintf("`%s` must be non-decreasing.", nm), call. = FALSE)
    }
    if (!is.null(lims) && is_len2_num(lims)) {
      if (any(bk < lims[1] | bk > lims[2], na.rm = TRUE)) {
        warning(sprintf("Some `%s` are outside plotting range [%g, %g].", nm, lims[1], lims[2]), call. = FALSE)
      }
    }
    invisible()
  }
  check_breaks(breaks.x, "breaks.x", limits.x)
  check_breaks(breaks.y, "breaks.y", limits.y)

  if (is.null(label.y)) {
    auto_label <- plot_default_y_label(survfit_object$type, type.y)
    if (!is.null(auto_label)) label.y <- auto_label
  }

  coerce_conf <- function(survfit_object, conf.type) {
    if (!is.null(survfit_object$lower) && !is.null(survfit_object$upper)) return(survfit_object)
    if (conf.type %in% c("none", "n") || length(survfit_object$strata) > 2) {
      x <- survfit_object
      x$lower <- x$surv
      x$upper <- x$surv
      return(x)
    }
    survfit_object
  }
  survfit_object <- coerce_conf(survfit_object, conf.type)


  type.y <- util_check_type_y(type.y)
#  type.y <- plot_normalize_type_y(type.y)
  target_type <- switch(
    survfit_object$type,
    "kaplan-meier"   = if (identical(type.y, "risk")) "risk" else "surv",
    "aalen-johansen" = if (identical(type.y, "surv")) "surv" else "risk",
    if (identical(type.y, "risk")) "risk" else "surv"
  )
  type.y <- if (identical(target_type, "risk")) "risk" else "surv"

  out_plot <- ggsurvfit(
    survfit_object,
    type         = target_type,
    linewidth    = linewidth,
    linetype_aes = linetype
  )

  list(
    out_survfit_object = out_plot,
    label.y            = label.y,
    type.y             = type.y
  )
}

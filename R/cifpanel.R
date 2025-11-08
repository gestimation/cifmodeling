#' @title Arrange multiple survival / CIF plots in a panel display
#'
#' @description
#' \code{cifpanel()} is the panel-building counterpart of \code{cifplot()}.
#' It takes one or more model formulas (or, alternatively, one formula and several
#' event-coding specifications) and returns a multi-panel figure, typically as a
#' \pkg{patchwork} object. Most display options (axis labels, marks, style, ggsave options)
#' are shared with \code{cifplot()}, but per-panel legends and risk tables are
#' suppressed to avoid duplicated display.
#'
#' Panel layout is specified by length-2 vector \code{rows.columns.panel}.
#' This function can also automatically determine the panel count in the following order:
#' (1) if \code{plots} is supplied, its length defines the number of plots,
#' (2) else if \code{formulas} is supplied, its length defines the number of plots,
#' (3) else if \code{code.events} is supplied, its length defines the number of plots
#' together with formula, and (4) otherwise \code{rows.columns.panel=c(1,1)}.
#'
#' -   `formula` or `formulas` — one formula or a list of formulas; each entry creates a panel.
#' -   `data`, `outcome.type`, `code.events`, `type.y` — recycled across panels unless a list is supplied for per-panel control.
#' -   `rows.columns.panel` — selects grid layout by c(rows, cols).
#' -   `inset.panel` — selects inset layout.
#' -   `title.panel`, `subtitle.panel`, `caption.panel`, `title.plot` — overall titles and captions.
#' -   `tag.panel` — panel tag style (e.g., "A", "a", "1").
#' -   `label.x`, `label.y`, `limits.x`, `limits.y`, `breaks.x`, `breaks.y` — shared axis control unless a list is supplied for per-panel control.
#'
#' @inheritParams cif-stat-arguments
#' @inheritParams cif-visual-arguments
#'
#' @param plots Optional list of already-built \code{ggplot} objects to be arranged.
#'   If supplied, these are used as-is (no fitting is done).
#' @param formula A model formula specifying the time-to-event outcome on the
#'   left-hand side (typically \code{Event(time, status)} or \code{Surv(time, status)})
#'   and, optionally, a stratification variable on the right-hand side.
#'   Unlike \code{\link{cifplot}}, this function does not accept a fitted
#'   \code{survfit} object.
#' @param formulas Optional list of formulas. When given, each formula defines
#'   **one panel**. This is the most common way to create “one variable per plot”
#'   panels.
#' @param code.events Optional numeric length-3 vector \code{c(event1, event2, censoring)}.
#'   When supplied, it overrides \code{code.event1}, \code{code.event2}, and \code{code.censoring}
#'   (primarily used when \code{cifpanel()} is called or when \code{printEachEvent = TRUE}).
#' @param legend.collect Logical; if \code{TRUE}, try to collect a single legend
#'   for all panels (passed to \pkg{patchwork}). Default \code{TRUE}.
#' @param inset.panel Logical. If \code{FALSE} (default), all panels are arranged
#'   in a regular grid using \code{patchwork::wrap_plots()} and \code{plot_layout()}.
#'   If \code{TRUE}, the function switches to “inset mode”: the **first** plot becomes
#'   the main plot and the **second** plot (only the second) is drawn on top of it
#'   as an inset. Additional plots beyond the second are ignored in inset mode.
#'   Use grid mode to display more than two panels (inset.panel = FALSE).
#' @param inset.left,inset.bottom,inset.right,inset.top Numeric values in the range
#'   \code{[0, 1]} that define the inset box as fractions of the reference area.
#'   \code{inset.left} / \code{inset.right} control the horizontal position,
#'   \code{inset.bottom} / \code{inset.top} control the vertical position.
#'   Values are interpreted as “from the left/bottom” of the reference.
#'   For example, \code{inset.left = 0.4}, \code{inset.right = 1.0} draws the inset
#'   over the right 60% of the reference area.
#' @param inset.align.to Character string specifying the coordinate system for the
#'   inset box. One of:
#'   \itemize{
#'     \item \code{"panel"} (default): the box is placed relative to the **panel area**
#'       (i.e. the plotting region, excluding outer titles/margins);
#'     \item \code{"plot"}: the box is placed relative to the **entire plot** area,
#'       including axes and titles of the main plot;
#'     \item \code{"full"}: the box is placed relative to the **full patchwork canvas**.
#'   }
#'   This argument is passed to \code{patchwork::inset_element()}.
#' @param inset.legend.position Optional legend position **for the inset plot only**.
#'   If \code{NULL} (default), the inset plot keeps whatever legend position was
#'   defined for it (often this means a legend will also be inset).
#'   Set, for example, \code{"none"} to hide the legend inside the inset,
#'   while still showing the main plot's legend.
#' @param title.panel,subtitle.panel,caption.panel Character annotations applied to the
#'   **whole** panel layout (not to individual plots). These are passed to
#'   \code{patchwork::plot_annotation()} and are useful for creating figure-like
#'   outputs (title + subfigures + caption).
#' @param tag.panel Passed to \code{patchwork::plot_annotation()} to auto-label
#'   individual panels (e.g. \code{"A"}, \code{"B"}, \code{"C"}). Typical values are
#'   \code{"A"}, \code{"1"}, or \code{"a"}. See \code{?patchwork::plot_annotation}.
#' @param title.plot Character vector of titles for **each panel** in the order they
#'   are drawn. Length-1 values are recycled to all panels. In inset mode, the first
#'   element refers to the main plot and the second (if present) to the inset.
#' @param print.panel Logical. If \code{TRUE}, the composed patchwork object is
#'   printed immediately (for interactive use). If \code{FALSE}, the object is
#'   returned invisibly so that it can be assigned, modified, or saved. Kept for
#'   backward compatibility.
#'
#' @param ... Additional arguments forwarded to the internal \code{cifplot_single()}
#'   calls that build each panel. Use this to pass low-level options such as
#'   \code{competing.risk.time}, \code{intercurrent.event.time}, or styling overrides.
#'
#' @details
#'
#' ### Overview
#' `cifpanel()` composes multiple survival/CIF plots into a single figure.
#' For each panel, it estimates curves via `cifcurve()` and renders them with
#' `cifplot()`. You can supply a single `formula` reused across panels or a
#' list in `formulas` (one per panel). When both are provided, `formulas` wins.
#'
#' ### Outcome type & event coding
#'
#' - Use `outcome.type` to set per-panel estimator (`"SURVIVAL"`=KM, `"COMPETING-RISK"`=AJ).
#' - Alternatively, pass `code.events` per panel to infer the type:
#'   - length 2 = SURVIVAL: `c(event1, censor)`
#'   - length 3 = COMPETING-RISK: `c(event1, event2, censor)`
#' - If `outcome.type` is `NULL`, the function infers each panel from its
#'   `code.events[[i]]` length. When both are given, `outcome.type` takes precedence.
#'
#' ### Panel-wise vs shared arguments
#'
#' Many arguments accept a **scalar** (recycled to all panels) or a **list/vector**
#' (one entry per panel). Precedence: **panel-wise explicit values** >
#' **shared scalar** > **internal defaults**. Length-1 inputs are recycled.
#'
#' ### Grid vs inset composition
#'
#' - **Grid mode** (`inset.panel = FALSE`, default): plots are arranged with
#'   `patchwork::wrap_plots()` and `plot_layout()`. If `legend.collect = TRUE`,
#'   legends are collected across panels where possible.
#' - **Inset mode** (`inset.panel = TRUE`): the **second** plot is overlaid
#'   into the **first** using `patchwork::inset_element()`. Only the first two
#'   plots are used; extra plots are ignored. Control the inset box with
#'   `inset.left`, `inset.bottom`, `inset.right`, `inset.top`, and its
#'   reference frame via `inset.align.to` (`"panel"`, `"plot"`, or `"full"`).
#'
#' ### Advanced panel controls (forwarded to `cifplot()`)
#'
#' The following arguments allow **per-panel** control by supplying vectors/lists,
#' or **shared** control by supplying scalars. They are forwarded to `cifplot()`.
#'
#' #### Scale & labels
#'
#' | Argument | Meaning | Default |
#' |---|---|---|
#' | `type.y` | `"risk"` (CIF y-axis) or `NULL` (survival). | inferred |
#' | `label.x`, `label.y` | Axis labels per panel. | auto |
#' | `label.strata` | Legend labels per panel. | from data |
#' | `limits.x`, `limits.y` | Axis limits `c(min, max)`. | auto |
#' | `breaks.x`, `breaks.y` | Axis breaks (forwarded to `breaks.x`/`breaks.y`). | auto |
#'
#' #### Plot layers (toggles)
#'
#' | Argument | Effect | Default |
#' |---|---|---|
#' | `addConfidenceInterval` | CI ribbon. | `TRUE` |
#' | `addCensorMark` | Censor marks. | `TRUE` |
#' | `addCompetingRiskMark` | Marks for event2 at supplied times. | `FALSE` |
#' | `addIntercurrentEventMark` | User-specified intercurrent marks. | `FALSE` |
#' | `addQuantileLine` | Quantile line(s). | `FALSE` |
#'
#' *(Time marks inputs such as `competing.risk.time` / `intercurrent.event.time`
#' can be given via `...` if needed; names must match strata labels.)*
#'
#' ### Legend & annotations
#'
#' - `legend.position`: `"top"`, `"right"`, `"bottom"`, `"left"`, or `"none"` (applies to all panels).
#' - Grid mode: `legend.collect = TRUE` attempts a shared legend.
#' - Panel annotations: `title.panel`, `subtitle.panel`, `caption.panel`.
#' - Tagging: `tag.panel` is passed to `patchwork::plot_annotation()`.
#' - In inset mode, `title.plot = c(title_base, title_inset)` labels the two plots.
#'
#' ### Export (optional)
#'
#' If `filename.ggsave` is non-`NULL`, the composed panel is saved with
#' `ggsave()` using `width.ggsave`, `height.ggsave`, and `dpi.ggsave`.
#' Otherwise, the function returns objects without saving.
#'
#' ### Value
#'
#' Returns **invisibly**:
#' `list(plots = <list of ggplot objects>, out_patchwork = <patchwork object>)`.
#' Print the latter to display the composed panel. If `print.panel = TRUE`,
#' printing is done automatically.
#'
#' ### Notes & tips
#'
#' - Mixed panel types are supported (e.g., AJ in panel 1; KM in panel 2).
#' - If `formulas` is shorter than the grid capacity, empty slots are ignored.
#' - When supplying vectors/lists per panel, their lengths must match the number
#'   of panels; length-1 inputs are recycled; otherwise an error is thrown.
#' - For CIF displays, set `type.y = "risk"` in the relevant panels.
#' - ADaM-style coding can be expressed via `code.events` (e.g., `c(0,1)` for KM:
#'   `event1=0`, `censor=1`).
#' - Additional graphical options (e.g., theme) can be added post-hoc to each
#'   element of `plots` or to the composed `out_patchwork`.

#' @importFrom patchwork wrap_plots plot_layout inset_element plot_annotation
#' @return An invisible list: \code{list(plots = <list of ggplot objects>, out_patchwork = <patchwork object>)}.
#' Print the returned object to display the panel, or access individual panels via
#' `out_patchwork$plots[[1]]`, `out_patchwork$plots[[2]]`, ...
#'
#' @keywords internal
#' @param survfit.info,axis.info,visual.info,panel.info,style.info,print.info,ggsave.info,inset.info
#'   Internal lists used for programmatic control. Not intended for direct user input.
#'
#' @examples
#' data(diabetes.complications)
#' cifpanel(
#'   title.panel = "A comparison of cumulative incidence of competing events",
#'   rows.columns.panel = c(1,2),
#'   formula = Event(t, epsilon) ~ fruitq,
#'   data = diabetes.complications,
#'   outcome.type = "COMPETING-RISK",
#'   code.events = list(c(1,2,0), c(2,1,0)),
#'   label.y = c("Diabetic retinopathy", "Macrovascular complications"),
#'   label.x = "Years from registration",
#'   subtitle.panel = "Stratified by fruit intake",
#'   caption.panel  = "Data: diabetes.complications",
#'   title.plot = c("Diabetic retinopathy", "Macrovascular complications"),
#'   legend.position = "bottom",
#'   legend.collect=TRUE
#' )
#'
#' cifpanel(
#'   title.plot = c("Associations between fruit intake and macrovascular complications", "Details"),
#'   inset.panel = TRUE,
#'   formula = Event(t, epsilon) ~ fruitq,
#'   data = diabetes.complications,
#'   outcome.type = "COMPETING-RISK",
#'   code.events = list(c(2,1,0), c(2,1,0)),
#'   label.y = c("CIF of macrovascular complications", ""),
#'   label.x = c("Years from registration", ""),
#'   limits.y     = list(c(0,1), c(0,0.15)),
#'   inset.left   = 0.40, inset.bottom = 0.45,
#'   inset.right  = 1.00, inset.top    = 0.95,
#'   inset.align.to = "plot",
#'   inset.legend.position = "none",
#'   legend.position = "bottom",
#'   addConfidenceInterval = FALSE
#' )
#'
#' output1 <- cifplot(Event(t,epsilon) ~ fruitq,
#'                    data = diabetes.complications,
#'                    outcome.type="COMPETING-RISK",
#'                    code.event1=2,
#'                    code.event2=1,
#'                    addConfidenceInterval = FALSE,
#'                    addRiskTable = FALSE,
#'                    label.y='CIF of macrovascular complications',
#'                    label.x='Years from registration')
#' output2 <- cifplot(Event(t,epsilon) ~ fruitq,
#'                    data = diabetes.complications,
#'                    outcome.type="COMPETING-RISK",
#'                    code.event1=2,
#'                    code.event2=1,
#'                    addConfidenceInterval = FALSE,
#'                    addRiskTable = FALSE,
#'                    label.y='CIF of macrovascular complications',
#'                    label.x='Years from registration',
#'                    limits.y=c(0,0.15))
#' output3 <- list(a=output1, b=output2)
#' cifpanel(plots = output3,
#'          inset.panel = TRUE,
#'          inset.left = 0.40, inset.bottom = 0.45,
#'          inset.right = 1.00, inset.top = 0.95,
#'          inset.align.to = "plot",
#'          inset.legend.position = "none",
#'          legend.position = "bottom")
#'
#' @importFrom ggplot2 ggplot theme_void ggsave theme element_text labs
#' @importFrom patchwork wrap_plots plot_layout inset_element plot_annotation

#' @name cifpanel
#' @section Lifecycle:
#' \lifecycle{experimental}
#'
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifcurve()] for KM/AJ estimators; [cifplot()] for display of a CIF; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
cifpanel <- function(
    plots                         = NULL,
    formula                       = NULL,
    formulas                      = NULL,
    data                          = NULL,
    weights                       = NULL,
    subset.condition              = NULL,
    na.action                     = na.omit,
    outcome.type                  = NULL,
    code.events                   = NULL,
    error                         = NULL,
    conf.type                     = NULL,
    conf.int                      = NULL,
    type.y                        = NULL,
    label.x                       = NULL,
    label.y                       = NULL,
    label.strata                  = NULL,
    order.strata                  = NULL,
    level.strata                  = NULL,
    limits.x                      = NULL,
    limits.y                      = NULL,
    breaks.x                      = NULL,
    breaks.y                      = NULL,
    addConfidenceInterval         = NULL,
    addRiskTable                  = NULL,
    addEstimateTable              = NULL,
    symbol.risktable              = NULL,
    font.size.risktable           = NULL,
    addCensorMark                 = NULL,
    shape.censor.mark             = NULL,
    size.censor.mark              = NULL,
    addCompetingRiskMark          = NULL,
    competing.risk.time           = NULL,
    shape.competing.risk.mark     = NULL,
    size.competing.risk.mark      = NULL,
    addIntercurrentEventMark      = NULL,
    intercurrent.event.time       = NULL,
    shape.intercurrent.event.mark = NULL,
    size.intercurrent.event.mark  = NULL,
    addQuantileLine               = NULL,
    quantile                      = NULL,
    rows.columns.panel            = c(1, 1),
    inset.panel                   = FALSE,
    title.panel                   = NULL,
    subtitle.panel                = NULL,
    caption.panel                 = NULL,
    tag.panel                     = NULL,
    title.plot                    = NULL,
    style                         = "CLASSIC",
    palette                       = NULL,
    font.family                   = "sans",
    font.size                     = 8,
    legend.position               = "top",
    legend.collect                = TRUE,
    inset.left                    = 0.60,
    inset.bottom                  = 0.05,
    inset.right                   = 0.98,
    inset.top                     = 0.45,
    inset.align.to                = c("panel","plot","full"),
    inset.legend.position         = NULL,
    print.panel                   = TRUE,
    filename.ggsave               = NULL,
    width.ggsave                  = NULL,
    height.ggsave                 = NULL,
    dpi.ggsave                    = 300,
    survfit.info                  = NULL,
    axis.info                     = NULL,
    visual.info                   = NULL,
    panel.info                    = NULL,
    style.info                    = NULL,
    inset.info                    = NULL,
    print.info                    = NULL,
    ggsave.info                   = NULL,
    engine                        = "cifplot",
    ...
){
#  if (!is.null(label.strata)) {
#    .warn("panel_disables_labelstrata")
#  }
#  if (isTRUE(addRiskTable) || isTRUE(addEstimateTable)) {
#    .warn("panel_disables_tables")
#  }
  legend.position  <- "none"
  addRiskTable     <- FALSE
  addEstimateTable <- FALSE
  inset.align.to <- match.arg(inset.align.to)

  dots <- list(...)

  survfit.info.user <- survfit.info
  axis.info.user    <- axis.info
  visual.info.user  <- visual.info
  panel.info.user   <- panel.info
  style.info.user   <- style.info
  inset.info.user   <- inset.info
  print.info.user   <- print.info
  ggsave.info.user  <- ggsave.info

  survfit.info <- modifyList(list(
    error     = error,
    conf.type = conf.type,
    conf.int  = conf.int
  ), survfit.info %||% list())

  axis.info <- modifyList(list(
    type.y            = type.y,
    label.x           = label.x,
    label.y           = label.y,
    level.strata      = level.strata,
    order.strata      = order.strata,
    label.strata      = label.strata,
    limits.x          = limits.x,
    limits.y          = limits.y,
    breaks.x          = breaks.x,
    breaks.y          = breaks.y,
    use_coord_cartesian = get0("use_coord_cartesian", ifnotfound = NULL)
  ), axis.info %||% list())

  visual.info <- modifyList(list(
    addConfidenceInterval    = addConfidenceInterval,
    ci.alpha                 = 0.25,
    addRiskTable             = FALSE,
    addEstimateTable         = FALSE,
    addCensorMark            = addCensorMark,
    shape.censor.mark        = 3,
    size.censor.mark         = 2,
    addCompetingRiskMark     = addCompetingRiskMark,
    competing.risk.time      = list(),
    shape.competing.risk.mark= 16,
    size.competing.risk.mark = 2,
    addIntercurrentEventMark = addIntercurrentEventMark,
    intercurrent.event.time  = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark  = 2,
    addQuantileLine          = addQuantileLine,
    quantile                 = 0.5,
    line.size                = 0.9,
    symbol.risktable         = NULL,
    font.size.risktable      = NULL
  ), visual.info %||% list())

  panel.info <- modifyList(list(
    printEachEvent     = FALSE,
    printEachVar       = FALSE,
    rows.columns.panel = rows.columns.panel,
    title.panel        = title.panel,
    subtitle.panel     = subtitle.panel,
    caption.panel      = caption.panel,
    tag.panel          = tag.panel,
    title.plot         = title.plot
  ), panel.info %||% list())

  style.info <- style.info %||% list()
  style.info$style           <- style.info$style           %||% style
  style.info$palette         <- style.info$palette         %||% palette
  style.info$font.family     <- style.info$font.family     %||% font.family
  style.info$font.size       <- style.info$font.size       %||% font.size
  style.info$legend.position <- style.info$legend.position %||% legend.position
  style.info$legend.collect  <- style.info$legend.collect  %||% legend.collect

  style.info <- modifyList(list(
    style           = "CLASSIC",
    palette         = NULL,
    font.family     = "sans",
    font.size       = 12,
    legend.position = "top",
    legend.collect  = FALSE
  ), style.info)

  inset.info <- modifyList(list(
    inset.panel           = inset.panel,
    inset.align.to        = inset.align.to,
    inset.left            = inset.left,
    inset.bottom          = inset.bottom,
    inset.right           = inset.right,
    inset.top             = inset.top,
    inset.legend.position = inset.legend.position
  ), inset.info %||% list())

  print.info <- modifyList(list(
    print.panel = print.panel
  ), print.info %||% list())

  ggsave.info <- modifyList(list(
    filename.ggsave = filename.ggsave,
    width.ggsave    = width.ggsave,
    height.ggsave   = height.ggsave,
    dpi.ggsave      = dpi.ggsave,
    units           = "in"
  ), ggsave.info %||% list())

  inset.info$inset.align.to <- match.arg(inset.info$inset.align.to, c("panel","plot","full"))

  rows.columns.panel <- panel.info$rows.columns.panel
  title.panel        <- panel.info$title.panel
  subtitle.panel     <- panel.info$subtitle.panel
  caption.panel      <- panel.info$caption.panel
  tag.panel          <- panel.info$tag.panel
  title.plot         <- panel.info$title.plot

  legend.position    <- style.info$legend.position
  legend.collect     <- isTRUE(style.info$legend.collect)

  inset.panel        <- isTRUE(inset.info$inset.panel)
  inset.align.to     <- inset.info$inset.align.to
  inset.left         <- inset.info$inset.left
  inset.bottom       <- inset.info$inset.bottom
  inset.right        <- inset.info$inset.right
  inset.top          <- inset.info$inset.top
  inset.legend.position <- inset.info$inset.legend.position

  print.panel        <- isTRUE(print.info$print.panel)

  filename.ggsave    <- ggsave.info$filename.ggsave
  width.ggsave       <- ggsave.info$width.ggsave
  height.ggsave      <- ggsave.info$height.ggsave
  dpi.ggsave         <- ggsave.info$dpi.ggsave
  ggsave.units       <- ggsave.info$units %||% "in"

  fonts <- panel_extract_fonts(c(style.info, dots))
  style.info$font.family <- fonts$family
  style.info$font.size   <- fonts$size
  theme.panel.unified    <- panel_build_theme(font.family = fonts$family, font.size = fonts$size)

  nrow <- as.integer(rows.columns.panel[1]); ncol <- as.integer(rows.columns.panel[2])
  n_slots <- nrow * ncol

  level_input <- axis.info$level.strata
  order_input <- axis.info$order.strata

  norm <- normalize_strata_info(
    level.strata = axis.info$level.strata,
    order.strata = axis.info$order.strata,
    label.strata = axis.info$label.strata
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

  type.y              <- axis.info$type.y
  label.x             <- axis.info$label.x
  label.y             <- axis.info$label.y
  label.strata        <- axis.info$label.strata
  order.strata        <- axis.info$order.strata
  limits.x            <- axis.info$limits.x
  limits.y            <- axis.info$limits.y
  breaks.x            <- axis.info$breaks.x
  breaks.y            <- axis.info$breaks.y
  use_coord_cartesian <- axis.info$use_coord_cartesian

  addConfidenceInterval    <- visual.info$addConfidenceInterval
  addCensorMark            <- visual.info$addCensorMark
  addCompetingRiskMark     <- visual.info$addCompetingRiskMark
  addIntercurrentEventMark <- visual.info$addIntercurrentEventMark
  addQuantileLine          <- visual.info$addQuantileLine

  # ------------------------------------------------------------
  # 1) plots が指定されているときは「並べるだけ」モード
  # ------------------------------------------------------------
  if (!is.null(plots)) {
    if (!is.list(plots)) {
      stop("`plots` must be a list of ggplot objects.")
    }
    if (length(plots) && !all(vapply(plots, function(p) inherits(p, "ggplot"), logical(1)))) {
      stop("All elements of `plots` must inherit from 'ggplot'.")
    }

    # ここは ggsurvfit の warning を避けるために touch_colour = FALSE にしてもよい
    plots <- apply_strata_to_plots(
      plots,
      order_data   = axis.info$order.strata,
      label_map    = axis.info$label.strata,
      touch_colour = TRUE
    )

    plots_out <- plots
    if (isTRUE(inset.panel)) {
      if (length(plots) < 2L) .err("inset_need_two")
      if (length(plots) > 2L) .warn("inset_extra_drop")
      p_base  <- plots[[1]] + ggplot2::theme(legend.position = legend.position)
      if (!is.null(title.plot) && length(title.plot) >= 1L) {
        p_base <- p_base + ggplot2::labs(title = title.plot[[1]])
      }
      p_inset <- plots[[2]] + ggplot2::theme(legend.position = inset.legend.position)
      if (!is.null(title.plot) && length(title.plot) >= 2L) {
        p_inset <- p_inset + ggplot2::labs(title = title.plot[[2]])
      }
      out_patchwork <- p_base +
        patchwork::inset_element(
          p_inset,
          left = inset.left,
          bottom = inset.bottom,
          right = inset.right,
          top = inset.top,
          align_to = inset.align.to
        )
    } else {
      plots2 <- plots
      if (length(plots2) < n_slots) {
        plots2 <- c(
          plots2,
          rep(list(ggplot2::ggplot() + ggplot2::theme_void()), n_slots - length(plots2))
        )
      } else if (length(plots2) > n_slots) {
        .warn("plots_extra_dropped", n_plots = length(plots2), n_slots = n_slots)
        plots2 <- plots2[seq_len(n_slots)]
      }
      out_patchwork <- patchwork::wrap_plots(plots2, nrow = nrow, ncol = ncol)
      if (isTRUE(legend.collect)) {
        out_patchwork <- out_patchwork +
          patchwork::plot_layout(guides = "collect") &
          ggplot2::theme(legend.position = legend.position)
      } else {
        out_patchwork <- out_patchwork &
          ggplot2::theme(legend.position = legend.position)
      }
      plots_out <- plots2
    }

    out_patchwork <- out_patchwork + patchwork::plot_annotation(
      title      = title.panel,
      subtitle   = subtitle.panel,
      caption    = caption.panel,
      tag_levels = tag.panel,
      theme      = theme.panel.unified
    )

    if (isTRUE(print.panel)) print(out_patchwork)
    if (!is.null(filename.ggsave)) {
      if (is.null(width.ggsave))  width.ggsave  <- if (isTRUE(inset.panel)) 6 else max(6, 5 * rows.columns.panel[2])
      if (is.null(height.ggsave)) height.ggsave <- if (isTRUE(inset.panel)) 6 else max(6, 5 * rows.columns.panel[1])
      ggplot2::ggsave(filename.ggsave, plot = out_patchwork,
                      width = width.ggsave, height = height.ggsave,
                      dpi = dpi.ggsave, units = ggsave.units)
    }

    return(invisible(list(
      plots        = plots_out,
      out_patchwork= out_patchwork,
      axis.info    = axis.info,
      survfit.info = survfit.info,
      visual.info  = visual.info,
      panel.info   = panel.info,
      style.info   = style.info,
      inset.info   = inset.info,
      print.info   = print.info,
      ggsave.info  = ggsave.info
    )))
  }

  rows.columns.panel <- panel.info$rows.columns.panel
  nrow   <- as.integer(rows.columns.panel[1])
  ncol   <- as.integer(rows.columns.panel[2])
  n_slots <- nrow * ncol

  # ------------------------------------------------------------
  # 2) ここから「推定して描く」通常モード
  # ------------------------------------------------------------
  if (is.null(data)) stop("data must be provided.")
  if (is.null(code.events) || !is.list(code.events) || length(code.events) == 0)
    .err("need_code_events")
  if (!is.null(formulas) && !is.null(formula))
    .warn("both_formula_forms")
  if (is.null(formulas) && is.null(formula))
    .err("need_formula_or_formulas")

  layout_info <- panel_update_rows.columns.panel(
    formulas   = formulas,
    code.events = code.events,
    panel.info  = panel.info,
    n_slots     = n_slots
  )
  K                  <- layout_info$K
  code.events        <- layout_info$code.events
  panel.info$rows.columns.panel <- layout_info$rows.columns.panel
  nrow               <- layout_info$nrow
  ncol               <- layout_info$ncol
  n_slots            <- layout_info$n_slots

  use_formula_list <- !is.null(formulas)
  if (use_formula_list) {
    stopifnot(is.list(formulas))
    formulas <- lapply(formulas, panel_as_formula_global)
    formulas <- panel_recycle_to(formulas, K)
  } else {
    formula  <- panel_as_formula_global(formula)
    formulas <- rep(list(formula), K)
  }

  toL <- panel_to_list; rec <- panel_recycle_to

  outcome.list <- toL(outcome.type);      if (!is.null(outcome.list)) outcome.list <- rec(outcome.list, K)
  typey.list   <- toL(type.y);            if (!is.null(typey.list))   typey.list   <- rec(typey.list, K)
  labely.list  <- toL(label.y);           if (!is.null(labely.list))  labely.list  <- rec(labely.list, K)
  labelx.list  <- toL(label.x);           if (!is.null(labelx.list))  labelx.list  <- rec(labelx.list, K)

  make_panel_list_preserve_vector <- function(x, K) {
    if (is.list(x)) return(panel_recycle_to(x, K))
    rep(list(x), K)
  }
  labelstrata.list  <- make_panel_list_preserve_vector(label.strata,  K)
  orderstrata.list  <- make_panel_list_preserve_vector(order.strata,  K)

  # limits/breaksをパネルごとにしておく
  limsx.list <- NULL
  if (!is.null(limits.x)) {
    limsx.list <- if (is.list(limits.x)) limits.x else list(limits.x)
    limsx.list <- rec(limsx.list, K)
  }
  limsy.list <- NULL
  if (!is.null(limits.y)) {
    limsy.list <- if (is.list(limits.y)) limits.y else list(limits.y)
    limsy.list <- rec(limsy.list, K)
  }

  breakx.list <- toL(breaks.x); if (!is.null(breakx.list)) breakx.list <- rec(breakx.list, K)
  breaky.list <- toL(breaks.y); if (!is.null(breaky.list)) breaky.list <- rec(breaky.list, K)

  addCI.list   <- toL(addConfidenceInterval);    if (!is.null(addCI.list))   addCI.list   <- rec(addCI.list, K)
  addCen.list  <- toL(addCensorMark);            if (!is.null(addCen.list))  addCen.list  <- rec(addCen.list, K)
  addCR.list   <- toL(addCompetingRiskMark);     if (!is.null(addCR.list))   addCR.list   <- rec(addCR.list, K)
  addIC.list   <- toL(addIntercurrentEventMark); if (!is.null(addIC.list))   addIC.list   <- rec(addIC.list, K)
  addQ.list    <- toL(addQuantileLine);          if (!is.null(addQ.list))    addQ.list    <- rec(addQ.list, K)

  if (!is.null(addCI.list))   visual.info$addConfidenceInterval    <- NULL
  if (!is.null(addCen.list))  visual.info$addCensorMark            <- NULL
  if (!is.null(addCR.list))   visual.info$addCompetingRiskMark     <- NULL
  if (!is.null(addIC.list))   visual.info$addIntercurrentEventMark <- NULL
  if (!is.null(addQ.list))    visual.info$addQuantileLine          <- NULL

  # outcome.flag 判定
  infer_flag_by_codes <- function(v) if (length(v) == 2L) "S" else if (length(v) == 3L) "C" else NA_character_
  if (!is.null(outcome.list)) {
    outcome.flags <- vapply(outcome.list, panel_norm_outcome, character(1))
  } else {
    outcome.flags <- vapply(code.events, infer_flag_by_codes, character(1))
    if (anyNA(outcome.flags)) .err("infer_outcome_fail")
  }
  panel_validate_code_events(code.events, outcome.flags)

  # dots の中にある「パネルで決めた値と衝突するやつ」を抜く
  kill_names <- c()
  if (!is.null(outcome.list))     kill_names <- c(kill_names, "outcome.type")
  if (!is.null(typey.list))       kill_names <- c(kill_names, "type.y")
  if (!is.null(labely.list))      kill_names <- c(kill_names, "label.y")
  if (!is.null(limsy.list))       kill_names <- c(kill_names, "limits.y")
  if (!is.null(labelx.list))      kill_names <- c(kill_names, "label.x")
  if (!is.null(limsx.list))       kill_names <- c(kill_names, "limits.x")
  if (!is.null(labelstrata.list)) kill_names <- c(kill_names, "label.strata")
  if (!is.null(orderstrata.list)) kill_names <- c(kill_names, "order.strata")
  if (!is.null(breakx.list))      kill_names <- c(kill_names, "breaks.x","breaks.x")
  if (!is.null(breaky.list))      kill_names <- c(kill_names, "breaks.y","breaks.y")
  if (!is.null(addCI.list))       kill_names <- c(kill_names, "addConfidenceInterval")
  if (!is.null(addCen.list))      kill_names <- c(kill_names, "addCensorMark")
  if (!is.null(addCR.list))       kill_names <- c(kill_names, "addCompetingRiskMark")
  if (!is.null(addIC.list))       kill_names <- c(kill_names, "addIntercurrentEventMark")
  if (!is.null(addQ.list))        kill_names <- c(kill_names, "addQuantileLine")

  dots <- panel_strip_overrides_from_dots(dots, unique(kill_names))

  # engine をパネル数にそろえる（★ここが今回の肝）
  engine.list <- panel_to_list(engine)
  engine.list <- panel_recycle_to(engine.list, K)

  # パネル内部は panel_prepare でまとめて生成させる（いままで通り）
  prep <- panel_prepare(
    K               = K,
    formulas        = formulas,
    data            = data,
    code.events     = code.events,
    outcome.flags   = outcome.flags,
    outcome.list    = outcome.list,
    typey.list      = typey.list,
    labely.list     = labely.list,
    typex.list      = typex.list,
    labelx.list     = labelx.list,
    limsx.list      = limsx.list,
    limsy.list      = limsy.list,
    breakx.list     = breakx.list,
    breaky.list     = breaky.list,
    addCI.list      = addCI.list,
    addCen.list     = addCen.list,
    addCR.list      = addCR.list,
    addIC.list      = addIC.list,
    addQ.list       = addQ.list,
    strata.list     = make_panel_list_preserve_vector(label.strata, K),
    legend.position = legend.position,
    survfit.info    = survfit.info,
    style.info      = style.info,
    dots            = dots,
    fonts           = fonts
  )

  plots <- lapply(seq_len(prep$K), function(i) {
    pa <- prep$plot_args[[i]]

    # ① 親を入れる（順番は親→子）
    pa$axis.info    <- modifyList(axis.info,    pa$axis.info %||% list())
    pa$visual.info  <- modifyList(visual.info,  pa$visual.info %||% list())
    pa$style.info   <- modifyList(style.info,   pa$style.info %||% list())
    pa$survfit.info <- modifyList(survfit.info, pa$survfit.info %||% list())

    # ② ★ここでパネルのものを全部たたき込む（これが最優先）
    pa <- panel_force_apply(
      pa,
      i,
      labelx.list = labelx.list,
      labely.list = labely.list,
      limsx.list  = limsx.list,
      limsy.list  = limsy.list,
      breakx.list = breakx.list,
      breaky.list = breaky.list,
      addCI.list  = addCI.list,
      addCen.list = addCen.list,
      addCR.list  = addCR.list,
      addIC.list  = addIC.list,
      addQ.list   = addQ.list
    )

    # --- ここで competing.risk.time を自動生成 ---
    if (isTRUE(pa$addCompetingRiskMark)) {
      ce <- code.events[[i]]
      has_event2 <- !is.null(ce) && length(ce) >= 3L && !is.na(ce[2])
      has_time   <- !is.null(pa$competing.risk.time) && length(pa$competing.risk.time) > 0

      if (has_event2 && !has_time) {
        pa$visual.info$competing.risk.time <- extract_time_to_event(
          formula          = formulas[[i]],
          data             = data,
          subset.condition = subset.condition,
          na.action        = na.action,
          which_event      = "event2",
          code.event1      = ce[1],
          code.event2      = ce[2],
          code.censoring   = ce[3]
        )
      }
    }
    eng_i <- engine.list[[i]]
    p_i <- if (identical(eng_i, "ggsurvfit")) {
      sf_i <- prep$curves[[i]]

      axis_i <- list(
        type.y            = pa$type.y,
        label.x           = pa$label.x,
        label.y           = pa$label.y,
        label.strata      = pa$label.strata,
        level.strata      = pa$level.strata,
        order.strata      = pa$order.strata,
        limits.x          = pa$limits.x,
        limits.y          = pa$limits.y,
        breaks.x          = pa$breaks.x,
        breaks.y          = pa$breaks.y,
        use_coord_cartesian = pa$use_coord_cartesian
      )

      visual_i <- modifyList(visual.info, list(
        addConfidenceInterval    = pa$addConfidenceInterval,
        addRiskTable             = pa$addRiskTable,
        addEstimateTable         = pa$addEstimateTable,
        addCensorMark            = pa$addCensorMark,
        shape.censor.mark        = pa$shape.censor.mark,
        size.censor.mark         = pa$size.censor.mark,
        addCompetingRiskMark     = pa$addCompetingRiskMark,
        competing.risk.time      = pa$competing.risk.time,
        shape.competing.risk.mark= pa$shape.competing.risk.mark,
        size.competing.risk.mark = pa$size.competing.risk.mark,
        addIntercurrentEventMark = pa$addIntercurrentEventMark,
        intercurrent.event.time  = pa$intercurrent.event.time,
        shape.intercurrent.event.mark = pa$shape.intercurrent.event.mark,
        size.intercurrent.event.mark  = pa$size.intercurrent.event.mark,
        addQuantileLine          = pa$addQuantileLine,
        quantile                 = pa$quantile
      ))

      panel_i <- list(
        printEachEvent     = FALSE,
        printEachVar       = FALSE,
        rows.columns.panel = NULL
      )

      call_ggsurvfit(
        survfit_object   = sf_i,
        survfit.info     = survfit.info,
        axis.info        = axis_i,
        visual.info      = visual_i,
        panel.info       = panel_i,
        style.info       = style.info,
        ggsave.info      = ggsave.info
      )
    } else {
      allowed <- setdiff(names(formals(cifplot_single)), "...")
      if (!is.null(names(pa))) {
        pa <- pa[intersect(names(pa), allowed)]
      }
      if (!"formula_or_fit" %in% names(pa)) {
        pa <- c(list(formula_or_fit = prep$curves[[i]]), pa)
      }
      do.call(cifplot_single, pa)
    }
    if (!is.null(title.plot)) {
      if (is.list(title.plot)) {
        title_k <- title.plot[[ min(i, length(title.plot)) ]]
      } else {
        title_k <- title.plot[ min(i, length(title.plot)) ]
      }
      if (!is.null(title_k) && !identical(title_k, "")) {
        p_i <- p_i + ggplot2::labs(title = title_k)
      }
    }
    p_i
  })

  has_ggsurvfit <- any(vapply(engine.list, identical, logical(1), y = "ggsurvfit"))

  plots <- apply_strata_to_plots(
    plots,
    order_data   = axis.info$order.strata,
    label_map    = axis.info$label.strata,
    touch_colour = !has_ggsurvfit
  )

  if (!is.null(title.plot)) {
    for (i in seq_along(plots)) {
      if (is.list(title.plot)) {
        title_k <- title.plot[[ min(i, length(title.plot)) ]]
      } else {
        title_k <- title.plot[ min(i, length(title.plot)) ]
      }
      if (!is.null(title_k) && !identical(title_k, "")) {
        plots[[i]] <- plots[[i]] + ggplot2::labs(title = title_k)
      }
    }
  }

  if (isTRUE(inset.panel)) {
    if (length(plots) < 2L) .err("inset_need_two")
    if (length(plots) > 2L) .warn("inset_extra_drop")
    p_base  <- plots[[1]] + ggplot2::theme(legend.position = legend.position)
    p_inset <- plots[[2]] + ggplot2::theme(legend.position = inset.legend.position)
    if (!is.null(title.plot)) {
      if (length(title.plot) >= 1L) p_base  <- p_base  + ggplot2::labs(title = title.plot[[1]])
      if (length(title.plot) >= 2L) p_inset <- p_inset + ggplot2::labs(title = title.plot[[2]])
    }
    out_patchwork <- p_base +
      patchwork::inset_element(
        p_inset, left = inset.left, bottom = inset.bottom,
        right = inset.right, top = inset.top, align_to = inset.align.to
      )
  } else {
    if (length(plots) < n_slots) {
      plots <- c(plots, rep(list(ggplot2::ggplot() + ggplot2::theme_void()), n_slots - length(plots)))
    } else if (length(plots) > n_slots) {
      .warn("plots_extra_dropped", n_plots = length(plots), n_slots = n_slots)
      plots <- plots[seq_len(n_slots)]
    }
    out_patchwork <- patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol)
    if (isTRUE(legend.collect)) {
      out_patchwork <- out_patchwork +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = legend.position)
    } else {
      out_patchwork <- out_patchwork &
        ggplot2::theme(legend.position = legend.position)
    }
  }

  out_patchwork <- out_patchwork + patchwork::plot_annotation(
    title      = title.panel,
    subtitle   = subtitle.panel,
    caption    = caption.panel,
    tag_levels = tag.panel,
    theme      = theme.panel.unified
  )

  if (isTRUE(print.panel)) print(out_patchwork)
  if (!is.null(filename.ggsave)) {
    if (is.null(width.ggsave))  width.ggsave  <- if (isTRUE(inset.panel)) 6 else max(6, 5 * rows.columns.panel[2])
    if (is.null(height.ggsave)) height.ggsave <- if (isTRUE(inset.panel)) 6 else max(6, 5 * rows.columns.panel[1])
    ggplot2::ggsave(filename.ggsave, plot = out_patchwork,
                    width = width.ggsave, height = height.ggsave,
                    dpi = dpi.ggsave, units = ggsave.units)
  }

  invisible(list(
    plots        = plots,
    out_patchwork= out_patchwork,
    axis.info    = axis.info,
    survfit.info = survfit.info,
    visual.info  = visual.info,
    panel.info   = panel.info,
    style.info   = style.info,
    inset.info   = inset.info,
    print.info   = print.info,
    ggsave.info  = ggsave.info
  ))
}

panel_force_apply <- function(
    pa,
    i,
    labelx.list = NULL,
    labely.list = NULL,
    limsx.list  = NULL,
    limsy.list  = NULL,
    breakx.list = NULL,
    breaky.list = NULL,
    addCI.list  = NULL,
    addCen.list = NULL,
    addCR.list  = NULL,
    addIC.list  = NULL,
    addQ.list   = NULL
) {
  if (!is.null(labelx.list)) {
    pa$label.x <- labelx.list[[i]]
    pa$axis.info$label.x <- labelx.list[[i]]
  }
  if (!is.null(labely.list)) {
    pa$label.y <- labely.list[[i]]
    pa$axis.info$label.y <- labely.list[[i]]
  }

  if (!is.null(limsx.list)) {
    pa$limits.x <- limsx.list[[i]]
    pa$axis.info$limits.x <- limsx.list[[i]]
  }
  if (!is.null(limsy.list)) {
    pa$limits.y <- limsy.list[[i]]
    pa$axis.info$limits.y <- limsy.list[[i]]
  }
  if (!is.null(breakx.list)) {
    pa$breaks.x <- breakx.list[[i]]
    pa$axis.info$breaks.x <- breakx.list[[i]]
  }
  if (!is.null(breaky.list)) {
    pa$breaks.y <- breaky.list[[i]]
    pa$axis.info$breaks.y <- breaky.list[[i]]
  }

  if (!is.null(addCI.list)) {
    v <- isTRUE(addCI.list[[i]])
    pa$addConfidenceInterval <- v
    pa$visual.info$addConfidenceInterval <- v
  }
  if (!is.null(addCen.list)) {
    v <- isTRUE(addCen.list[[i]])
    pa$addCensorMark <- v
    pa$visual.info$addCensorMark <- v
  }
  if (!is.null(addCR.list)) {
    v <- isTRUE(addCR.list[[i]])
    pa$addCompetingRiskMark <- v
    pa$visual.info$addCompetingRiskMark <- v
  }
  if (!is.null(addIC.list)) {
    v <- isTRUE(addIC.list[[i]])
    pa$addIntercurrentEventMark <- v
    pa$visual.info$addIntercurrentEventMark <- v
  }
  if (!is.null(addQ.list)) {
    v <- isTRUE(addQ.list[[i]])
    pa$addQuantileLine <- v
    pa$visual.info$addQuantileLine <- v
  }

  pa
}




normalize_strata_info <- function(level.strata = NULL,
                                  order.strata = NULL,
                                  label.strata = NULL) {
  ch <- function(x) if (is.null(x)) NULL else as.character(x)

  level <- ch(level.strata)
  if (is.null(level) || !length(level)) {
    return(list(level = NULL, order_data = NULL, label_map = NULL))
  }

  if (!is.null(label.strata)) {
    if (is.null(names(label.strata)) || !any(nzchar(names(label.strata)))) {
      if (length(label.strata) == length(level)) {
        label_map <- stats::setNames(as.character(label.strata), level)
      } else {
        label_map <- stats::setNames(level, level)
      }
    } else {
      lab_names <- names(label.strata)
      lab_vals  <- as.character(unname(label.strata))

      keep <- intersect(lab_names, level)

      label_map <- stats::setNames(label.strata[keep], keep)

      miss <- setdiff(level, names(label_map))
      if (length(miss)) {
        label_map <- c(label_map, stats::setNames(miss, miss))
      }

      label_map <- label_map[level]
    }
  } else {
    label_map <- stats::setNames(level, level)
  }

  ord <- ch(order.strata)
  if (is.null(ord)) ord <- level

  if (!all(ord %in% level)) {
    return(list(level = level, order_data = NULL, label_map = NULL))
  }

  list(level = level, order_data = ord, label_map = label_map)
}

apply_strata_to_plots <- function(plots, order_data, label_map, touch_colour = TRUE) {
  if (!length(plots) || is.null(order_data) || is.null(label_map)) return(plots)
  brks <- order_data
  labs <- unname(label_map[brks])

  relabel_one <- function(p) {
    layers <- list(
      ggplot2::scale_linetype_discrete(breaks = brks, labels = labs),
      ggplot2::scale_shape_discrete(breaks = brks, labels = labs),
      ggplot2::scale_fill_discrete(breaks = brks, labels = labs)
    )
    if (isTRUE(touch_colour)) {
      layers <- c(layers, list(ggplot2::scale_color_discrete(breaks = brks, labels = labs)))
    }
    suppressMessages(Reduce(`+`, layers, init = p))
  }

  lapply(plots, relabel_one)
}

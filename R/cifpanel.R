#' @title Generates a multi-panel figure for survival or cumulative incidence curves,
#' arranged either in a grid layout or as an inset overlay
#' @description
#' Build a panel of publication-ready time-to-event plots by combining
#' \code{{cifcurve}} (estimation) and \code{{cifplot}} (rendering).
#' Supports grid layout or an inset plot overlay. Panel-wise overrides for
#' y-scale type, labels, limits, and various plot decorations are supported.
#'
#' @details
#' \strong{How it works:}
#' For each panel slot, the function first calls \code{{cifcurve}} to obtain
#' a \code{survfit}-like object, and then calls \code{{cifplot}} to render it.
#' You can pass \code{code.events = list(c(e1, e2, c), ...)} for competing risks
#' (Aalen–Johansen), or \code{code.events = list(c(e1, c), ...)} for standard
#' survival (Kaplan–Meier). If \code{outcome.type} is not specified, each panel's
#' outcome is inferred from the length of its \code{code.events[[i]]}.
#'
#' \strong{Grid vs inset:}
#' - \code{use_inset_element = FALSE} (default): arrange plots in a grid using
#'   \pkg{patchwork} \code{wrap_plots()}. With \code{legend.collect = TRUE}, legends are
#'   collected across panels.
#' - \code{use_inset_element = TRUE}: overlay the second plot into the first one using
#'   \code{patchwork::inset_element()}. Only the first two plots are used; extra plots are ignored.
#'
#' \strong{Notes:}
#' @param plots list of ggplot objects. If supplied, \code{cifpanel()} will skip
#'   estimation/plotting and only compose the provided plots.
#' @param formula A single model formula evaluated in \code{data}, used for all panels
#'   when \code{formulas} is not provided.
#' @param formulas A list of formulas (one per panel). If provided, overrides \code{formula}.
#' @param data A data.frame containing variables used in the formula(s).
#' @param outcome.type Optional vector/list of outcome types per panel.
#' One of \code{"SURVIVAL"} (Kaplan–Meier type) or \code{"COMPETING-RISK"} (Aalen–Johansen type).
#' If \code{NULL} (default), the function automatically infers the outcome type
#' from the data: if the event variable has more than two unique levels,
#' \code{"COMPETING-RISK"} is assumed; otherwise, \code{"SURVIVAL"} is used.
#' You can also use abbreviations such as \code{"S"} or \code{"C"}.
#' Mixed or ambiguous inputs (e.g., \code{c("S", "C")}) trigger automatic
#' detection based on the event coding in \code{data}.

#' @param code.events A list of numeric vectors per panel.
#'   For SURVIVAL: \code{c(code.event1, code.censoring)};
#'   for COMPETING-RISK: \code{c(code.event1, code.event2, code.censoring)}.
#' @param type.x Reserved for future x-scale control (currently ignored).
#' @param type.y Optional vector/list per panel: \code{NULL} (survival) or \code{"risk"} (1-survival).
#' @param label.x,label.y Optional vectors/lists of axis labels per panel.
#' @param label.strata Optional list of character vectors for legend labels per panel
#'   (passed to [cifplot()]).
#' @param order.strata Optional list of character vectors for ordering labels per panel
#'   (passed to [cifplot()]).
#' @param level.strata Optional character vector describing the full set of strata
#'   levels expected in the plots. When supplied, `order.strata` and
#'   `label.strata` are validated against these levels before being applied.
#' @param limits.x,limits.y Optional vectors/lists of numeric length-2 axis limits per panel.
#' @param breaks.x,breaks.y Optional vectors/lists of axis breaks per panel (forwarded to
#'   \code{breaks.x} / \code{breaks.y} in [cifplot()]).
#' @param addConfidenceInterval,addCensorMark,addCompetingRiskMark,addIntercurrentEventMark,addQuantileLine
#'   Optional logical vectors/lists per panel to toggle features in [cifplot()].
#'   If \code{NULL}, sensible defaults are used (CI/Censor on; others off).
#' @param rows.columns.panel Integer vector \code{c(nrow, ncol)} specifying the grid size.
#' @param title.panel,subtitle.panel,caption.panel Optional strings for panel annotation.
#' @param tag_levels.panel Passed to \code{patchwork::plot_annotation(tag_levels = ...)}.
#' @param title.plot Optional length-2 character vector, titles for base/inset plots when
#'   \code{use_inset_element = TRUE}.
#' @param legend.position Position of legends: \code{"top"}, \code{"right"}, \code{"bottom"},
#'   \code{"left"}, or \code{"none"}.
#' @param legend.collect If \code{TRUE} (grid mode), collect legends across subplots.
#' @param use_inset_element If \code{TRUE}, place the second plot as an inset over the first.
#' @param inset.left,inset.bottom,inset.right,inset.top Numeric positions (0–1) of the inset box.
#' @param inset.align_to One of \code{"panel"}, \code{"plot"}, or \code{"full"}.
#' @param inset.legend.position Legend position for the inset plot (e.g., \code{"none"}).
#' @param print.panel If \code{TRUE}, print the composed panel.
#' @param filename.ggsave Character save the composed panel with the path and name specified.
#' @param width.ggsave Numeric specify width of the composed panel.
#' @param height.ggsave Numeric specify height of the composed panel.
#' @param dpi.ggsave Numeric specify dpi of the composed panel.
#' @param survfit.info,axis.info,visual.info,panel.info,style.info,inset.info,print.info,ggsave.info
#'   Optional lists of settings. Each list is merged with the corresponding scalar
#'   arguments (e.g., `axis.info$list` overrides `label.x`, `limits.x`, etc.) so that
#'   existing code using scalar inputs continues to work while bulk updates can be
#'   provided via a single structure.
#' @param ... Additional arguments forwarded to \code{{cifplot}} (e.g., \code{style},
#'   \code{font.family}, \code{font.size}, etc.). Panel-wise overrides provided via explicit
#'   arguments take precedence over \code{...}.
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
#' - **Grid mode** (`use_inset_element = FALSE`, default): plots are arranged with
#'   `patchwork::wrap_plots()` and `plot_layout()`. If `legend.collect = TRUE`,
#'   legends are collected across panels where possible.
#' - **Inset mode** (`use_inset_element = TRUE`): the **second** plot is overlaid
#'   into the **first** using `patchwork::inset_element()`. Only the first two
#'   plots are used; extra plots are ignored. Control the inset box with
#'   `inset.left`, `inset.bottom`, `inset.right`, `inset.top`, and its
#'   reference frame via `inset.align_to` (`"panel"`, `"plot"`, or `"full"`).
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
#' - Tagging: `tag_levels.panel` is passed to `patchwork::plot_annotation()`.
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
#'   use_inset_element = TRUE,
#'   formula = Event(t, epsilon) ~ fruitq,
#'   data = diabetes.complications,
#'   outcome.type = "COMPETING-RISK",
#'   code.events = list(c(2,1,0), c(2,1,0)),
#'   label.y = c("CIF of macrovascular complications", ""),
#'   label.x = c("Years from registration", ""),
#'   limits.y     = list(c(0,1), c(0,0.15)),
#'   inset.left   = 0.40, inset.bottom = 0.45,
#'   inset.right  = 1.00, inset.top    = 0.95,
#'   inset.align_to = "plot",
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
#'          use_inset_element = TRUE,
#'          inset.left   = 0.40, inset.bottom = 0.45,
#'          inset.right  = 1.00, inset.top    = 0.95,
#'          inset.align_to = "plot",
#'          inset.legend.position = "none",
#'          legend.position = "bottom")
#'
#' @importFrom ggplot2 ggplot theme_void ggsave theme element_text labs
#' @importFrom patchwork wrap_plots plot_layout inset_element plot_annotation

#' @name cifpanel
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifcurve()] for KM/AJ estimators; [cifplot()] for display of a CIF; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
cifpanel <- function(
    plots                   = NULL,
    formula                 = NULL,
    formulas                = NULL,
    data                    = NULL,
    subset.condition        = NULL,
    na.action               = na.omit,
    outcome.type            = NULL,
    code.events             = NULL,
    error                   = NULL,
    conf.type               = NULL,
    conf.int                = NULL,
    type.y                  = NULL,
    label.x                 = NULL,
    label.y                 = NULL,
    label.strata            = NULL,
    order.strata            = NULL,
    level.strata            = NULL,
    limits.x                = NULL,
    limits.y                = NULL,
    breaks.x                = NULL,
    breaks.y                = NULL,
    addConfidenceInterval   = NULL,
    addCensorMark           = NULL,
    addCompetingRiskMark    = NULL,
    addIntercurrentEventMark= NULL,
    addQuantileLine         = NULL,
    rows.columns.panel      = c(1, 1),
    title.panel             = NULL,
    subtitle.panel          = NULL,
    caption.panel           = NULL,
    tag_levels.panel        = NULL,
    title.plot              = NULL,
    style                   = "CLASSIC",
    palette                 = NULL,
    font.family             = "sans",
    font.size               = 12,
    legend.position         = "top",
    legend.collect          = FALSE,
    use_inset_element       = FALSE,
    inset.left              = 0.60,
    inset.bottom            = 0.05,
    inset.right             = 0.98,
    inset.top               = 0.45,
    inset.align_to          = c("panel","plot","full"),
    inset.legend.position   = NULL,
    print.panel             = TRUE,
    filename.ggsave         = NULL,
    width.ggsave            = NULL,
    height.ggsave           = NULL,
    dpi.ggsave              = 300,
    survfit.info            = NULL,
    axis.info               = NULL,
    visual.info             = NULL,
    panel.info              = NULL,
    style.info              = NULL,
    inset.info              = NULL,
    print.info              = NULL,
    ggsave.info             = NULL,
    engine                  = "cifplot",
    ...
){
  inset.align_to <- match.arg(inset.align_to)
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
    addConfidenceInterval    = if (is.null(addConfidenceInterval)) NULL else isTRUE(addConfidenceInterval),
    ci.alpha                 = 0.25,
    addRiskTable             = FALSE,
    addEstimateTable         = FALSE,
    addCensorMark            = if (is.null(addCensorMark)) NULL else isTRUE(addCensorMark),
    shape.censor.mark        = 3,
    size.censor.mark         = 2,
    addCompetingRiskMark     = isTRUE(addCompetingRiskMark),
    competing.risk.time      = list(),
    shape.competing.risk.mark= 16,
    size.competing.risk.mark = 2,
    addIntercurrentEventMark = isTRUE(addIntercurrentEventMark),
    intercurrent.event.time  = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark  = 2,
    addQuantileLine          = isTRUE(addQuantileLine),
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
    tag_levels.panel   = tag_levels.panel,
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
    use_inset_element     = use_inset_element,
    inset.align_to        = inset.align_to,
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

  inset.info$inset.align_to <- match.arg(inset.info$inset.align_to, c("panel","plot","full"))

  rows.columns.panel <- panel.info$rows.columns.panel
  title.panel        <- panel.info$title.panel
  subtitle.panel     <- panel.info$subtitle.panel
  caption.panel      <- panel.info$caption.panel
  tag_levels.panel   <- panel.info$tag_levels.panel
  title.plot         <- panel.info$title.plot

  legend.position    <- style.info$legend.position
  legend.collect     <- isTRUE(style.info$legend.collect)

  use_inset_element  <- isTRUE(inset.info$use_inset_element)
  inset.align_to     <- inset.info$inset.align_to
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
    if (isTRUE(use_inset_element)) {
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
          align_to = inset.align_to
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
      tag_levels = tag_levels.panel,
      theme      = theme.panel.unified
    )

    if (isTRUE(print.panel)) print(out_patchwork)
    if (!is.null(filename.ggsave)) {
      if (is.null(width.ggsave))  width.ggsave  <- if (isTRUE(use_inset_element)) 6 else max(6, 5 * rows.columns.panel[2])
      if (is.null(height.ggsave)) height.ggsave <- if (isTRUE(use_inset_element)) 6 else max(6, 5 * rows.columns.panel[1])
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

  K <- max(n_slots, length(code.events))
  code.events <- panel_recycle_to(code.events, K)

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

    # panel 全体の設定をパネルに上書きする
    pa$axis.info    <- modifyList(pa$axis.info %||% list(), axis.info)
    pa$visual.info  <- modifyList(pa$visual.info %||% list(), visual.info)
    pa$style.info   <- modifyList(pa$style.info %||% list(), style.info)
    pa$survfit.info <- modifyList(pa$survfit.info %||% list(), survfit.info)

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

    if (identical(eng_i, "ggsurvfit")) {
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

      p_i <- call_ggsurvfit(
        survfit_object   = sf_i,
        survfit.info     = survfit.info,
        axis.info        = axis_i,
        visual.info      = visual_i,
        panel.info       = panel_i,
        style.info       = style.info,
        ggsave.info      = ggsave.info
      )
      return(p_i)
    } else {
      allowed <- setdiff(names(formals(cifplot_single)), "...")
      if (!is.null(names(pa))) {
        pa <- pa[intersect(names(pa), allowed)]
      }
      if (!"formula_or_fit" %in% names(pa)) {
        pa <- c(list(formula_or_fit = prep$curves[[i]]), pa)
      }
      return(do.call(cifplot_single, pa))
    }
  })

  has_ggsurvfit <- any(vapply(engine.list, identical, logical(1), y = "ggsurvfit"))

  plots <- apply_strata_to_plots(
    plots,
    order_data   = axis.info$order.strata,
    label_map    = axis.info$label.strata,
    touch_colour = !has_ggsurvfit
  )

  if (isTRUE(use_inset_element)) {
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
        right = inset.right, top = inset.top, align_to = inset.align_to
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
    tag_levels = tag_levels.panel,
    theme      = theme.panel.unified
  )

  if (isTRUE(print.panel)) print(out_patchwork)
  if (!is.null(filename.ggsave)) {
    if (is.null(width.ggsave))  width.ggsave  <- if (isTRUE(use_inset_element)) 6 else max(6, 5 * rows.columns.panel[2])
    if (is.null(height.ggsave)) height.ggsave <- if (isTRUE(use_inset_element)) 6 else max(6, 5 * rows.columns.panel[1])
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

call_ggsurvfit <- function(
    survfit_object,
    out_readSurv = NULL,
    survfit.info = NULL,
    axis.info    = NULL,
    visual.info  = NULL,
    panel.info   = NULL,
    style.info   = NULL,
    ggsave.info  = NULL
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
  use_coord_cartesian <- isTRUE(axis.info$use_coord_cartesian)

  addConfidenceInterval         <- visual.info$addConfidenceInterval
  addRiskTable                  <- visual.info$addRiskTable
  addEstimateTable              <- visual.info$addEstimateTable
  symbol.risktable              <- visual.info$symbol.risktable
  font.size.risktable           <- visual.info$font.size.risktable
  addCensorMark                 <- visual.info$addCensorMark
  shape.censor.mark             <- visual.info$shape.censor.mark
  size.censor.mark              <- visual.info$size.censor.mark
  addCompetingRiskMark          <- visual.info$addCompetingRiskMark
  competing.risk.time           <- visual.info$competing.risk.time
  shape.competing.risk.mark     <- visual.info$shape.competing.risk.mark
  size.competing.risk.mark      <- visual.info$size.competing.risk.mark
  addIntercurrentEventMark      <- visual.info$addIntercurrentEventMark
  intercurrent.event.time       <- visual.info$intercurrent.event.time
  shape.intercurrent.event.mark <- visual.info$shape.intercurrent.event.mark
  size.intercurrent.event.mark  <- visual.info$size.intercurrent.event.mark
  addQuantileLine               <- visual.info$addQuantileLine
  quantile                      <- visual.info$quantile

  printEachEvent     <- isTRUE(panel.info$printEachEvent)
  printEachVar       <- isTRUE(panel.info$printEachVar)

  style              <- style.info$style
  palette            <- style.info$palette
  font.family        <- style.info$font.family
  font.size          <- style.info$font.size
  legend.position    <- style.info$legend.position

  filename.ggsave    <- ggsave.info$filename.ggsave
  width.ggsave       <- ggsave.info$width.ggsave
  height.ggsave      <- ggsave.info$height.ggsave
  dpi.ggsave         <- ggsave.info$dpi.ggsave
  ggsave.units       <- ggsave.info$units %||% "in"

  label.strata.map   <- plot_make_label.strata.map(
    survfit_object   = survfit_object,
    label.strata     = label.strata,
    level.strata     = level.strata
  )

  out_cg <- check_ggsurvfit(
    survfit_object   = survfit_object,
    survfit.info     = survfit.info,
    axis.info        = axis.info,
    visual.info      = visual.info,
    style.info       = style.info,
    out_readSurv     = out_readSurv
  )

  res <- plot_reconcile_order_and_labels(
    survfit_object   = survfit_object,
    label.strata.map = label.strata.map,
    level.strata     = level.strata,
    order.strata     = order.strata
  )

  limits_arg                 <- res$limits_arg
  label.strata.map           <- res$label.strata.map
  strata_levels_final        <- res$strata_levels_final
  strata_labels_final        <- res$strata_labels_final
  n_strata_effective         <- length(limits_arg)

  p <- out_cg$out_survfit_object +
    ggplot2::labs(
      x = label.x.user %||% "Time",
      y = label.y.user %||% out_cg$label.y
    )

  if (isTRUE(addConfidenceInterval)) {
    p <- p + add_confidence_interval()
  }
  if (isTRUE(addCensorMark)) {
    p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)
  }
  if (isTRUE(addQuantileLine)) {
    p <- p + ggsurvfit::add_quantile(y_value = quantile)
  }

  apply_add_risktable_strata_symbol <- function (p, symbol.risktable) {
    if (symbol.risktable=="square") {
      p <- p + add_risktable_strata_symbol(symbol = "\U25A0", size = 14)
    } else if (symbol.risktable=="circle") {
      p <- p + add_risktable_strata_symbol(symbol = "\U25CF", size = 14)
    } else if (symbol.risktable=="triangle") {
      p <- p + add_risktable_strata_symbol(symbol = "\U25B2", size = 14)
    }
    p
  }

  if (isTRUE(addEstimateTable) && isTRUE(addRiskTable)) {
    p <- p + add_risktable(
      risktable_stats = c("n.risk", "{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"),
      stats_label     = c("No. at risk", "Estimate (95% CI)"),
      theme           = plot_theme_risktable_font(font.family = font.family, plot.title.size = font.size),
      risktable_group = "risktable_stats",
      size            = font.size.risktable
    )
    p <- apply_add_risktable_strata_symbol(p, symbol.risktable)
  } else if (isTRUE(addEstimateTable)) {
    p <- p + add_risktable(
      risktable_stats = c("{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"),
      stats_label     = c("Estimate (95% CI)"),
      theme           = plot_theme_risktable_font(font.family = font.family, plot.title.size = font.size),
      size            = font.size.risktable
    )
    p <- apply_add_risktable_strata_symbol(p, symbol.risktable)
  } else if (isTRUE(addRiskTable)) {
    p <- p + add_risktable(
      risktable_stats = c("n.risk"),
      stats_label     = c("No. at risk"),
      theme           = plot_theme_risktable_font(font.family = font.family, plot.title.size = font.size),
      size            = font.size.risktable
    )
    p <- apply_add_risktable_strata_symbol(p, symbol.risktable)
  }

  if (isTRUE(addCompetingRiskMark) && length(competing.risk.time)>0) {
    p <- plot_draw_marks(
      p, survfit_object,
      competing.risk.time, out_cg$type.y,
      shape = shape.competing.risk.mark,
      size  = size.competing.risk.mark
    )
  }
  if (isTRUE(addIntercurrentEventMark) && length(intercurrent.event.time)>0) {
    p <- plot_draw_marks(
      p, survfit_object,
      intercurrent.event.time, out_cg$type.y,
      shape = shape.intercurrent.event.mark,
      size  = size.intercurrent.event.mark
    )
  }

  x_max <- plot_make_x_max(survfit_object)
  if (isTRUE(use_coord_cartesian)) {
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

  if (!identical(style, "GGSURVFIT")) {
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
  p
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

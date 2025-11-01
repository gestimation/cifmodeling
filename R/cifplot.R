#' @title Generate a survival or cumulative incidence curve with marks that represent
#' censoring, competing risks and intermediate events
#' @description
#' This function produces the Kaplan–Meier survival or Aalen–Johansen cumulative
#' incidence curve from a unified formula + data interface (\code{Event()} or \code{Surv()} on
#' the left-hand side). It auto-labels axes based on `\code{outcome.type} and \code{type.y}, can
#' add censoring/competing-risk/intercurrent-event marks, and returns a regular \code{ggplot}
#' object (compatible with \code{+} and \code{%+%}). You may also pass a survfit-compatible object directly.
#'
#' @param formula_or_fit A model formula or a survfit object.
#' @param data A data frame containing variables in \code{formula}.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be positive.
#' @param subset.condition Optional character expression to subset \code{data} before analysis.
#' @param na.action Function to handle missing values (default: \code{na.omit} in \pkg{stats}).
#' @param outcome.type
#' Character string specifying the type of time-to-event outcome.
#' One of \code{"SURVIVAL"} (Kaplan–Meier type) or \code{"COMPETING-RISK"} (Aalen–Johansen type).
#' If \code{NULL} (default), the function automatically infers the outcome type
#' from the data: if the event variable has more than two unique levels,
#' \code{"COMPETING-RISK"} is assumed; otherwise, \code{"SURVIVAL"} is used.
#' You can also use abbreviations such as \code{"S"} or \code{"C"}.
#' Mixed or ambiguous inputs (e.g., \code{c("S", "C")}) trigger automatic
#' detection based on the event coding in \code{data}.
#' @param code.event1 Integer code of the event of interest (default \code{1}).
#' @param code.event2 Integer code of the competing event (default \code{2}).
#' @param code.censoring Integer code of censoring (default \code{0}).
#' @param code.events Optional numeric length-3 vector \code{c(event1, event2, censoring)}.
#'   When supplied, it overrides \code{code.event1}, \code{code.event2}, and \code{code.censoring}
#'   (primarily used when \code{printEachEvent = TRUE}).
#' @param error Character specifying variance type used internally. For \code{"SURVIVAL"} typically \code{"greenwood"}.
#'   For \code{"COMPETING-RISK"} pass options supported by \code{calculateAalenDeltaSE()}
#'   (\code{"aalen"}, \code{"delta"}, \code{"none"}).
#' @param conf.type Character transformation for CI on the probability scale (default \code{"arcsine-square root"}).
#' @param conf.int numeric two-sided confidence level (default \code{0.95}).
#' @param type.y \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param label.x Character x-axis labels (default \code{"Time"}).
#' @param label.y Character y-axis labels (default internally set to \code{"Survival"} or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata. When \code{printEachVar = TRUE}, you may
#'   also supply a named list \code{list(var = c("L1","L2", ...))}.
#' @param order.strata Optional ordering of strata levels.
#'   - When \code{printEachVar = TRUE}, supply a named list
#'     \code{list(var = c("L1","L2",...))} for each RHS variable; unmatched levels are dropped.
#'   - When \code{printEachVar = FALSE}, supply a character vector \code{c("L1","L2",...)}
#'     that specifies the display order (legend/risktable) of the single stratification factor.
#'     Levels not listed are dropped.
#'   If \code{label.strata} is a named vector, its names must match the (re-ordered) levels.
#' @param level.strata Optional character vector giving the full set of strata levels
#'   expected in the data. When provided, \code{order.strata} and \code{label.strata} are
#'   validated against these levels before being applied.
#' @param limits.x Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,max(out_readSurv$t))}.
#' @param limits.y Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,1)}.
#' @param breaks.x Numeric vectors for axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vectors for axis breaks (default \code{NULL}).
#' @param use_coord_cartesian Logical specify use of coord_cartesian() (default \code{FALSE}).
#'
#' @param addConfidenceInterval Logical add \code{add_confidence_interval()} to plot.
#' It calls geom_ribbon() (default \code{TRUE}).
#' @param addRiskTable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param addEstimateTable Logical add \code{add_risktable(risktable_stats="estimate (conf.low, conf.high)")}
#' to plot (default \code{FALSE}).
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot. It calls geom_point() (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addCompetingRiskMark Logical add time marks to describe event2 specified by Event(), usually the competing events.
#' It calls geom_point() (default \code{TRUE}).
#' @param competing.risk.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param addIntercurrentEventMark Logical overlay user-specified time marks per strata.
#' It calls geom_point() (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param addQuantileLine Logical add \code{add_quantile()} to plot. It calls geom_segment() (default \code{TRUE}).
#' @param quantile Numeric specify quantile for \code{add_quantile()} (default \code{0.5}).

#' @param printEachEvent Logical. If \code{TRUE} and \code{outcome.type == "COMPETING-RISK"},
#'   \code{cifplot()} internally calls \code{cifpanel} to display both event-specific
#'   cumulative incidence curves side-by-side (event 1 and event 2). Defaults to \code{FALSE}.
#'   Ignored for non-competing-risk outcomes.
#' @param printEachVar Logical. If \code{TRUE}, when multiple covariates are listed
#'   on the right-hand side (e.g., \code{~ a + b + c}), the function produces
#'   a panel of CIF plots, each stratified by one variable at a time.
#' @param rows.columns.panel Optional integer vector \code{c(nrow, ncol)} controlling
#'   the panel layout. If \code{NULL}, an automatic layout is used.

#' @param style Character plot theme controls (default \code{"CLASSIC"}).
#' @param palette Optional character vector specify color palette, e.g. palette=c("blue", "cyan", "navy", "green")
#'  (default \code{NULL}).
#' @param font.family Character plot theme controls (e.g. "sans", "serif", and "mono". default \code{"sans"}).
#' @param font.size Integer plot theme controls (default \code{12}).
#' @param legend.position Character specify position of legend: \code{"top"}, \code{"right"},
#' \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).
#' @param filename.ggsave Character save the \pkg{ggsurvfit} plot with the path and name specified.
#' @param width.ggsave Numeric specify width of the \pkg{ggsurvfit} plot.
#' @param height.ggsave Numeric specify height of the \pkg{ggsurvfit} plot.
#' @param dpi.ggsave Numeric specify dpi of the \pkg{ggsurvfit} plot.
#' @param survfit.info,axis.info,visual.info,panel.info,style.info,print.info,ggsave.info
#'   Optional lists providing batched overrides for the corresponding scalar arguments
#'   (e.g., \code{axis.info$limits.x} overrides \code{limits.x}). Scalar arguments are
#'   merged with these lists to preserve backwards compatibility. \code{inset.info} is
#'   accepted for parity with [cifpanel()] but currently ignored.
#' @param ... Additional arguments forwarded to \code{cifpanel()} when \code{printEachEvent = TRUE}.

#' @details
#' This function calls an internal helper \code{call_ggsurvfit()} which adds confidence intervals,
#' risk table, censoring marks, and optional competing-risk and intercurrent-event marks.
#'
#' When \code{printEachEvent = TRUE}, two panels are created with
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
#' ### Standard error and confidence intervals
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `error` | Standard error for KM: `"greenwood"`, `"tsiatis"`. For CIF: `"aalen"`, `"delta"`, `"none"`. | Automatically chosen `"greenwood"` or `"delta"` |
#' | `conf.type` | Transformation for confidence intervals: `"plain"`, `"log"`, `"log-log"`, `"arcsin"`, `"logit"`, or `"none"`. | `"arcsin"` |
#' | `conf.int` | Two-sided confidence intervals level. | `0.95` |
#'
#' #### Graphical layers
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `addConfidenceInterval` | Add confidence interval ribbon. | `TRUE` |
#' | `addRiskTable` | Add numbers-at-risk table below the plot. | `TRUE` |
#' | `addEstimateTable` | Add estimates & CIs table. | `FALSE` |
#' | `addCensorMark` | Add censoring marks. | `TRUE` |
#' | `addCompetingRiskMark` | Add marks for event2 of "COMPETING-RISK" outcome. | `FALSE` |
#' | `addIntercurrentEventMark` | Add intercurrent event marks at user-specified times. | `FALSE` |
#' | `addQuantileLine` | Add quantile lines. | `FALSE` |
#' | `quantile` | Quantile for `addQuantileLine`. | `0.5` |
#'
#' #### Time for marks
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `competing.risk.time` | **Named list** of numeric vectors that contains times of competing risks. Names must match strata labels. Typically created internally. | `list()` |
#' | `intercurrent.event.time` | **Named list** of numeric vectors that contains times of intercurrent events. Names must match strata labels. Typically created by `extract_time_to_event()`. | `list()` |
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
#' #### Axes and legend
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `limits.x`, `limits.y` | Axis limits (`c(min, max)`). | Auto |
#' | `breaks.x`, `breaks.y` | Tick breaks for x and y axes. | Auto |
#' | `use_coord_cartesian` | For zooming use `coord_cartesian()`. | `FALSE` |
#' | `legend.position` |`"top"`, `"right"`, `"bottom"`, `"left"`, `"none"`. | `"top"` |
#'
#' #### Export (optional convenience)
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `filename.ggsave` | If non-`NULL`, save the plot using `ggsave()`. | `NULL` |
#' | `width.ggsave` | Size passed to `ggsave()`. | `6`|
#' | `height.ggsave` | Size passed to `ggsave()`. | `6`|
#' | `dpi.ggsave` | DPI passed to `ggsave()`. | `300` |
#'
#' **Notes.**
#' - For CIF displays, set `type.y = "risk"`. For survival scale, use `type.y = NULL` or `= "surv"`.
#' - Event coding can be controlled via `code.event1`, `code.event2`, `code.censoring`.
#'   For ADaM-style data, use `code.event1 = 0`, `code.censoring = 1`.
#' - Per-stratum time lists should have names identical to plotted strata labels.

#' @return A \code{ggplot} object. When \code{printEachVar = TRUE}, a \pkg{patchwork}
#'   object is returned with an attribute \code{attr(x, "plots")} containing the
#'   individual \code{ggplot} panels.
#'
#' @examples
#' data(diabetes.complications)
#' cifplot(Event(t,epsilon) ~ fruitq,
#'         data = diabetes.complications,
#'         outcome.type="COMPETING-RISK",
#'         addRiskTable = FALSE,
#'         label.y='CIF of diabetic retinopathy',
#'         label.x='Years from registration')

#' @importFrom ggsurvfit ggsurvfit add_confidence_interval add_risktable add_risktable_strata_symbol add_censor_mark add_quantile
#' @importFrom ggplot2 theme_classic theme_bw element_text labs lims geom_point aes ggsave guides scale_color_discrete scale_fill_discrete element_text element_rect element_blank scale_color_manual scale_fill_manual scale_linetype_manual scale_shape_manual scale_linetype_discrete scale_shape_discrete
#' @importFrom grDevices gray
#' @importFrom patchwork wrap_plots

#' @name cifplot
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifcurve()] for KM/AJ estimators; [cifpanel()] for display of multiple CIFs; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
cifplot <- function(
    formula_or_fit,
    data = NULL,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = c("COMPETING-RISK", "SURVIVAL"),
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.events = NULL,
    error = NULL,
    conf.type = "arcsine-square root",
    conf.int = 0.95,
    type.y = NULL,
    label.x = "Time",
    label.y = NULL,
    level.strata = NULL,
    label.strata = NULL,
    order.strata = NULL,
    limits.x = NULL,
    limits.y = NULL,
    breaks.x = NULL,
    breaks.y = NULL,
    use_coord_cartesian = FALSE,
    addConfidenceInterval = TRUE,
    addRiskTable = TRUE,
    addEstimateTable = FALSE,
    symbol.risktable = "square",
    font.size.risktable = 3,
    addCensorMark = TRUE,
    shape.censor.mark = 3,
    size.censor.mark = 2,
    addCompetingRiskMark = FALSE,
    competing.risk.time = list(),
    shape.competing.risk.mark = 16,
    size.competing.risk.mark = 2,
    addIntercurrentEventMark = FALSE,
    intercurrent.event.time = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark = 2,
    addQuantileLine = FALSE,
    quantile = 0.5,
    printEachEvent = FALSE,
    printEachVar = FALSE,
    rows.columns.panel = NULL,
    style = "CLASSIC",
    palette = NULL,
    font.family = "sans",
    font.size = 12,
    legend.position = "top",
    filename.ggsave = NULL,
    width.ggsave = 6,
    height.ggsave = 6,
    dpi.ggsave = 300,
    survfit.info = NULL,
    axis.info = NULL,
    visual.info = NULL,
    panel.info = NULL,
    style.info = NULL,
    inset.info = NULL,
    print.info = NULL,
    ggsave.info = NULL,
    ...
) {
  dots <- list(...)

  print_panel <- isTRUE(printEachVar) || isTRUE(printEachEvent)
  if (print_panel) {
    legend.position <- "none"
    addRiskTable     <- FALSE
    addEstimateTable <- FALSE
    if (!is.null(label.strata)) {
      .warn("panel_disables_labelstrata")
    }
    if (isTRUE(addRiskTable) || isTRUE(addEstimateTable)) {
      .warn("panel_disables_tables")
    }
  }


  # 1) まずはすべて info にまとめる（ユーザ直指定 > info 引数 > デフォルト）
  infos <- cifplot_build_info(
    error     = error,
    conf.type = conf.type,
    conf.int  = conf.int,

    type.y    = type.y,
    label.x   = label.x,
    label.y   = label.y,
    level.strata = level.strata,
    label.strata = label.strata,
    order.strata = order.strata,
    limits.x  = limits.x,
    limits.y  = limits.y,
    breaks.x  = breaks.x,
    breaks.y  = breaks.y,
    use_coord_cartesian = use_coord_cartesian,

    addConfidenceInterval     = addConfidenceInterval,
    addRiskTable              = addRiskTable,
    symbol.risktable          = symbol.risktable,
    addEstimateTable          = addEstimateTable,
    font.size.risktable       = font.size.risktable,
    addCensorMark             = addCensorMark,
    shape.censor.mark         = shape.censor.mark,
    size.censor.mark          = size.censor.mark,
    addCompetingRiskMark      = addCompetingRiskMark,
    competing.risk.time       = competing.risk.time,
    shape.competing.risk.mark = shape.competing.risk.mark,
    size.competing.risk.mark  = size.competing.risk.mark,
    addIntercurrentEventMark  = addIntercurrentEventMark,
    intercurrent.event.time   = intercurrent.event.time,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark,
    size.intercurrent.event.mark  = size.intercurrent.event.mark,
    addQuantileLine               = addQuantileLine,
    quantile                      = quantile,

    printEachEvent                = printEachEvent,
    printEachVar                  = printEachVar,
    rows.columns.panel            = rows.columns.panel,

    style           = style,
    palette         = palette,
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

  # 2) 最低保証
  style.info$font.family <- style.info$font.family %||% "sans"
  style.info$font.size   <- style.info$font.size   %||% 12

  # 3) strata 情報を正規化（ここで level/order/label の対応を1回で確定させる）
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

  # 4) printEachVarとprintEachEventの排他チェック
  .assert(!(isTRUE(panel.info$printEachVar) && isTRUE(panel.info$printEachEvent)),
          "incompatible_flags",
          which = "printEachVar and printEachEvent")

  # 5) outcome.type をここで決める
  outcome.type <- util_check_outcome_type(
    outcome.type, formula = if (inherits(formula_or_fit, "survfit")) NULL else formula_or_fit,
    data = data
  )

  # 6) code.events を code.event1.. に落とし込む（単一描画のときも EachVar のときも使うのでここで）
  if (!is.null(code.events)) {
    ce <- plot_check_code_events(code.events)
    .assert(length(ce) == 3L, "code_events_len_vec")
    .assert(all(ce == as.integer(ce)), "code_events_integer")
    .assert(ce[1L] != ce[2L], "code_events_distinct")
    code.event1    <- ce[1L]
    code.event2    <- ce[2L]
    code.censoring <- ce[3L]
  }

  # 7) printEachVar = TRUE のときは info をそのまま渡して終了
  if (isTRUE(panel.info$printEachVar)) {
    .assert(!inherits(formula_or_fit, "survfit"), "need_formula_for_printEachVar")
    .assert(inherits(formula_or_fit, "formula"),  "formula_must_be")
    .assert(!is.null(data),                       "need_data")

    return(
      cifplot_printEachVar(
        formula         = formula_or_fit,
        data            = data,
        weights         = weights,
        subset.condition= subset.condition,
        na.action       = na.action,
        outcome.type    = outcome.type,
        code.event1     = code.event1,
        code.event2     = code.event2,
        code.censoring  = code.censoring,
        code.events     = NULL,     # ここは各パネル側で使わないのでNULLでよい
        survfit.info    = survfit.info,
        axis.info       = axis.info,
        visual.info     = visual.info,
        panel.info      = panel.info,
        style.info      = style.info,
        ggsave.info     = ggsave.info,
        !!!dots         # dotsも渡しておく
      )
    )
  }

  # 8) ここからは「単一プロット」のときだけ
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

  # strata の order/label を最後に一応適用
  out_plot <- apply_strata_to_plots(
    list(out_plot),
    order_data = axis.info$order.strata,
    label_map  = axis.info$label.strata
  )[[1]]

  out_plot
}

cifplot_single <- function(
    formula_or_fit,
    data = NULL,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = c("COMPETING-RISK", "SURVIVAL"),
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.events = NULL,
    # info 系
    survfit.info = NULL,
    axis.info    = NULL,
    visual.info  = NULL,
    panel.info   = NULL,
    style.info   = NULL,
    ggsave.info  = NULL,
    ...
) {
  dots <- list(...)

  ## 0) まず全部を list にしておく
  survfit.info <- survfit.info %||% list()
  axis.info    <- axis.info    %||% list()
  visual.info  <- visual.info  %||% list()
  panel.info   <- panel.info   %||% list()
  style.info   <- style.info   %||% list()
  ggsave.info  <- ggsave.info  %||% list()

  # 0.1) 万一listじゃなかったときに包む
  if (!is.list(survfit.info)) survfit.info <- list(value = survfit.info)
  if (!is.list(axis.info))    axis.info    <- list(value = axis.info)
  if (!is.list(visual.info))  visual.info  <- list(value = visual.info)
  if (!is.list(panel.info))   panel.info   <- list(value = panel.info)
  if (!is.list(style.info))   style.info   <- list(style = style.info)
  if (!is.list(ggsave.info))  ggsave.info  <- list(value = ggsave.info)

  # 0.2) ... に style / palette / font.* / legend.position が来てたら style.info に吸収
  if (!is.null(dots$style)) {
    style.info$style <- dots$style
    dots$style <- NULL
  }
  if (!is.null(dots$palette)) {
    style.info$palette <- dots$palette
    dots$palette <- NULL
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

  ## 0.5) ここで必ずデフォルトをマージする（←これがなかった）
  survfit.info <- modifyList(list(
    error     = NULL,
    conf.type = "arcsine-square root",
    conf.int  = 0.95
  ), survfit.info)

  axis.info <- modifyList(list(
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
    use_coord_cartesian = FALSE
  ), axis.info)

  visual.info <- modifyList(list(
    addConfidenceInterval     = TRUE,
    addRiskTable              = FALSE,
    addEstimateTable          = FALSE,
    symbol.risktable          = "square",
    font.size.risktable       = 3,
    addCensorMark             = TRUE,
    shape.censor.mark         = 3,
    size.censor.mark          = 2,
    addCompetingRiskMark      = FALSE,
    competing.risk.time       = list(),
    shape.competing.risk.mark = 16,
    size.competing.risk.mark  = 2,
    addIntercurrentEventMark  = FALSE,
    intercurrent.event.time   = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark  = 2,
    addQuantileLine           = FALSE,
    quantile                  = 0.5
  ), visual.info)

  panel.info <- modifyList(list(
    printEachEvent     = FALSE,
    printEachVar       = FALSE,
    rows.columns.panel = NULL
  ), panel.info)

  style.info <- modifyList(list(
    style           = "CLASSIC",
    palette         = NULL,
    font.family     = "sans",
    font.size       = 12,
    legend.position = "top"
  ), style.info)

  ggsave.info <- modifyList(list(
    filename.ggsave = NULL,
    width.ggsave    = 6,
    height.ggsave   = 6,
    dpi.ggsave      = 300,
    units           = "in"
  ), ggsave.info)

  error     <- survfit.info$error
  conf.type <- survfit.info$conf.type
  conf.int  <- survfit.info$conf.int

  type.y            <- axis.info$type.y
  label.x           <- axis.info$label.x
  label.y           <- axis.info$label.y
  label.strata      <- axis.info$label.strata
  level.strata      <- axis.info$level.strata
  order.strata      <- axis.info$order.strata
  limits.x          <- axis.info$limits.x
  limits.y          <- axis.info$limits.y
  breaks.x          <- axis.info$breaks.x
  breaks.y          <- axis.info$breaks.y
  use_coord_cartesian <- isTRUE(axis.info$use_coord_cartesian)

  addConfidenceInterval        <- visual.info$addConfidenceInterval
  addRiskTable                 <- visual.info$addRiskTable
  addEstimateTable             <- visual.info$addEstimateTable
  symbol.risktable             <- visual.info$symbol.risktable
  font.size.risktable          <- visual.info$font.size.risktable
  addCensorMark                <- visual.info$addCensorMark
  shape.censor.mark            <- visual.info$shape.censor.mark
  size.censor.mark             <- visual.info$size.censor.mark
  addCompetingRiskMark         <- visual.info$addCompetingRiskMark
  competing.risk.time          <- visual.info$competing.risk.time
  shape.competing.risk.mark    <- visual.info$shape.competing.risk.mark
  size.competing.risk.mark     <- visual.info$size.competing.risk.mark
  addIntercurrentEventMark     <- visual.info$addIntercurrentEventMark
  intercurrent.event.time      <- visual.info$intercurrent.event.time
  shape.intercurrent.event.mark<- visual.info$shape.intercurrent.event.mark
  size.intercurrent.event.mark <- visual.info$size.intercurrent.event.mark
  addQuantileLine              <- visual.info$addQuantileLine
  quantile                     <- visual.info$quantile

  printEachEvent     <- isTRUE(panel.info$printEachEvent)
  printEachVar       <- isTRUE(panel.info$printEachVar)
  rows.columns.panel <- panel.info$rows.columns.panel

  style           <- style.info$style
  palette         <- style.info$palette
  font.family     <- style.info$font.family
  font.size       <- style.info$font.size
  legend.position <- style.info$legend.position

  filename.ggsave <- ggsave.info$filename.ggsave
  width.ggsave    <- ggsave.info$width.ggsave
  height.ggsave   <- ggsave.info$height.ggsave
  dpi.ggsave      <- ggsave.info$dpi.ggsave
  ggsave.units    <- ggsave.info$units %||% "in"

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
  level.strata <- axis.info$level.strata
  label.strata <- axis.info$label.strata
  order.strata <- axis.info$order.strata

  if (is.null(outcome.type)) {
    outcome.type <- util_check_outcome_type(formula = if (inherits(formula_or_fit,"survfit")) NULL
                                            else formula_or_fit, data = data, na.action = na.action, auto_message = FALSE)
  } else {
    outcome.type <- match.arg(outcome.type, c("COMPETING-RISK","SURVIVAL"))
  }

  if (isTRUE(printEachEvent)) {
    if (inherits(formula_or_fit, "survfit")) {
      warning("printEachEvent=TRUE requires a formula interface; falling back to single-plot.")
    } else {
      outcome.type_corrected <- outcome.type
      if (is.null(outcome.type_corrected)) {
        outcome.type_corrected <- util_check_outcome_type(
          formula = formula_or_fit, data = data, na.action = na.action, auto_message = FALSE
        )
      }
      if (outcome.type_corrected != "COMPETING-RISK") {
        warning("printEachEvent=TRUE is only for COMPETING-RISK; falling back to single-plot.")
      } else {
        ce_panel <- plot_check_code_events(c(code.event1, code.event2, code.censoring))

        if (!is.null(dots$title.plot)) dots$title.plot <- NULL

        ylabs_vec <- label.y
        if (is.null(ylabs_vec) && !is.null(dots$label.y)) ylabs_vec <- dots$label.y
        if (is.null(ylabs_vec)) {
          ylabs_vec <- plot_default_event_y_labels()
        } else {
          if (length(ylabs_vec) == 1L) ylabs_vec <- rep(ylabs_vec, 2L)
          if (length(ylabs_vec)  > 2L) ylabs_vec <- ylabs_vec[1:2]
        }
        if (!is.null(dots$label.y)) dots$label.y <- NULL

        axis.info.panel <- modifyList(axis.info, list(
          type.y        = type.y,
          label.x       = label.x,
          label.y       = ylabs_vec,
          label.strata  = label.strata,
          order.strata  = order.strata,
          level.strata  = level.strata,
          limits.x      = limits.x,
          limits.y      = limits.y,
          breaks.x      = breaks.x,
          breaks.y      = breaks.y
        ))
        visual.info.panel <- modifyList(visual.info, list(
          addConfidenceInterval         = addConfidenceInterval,
          addRiskTable                  = addRiskTable,
          addEstimateTable              = addEstimateTable,
          addCensorMark                 = addCensorMark,
          shape.censor.mark             = shape.censor.mark,
          size.censor.mark              = size.censor.mark,
          addCompetingRiskMark          = addCompetingRiskMark,
          competing.risk.time           = competing.risk.time,
          shape.competing.risk.mark     = shape.competing.risk.mark,
          size.competing.risk.mark      = size.competing.risk.mark,
          addIntercurrentEventMark      = addIntercurrentEventMark,
          intercurrent.event.time       = intercurrent.event.time,
          shape.intercurrent.event.mark = shape.intercurrent.event.mark,
          size.intercurrent.event.mark  = size.intercurrent.event.mark,
          addQuantileLine               = addQuantileLine,
          quantile                      = quantile
        ))
        panel.info.panel <- modifyList(panel.info, list(
          rows.columns.panel = if (is.null(rows.columns.panel)) c(1L, 2L) else rows.columns.panel
        ))
        ggsave.info.panel <- modifyList(ggsave.info, list(
          filename.ggsave = filename.ggsave,
          width.ggsave    = width.ggsave,
          height.ggsave   = height.ggsave,
          dpi.ggsave      = dpi.ggsave,
          units           = ggsave.units
        ))

        style_cur <- style.info$style
        palette_cur <- style.info$palette
        ff_cur <- style.info$font.family
        fs_cur <- style.info$font.size
        lg_cur <- style.info$legend.position

        panel_args <- list(
          formula      = formula_or_fit,
          data         = data,
          outcome.type = "COMPETING-RISK",
          code.events  = list(ce_panel, c(ce_panel[2L], ce_panel[1L], ce_panel[3L])),
          axis.info    = axis.info.panel,
          visual.info  = visual.info.panel,
          panel.info   = panel.info.panel,
          style.info   = list(
            style           = style_cur,
            palette         = palette_cur,
            font.family     = ff_cur,
            font.size       = fs_cur,
            legend.position = lg_cur
          ),
          ggsave.info  = ggsave.info.panel,
          survfit.info = survfit.info
        )

        if (is.null(dots$rows.columns.panel)) dots$rows.columns.panel <- panel.info.panel$rows.columns.panel
        if (is.null(dots$legend.collect))     dots$legend.collect     <- TRUE
        if (is.null(dots$print.panel))        dots$print.panel        <- FALSE

        panel_out <- do.call(cifpanel, c(panel_args, dots))

        if (is.list(panel_out) && !is.null(panel_out$out_patchwork)) {
          attr(panel_out$out_patchwork, "plots") <- panel_out$plots
          return(panel_out$out_patchwork)
        }
        return(panel_out)
      }
    }
  }

  if (!inherits(formula_or_fit, "survfit")) {
    if (is.null(data)) stop("When `formula` is a formula, `data` must be provided.")
    norm_inputs <- plot_normalize_formula_data(formula_or_fit, data)
    data_working <- norm_inputs$data
    if (!isTRUE(printEachEvent) &&
        isTRUE(addCompetingRiskMark) &&
        length(competing.risk.time) == 0) {
      competing.risk.time <- extract_time_to_event(
        formula_or_fit, data = data_working, which_event = "event2",
        code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring)
    }
    formula_or_fit <- cifcurve(formula_or_fit, data = data_working, weights = weights, subset.condition = subset.condition, na.action = na.action,
                               outcome.type = outcome.type, code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
                               error = error, conf.type = conf.type, conf.int = conf.int)
    formula_or_fit <- plot_survfit_short_strata_names(formula_or_fit)
#    if (!is.null(axis.info$label.strata) && !is.null(formula_or_fit$strata)) {
#      formula_or_fit <- plot_survfit_strata_labels(
#        formula_or_fit, axis.info$label.strata
#      )
#    }
  }
  p <- call_ggsurvfit(
    survfit_object = if (inherits(formula_or_fit, "survfit")) formula_or_fit else stop("..."),
    out_readSurv   = NULL,
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

cifplot_printEachVar <- function(
    rows.columns.panel = NULL,
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = c("COMPETING-RISK", "SURVIVAL"),
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.events = NULL,
    error = NULL,
    conf.type = "arcsine-square root",
    conf.int = 0.95,
    type.y = NULL,
    label.x = "Time",
    label.y = NULL,
    label.strata = NULL,
    level.strata = NULL,
    order.strata = NULL,
    limits.x = NULL,
    limits.y = NULL,
    breaks.x = NULL,
    breaks.y = NULL,
    use_coord_cartesian = FALSE,
    addConfidenceInterval = TRUE,
    addRiskTable = TRUE,
    addEstimateTable = FALSE,
    addCensorMark = TRUE,
    shape.censor.mark = 3,
    size.censor.mark = 2,
    addCompetingRiskMark = FALSE,
    competing.risk.time = list(),
    shape.competing.risk.mark = 16,
    size.competing.risk.mark = 2,
    addIntercurrentEventMark = FALSE,
    intercurrent.event.time = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark = 2,
    addQuantileLine = FALSE,
    quantile = 0.5,
    style = "CLASSIC",
    palette = NULL,
    font.family = "sans",
    font.size = 12,
    legend.position = "top",
    filename.ggsave = NULL,
    width.ggsave = 6,
    height.ggsave = 6,
    dpi.ggsave = 300,
    printEachEvent = FALSE,
    printEachVar = TRUE,
    survfit.info = NULL,
    axis.info    = NULL,
    visual.info  = NULL,
    panel.info   = NULL,
    style.info   = NULL,
    ggsave.info  = NULL,
    ...
){
  dots <- list(...)
  survfit.info.user <- survfit.info
  axis.info.user    <- axis.info
  visual.info.user  <- visual.info
  panel.info.user   <- panel.info
  style.info.user   <- style.info
  ggsave.info.user  <- ggsave.info

  survfit.info <- modifyList(list(
    error     = error,
    conf.type = conf.type,
    conf.int  = conf.int
  ), survfit.info %||% list())

  axis.info <- modifyList(list(
    type.x            = NULL,
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
    use_coord_cartesian = use_coord_cartesian
  ), axis.info %||% list())

  visual.info <- modifyList(list(
    addConfidenceInterval    = addConfidenceInterval,
    addRiskTable             = addRiskTable,
    addEstimateTable         = addEstimateTable,
    addCensorMark            = addCensorMark,
    shape.censor.mark        = shape.censor.mark,
    size.censor.mark         = size.censor.mark,
    addCompetingRiskMark     = addCompetingRiskMark,
    competing.risk.time      = competing.risk.time,
    shape.competing.risk.mark= shape.competing.risk.mark,
    size.competing.risk.mark = size.competing.risk.mark,
    addIntercurrentEventMark = addIntercurrentEventMark,
    intercurrent.event.time  = intercurrent.event.time,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark,
    size.intercurrent.event.mark  = size.intercurrent.event.mark,
    addQuantileLine          = addQuantileLine,
    quantile                 = quantile
  ), visual.info %||% list())

  panel.info <- modifyList(list(
    printEachEvent     = printEachEvent,
    printEachVar       = printEachVar,
    rows.columns.panel = rows.columns.panel
  ), panel.info %||% list())

  style.info <- modifyList(list(
    style           = style,
    palette         = palette,
    font.family     = font.family,
    font.size       = font.size,
    legend.position = legend.position
  ), style.info %||% list())

  ggsave.info <- modifyList(list(
    filename.ggsave = filename.ggsave,
    width.ggsave    = width.ggsave,
    height.ggsave   = height.ggsave,
    dpi.ggsave      = dpi.ggsave,
    units           = "in"
  ), ggsave.info %||% list())

  style.info$font.family <- style.info$font.family %||% "sans"
  style.info$font.size   <- style.info$font.size   %||% 12

  legend.position <- style.info$legend.position
  style           <- style.info$style
  palette         <- style.info$palette
  font.family     <- style.info$font.family
  font.size       <- style.info$font.size

  filename.ggsave <- ggsave.info$filename.ggsave
  width.ggsave    <- ggsave.info$width.ggsave
  height.ggsave   <- ggsave.info$height.ggsave
  dpi.ggsave      <- ggsave.info$dpi.ggsave
  ggsave.units    <- ggsave.info$units %||% "in"

  printEachEvent <- isTRUE(panel.info$printEachEvent)
  printEachVar   <- isTRUE(panel.info$printEachVar)
  rows.columns.panel <- panel.info$rows.columns.panel

  type.y <- axis.info$type.y
  label.x <- axis.info$label.x
  label.y <- axis.info$label.y
  label.strata <- axis.info$label.strata
  order.strata <- axis.info$order.strata
  level.strata <- axis.info$level.strata
  limits.x <- axis.info$limits.x
  limits.y <- axis.info$limits.y
  breaks.x <- axis.info$breaks.x
  breaks.y <- axis.info$breaks.y
  use_coord_cartesian <- isTRUE(axis.info$use_coord_cartesian)

  addConfidenceInterval <- visual.info$addConfidenceInterval
  addRiskTable          <- visual.info$addRiskTable
  addEstimateTable      <- visual.info$addEstimateTable
  addCensorMark         <- visual.info$addCensorMark
  shape.censor.mark     <- visual.info$shape.censor.mark
  size.censor.mark      <- visual.info$size.censor.mark
  addCompetingRiskMark  <- visual.info$addCompetingRiskMark
  competing.risk.time   <- visual.info$competing.risk.time
  shape.competing.risk.mark <- visual.info$shape.competing.risk.mark
  size.competing.risk.mark  <- visual.info$size.competing.risk.mark
  addIntercurrentEventMark  <- visual.info$addIntercurrentEventMark
  intercurrent.event.time   <- visual.info$intercurrent.event.time
  shape.intercurrent.event.mark <- visual.info$shape.intercurrent.event.mark
  size.intercurrent.event.mark  <- visual.info$size.intercurrent.event.mark
  addQuantileLine        <- visual.info$addQuantileLine
  quantile               <- visual.info$quantile

  error     <- survfit.info$error     %||% error
  conf.type <- survfit.info$conf.type %||% conf.type
  conf.int  <- survfit.info$conf.int  %||% conf.int

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
  label.strata <- axis.info$label.strata
  order.strata <- axis.info$order.strata
  level.strata <- axis.info$level.strata

  Terms <- stats::terms(formula, data = data)
  rhs_vars <- attr(Terms, "term.labels")

  .assert(length(rhs_vars) > 0L, "need_rhs_vars_for_printEachVar")
  invalid_terms <- rhs_vars[grepl("[:*()]", rhs_vars)]
  .assert(length(invalid_terms) == 0L, "no_transform_for_printEachVar", which = paste(invalid_terms, collapse = ", "))

  response_str <- deparse(formula[[2L]])
  lab_is_list <- is.list(label.strata)
  ord_is_list <- is.list(order.strata)

  plots <- lapply(rhs_vars, function(var_name) {
    .assert(var_name %in% names(data), "var_not_found", var = var_name)
    dv <- data
    original_vec <- dv[[var_name]]
    strata_var_name <- paste0(var_name, "_strata")

    norm_res <- plot_normalize_strata(original_vec)
    strata_vec <- norm_res$values
    dv[[strata_var_name]] <- strata_vec

    ord_vec <- NULL
    if (!is.null(order.strata)) {
      if (ord_is_list) {
        ord_vec <- order.strata[[var_name]]
      } else if (length(rhs_vars) == 1L) {
        ord_vec <- order.strata
      }
    }
    if (!is.null(ord_vec)) {
      .assert(is.character(ord_vec), "order_strata_char", var = var_name)
      current_levels <- levels(dv[[strata_var_name]])
      keep_levels <- ord_vec[ord_vec %in% current_levels]
      .assert(length(keep_levels) > 0L, "order_strata_overlap_none", var = var_name,
              levels = paste(current_levels, collapse = ", "))
      strata_factor_tmp <- factor(as.character(dv[[strata_var_name]]), levels = keep_levels)
      keep_idx <- !is.na(strata_factor_tmp)
      dv <- dv[keep_idx, , drop = FALSE]
      strata_vec <- droplevels(strata_vec[keep_idx])
      dv[[strata_var_name]] <- factor(as.character(dv[[strata_var_name]]), levels = keep_levels)
      ord_vec <- NULL
    }

    dv <- dv[!is.na(dv[[strata_var_name]]), , drop = FALSE]
    dv[[strata_var_name]] <- droplevels(dv[[strata_var_name]])

    lab_vec <- NULL
    if (!is.null(label.strata)) {
      if (lab_is_list) {
        lab_vec <- label.strata[[var_name]]
      } else if (length(rhs_vars) == 1L) {
        lab_vec <- label.strata
      }
    }
    if (!is.null(lab_vec)) {
      lv <- levels(dv[[strata_var_name]])
      if (!is.null(names(lab_vec)) && any(nzchar(names(lab_vec)))) {
        .assert(all(lv %in% names(lab_vec)), "label_strata_named_mismatch",
                var = var_name, levels = paste(lv, collapse = ", "))
        lab_vec <- lab_vec[lv]
        names(lab_vec) <- lv
      } else {
        .assert(length(lab_vec) == length(lv), "label_strata_len_mismatch",
                var = var_name, nlab = length(lab_vec), nlv = length(lv))
      }
    }
    if (is.null(lab_vec)) lab_vec <- levels(dv[[strata_var_name]])

    single_formula <- stats::as.formula(sprintf("%s ~ %s", response_str, strata_var_name))

    args_var <- c(
      list(
        formula_or_fit             = single_formula,
        data                       = dv,
        weights                    = weights,
        subset.condition           = subset.condition,
        na.action                  = na.action,
        outcome.type               = outcome.type,
        code.event1                = code.event1,
        code.event2                = code.event2,
        code.censoring             = code.censoring,
        code.events                = code.events,
        error                      = error,
        conf.type                  = conf.type,
        conf.int                   = conf.int,
        type.y                     = type.y,
        printEachEvent             = FALSE,
        label.x                    = label.x,
        label.y                    = label.y,
        label.strata               = lab_vec,
        order.strata               = NULL,
        limits.x                   = limits.x,
        limits.y                   = limits.y,
        breaks.x                   = breaks.x,
        breaks.y                   = breaks.y,
        use_coord_cartesian        = use_coord_cartesian,
        addConfidenceInterval      = addConfidenceInterval,
        addRiskTable               = addRiskTable,
        addEstimateTable           = addEstimateTable,
        addCensorMark              = addCensorMark,
        shape.censor.mark          = shape.censor.mark,
        size.censor.mark           = size.censor.mark,
        addCompetingRiskMark       = addCompetingRiskMark,
        competing.risk.time        = competing.risk.time,
        shape.competing.risk.mark  = shape.competing.risk.mark,
        size.competing.risk.mark   = size.competing.risk.mark,
        addIntercurrentEventMark   = addIntercurrentEventMark,
        intercurrent.event.time    = intercurrent.event.time,
        shape.intercurrent.event.mark = shape.intercurrent.event.mark,
        size.intercurrent.event.mark  = size.intercurrent.event.mark,
        addQuantileLine            = addQuantileLine,
        quantile                   = quantile,
        style                      = style,
        palette                    = palette,
        font.family                = font.family,
        font.size                  = font.size,
        legend.position            = legend.position,
        filename.ggsave            = NULL,
        width.ggsave               = width.ggsave,
        height.ggsave              = height.ggsave,
        dpi.ggsave                 = dpi.ggsave
      ),
      dots
    )
    plot_i <- do.call(cifplot_single, args_var)
    plot_i + ggplot2::ggtitle(var_name) +
      ggplot2::labs(color = var_name, fill = var_name, linetype = var_name, shape = var_name)
  })

  nrow <- ncol <- NULL
  if (!is.null(rows.columns.panel)) {
    .assert(is.numeric(rows.columns.panel) && length(rows.columns.panel) == 2L, "rows_columns_panel_len2")
    nrow <- rows.columns.panel[1L]
    ncol <- rows.columns.panel[2L]
  }
  patch <- patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol, guides = "keep")
  attr(patch, "plots") <- plots
  return(patch)
}

#' Plot survival or cumulative incidence curves with ggsurvfit
#'
#' @param survfit_object A \code{survfit} object.
#' @param out_readSurv (optional) List returned by your \code{util_read_surv()} to auto-set x limits.
#' @param conf.type Character transformation for CI on the probability scale (default \code{"arcsine-square root"}).

#' @param addConfidenceInterval Logical add \code{add_confidence_interval()} to plot. It calls geom_ribbon() (default \code{TRUE}).
#' @param addRiskTable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param addEstimateTable Logical add \code{add_risktable(risktable_stats="estimate (conf.low, conf.high)")} to plot (default \code{FALSE}).
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot. It calls geom_point() (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addCompetingRiskMark Logical add time marks to describe event2 specified by Event(), usually the competing events. It calls geom_point() (default \code{TRUE}).
#' @param competing.risk.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param addIntercurrentEventMark Logical overlay user-specified time marks per strata calls geom_point() (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param addQuantileLine Logical add \code{add_quantile()} to plot. It calls geom_segment() (default \code{TRUE}).
#' @param quantile Numeric specify quantile for \code{add_quantile()} (default \code{0.5}).

#' @param type.y \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param label.x Character x-axis labels (default \code{"Time"}).
#' @param label.y Character y-axis labels (default internally set to \code{"Survival"} or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata.

#' @param limits.x Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,max(out_readSurv$t))}.
#' @param limits.y Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,1)}.
#' @param breaks.x Numeric vectors for axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vectors for axis breaks (default \code{NULL}).
#' @param use_coord_cartesian Logical specify use of coord_cartesian() (default \code{FALSE}).

#' @param style Character plot theme controls (default \code{"CLASSIC"}).
#' @param font.family Character plot theme controls (default \code{"sans"}).
#' @param font.size Integer plot theme controls (default \code{14}).
#' @param legend.position Character specify position of legend: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).

#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
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
  survfit.info <- survfit.info %||% list()
  axis.info    <- axis.info    %||% list()
  visual.info  <- visual.info  %||% list()
  panel.info   <- panel.info   %||% list()
  style.info   <- style.info   %||% list()
  ggsave.info  <- ggsave.info  %||% list()

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
  rows.columns.panel <- panel.info$rows.columns.panel

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
  forbid_limits_due_to_order <- res$forbid_limits_due_to_order
  n_strata_effective         <- length(limits_arg)

  p <- out_cg$out_survfit_object +
    ggplot2::labs(x = label.x, y = out_cg$label.y)

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
      p <- p + add_risktable_strata_symbol(symbol = "\U25CF", size = 14)
    }
    return(p)
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

  if (isTRUE(addCompetingRiskMark) && length(competing.risk.time)) {
    p <- plot_draw_marks(
      p, survfit_object,
      competing.risk.time, out_cg$type.y,
      shape = shape.competing.risk.mark,
      size  = size.competing.risk.mark
    )
  }
  if (isTRUE(addIntercurrentEventMark) && length(intercurrent.event.time)) {
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

#  built <- ggplot2::ggplot_build(p)
#  print("ggplot2::ggplot_build(p)$data")
#  print(built$data[[1]]$group)
#  print(strata_labels_final)
#  print(strata_levels_final)
#  print(limits_arg)

#  strata_labels_final <- c("Hi", "China")
#  strata_levels_final <- c("0", "1")
#  limits_arg <- c("0", "1")

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
    survfit.info = NULL,
    axis.info    = NULL,
    visual.info  = NULL,
    style.info   = NULL,
    out_readSurv = NULL
){
  # 0) デフォルト化
  survfit.info <- survfit.info %||% list()
  axis.info    <- axis.info    %||% list()
  visual.info  <- visual.info  %||% list()
  style.info   <- style.info   %||% list()

  #############################################################
  # 1) survfit.info
  #############################################################
  conf.type <- survfit.info$conf.type
  # error, conf.int もあっていいけどここでは conf.type だけ使ってる
  # error     <- survfit.info$error
  # conf.int  <- survfit.info$conf.int

  #############################################################
  # 2) axis.info
  #############################################################
  type.y              <- axis.info$type.y
  label.y             <- axis.info$label.y
  limits.x            <- axis.info$limits.x
  limits.y            <- axis.info$limits.y
  breaks.x            <- axis.info$breaks.x
  breaks.y            <- axis.info$breaks.y
  use_coord_cartesian <- isTRUE(axis.info$use_coord_cartesian)

  #############################################################
  # 3) visual.info
  #############################################################
  addConfidenceInterval         <- visual.info$addConfidenceInterval
  addCensorMark                 <- visual.info$addCensorMark
  addCompetingRiskMark          <- visual.info$addCompetingRiskMark
  addIntercurrentEventMark      <- visual.info$addIntercurrentEventMark
  shape.censor.mark             <- visual.info$shape.censor.mark
  shape.competing.risk.mark     <- visual.info$shape.competing.risk.mark
  shape.intercurrent.event.mark <- visual.info$shape.intercurrent.event.mark

  #############################################################
  # 4) style.info
  #############################################################
  style   <- style.info$style
  palette <- style.info$palette
  #############################################################

  # ===== ここから下は元のロジックそのまま =====

  if (isTRUE(addCensorMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.censor.mark, shape.intercurrent.event.mark)) {
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.intercurrent.event.mark")
  }
  if (isTRUE(addCensorMark) && isTRUE(addCompetingRiskMark) &&
      identical(shape.censor.mark, shape.competing.risk.mark)) {
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.competing.risk.mark")
  }
  if (isTRUE(addCompetingRiskMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.competing.risk.mark, shape.intercurrent.event.mark)) {
    .warn("shape_identical", a = "shape.competing.risk.mark", b = "shape.intercurrent.event.mark")
  }

  is_len2_num <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  is_nondec   <- function(x) all(diff(x) >= 0, na.rm = TRUE)

  # ---- limits.x ----
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
  } else if (!is.null(out_readSurv) && !is.null(out_readSurv$t)) {
    tmax <- suppressWarnings(max(out_readSurv$t, na.rm = TRUE))
    if (!is.finite(tmax) || tmax <= 0) .warn("ors_tmax_bad")
  }

  # ---- limits.y ----
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

    if (isTRUE(addConfidenceInterval)) {
      if (!is.null(upper) && any(upper < limits.y[1] | upper > limits.y[2], na.rm = TRUE))
        .warn("upper_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])
      if (!is.null(lower) && any(lower < limits.y[1] | lower > limits.y[2], na.rm = TRUE))
        .warn("lower_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])
    }
  }

  # ---- breaks check ----
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

  # ---- label.y の自動決定 ----
  if (is.null(label.y)) {
    auto_label <- plot_default_y_label(survfit_object$type, type.y)
    if (!is.null(auto_label)) label.y <- auto_label
  }

  # ---- 信頼区間を補う ----
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

  # ---- type.y を確定 ----
  type.y <- plot_normalize_type_y(type.y)
  target_type <- switch(
    survfit_object$type,
    "kaplan-meier"   = if (identical(type.y, "risk")) "risk" else "survival",
    "aalen-johansen" = if (identical(type.y, "survival")) "survival" else "risk",
    if (identical(type.y, "risk")) "risk" else "survival"
  )
  type.y <- if (identical(target_type, "risk")) "risk" else "survival"

  # ---- 色が1色しかないときは線種を立てる ----
  decide_linetype_flag <- function(style, palette) {
    if (identical(style, "MONOCHROME")) return(TRUE)
    if (is.null(palette)) return(FALSE)
    pal <- palette
    pal <- ifelse(grepl("^#", pal), pal, plot_validate_fix_color(pal))
    length(unique(tolower(pal))) == 1L
  }
  linetype_aes_flag <- decide_linetype_flag(style, palette)

  old_opt <- getOption("ggsurvfit.switch-color-linetype", FALSE)
  out_plot <- ggsurvfit::ggsurvfit(
    survfit_object,
    type        = target_type,
    linetype_aes = linetype_aes_flag
  )

  list(
    out_survfit_object = out_plot,
    label.y            = label.y,
    type.y             = type.y
  )
}

create_rr_text <- function(coefficient, cov, index, omit.conf.int=TRUE, conf.int=0.95) {
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  coef <- coefficient[index]
  coef_se <- sqrt(diag(cov)[index])
  conf_low <- coef - critical_value * coef_se
  conf_high <- coef + critical_value * coef_se
  p_value <- floor(2 * (1 - pnorm(abs(coef) / coef_se)))
  if (omit.conf.int==TRUE) {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digits=2), ", p<0.01")
    else text <- paste0("RR=", round(exp(coef), digits=2), ", p=", p_value)
  } else {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digits=2), " (", round(exp(conf_low), digits=2), " to ", round(exp(conf_high), digits=2), ", p<0.01", ")")
    else text <- paste0("RR=", round(exp(coef), digits=2), " (", round(exp(conf_low), digits=2), " to ", round(exp(conf_high), digits=2), ", p=", p_value, ")")
  }
  return(text)
}


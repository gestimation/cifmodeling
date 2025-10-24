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
#' - \code{type.x} is reserved for future x-scale control (currently not applied).
#' - If both \code{formula} and \code{formulas} are provided, the latter takes precedence.
#'
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
#' @param palette Optional palette specification forwarded to each [cifplot()] call.
#' @param color Single-stratum line color forwarded to [cifplot()] when applicable.
#' @param fill  Single-stratum ribbon fill forwarded to [cifplot()] when applicable.
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
#' `ggplot2::ggsave()` using `width.ggsave`, `height.ggsave`, and `dpi.ggsave`.
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
    plots        = NULL,
    formula      = NULL,
    formulas     = NULL,
    data         = NULL,
    outcome.type = NULL,
    code.events  = NULL,
    type.x       = NULL,
    type.y       = NULL,
    label.x      = NULL,
    label.y      = NULL,
    label.strata = NULL,
    palette = NULL,
    color = NULL,
    fill = NULL,
    limits.x     = NULL,
    limits.y     = NULL,
    breaks.x      = NULL,
    breaks.y      = NULL,
    addConfidenceInterval   = NULL,
    addCensorMark           = NULL,
    addCompetingRiskMark    = NULL,
    addIntercurrentEventMark= NULL,
    addQuantileLine         = NULL,
    rows.columns.panel = c(1, 1),
    title.panel = NULL,
    subtitle.panel = NULL,
    caption.panel = NULL,
    tag_levels.panel = NULL,
    title.plot = NULL,
    legend.position = "top",
    legend.collect = FALSE,
    use_inset_element = FALSE,
    inset.left   = 0.60,
    inset.bottom = 0.05,
    inset.right  = 0.98,
    inset.top    = 0.45,
    inset.align_to = c("panel","plot","full"),
    inset.legend.position = NULL,
    print.panel = TRUE,
    filename.ggsave = NULL,
    width.ggsave = NULL,
    height.ggsave = NULL,
    dpi.ggsave = 300,
    ...
){
  inset.align_to <- match.arg(inset.align_to)
  dots <- list(...)
  if (is.null(dots$palette)) dots$palette <- palette
  if (is.null(dots$color))   dots$color   <- color
  if (is.null(dots$fill))    dots$fill    <- fill
  nrow <- as.integer(rows.columns.panel[1]); ncol <- as.integer(rows.columns.panel[2])
  n_slots <- nrow * ncol

  if (!is.null(plots)) {
    fonts <- panel_extract_fonts(dots)
    theme.panel.unified <- panel_build_theme(font.family = fonts$family, font.size = fonts$size)
    if (!is.list(plots)) {
      stop("`plots` must be a list of ggplot objects.")
    }
    if (length(plots) && !all(vapply(plots, function(p) inherits(p, "ggplot"), logical(1)))) {
      stop("All elements of `plots` must inherit from 'ggplot'.")
    }

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
      ggplot2::ggsave(filename.ggsave, plot = out_patchwork, width = width.ggsave, height = height.ggsave, dpi = dpi.ggsave)
    }

    return(invisible(list(plots = plots_out, out_patchwork = out_patchwork)))
  }

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

  outcome.list <- toL(outcome.type); if (!is.null(outcome.list)) outcome.list <- rec(outcome.list, K)
  typey.list   <- toL(type.y);       if (!is.null(typey.list))   typey.list   <- rec(typey.list, K)
  labely.list  <- toL(label.y);      if (!is.null(labely.list))  labely.list  <- rec(labely.list, K)

  typex.list   <- toL(type.x);       if (!is.null(typex.list))   typex.list   <- rec(typex.list, K)
  labelx.list  <- toL(label.x);      if (!is.null(labelx.list))  labelx.list  <- rec(labelx.list, K)

  limsx.list <- NULL
  if (!is.null(limits.x)) {
    limsx.list <- if (is.list(limits.x)) limits.x else list(limits.x)
    limsx.list <- rec(limsx.list, K)
    for (i in seq_len(K)) {
      li <- limsx.list[[i]]
      if (!is.null(li) && (!is.numeric(li) || length(li) != 2))
        stop(sprintf("limits.x[[%d]] must be numeric length-2 or NULL.", i))
    }
  }
  limsy.list <- NULL
  if (!is.null(limits.y)) {
    limsy.list <- if (is.list(limits.y)) limits.y else list(limits.y)
    limsy.list <- rec(limsy.list, K)
    for (i in seq_len(K)) {
      li <- limsy.list[[i]]
      if (!is.null(li) && (!is.numeric(li) || length(li) != 2))
        stop(sprintf("limits.y[[%d]] must be numeric length-2 or NULL.", i))
    }
  }

  breakx.list <- toL(breaks.x); if (!is.null(breakx.list)) breakx.list <- rec(breakx.list, K)
  breaky.list <- toL(breaks.y); if (!is.null(breaky.list)) breaky.list <- rec(breaky.list, K)

  addCI.list   <- toL(addConfidenceInterval);    if (!is.null(addCI.list))   addCI.list   <- rec(addCI.list, K)
  addCen.list  <- toL(addCensorMark);            if (!is.null(addCen.list))  addCen.list  <- rec(addCen.list, K)
  addCR.list   <- toL(addCompetingRiskMark);     if (!is.null(addCR.list))   addCR.list   <- rec(addCR.list, K)
  addIC.list   <- toL(addIntercurrentEventMark); if (!is.null(addIC.list))   addIC.list   <- rec(addIC.list, K)
  addQ.list    <- toL(addQuantileLine);          if (!is.null(addQ.list))    addQ.list    <- rec(addQ.list, K)

  strata.list  <- toL(label.strata);             if (!is.null(strata.list))  strata.list  <- rec(strata.list, K)

  infer_flag_by_codes <- function(v) if (length(v) == 2L) "S" else if (length(v) == 3L) "C" else NA_character_
  if (!is.null(outcome.list)) {
    outcome.flags <- vapply(outcome.list, panel_norm_outcome, character(1))
  } else {
    outcome.flags <- vapply(code.events, infer_flag_by_codes, character(1))
    if (anyNA(outcome.flags)) .err("infer_outcome_fail")
  }
  panel_validate_code_events(code.events, outcome.flags)

  kill_names <- c()
  if (!is.null(outcome.list)) kill_names <- c(kill_names, "outcome.type")
  if (!is.null(typey.list))   kill_names <- c(kill_names, "type.y")
  if (!is.null(labely.list))  kill_names <- c(kill_names, "label.y")
  if (!is.null(limsy.list))   kill_names <- c(kill_names, "limits.y")
  if (!is.null(typex.list))   kill_names <- c(kill_names, "type.x")
  if (!is.null(labelx.list))  kill_names <- c(kill_names, "label.x")
  if (!is.null(limsx.list))   kill_names <- c(kill_names, "limits.x")
  if (!is.null(breakx.list))  kill_names <- c(kill_names, "breaks.x","breaks.x")
  if (!is.null(breaky.list))  kill_names <- c(kill_names, "breaks.y","breaks.y")
  if (!is.null(addCI.list))   kill_names <- c(kill_names, "addConfidenceInterval")
  if (!is.null(addCen.list))  kill_names <- c(kill_names, "addCensorMark")
  if (!is.null(addCR.list))   kill_names <- c(kill_names, "addCompetingRiskMark")
  if (!is.null(addIC.list))   kill_names <- c(kill_names, "addIntercurrentEventMark")
  if (!is.null(addQ.list))    kill_names <- c(kill_names, "addQuantileLine")
  if (!is.null(strata.list))  kill_names <- c(kill_names, "label.strata")
  dots <- panel_strip_overrides_from_dots(dots, unique(kill_names))

  fonts <- panel_extract_fonts(dots)
  theme.panel.unified <- panel_build_theme(font.family = fonts$family, font.size = fonts$size)

  prep <- panel_prepare(
    K = K,
    formulas = formulas,
    data = data,
    code.events = code.events,
    outcome.flags = outcome.flags,
    outcome.list = outcome.list,
    typey.list = typey.list,
    labely.list = labely.list,
    typex.list = typex.list,
    labelx.list = labelx.list,
    limsx.list = limsx.list,
    limsy.list = limsy.list,
    breakx.list = breakx.list,
    breaky.list = breaky.list,
    addCI.list = addCI.list,
    addCen.list = addCen.list,
    addCR.list = addCR.list,
    addIC.list = addIC.list,
    addQ.list = addQ.list,
    strata.list = strata.list,
    legend.position = legend.position,
    dots = dots,
    fonts = fonts
  )
  plots <- lapply(seq_len(prep$K), function(i) {
    pa <- prep$plot_args[[i]]

    # 念のため: cifplot_single が受け付けない名前は落とす
    allowed <- setdiff(names(formals(cifplot_single)), "...")
    if (!is.null(names(pa))) {
      pa <- pa[intersect(names(pa), allowed)]
    }

    # サバイバル曲線を直接渡す
    args_i <- c(list(formula_or_fit = prep$curves[[i]]), pa)

    do.call(cifplot_single, args_i)
  })

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
    ggplot2::ggsave(filename.ggsave, plot = out_patchwork, width = width.ggsave, height = height.ggsave, dpi = dpi.ggsave)
  }
  invisible(list(plots = plots, out_patchwork = out_patchwork))
}

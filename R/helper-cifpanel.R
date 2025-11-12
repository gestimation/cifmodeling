#' Prepare curve fits and plotting arguments for cifpanel().
#'
#' This internal helper extracts per-panel arguments, evaluates
#' [cifcurve()] for each panel, and prepares the argument lists needed by
#' [cifplot()]. It keeps the logic centralized so that `cifpanel()` only
#' focuses on layout/patchwork composition.
#'
#' @param K Number of panels.
#' @param formulas List of formulas (one per panel).
#' @param data Data frame passed to [cifcurve()].
#' @param code.events List of numeric codes describing events per panel.
#' @param outcome.flags Character vector of normalized outcome flags
#'   (`"S"` or `"C"`).
#' @param outcome.list Optional list of outcome types (passed through to
#'   [cifcurve()]).
#' @param typey.list,labely.list,labelx.list Panel-wise plot
#'   controls forwarded to [cifplot()].
#' @param limsx.list,limsy.list Panel-wise axis limits for [cifplot()].
#' @param breakx.list,breaky.list Panel-wise axis breaks for [cifplot()].
#' @param addCI.list,addCen.list,addCR.list,addIC.list,addQ.list Panel-wise
#'   logical toggles forwarded to [cifplot()].
#' @param strata.list Optional legend labels per panel.
#' @param legend.position Legend position forwarded to [cifplot()].
#' @param dots List of additional arguments (e.g., style, fonts).
#' @param fonts Optional list with `family` and `size`; if `NULL`, extracted
#'   from `dots` via [cifpanel()] helpers.
#'
#' @return A list with elements `curves`, `plot_args`, and `K`.
#' @keywords internal
#' @noRd
panel_prepare <- function(
    K, formulas, data, code.events, outcome.flags,
    outcome.list = NULL,
    typey.list = NULL, labely.list = NULL, labelx.list = NULL,
    limsx.list = NULL, limsy.list = NULL,
    breakx.list = NULL, breaky.list = NULL,
    addCI.list = NULL, addCen.list = NULL, addCR.list = NULL, addIC.list = NULL, addQ.list = NULL,
    strata.list = NULL, legend.position = NULL,
    survfit.info = list(), style.info = list(), dots = list(), fonts = NULL,
    na.action = stats::na.omit                 # ★ 追加（既定は stats::na.omit）
){
  if (is.null(fonts)) fonts <- panel_extract_fonts(dots)

  curves    <- vector("list", length = K)
  plot_args <- vector("list", length = K)

  for (i in seq_len(K)) {
    pair <- code.events[[i]]

    # ★ flags は小文字化して判定
    flag <- tolower(outcome.flags[i])
    if (flag == "s") {
      ce1 <- pair[1]; ce2 <- NULL; cc <- pair[2]
    } else {
      ce1 <- pair[1]; ce2 <- pair[2]; cc <- pair[3]
    }

    norm_inputs <- plot_normalize_formula_data(formulas[[i]], data)
    cur_formula <- norm_inputs$formula %||% formulas[[i]]   # ★ ここで cur_formula を定義
    data_i      <- norm_inputs$data

    args_est <- panel_drop_nulls(c(
      list(
        formula        = cur_formula,            # ★ 正規化済み formula を使用
        data           = data_i,
        outcome.type   = if (!is.null(outcome.list)) outcome.list[[i]] else NULL,
        code.event1    = ce1,
        code.event2    = ce2,
        code.censoring = cc,
        na.action      = na.action               # ★ cifcurve にも渡す
      ),
      survfit.info
    ))

    ot <- args_est$outcome.type
    # ★ 正規化判定にはこのパネルの pair をそのまま渡す
    args_est$outcome.type <- panel_normalize_outcome_type(
      outcome.type = ot,
      code.events  = pair,
      formula      = cur_formula,
      data         = data_i,
      na.action    = na.action
    )

    fit_i <- do.call(cifcurve, args_est)
    curves[[i]] <- fit_i

    args_plot <- list(
      formula_or_fit             = fit_i,
      type.y                     = if (!is.null(typey.list))   typey.list[[i]]   else NULL,
      label.y                    = if (!is.null(labely.list))  labely.list[[i]]  else NULL,
      label.x                    = if (!is.null(labelx.list))  labelx.list[[i]]  else NULL,
      limits.y                   = if (!is.null(limsy.list))   limsy.list[[i]]   else NULL,
      limits.x                   = if (!is.null(limsx.list))   limsx.list[[i]]   else NULL,
      breaks.x                   = if (!is.null(breakx.list))  breakx.list[[i]]  else NULL,
      breaks.y                   = if (!is.null(breaky.list))  breaky.list[[i]]  else NULL,
      add.conf                   = if (!is.null(addCI.list))   addCI.list[[i]]   else TRUE,
      add.censor.mark            = if (!is.null(addCen.list))  addCen.list[[i]]  else TRUE,
      add.competing.risk.mark    = if (!is.null(addCR.list))   addCR.list[[i]]   else FALSE,
      add.intercurrent.event.mark= if (!is.null(addIC.list))   addIC.list[[i]]   else FALSE,
      add.quantile               = if (!is.null(addQ.list))    addQ.list[[i]]    else FALSE,
      label.strata               = if (!is.null(strata.list))  strata.list[[i]]  else NULL,
      font.family                = fonts$family,
      font.size                  = fonts$size
    )
    if (!is.null(dots$style))     args_plot$style     <- dots$style
    if (!is.null(dots$palette))   args_plot$palette   <- dots$palette
    if (!is.null(dots$linewidth)) args_plot$linewidth <- dots$linewidth
    if (!is.null(dots$linetype))  args_plot$linetype  <- dots$linetype
    if (!is.null(legend.position)) args_plot$legend.position <- legend.position

    plot_args[[i]] <- panel_drop_nulls(args_plot)
  }

  list(curves = curves, plot_args = plot_args, K = K)
}

panel_as_formula_global <- function(f) {
  if (is.character(f))   return(stats::as.formula(f, env = .GlobalEnv))
  if (inherits(f, "formula")) return(stats::as.formula(f, env = .GlobalEnv))
  .err("formula_must_be")
}

panel_recycle_to <- function(x, n) {
  if (length(x) == n) return(x)
  if (length(x) == 0L) .err("recycle_empty")
  rep_len(x, n)
}
panel_to_list <- function(x) if (is.null(x)) NULL else if (is.list(x)) x else as.list(x)
panel_drop_nulls <- function(lst) lst[!vapply(lst, is.null, logical(1))]
panel_strip_overrides_from_dots <- function(dots, override_names) {
  if (length(override_names) == 0L) return(dots)
  dots[setdiff(names(dots), override_names)]
}

panel_is_surv <- function(x) {
  x <- tolower(as.character(x %||% ""))
  x %in% c("s", "survival")
}

panel_is_comp <- function(x) {
  x <- tolower(as.character(x %||% ""))
  x %in% c("c", "competing-risk", "competing_risk", "competingrisk")
}

panel_norm_outcome <- function(x) {
  if (panel_is_surv(x)) return("s")
  if (panel_is_comp(x)) return("c")
  stop("Unknown outcome.type: ", x, " (use 's'/'survival' or 'c'/'competing-risk').")
}

panel_extract_fonts <- function(dots) {
  list(
    family = dots$font.family %||% "sans",
    size   = dots$font.size   %||% 12
  )
}

panel_update_rows.columns.panel <- function(
    formulas = NULL,
    code.events = NULL,
    panel.info = list(rows.columns.panel = c(1, 1)),
    n_slots = 1
) {
  len_formulas <- if (!is.null(formulas)) length(formulas) else 0L

  K <- max(n_slots, length(code.events), len_formulas)

  rc <- panel.info$rows.columns.panel
  if (length(rc) == 2L && all(is.finite(rc)) && all(rc == 1) && K > 1L) {
    if (K %% 2L == 0L) {
      rc <- c(K %/% 2L, 2L)
    } else {
      rc <- c((K + 1L) %/% 2L, 2L)
    }
  }

  nrow <- as.integer(rc[1])
  ncol <- as.integer(rc[2])
  n_slots <- nrow * ncol

  list(
    K = K,
    code.events = panel_recycle_to(code.events, K),
    rows.columns.panel = rc,
    nrow = nrow,
    ncol = ncol,
    n_slots = n_slots
  )
}

panel_build_theme <- function(font.family = "sans", font.size = 12) {
  base  <- font.size
  big   <- base * 1.20
  mid   <- base * 1.00
  small <- base * 0.85
  ggplot2::theme(
    text          = ggplot2::element_text(family = font.family, size = base),
    plot.title    = ggplot2::element_text(family = font.family, size = big, face = "bold"),
    plot.subtitle = ggplot2::element_text(family = font.family, size = mid),
    plot.caption  = ggplot2::element_text(family = font.family, size = small),
    axis.title.x  = ggplot2::element_text(family = font.family, size = mid),
    axis.title.y  = ggplot2::element_text(family = font.family, size = mid),
    axis.text.x   = ggplot2::element_text(family = font.family, size = small),
    axis.text.y   = ggplot2::element_text(family = font.family, size = small),
    legend.title  = ggplot2::element_text(family = font.family, size = mid),
    legend.text   = ggplot2::element_text(family = font.family, size = small),
    strip.text    = ggplot2::element_text(family = font.family, size = mid)
  )
}

panel_modify_list <- function(x, y, keep.null = FALSE) {
  if (is.null(y) || !length(y)) return(x)
  for (nm in names(y)) {
    val <- y[[nm]]
    if (is.null(val) && !keep.null) x[[nm]] <- NULL else x[[nm]] <- val

      }
  x
}

panel_normalize_outcome_type <- function(outcome.type, code.events, formula, data, na.action = na.omit) {
  if (!is.null(outcome.type)) {
    return(match.arg(outcome.type, c("competing-risk", "survival")))
  }
  if (!is.null(code.events)) {
    le <- if (is.list(code.events)) lengths(code.events) else length(code.events)
    if (any(le == 3L)) return("competing-risk")
    if (all(le == 2L)) return("survival")
  }
  util_check_outcome_type(formula = formula, data = data, na.action = na.action, auto_message = FALSE)
}

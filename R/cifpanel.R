# ---- helpers (既存) ----
`%||%` <- function(a, b) if (is.null(a)) b else a
.cif_as_formula_global <- function(f) {
  if (is.character(f))   return(stats::as.formula(f, env = .GlobalEnv))
  if (inherits(f, "formula")) return(stats::as.formula(f, env = .GlobalEnv))
  stop("Each formula must be a character string or a formula object.")
}
.cif_recycle_to <- function(x, n) {
  if (length(x) == n) return(x)
  if (length(x) == 0L) stop("Cannot recycle an empty object to nonzero length.")
  rep_len(x, n)
}
.cif_to_list <- function(x) if (is.null(x)) NULL else if (is.list(x)) x else as.list(x)
.cif_drop_nulls <- function(lst) lst[!vapply(lst, is.null, logical(1))]
.cif_strip_overrides_from_dots <- function(dots, override_names) {
  if (length(override_names) == 0L) return(dots)
  dots[setdiff(names(dots), override_names)]
}
.cif_is_surv <- function(x) {
  x <- toupper(as.character(x %||% ""))
  x %in% c("S", "SURVIVAL")
}
.cif_is_comp <- function(x) {
  x <- toupper(as.character(x %||% ""))
  x %in% c("C", "COMPETING-RISK", "COMPETING_RISK", "COMPETINGRISK")
}
.cif_norm_outcome <- function(x) {
  if (.cif_is_surv(x)) return("S")
  if (.cif_is_comp(x)) return("C")
  stop("Unknown outcome.type: ", x, " (use 'S'/'SURVIVAL' or 'C'/'COMPETING-RISK').")
}
.cif_validate_code_events <- function(code_events_list, outcome_flags) {
  stopifnot(length(code_events_list) == length(outcome_flags))
  for (i in seq_along(code_events_list)) {
    pair <- code_events_list[[i]]
    if (outcome_flags[i] == "S") {
      if (!(is.numeric(pair) && length(pair) == 2L))
        stop(sprintf("code.events[[%d]] must be c(event.code1, censoring) for SURVIVAL.", i))
    } else {
      if (!(is.numeric(pair) && length(pair) == 3L))
        stop(sprintf("code.events[[%d]] must be c(event.code1, event.code2, censoring) for COMPETING-RISK.", i))
    }
  }
  invisible(TRUE)
}
.cif_extract_fonts <- function(dots) {
  list(
    family = dots$font.family %||% "sans",
    size   = dots$font.size   %||% 12
  )
}
.cif_build_theme <- function(font.family = "sans", font.size = 12) {
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

# ---- main ----
cifpanel <- function(
    formula = NULL,
    formulas = NULL,
    data = NULL,
    code.events = NULL,
    outcome.type = NULL,
    # y軸まわり（既存）
    type.y       = NULL,
    label.y      = NULL,
    limits.y     = NULL,
    # ★ 新規: x軸/共通オプション（per-panel可）
    type.x       = NULL,
    label.x      = NULL,
    limits.x     = NULL,
    break.x      = NULL,
    break.y      = NULL,
    # ★ 新規: マークやCIなどのフラグ（per-panel可）
    addConfidenceInterval   = NULL,
    addCensorMark           = NULL,
    addCompetingRiskMark    = NULL,
    addIntercurrentEventMark= NULL,
    addQuantileLine         = NULL,
    # ★ 新規: strataラベル（per-panel可）
    label.strata            = NULL,
    # パネル設定
    rows.columns.panel = c(1, 1),
    title.panel = NULL,
    subtitle.panel = NULL,
    caption.panel = NULL,
    tag_levels.panel = NULL,
    title.plot = NULL,
    legend.position = "top",
    legend.collect = FALSE,
    # インセット
    use_inset_element = FALSE,
    inset.left   = 0.60,
    inset.bottom = 0.05,
    inset.right  = 0.98,
    inset.top    = 0.45,
    inset.align_to = c("panel","plot","full"),
    inset.legend.position = NULL,
    # 出力
    print.panel = TRUE,
    filename.panel = NULL,
    width.panel = NULL,
    height.panel = NULL,
    ...
){
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required. install.packages('patchwork')")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  # Grid
  stopifnot(is.numeric(rows.columns.panel), length(rows.columns.panel) == 2, all(rows.columns.panel >= 1))
  nrow <- as.integer(rows.columns.panel[1]); ncol <- as.integer(rows.columns.panel[2])
  n_slots <- nrow * ncol

  # Inputs
  if (is.null(data)) stop("`data` must be provided.")
  if (is.null(code.events) || !is.list(code.events) || length(code.events) == 0)
    stop("Provide non-empty `code.events` as a list of c(event.code1, event.code2[, censor]) per panel.")
  inset.align_to <- match.arg(inset.align_to)

  # Panel count
  K <- max(n_slots, length(code.events))
  code.events <- .cif_recycle_to(code.events, K)

  # Formulas
  use_formula_list <- !is.null(formulas)
  if (use_formula_list && !is.null(formula)) {
    warning("Both `formula` and `formulas` provided. Using `formulas` only.")
  }
  if (!use_formula_list && is.null(formula)) {
    stop("Provide either `formula` (single) or `formulas` (list).")
  }
  if (use_formula_list) {
    stopifnot(is.list(formulas))
    formulas <- lapply(formulas, .cif_as_formula_global)
    formulas <- .cif_recycle_to(formulas, K)
  } else {
    formula  <- .cif_as_formula_global(formula)
    formulas <- rep(list(formula), K)
  }

  # ----- per-panel overrides: to_list + recycle -----
  toL <- .cif_to_list; rec <- .cif_recycle_to

  outcome.list <- toL(outcome.type); if (!is.null(outcome.list)) outcome.list <- rec(outcome.list, K)
  typey.list   <- toL(type.y);       if (!is.null(typey.list))   typey.list   <- rec(typey.list, K)
  labely.list  <- toL(label.y);      if (!is.null(labely.list))  labely.list  <- rec(labely.list, K)

  typex.list   <- toL(type.x);       if (!is.null(typex.list))   typex.list   <- rec(typex.list, K)
  labelx.list  <- toL(label.x);      if (!is.null(labelx.list))  labelx.list  <- rec(labelx.list, K)

  # limits.x / limits.y のバリデーション
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

  # breaks は型が色々ありうるので検証は最小限（そのまま渡す）
  breakx.list <- toL(break.x); if (!is.null(breakx.list)) breakx.list <- rec(breakx.list, K)
  breaky.list <- toL(break.y); if (!is.null(breaky.list)) breaky.list <- rec(breaky.list, K)

  # flags
  addCI.list   <- toL(addConfidenceInterval);    if (!is.null(addCI.list))   addCI.list   <- rec(addCI.list, K)
  addCen.list  <- toL(addCensorMark);            if (!is.null(addCen.list))  addCen.list  <- rec(addCen.list, K)
  addCR.list   <- toL(addCompetingRiskMark);     if (!is.null(addCR.list))   addCR.list   <- rec(addCR.list, K)
  addIC.list   <- toL(addIntercurrentEventMark); if (!is.null(addIC.list))   addIC.list   <- rec(addIC.list, K)
  addQ.list    <- toL(addQuantileLine);          if (!is.null(addQ.list))    addQ.list    <- rec(addQ.list, K)

  # strata labels
  strata.list  <- toL(label.strata);             if (!is.null(strata.list))  strata.list  <- rec(strata.list, K)

  # Outcome flags (default "C") + code.events validation
  outcome.flags <- if (!is.null(outcome.list)) vapply(outcome.list, .cif_norm_outcome, character(1)) else rep("C", K)
  .cif_validate_code_events(code.events, outcome.flags)

  # Dots & fonts（... と衝突するキーを削除）
  dots <- list(...)
  kill_names <- c()
  if (!is.null(outcome.list)) kill_names <- c(kill_names, "outcome.type")
  if (!is.null(typey.list))   kill_names <- c(kill_names, "type.y")
  if (!is.null(labely.list))  kill_names <- c(kill_names, "label.y")
  if (!is.null(limsy.list))   kill_names <- c(kill_names, "limits.y")
  if (!is.null(typex.list))   kill_names <- c(kill_names, "type.x")
  if (!is.null(labelx.list))  kill_names <- c(kill_names, "label.x")
  if (!is.null(limsx.list))   kill_names <- c(kill_names, "limits.x")
  if (!is.null(breakx.list))  kill_names <- c(kill_names, "breaks.x","break.x")
  if (!is.null(breaky.list))  kill_names <- c(kill_names, "breaks.y","break.y")
  if (!is.null(addCI.list))   kill_names <- c(kill_names, "addConfidenceInterval")
  if (!is.null(addCen.list))  kill_names <- c(kill_names, "addCensorMark")
  if (!is.null(addCR.list))   kill_names <- c(kill_names, "addCompetingRiskMark")
  if (!is.null(addIC.list))   kill_names <- c(kill_names, "addIntercurrentEventMark")
  if (!is.null(addQ.list))    kill_names <- c(kill_names, "addQuantileLine")
  if (!is.null(strata.list))  kill_names <- c(kill_names, "label.strata")
  dots <- .cif_strip_overrides_from_dots(dots, unique(kill_names))

  fonts <- .cif_extract_fonts(dots)
  theme.panel.unified <- .cif_build_theme(font.family = fonts$family, font.size = fonts$size)

  plots <- lapply(seq_len(K), function(i) {
    pair <- code.events[[i]]
    if (outcome.flags[i] == "S") {
      ce1 <- pair[1]; ce2 <- NULL; cc <- pair[2]
    } else {
      ce1 <- pair[1]; ce2 <- pair[2]; cc <- pair[3]
    }

    ## --- 1) 推定引数だけで cifcurve() を呼ぶ（描画引数は入れない）---
    call_args_est <- .cif_drop_nulls(list(
      formula        = formulas[[i]],
      data           = data,
      outcome.type   = if (!is.null(outcome.list)) outcome.list[[i]] else NULL,
      code.event1    = ce1,
      code.event2    = ce2,
      code.censoring = cc
      # 必要に応じて：weights, subset.condition, na.action, error, conf.type, conf.int, report.survfit.std.err, label.strata など
      # ただし描画系(addRiskTable等, label.x/y, limits/breaks, print/return…)は絶対に入れない
    ))
    fit_i <- do.call(cifcurve, call_args_est)

    ## --- 2) 描画は cifplot() に渡す（ここに描画引数を集約）---
    call_args_plot <- .cif_drop_nulls(list(
      x = fit_i,
      type.y         = if (!is.null(typey.list))   typey.list[[i]]   else NULL,
      label.y        = if (!is.null(labely.list))  labely.list[[i]]  else NULL,
      label.x        = if (!is.null(labelx.list))  labelx.list[[i]]  else NULL,
      limits.y       = if (!is.null(limsy.list))   limsy.list[[i]]   else NULL,
      limits.x       = if (!is.null(limsx.list))   limsx.list[[i]]   else NULL,
      breaks.x       = if (!is.null(breakx.list))  breakx.list[[i]]  else NULL,
      breaks.y       = if (!is.null(breaky.list))  breaky.list[[i]]  else NULL,
      addConfidenceInterval    = if (!is.null(addCI.list))  addCI.list[[i]]   else TRUE,
      addCensorMark            = if (!is.null(addCen.list)) addCen.list[[i]]  else TRUE,
      addCompetingRiskMark     = if (!is.null(addCR.list))  addCR.list[[i]]   else FALSE,
      addIntercurrentEventMark = if (!is.null(addIC.list))  addIC.list[[i]]   else FALSE,
      addQuantileLine          = if (!is.null(addQ.list))   addQ.list[[i]]    else FALSE,
      label.strata             = if (!is.null(strata.list)) strata.list[[i]]  else NULL,
      # スタイル系は panel の共通設定をそのまま流用
      style           = dots$style        %||% "CLASSIC",
      font.family     = fonts$family,
      font.size       = fonts$size,
      legend.position = legend.position
    ))
    p_i <- do.call(cifplot, call_args_plot)

    p_i
  })

  obj_surv <- do.call(cifcurve, call_args_estimation_only)  # 描画引数は渡さない

  obj_plot <- cifplot(
    obj_surv,
    type.y         = if (!is.null(typey.list))   typey.list[[i]]       else NULL,
    label.y        = if (!is.null(labely.list))  labely.list[[i]]      else NULL,
    limits.y       = if (!is.null(limsy.list))   limsy.list[[i]]       else NULL,
    label.x        = if (!is.null(labelx.list))  labelx.list[[i]]      else NULL,
    limits.x       = if (!is.null(limsx.list))   limsx.list[[i]]       else NULL,
    breaks.x       = if (!is.null(breakx.list))  breakx.list[[i]]      else NULL,
    breaks.y       = if (!is.null(breaky.list))  breaky.list[[i]]      else NULL,
    addConfidenceInterval    = if (!is.null(addCI.list))  addCI.list[[i]]   else TRUE,
    addCensorMark            = if (!is.null(addCen.list)) addCen.list[[i]]  else TRUE,
    addCompetingRiskMark     = if (!is.null(addCR.list))  addCR.list[[i]]   else FALSE,
    addIntercurrentEventMark = if (!is.null(addIC.list))  addIC.list[[i]]   else FALSE,
    addQuantileLine          = if (!is.null(addQ.list))   addQ.list[[i]]    else FALSE,
    label.strata             = if (!is.null(strata.list)) strata.list[[i]]  else NULL,
    style                    = fonts$family %>% { NULL }  # ここは既存ロジックに合わせて
  )
  plots[[i]] <- obj_plot



  # ---- 合成部（inset or grid） ----
  if (isTRUE(use_inset_element)) {
    if (length(plots) < 2L) stop("use_inset_element=TRUE の場合は少なくとも2枚のプロットが必要です。")
    if (length(plots) > 2L) warning("use_inset_element=TRUE: 先頭2枚のみ使用します。")
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
      warning(sprintf("There are %d plots but grid holds %d. Extra plots are dropped.", length(plots), n_slots))
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

  # ---- アノテーション（統一テーマ） ----
  out_patchwork <- out_patchwork + patchwork::plot_annotation(
    title      = title.panel,
    subtitle   = subtitle.panel,
    caption    = caption.panel,
    tag_levels = tag_levels.panel,
    theme      = theme.panel.unified
  )

  if (isTRUE(print.panel)) print(out_patchwork)
  if (!is.null(filename.panel)) {
    if (is.null(width.panel))  width.panel  <- if (isTRUE(use_inset_element)) 6 else max(6, 5 * rows.columns.panel[2])
    if (is.null(height.panel)) height.panel <- if (isTRUE(use_inset_element)) 6 else max(6, 5 * rows.columns.panel[1])
    ggplot2::ggsave(filename.panel, plot = out_patchwork, width = width.panel, height = height.panel, dpi = 300)
  }
  invisible(list(plots = plots, out_patchwork = out_patchwork))
}

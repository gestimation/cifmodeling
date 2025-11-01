plot_fix_palette_vector_arg <- function(p) {
  scs <- p$scales$scales
  if (!length(scs)) return(p)
  for (k in seq_along(scs)) {
    sc <- scs[[k]]
    if (inherits(sc, "ScaleDiscrete") && is.function(sc$palette)) {
      orig <- sc$palette
      sc$palette <- function(n) {
        n1 <- if (length(n) > 1L) max(n, na.rm = TRUE) else n
        unname(orig(n1))
      }
      p$scales$scales[[k]] <- sc
    }
  }
  scs <- p$scales$scales
  if (!length(scs)) return(p)
  for (k in seq_along(scs)) {
    sc <- scs[[k]]
    if (inherits(sc, "ScaleDiscrete")) {
      is_manual <- identical(sc$scale_name, "manual")
      if (is_manual || isTRUE(grepl("manual", paste0(sc$name, collapse=""), ignore.case=TRUE))) {
        sc$scale_name <- "manual"
        if (!inherits(sc, "ScaleDiscreteManual")) {
          class(sc) <- c("ScaleDiscreteManual", class(sc))
        }
        scs[[k]] <- sc
      }
    }
  }
  p$scales$scales <- scs
  return(p)
}

plot_make_dots_clean <- function(dots) {
  plot_drop_panel_only_args <- function(dots) {
    panel_only <- c(
      "rows.columns.panel", "legend.collect", "title.panel", "subtitle.panel",
      "caption.panel", "print.panel", "title.plot", "zoom.position"
    )
    if (length(dots) && !is.null(names(dots))) {
      dots[setdiff(names(dots), panel_only)]
    } else dots
  }
  dots1 <- plot_drop_panel_only_args(dots)
  allowed <- setdiff(names(formals(cifplot_single)), "...")
  drop_extra <- c("printEachVar")
  dots2 <- dots1[setdiff(names(dots1), drop_extra)]

  dots_clean <- if (!is.null(names(dots2))) {
    dots2[intersect(names(dots2), allowed)]
  } else dots2
  return(dots_clean)
}

plot_check_code_events <- function(code_events) {
  if (!(is.numeric(code_events) && length(code_events) == 3L)) {
    .err("code_events_len_vec")
  }
  if (anyNA(code_events)) .err("na", arg = "code.events")
  if (any(!is.finite(code_events))) .err("finite", arg = "code.events")
  out <- as.integer(code_events)
  if (any(abs(code_events - out) > .Machine$double.eps^0.5)) {
    .err("code_events_integer")
  }
  if (out[1L] == out[2L]) .err("code_events_distinct")
  return(out)
}

plot_default_fallback_color <- function(k) {
  hue_fn <- scales::hue_pal(h = c(15, 375), c = 100, l = 65, h.start = 0)
  hue_fn(k)
}

plot_validate_fix_color <- function(x) {
  is_hex <- grepl("^#([0-9A-Fa-f]{3}|[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)
  known_names <- tolower(grDevices::colors())
  name_candidates <- tolower(x[!is_hex])

  bad_name_idx <- which(!(name_candidates %in% known_names))
  if (length(bad_name_idx)) {
    bad_vals <- unique(x[!is_hex][bad_name_idx])
    warning(
      "Unknown color name: ", paste(bad_vals, collapse = ", "),
      " (defaulting to 'black')"
    )
    x[!is_hex][bad_name_idx] <- "black"
  }
  x[!is_hex] <- tolower(x[!is_hex])
  return(x)
}

plot_resolve_palette_color <- function(levels_final, palette, n, fallback_colors = NULL) {
  if (is.null(fallback_colors)) fallback_colors <- plot_default_fallback_color(max(1, n))
  if (is.null(palette)) {
    cols <- rep_len(fallback_colors, max(1, n))
    return(list(values = unname(cols), supplied = FALSE))
  }

  stopifnot(is.character(palette))
  pal <- plot_validate_fix_color(palette)

  if (!is.null(levels_final) && length(levels_final)) {
    cols <- rep_len(fallback_colors, length(levels_final))
    if (is.null(names(palette))) {
      cols <- rep_len(pal, length(levels_final))
    } else {
      cols <- vapply(
        levels_final,
        function(s) pal[[s]] %||% NA_character_,
        FUN.VALUE = character(1)
      )
      if (anyNA(cols)) {
        miss <- which(is.na(cols))
        cols[miss] <- rep_len(fallback_colors, length(levels_final))[miss]
      }
    }
    return(list(values = unname(cols), supplied = TRUE))
  }

  cols <- rep_len(pal, max(1, n))
  list(values = unname(cols), supplied = TRUE)
}

plot_apply_style <- function(
    p,
    style = c("CLASSIC", "BOLD", "FRAMED", "GRID", "GRAY"),
    font.family = "sans",
    font.size = 14,
    legend.position = "top",
    n_strata = 6,
    palette_colors = NULL,
    strata_levels_final = NULL,
    strata_labels_final = NULL
) {
  if (style=="G") style <- "GRID"
  style <- match.arg(style)
  style_theme <- switch(
    style,
    CLASSIC    = plot_style_classic(font.family, font.size, legend.position),
    BOLD       = plot_style_bold(font.family, font.size, legend.position),
    FRAMED     = plot_style_framed(font.family, font.size, legend.position),
    GRID       = plot_style_grid(font.family, font.size, legend.position),
    GRAY       = plot_style_gray(font.family, font.size, legend.position)
  )

  # ★ここで「レジェンドを絶対に出す」方向に倒す
  p <- p + style_theme

  # strata が分かれてるなら、ここで上書き
  if (!is.null(strata_levels_final) && length(strata_levels_final) > 1L) {
    p <- p + ggplot2::theme(legend.position = legend.position %||% "top")
  }

  p
}

plot_apply_all_scales <- function(
    p,
    style,
    palette,
    n_strata,
    strata_levels_final,
    strata_labels_final,
    limits_arg = NULL   # ★ 追加
) {
  p <- plot_strip_mapped_scales(p)

  # ★ ここで使う limits は「上流で確定した limits_arg 」を“唯一の真実”として使う
  lvls <- limits_arg
  # no-overlap のとき、labels も固定しない（NULLにする）
  labs <- if (is.null(lvls) || identical(lvls, character())) NULL else strata_labels_final

  n_effective <- if (!is.null(lvls) && length(lvls)) length(lvls) else n_strata %||% 1L
  n_effective <- max(1L, n_effective)

  use_manual <- !is.null(palette)
  palette_info <- plot_resolve_palette_color(lvls, palette, n_effective)
  col_values   <- unname(rep_len(palette_info$values, n_effective))
  if (!is.null(lvls) && length(lvls)) {
    names(col_values) <- NULL
  }

  if (use_manual) {
    color_scale <- ggplot2::scale_color_manual(
      values = col_values,
      limits = lvls,   # ★ ここが肝
      labels = labs,
      drop   = FALSE
    )
    fill_scale <- ggplot2::scale_fill_manual(
      values = col_values,
      limits = lvls,   # ★
      labels = labs,
      drop   = FALSE,
      guide  = "none"
    )

    color_scale$scale_name <- "manual"
    fill_scale$scale_name  <- "manual"
    if (!inherits(color_scale, "ScaleDiscreteManual")) {
      class(color_scale) <- c("ScaleDiscreteManual", class(color_scale))
    }
    if (!inherits(fill_scale, "ScaleDiscreteManual")) {
      class(fill_scale) <- c("ScaleDiscreteManual", class(fill_scale))
    }
  } else {
    color_scale <- ggplot2::scale_color_discrete(
      limits = lvls, labels = labs, drop = FALSE
    )
    fill_scale <- ggplot2::scale_fill_discrete(
      limits = lvls, labels = labs, drop = FALSE, guide = "none"
    )
  }

  p +
    color_scale +
    fill_scale +
    ggplot2::scale_linetype_discrete(
      limits = lvls,   # ★ linetype/shape も同じ limits を共有
      labels = labs,
      drop = FALSE
    ) +
    ggplot2::scale_shape_discrete(
      limits = lvls,
      labels = labs,
      drop = FALSE,
      guide = "none"
    )
}

plot_apply_all_scales_old <- function(
    p,
    style,
    palette,
    n_strata,
    strata_levels_final,
    strata_labels_final
) {
  p <- plot_strip_mapped_scales(p)
  lvls <- strata_levels_final
  labs <- strata_labels_final
  n_effective <- if (!is.null(lvls) && length(lvls)) length(lvls) else n_strata %||% 1L
  n_effective <- max(1L, n_effective)

  use_manual <- !is.null(palette)
  palette_info <- plot_resolve_palette_color(lvls, palette, n_effective)
  col_values   <- unname(rep_len(palette_info$values, n_effective))
  if (!is.null(lvls) && length(lvls)) {
    names(col_values) <- NULL
  }

  if (use_manual) {
    color_scale <- ggplot2::scale_color_manual(
      values = col_values,
      limits = lvls,
      labels = labs,
      drop   = FALSE
    )
    fill_scale <- ggplot2::scale_fill_manual(
      values = col_values,
      limits = lvls,
      labels = labs,
      drop   = FALSE,
      guide  = "none"
    )

    color_scale$scale_name <- "manual"
    fill_scale$scale_name  <- "manual"
    if (!inherits(color_scale, "ScaleDiscreteManual")) {
      class(color_scale) <- c("ScaleDiscreteManual", class(color_scale))
    }
    if (!inherits(fill_scale, "ScaleDiscreteManual")) {
      class(fill_scale) <- c("ScaleDiscreteManual", class(fill_scale))
    }
  } else {
    color_scale <- ggplot2::scale_color_discrete(
      limits = lvls, labels = labs, drop = FALSE
    )
    fill_scale <- ggplot2::scale_fill_discrete(
      limits = lvls, labels = labs, drop = FALSE, guide = "none"
    )
  }
  p +
    color_scale +
    fill_scale +
    ggplot2::scale_linetype_discrete(
      limits = lvls,
      labels = labs,
      drop = FALSE
    ) +
    ggplot2::scale_shape_discrete(
      limits = lvls,
      labels = labs,
      drop = FALSE,
      guide = "none"
    )
}

plot_strip_mapped_scales <- function(p, aes = c("colour","fill","linetype","shape")) {
  scs <- p$scales$scales
  if (!length(scs)) return(p)
  keep <- vapply(
    scs,
    function(sc) {
      a <- tryCatch(sc$aesthetics, error = function(e) NULL)
      if (is.null(a)) return(TRUE)
      !any(a %in% aes)
    },
    logical(1)
  )
  p$scales$scales <- scs[keep]
  return(p)
}

plot_default_labels <- list(
  strata = list(
    below_median = "Below median",
    above_median = "Above median"
  ),
  y_axis = list(
    survival = "Survival",
    survival_risk = "Risk",
    cif = "Cumulative incidence",
    cif_survival = "1 - cumulative incidence"
  ),
  event_panels = list(
    interest = "Cumulative incidence of interest",
    competing = "Cumulative incidence of competing risk"
  )
)

plot_normalize_type_y <- function(type.y) {
  if (is.null(type.y) || length(type.y) == 0L) return(NULL)
  ty <- tolower(as.character(type.y[1L]))
  if (is.na(ty) || !nzchar(ty)) return(NULL)
  if (ty %in% c("risk", "r")) return("risk")
  if (ty %in% c("survival", "s")) return("survival")
  type.y
}

plot_default_y_label <- function(fit_type, type.y = NULL) {
  ft <- tolower(as.character(fit_type %||% ""))
  ty <- plot_normalize_type_y(type.y)
  if (ft %in% c("kaplan-meier", "kaplan_meier", "km")) {
    if (identical(ty, "risk")) return(plot_default_labels$y_axis$survival_risk)
    return(plot_default_labels$y_axis$survival)
  }
  if (ft %in% c("aalen-johansen", "aalen_johansen", "aj")) {
    if (identical(ty, "survival")) return(plot_default_labels$y_axis$cif_survival)
    return(plot_default_labels$y_axis$cif)
  }
  NULL
}

plot_default_event_y_labels <- function() {
  unname(unlist(plot_default_labels$event_panels, use.names = FALSE))
}

plot_normalize_strata <- function(x, median_threshold = 9L) {
  res <- list(strategy = "factor", threshold = median_threshold)
  if (is.factor(x)) {
    vals <- droplevels(x)
    res$values <- vals
    res$levels <- levels(vals)
    return(res)
  }

  is_datetime <- inherits(x, c("Date", "POSIXt", "POSIXct", "POSIXlt"))
  if (is.numeric(x) && !is_datetime) {
    uniq_vals <- unique(x[!is.na(x) & is.finite(x)])
    if (length(uniq_vals) >= median_threshold) {
      med <- stats::median(x, na.rm = TRUE)
      if (is.finite(med)) {
        below <- plot_default_labels$strata$below_median
        above <- plot_default_labels$strata$above_median
        vals <- ifelse(x > med, above, below)
        vals[is.na(x)] <- NA_character_
        fac <- factor(vals, levels = c(below, above))
        res$strategy <- "median"
        res$cutpoint <- med
        res$values <- fac
        res$levels <- levels(fac)
        return(res)
      }
    }
  }

  fac <- factor(x)
  res$values <- fac
  res$levels <- levels(fac)
  res
}

plot_normalize_formula_data <- function(formula, data, median_threshold = 9L) {
  if (!inherits(formula, "formula")) {
    return(list(data = data, info = list()))
  }

  Terms <- stats::terms(formula, data = data)
  rhs_vars <- attr(Terms, "term.labels")
  if (length(rhs_vars) == 0L) {
    return(list(data = data, info = list()))
  }

  out_data <- data
  info <- list()
  for (var_name in rhs_vars) {
    if (!nzchar(var_name) || grepl("\\(", var_name, fixed = FALSE)) next
    if (!var_name %in% names(out_data)) next
    norm <- plot_normalize_strata(out_data[[var_name]], median_threshold = median_threshold)
    out_data[[var_name]] <- norm$values
    info[[var_name]] <- norm
  }
  list(data = out_data, info = info)
}

plot_style_classic <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 4, family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 4, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

plot_style_bold <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

plot_style_framed <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA),
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 2)
    )
}

plot_style_grid <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 2)
    )
}

plot_style_gray <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = 'gray97'),
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 0.2)
    )
}

plot_draw_marks <- function(p, survfit_object, marks, type.y, shape, size) {
  time <- y <- strata <- NULL
  if (is.null(marks) || !length(marks)) return(p)
  mark_df <- plot_make_mark_data_frame(survfit_object, marks, type.y, extend = TRUE)
  if (is.null(mark_df) || !nrow(mark_df)) return(p)
  if (is.null(survfit_object$strata)) {
    p + geom_point(
      data = mark_df,
      aes(x = time, y = y),
      inherit.aes = FALSE,
      shape = shape,
      size  = size,
      show.legend = FALSE
    )
  } else {
    p + geom_point(
      data = mark_df,
      aes(x = time, y = y, group = strata, colour = strata),
      inherit.aes = FALSE,
      shape = shape,
      size  = size,
      show.legend = FALSE
    )
  }
}

plot_align_mark_keys <- function(fit, marks) {
  canon_str <- function(x) sub("^.*=", "", as.character(x))
  if (is.null(marks) || length(marks) == 0) return(marks)
  fit_names <- if (is.null(fit$strata)) "(all)" else names(fit$strata)
  fit_key   <- canon_str(fit_names)
  out <- list()
  for (k in names(marks)) {
    k_can <- canon_str(k)
    j <- match(k_can, fit_key)
    if (is.na(j)) next
    out[[ fit_names[j] ]] <- marks[[k]]
  }
  return(out)
}

plot_make_mark_data_frame <- function(fit, marks, type.y, extend = TRUE) {
  if (is.null(marks) || length(marks) == 0) return(NULL)
  out_alignMarkKeys <- plot_align_mark_keys(fit, marks)
  strata_names <- if (is.null(fit$strata)) "(all)" else names(fit$strata)

  out <- lapply(names(out_alignMarkKeys), function(st_lab) {
    tt <- out_alignMarkKeys[[st_lab]]
    if (length(tt) == 0) return(NULL)
    fit_st <- if (is.null(fit$strata)) {
      if (!identical(st_lab, "(all)")) return(NULL)
      fit
    } else {
      i <- match(st_lab, strata_names); if (is.na(i)) return(NULL)
      fit[i]
    }
    sm <- summary(fit_st, times = tt, extend = extend)
    y  <- if (identical(type.y, "risk")) 1 - sm$surv else sm$surv
    data.frame(strata = st_lab, time = sm$time, y = y, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

plot_theme_risktable_font <- function(
    axis.text.y.size = 10,
    plot.title.size = 10.75,
    font.family = "sans"
) {
  list(
    ggplot2::theme_bw(base_family = font.family),
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(
        size = 9, vjust = 1, hjust = 1, family = font.family
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(
        size = axis.text.y.size, colour = "black", face = "plain", family = font.family
      ),
      plot.margin = ggplot2::margin(t = 0, r = 5.5, b = 0, l = 5.5, unit = "pt"),
      plot.title  = ggplot2::element_text(
        hjust = 0, vjust = 0, size = plot.title.size, family = font.family
      ),
      legend.position = "none"
    ),
    ggplot2::xlab(NULL),
    ggplot2::ylab(NULL)
  )
}

plot_make_x_max <- function(sf) {
  if (!is.null(sf$time)) {
    xm <- suppressWarnings(max(sf$time, na.rm = TRUE))
    if (is.finite(xm)) return(xm)
  }
  sm <- try(suppressWarnings(summary(sf)), silent = TRUE)
  if (!inherits(sm, "try-error") && !is.null(sm$time)) {
    xm <- suppressWarnings(max(sm$time, na.rm = TRUE))
    if (is.finite(xm)) return(xm)
  }
  return(1)
}

plot_make_label.strata.map <- function(survfit_object,
                                       label.strata,
                                       level.strata = NULL) {
  if (is.null(label.strata)) return(NULL)

  cur_lvls_full  <- if (!is.null(survfit_object$strata)) unique(names(survfit_object$strata)) else NULL
  cur_lvls_short <- if (!is.null(cur_lvls_full)) sub(".*?=", "", cur_lvls_full) else NULL

  if ((is.null(names(label.strata)) || !any(nzchar(names(label.strata))))) {
    if (!is.null(level.strata)) {
      if (!is.null(cur_lvls_full) && length(label.strata) != length(cur_lvls_full)) {
        stop("`label.strata` length must match the number of strata.")
      }

      lbl_map <- stats::setNames(as.character(label.strata),
                                 as.character(level.strata))

      if (!is.null(cur_lvls_short) &&
          setequal(as.character(level.strata), cur_lvls_short)) {
        idx <- match(as.character(level.strata), cur_lvls_short)
        names(lbl_map) <- cur_lvls_full[idx]
      }

      return(lbl_map)
    }
    if (length(label.strata) != length(cur_lvls_full)) {
      stop("`label.strata` length must match the number of strata.")
    }
    lbl_map <- stats::setNames(label.strata, cur_lvls_full)
    return(lbl_map)
  }

  nm <- names(label.strata)

  if (!is.null(cur_lvls_full) && all(nm %in% cur_lvls_full)) {
    return(label.strata)
  }

  if (!is.null(cur_lvls_short) && all(nm %in% cur_lvls_short)) {
    idx <- match(nm, cur_lvls_short)
    names(label.strata) <- cur_lvls_full[idx]
    return(label.strata)
  }

  if (length(label.strata) == length(cur_lvls_full)) {
    lbl_map <- stats::setNames(unname(label.strata), cur_lvls_full)
    return(lbl_map)
  }

  stop("`label.strata` names must be strata levels (full or short).")
}

plot_survfit_strata_labels <- function(survfit_object,
                                       strata.label,
                                       recycle_ok = FALSE) {
  if (is.null(survfit_object)) return(survfit_object)

  nms <- names(survfit_object$strata)

  # strata に名前がないときは「文字列そのものを名前にする」しかできない
  if (is.null(nms) || !length(nms)) {
    new_nms <- as.character(strata.label)
    if (length(new_nms) != length(survfit_object$strata)) {
      if (recycle_ok) {
        new_nms <- rep_len(new_nms, length(survfit_object$strata))
      } else {
        stop("Length of `strata.label` does not match number of strata.", call. = FALSE)
      }
    }
    names(survfit_object$strata) <- new_nms
    return(survfit_object)
  }

  # ここからは「もともと名前がある」ケース
  # nms が "var=level" 形式かもしれないので左右に割る
  has_eq <- grepl("=", nms, fixed = TRUE)
  lhs    <- ifelse(has_eq, sub("^\\s*([^=]+)\\s*=.*$", "\\1", nms, perl = TRUE), NA_character_)
  rhs    <- ifelse(has_eq, sub("^.*?=\\s*", "", nms, perl = TRUE), nms)

  lab <- strata.label

  # 1) label が full 名できている（"fruitq1=0" など）→ そのまま採用
  if (!is.null(names(lab)) && length(intersect(names(lab), nms)) == length(nms)) {
    new_rhs <- unname(lab[nms])

    # 2) label が short 名できている（"0","1" など）→ 右側にマッチさせる
  } else if (!is.null(names(lab)) && length(intersect(names(lab), rhs)) == length(nms)) {
    new_rhs <- unname(lab[rhs])

    # 3) 名前なしベクトルなら長さで対応
  } else {
    lab <- as.character(lab)
    if (length(lab) != length(nms)) {
      if (recycle_ok) lab <- rep_len(lab, length(nms))
      else stop("Length of `strata.label` does not match number of strata.", call. = FALSE)
    }
    new_rhs <- lab
  }

  # ★★ ここが重要 ★★
  # "=" が元からあったものは、左側（変数名）を残して右側だけ差し替える
  new_nms <- character(length(nms))
  for (i in seq_along(nms)) {
    if (isTRUE(has_eq[i]) && !is.na(lhs[i])) {
      new_nms[i] <- paste0(lhs[i], "=", new_rhs[i])
    } else {
      # "=" がなかったやつは右側をそのまま名前にする
      new_nms[i] <- new_rhs[i]
    }
  }

  names(survfit_object$strata) <- new_nms
  return(survfit_object)
}

plot_reconcile_order_and_labels <- function(
    survfit_object,
    level.strata    = NULL,   # ユーザーが「元の層は0/1だよ」と言いたいとき
    order.strata    = NULL,   # 並べ順の指定（"0","1" でも "x=0","x=1" でもOK）
    label.strata.map = NULL   # 表示名の指定（名前あり/なしどっちでもOK）
) {
  # ------------------------------------------------------------
  # 0. survfit からいま使ってる層名をとる
  # ------------------------------------------------------------
  cur_lvls_full <- NULL
  if (!is.null(survfit_object$strata)) {
    cur_lvls_full <- unique(names(survfit_object$strata))
  }

  # "x=0" みたいな形なら右側を short にする
  cur_lvls_short <- NULL
  if (!is.null(cur_lvls_full)) {
    if (any(grepl("=", cur_lvls_full, fixed = TRUE))) {
      cur_lvls_short <- sub(".*?=", "", cur_lvls_full)
    } else {
      # すでに short になってるパターン
      cur_lvls_short <- cur_lvls_full
    }
  }

  # 実データが持ってる順をベースにする
  base_levels <- cur_lvls_full

  # ------------------------------------------------------------
  # 1. label.strata.map が「名前なし」で来たら名前をつける
  # ------------------------------------------------------------
  if (!is.null(label.strata.map)) {
    nm <- names(label.strata.map)
    if (is.null(nm) || !any(nzchar(nm))) {
      # 名前が無いなら、まず level.strata を優先
      if (!is.null(level.strata)) {
        label.strata.map <- stats::setNames(as.character(label.strata.map),
                                            as.character(level.strata))
      } else if (!is.null(cur_lvls_full)) {
        # それもなければ、実データのキーにつける
        label.strata.map <- stats::setNames(as.character(label.strata.map),
                                            cur_lvls_full)
      }
    }
  }

  # ------------------------------------------------------------
  # 2. ラベル指定がまったく無いときのデフォルト
  #    → 実データが "x=0","x=1" で short が "0","1" なら short を表示にする
  # ------------------------------------------------------------
  if (is.null(label.strata.map)) {
    if (!is.null(cur_lvls_full) &&
        !is.null(cur_lvls_short) &&
        length(cur_lvls_full) == length(cur_lvls_short) &&
        !identical(cur_lvls_full, cur_lvls_short)) {
      # 表示は short, key は full
      label.strata.map <- stats::setNames(cur_lvls_short, cur_lvls_full)
    } else if (!is.null(cur_lvls_full)) {
      # そのまま表示
      label.strata.map <- stats::setNames(cur_lvls_full, cur_lvls_full)
    } else {
      # strata 自体が無いケース（単一カーブ）
      label.strata.map <- NULL
    }
  }

  # ------------------------------------------------------------
  # 3. order.strata を “実データのキー” に寄せる
  #    ユーザーが "0","1" で渡しても "x=0","x=1" にくっつける
  # ------------------------------------------------------------
  map_to_full <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    out <- rep(NA_character_, length(x))

    # そのまま full に当ててみる
    if (!is.null(cur_lvls_full)) {
      hit_full <- x %in% cur_lvls_full
      out[hit_full] <- x[hit_full]
    }

    # full に無かったやつは short → full で当てる
    if (!is.null(cur_lvls_full) && !is.null(cur_lvls_short)) {
      hit_short <- x %in% cur_lvls_short
      out[hit_short] <- cur_lvls_full[match(x[hit_short], cur_lvls_short)]
    }

    out
  }
  order_full <- map_to_full(order.strata)

  # ------------------------------------------------------------
  # 4. ラベルを base_levels の順にそろえる（足りなければ埋める）
  # ------------------------------------------------------------
  if (!is.null(label.strata.map) && !is.null(base_levels)) {
    missing <- setdiff(base_levels, names(label.strata.map))
    if (length(missing)) {
      # 足りないところはキーそのまま
      label.strata.map <- c(label.strata.map,
                            stats::setNames(missing, missing))
    }
    # 実データの並びに並べ替える
    label.strata.map <- label.strata.map[base_levels]
  }

  # ------------------------------------------------------------
  # 5. limits を決める（これは ggplot の scale_* にそのまま渡すやつ）
  # ------------------------------------------------------------
  limits_arg <- base_levels
  used_order <- FALSE
  forbid_limits_due_to_order <- FALSE

  if (!is.null(order_full) && !all(is.na(order_full))) {
    ord_clean <- unique(order_full[!is.na(order_full)])
    known <- ord_clean[ord_clean %in% base_levels]
    rest  <- setdiff(base_levels, known)
    limits_arg <- c(known, rest)
    used_order <- TRUE
  }

  strata_levels_final <- if (length(limits_arg)) limits_arg else NULL
  strata_labels_final <- if (!is.null(label.strata.map)) unname(label.strata.map) else NULL

  list(
    limits_arg             = limits_arg,
    label.strata.map       = label.strata.map,
    strata_levels_final    = strata_levels_final,
    strata_labels_final    = strata_labels_final,
    forbid_limits_due_to_order = forbid_limits_due_to_order,
    used_order             = used_order
  )
}

plot_reconcile_order_and_labels_old <- function(
    cur_lvls_full,
    cur_lvls_short,
    level.strata = NULL,
    order.strata = NULL,
    label.strata.map = NULL
) {

  # 1) label.strata.map が「名前なしベクトル」で来たら level.strata で名前を付ける
  if (!is.null(label.strata.map)) {
    nm <- names(label.strata.map)
    if (is.null(nm) || !any(nzchar(nm))) {
      label.strata.map <- stats::setNames(as.character(label.strata.map),
                                          as.character(level.strata))
    }
  }

  # 1.5) ★ ラベル指定がまったく無いときのデフォルトを作る
  if (is.null(label.strata.map)) {
    if (!is.null(cur_lvls_short)) {
      # 例: c("A","B")
      label.strata.map <- stats::setNames(cur_lvls_short, cur_lvls_full)
    } else {
      # short がないなら full をそのまま表示
      label.strata.map <- stats::setNames(cur_lvls_full, cur_lvls_full)
    }
  }

  # 2) ラベルが short 名で来てたら full に直す
  .remap_to_full_if_needed <- function(lbl_map, full, short) {
    if (is.null(lbl_map) || is.null(full) || is.null(short)) return(lbl_map)
    nm <- names(lbl_map); if (!length(nm)) return(lbl_map)
    if (!all(nm %in% full) && all(nm %in% short)) {
      idx <- match(nm, short)
      names(lbl_map) <- full[idx]
    }
    lbl_map
  }
  label.strata.map <- .remap_to_full_if_needed(label.strata.map, cur_lvls_full, cur_lvls_short)

  # 3) order.strata も short → full に寄せる
  map_to_full <- function(x) {
    if (is.null(x)) return(NULL)
    out <- rep(NA_character_, length(x))
    if (!is.null(cur_lvls_full)) {
      hit_f <- x %in% cur_lvls_full
      out[hit_f] <- x[hit_f]
    }
    if (!is.null(cur_lvls_short)) {
      hit_s <- x %in% cur_lvls_short
      out[hit_s] <- cur_lvls_full[match(x[hit_s], cur_lvls_short)]
    }
    out
  }
  order_full_raw <- map_to_full(order.strata)

  # 4) 実データにある順をベースに
  base_levels <- cur_lvls_full %||% as.character(level.strata)

  # 5) ラベルを base_levels の順にそろえる（足りなければ追加）
  if (!is.null(label.strata.map)) {
    missing <- setdiff(base_levels, names(label.strata.map))
    if (length(missing)) {
      label.strata.map <- c(label.strata.map,
                            stats::setNames(missing, missing))
    }
    label.strata.map <- label.strata.map[base_levels]
  }

  limits_arg <- NULL
  used_order <- FALSE
  forbid_limits_due_to_order <- FALSE

  # 6) order.strata があればそれを最優先
  if (!is.null(order_full_raw) && !all(is.na(order_full_raw))) {
    ord_clean <- unique(order_full_raw[!is.na(order_full_raw)])
    known <- ord_clean[ord_clean %in% base_levels]
    rest  <- setdiff(base_levels, known)
    limits_arg <- c(known, rest)
    used_order <- TRUE

    if (!is.null(label.strata.map)) {
      label.strata.map <- label.strata.map[limits_arg]
    }
  } else {
    limits_arg <- base_levels
  }

  strata_levels_final <- if (length(limits_arg)) limits_arg else NULL
  strata_labels_final <- if (!is.null(label.strata.map)) unname(label.strata.map) else NULL

  list(
    limits_arg = limits_arg,
    label.strata.map = label.strata.map,
    strata_levels_final = strata_levels_final,
    strata_labels_final = strata_labels_final,
    forbid_limits_due_to_order = forbid_limits_due_to_order,
    used_order = used_order
  )
}

plot_ensure_factor_strata <- function(formula, data) {
  rhs <- all.vars(update(formula, . ~ .))
  rhs <- setdiff(rhs, all.vars(update(formula, . ~ 0)))
  for (v in rhs) {
    if (v %in% names(data)) {
      x <- data[[v]]
      if (is.numeric(x) || is.integer(x) || is.logical(x)) {
        data[[v]] <- factor(x)
      }
    }
  }
  list(formula = formula, data = data)
}

plot_needs_survfit_label_update <- function(
    survfit_object,
    label.strata = NULL,
    order.strata = NULL,
    level.strata = NULL
) {
  if (!inherits(survfit_object, "survfit")) {
    return(FALSE)
  }

  cur_lvls_full <- if (!is.null(survfit_object$strata)) unique(names(survfit_object$strata)) else NULL
  if (is.null(cur_lvls_full)) {
    return(FALSE)
  }

  cur_lvls_short <- sub(".*?=", "", cur_lvls_full)

  if (is.null(label.strata) && is.null(order.strata) && is.null(level.strata)) {
    return(FALSE)
  }

  cand <- character()

  if (!is.null(label.strata)) {
#    if (!is.null(names(label.strata)) && any(nzchar(names(label.strata)))) {
#      cand <- c(cand, names(label.strata))
#    } else {
      cand <- c(cand, as.character(label.strata))
#    }
  }

#  if (!is.null(order.strata)) {
#    cand <- c(cand, as.character(order.strata))
#  }

#  if (!is.null(level.strata)) {
#    cand <- c(cand, as.character(level.strata))
#  }

  cand <- cand[!is.na(cand) & nzchar(cand)]

  if (!length(cand)) {
    print("C")
    return(FALSE)
  }

  same_short <- setequal(cur_lvls_short, cand)

  if (same_short) {
    return(FALSE)
  }

  same_full <- setequal(cur_lvls_full, cand)
  if (same_full) {
    return(FALSE)
  }

  TRUE
}


# =========================================================
# 凡例を最後にゴリッと復活させるやつ
# =========================================================
plot_force_strata_legend <- function(p,
                                     strata_levels = NULL,
                                     strata_labels = NULL) {

  # 1) まず colour の scale を拾う
  scs <- p$scales$scales
  has_colour <- FALSE
  if (length(scs)) {
    for (i in seq_along(scs)) {
      sc <- scs[[i]]
      aes <- tryCatch(sc$aesthetics, error = function(e) NULL)
      if (!is.null(aes) && any(aes %in% c("colour", "color"))) {
        # ここでガイドを legend に戻す
        sc$guide <- "legend"

        # 上流から渡ってきた順番やラベルがあればここで上書き
        if (!is.null(strata_levels)) sc$limits <- strata_levels
        if (!is.null(strata_labels)) sc$labels <- strata_labels

        p$scales$scales[[i]] <- sc
        has_colour <- TRUE
      }
    }
  }

  # 2) それでも無いときは scale_color_discrete を1本足す
  if (!has_colour) {
    p <- p + ggplot2::scale_color_discrete(
      limits = strata_levels,
      labels = strata_labels,
      drop   = FALSE,
      guide  = "legend"
    )
  }

  # 3) guides() でもう一度強制
  p <- p + ggplot2::guides(
    colour  = ggplot2::guide_legend(override.aes = list(fill = NA)),
    linetype = ggplot2::guide_legend()
  )

  p
}

plot_survfit_short_strata_names <- function(survfit_object) {
  if (is.null(survfit_object) || is.null(survfit_object$strata)) return(survfit_object)

  nm <- names(survfit_object$strata)
  if (is.null(nm) || !length(nm)) return(survfit_object)

  # "=something" の形式だけを対象にする
  idx <- grepl("=", nm, fixed = TRUE)
  if (!any(idx)) return(survfit_object)

  short <- sub(".*?=", "", nm[idx])
  short <- trimws(short)

  nn <- nm
  nn[idx] <- short
  names(survfit_object$strata) <- nn
  survfit_object
}

cifplot_build_info <- function(
  error,
  conf.type,
  conf.int,

  type.y,
  label.x,
  label.y,
  level.strata,
  label.strata,
  order.strata,
  limits.x,
  limits.y,
  breaks.x,
  breaks.y,
  use_coord_cartesian,

  addConfidenceInterval,
  addRiskTable,
  symbol.risktable,
  addEstimateTable,
  symbol.estimatetable,
  addCensorMark,
  shape.censor.mark,
  size.censor.mark,
  addCompetingRiskMark,
  competing.risk.time,
  shape.competing.risk.mark,
  size.competing.risk.mark,
  addIntercurrentEventMark,
  intercurrent.event.time,
  shape.intercurrent.event.mark,
  size.intercurrent.event.mark,
  addQuantileLine,
  quantile,

  printEachEvent,
  printEachVar,
  rows.columns.panel,

  style,
  palette,
  font.family,
  font.size,
  legend.position,

  filename.ggsave,
  width.ggsave,
  height.ggsave,
  dpi.ggsave,

  # ユーザーがすでに info を渡してきた場合
  survfit.info = NULL,
  axis.info    = NULL,
  visual.info  = NULL,
  panel.info   = NULL,
  style.info   = NULL,
  ggsave.info  = NULL
) {

  survfit.info <- modifyList(list(
    error     = error,
    conf.type = conf.type,
    conf.int  = conf.int
  ), survfit.info %||% list())

  axis.info <- modifyList(list(
    type.y              = type.y,
    label.x             = label.x,
    label.y             = label.y,
    level.strata        = level.strata,
    label.strata        = label.strata,
    order.strata        = order.strata,
    limits.x            = limits.x,
    limits.y            = limits.y,
    breaks.x            = breaks.x,
    breaks.y            = breaks.y,
    use_coord_cartesian = use_coord_cartesian
  ), axis.info %||% list())

  visual.info <- modifyList(list(
    addConfidenceInterval        = addConfidenceInterval,
    addRiskTable                 = addRiskTable,
    symbol.risktable             = symbol.risktable,
    addEstimateTable             = addEstimateTable,
    symbol.estimatetable         = symbol.estimatetable,
    addCensorMark                = addCensorMark,
    shape.censor.mark            = shape.censor.mark,
    size.censor.mark             = size.censor.mark,
    addCompetingRiskMark         = addCompetingRiskMark,
    competing.risk.time          = competing.risk.time,
    shape.competing.risk.mark    = shape.competing.risk.mark,
    size.competing.risk.mark     = size.competing.risk.mark,
    addIntercurrentEventMark     = addIntercurrentEventMark,
    intercurrent.event.time      = intercurrent.event.time,
    shape.intercurrent.event.mark= shape.intercurrent.event.mark,
    size.intercurrent.event.mark = size.intercurrent.event.mark,
    addQuantileLine              = addQuantileLine,
    quantile                     = quantile
  ), visual.info %||% list())

  panel.info <- modifyList(list(
    printEachEvent     = printEachEvent,
    printEachVar       = printEachVar,
    rows.columns.panel = rows.columns.panel
  ), panel.info %||% list())

  style.info <- modifyList(list(
    style           = style,
    palette         = palette,
    font.family     = font.family %||% "sans",
    font.size       = font.size   %||% 12,
    legend.position = legend.position
  ), style.info %||% list())

  ggsave.info <- modifyList(list(
    filename.ggsave = filename.ggsave,
    width.ggsave    = width.ggsave,
    height.ggsave   = height.ggsave,
    dpi.ggsave      = dpi.ggsave,
    units           = "in"
  ), ggsave.info %||% list())

  list(
    survfit.info = survfit.info,
    axis.info    = axis.info,
    visual.info  = visual.info,
    panel.info   = panel.info,
    style.info   = style.info,
    ggsave.info  = ggsave.info
  )
}



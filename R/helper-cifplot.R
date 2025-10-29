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
  p + style_theme
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

plot_make_label.strata.map <- function(fit, label.strata) {
  .sf_strata_names <- function(fit) {
    if (is.null(fit$strata)) return(NULL)
    nm <- names(fit$strata)
    if (is.null(nm)) return(NULL)
    nm
  }
  .canon_str <- function(x) sub("^.*?=", "", as.character(x))  # "strata=level" > "level"

  if (is.null(label.strata)) return(NULL)

  fit_names <- .sf_strata_names(fit)
  if (is.null(fit_names)) return(NULL)
  if (!is.null(names(label.strata)) && any(nzchar(names(label.strata)))) {
    key_in  <- .canon_str(names(label.strata))
    key_fit <- .canon_str(fit_names)
    idx <- match(key_fit, key_in)
    if (all(!is.na(idx))) {
      out <- unname(label.strata[idx])
      names(out) <- fit_names
      return(out)
    } else {
      ok <- which(!is.na(idx))
      if (length(ok) > 0L) {
        out <- unname(label.strata[idx[ok]])
        names(out) <- fit_names[ok]
        warning("Some label.strata names did not match strata and were ignored: ",
                paste(setdiff(key_in, key_fit), collapse = ", "), call. = FALSE)
        return(out)
      } else {
        warning("No names in label.strata matched strata; falling back to order.", call. = FALSE)
      }
    }
  }
  if (length(label.strata) == length(fit_names)) {
    out <- label.strata
    names(out) <- fit_names
    return(out)
  }
  warning(sprintf(
    "Length of label.strata (%d) does not match number of strata (%d); labels ignored.",
    length(label.strata), length(fit_names)
  ), call. = FALSE)
  return(NULL)
}

plot_survfit_strata_labels <- function(survfit_object, strata.label, recycle_ok = FALSE) {
  if (is.null(survfit_object)) return(survfit_object)

  nms <- names(survfit_object$strata)
  if (is.null(nms) || !length(nms)) {
    new_nms <- as.character(strata.label)
    if (length(new_nms) != length(survfit_object$strata)) {
      if (recycle_ok) new_nms <- rep_len(new_nms, length(survfit_object$strata))
      else stop("Length of `strata.label` does not match number of strata.", call. = FALSE)
    }
    names(survfit_object$strata) <- new_nms
    return(survfit_object)
  }

  lhs <- sub("^\\s*([^=]+)\\s*=.*$", "\\1", nms, perl = TRUE)
  rhs <- sub("^.*?=\\s*", "", nms, perl = TRUE)

  lab <- strata.label

  if (!is.null(names(lab)) && length(intersect(names(lab), nms)) == length(nms)) {
    new_rhs <- unname(lab[nms])

  } else if (!is.null(names(lab)) && length(intersect(names(lab), rhs)) == length(nms)) {
    new_rhs <- unname(lab[rhs])
  } else {
    lab <- as.character(lab)
    if (length(lab) != length(nms)) {
      if (recycle_ok) lab <- rep_len(lab, length(nms))
      else stop("Length of `strata.label` does not match number of strata.", call. = FALSE)
    }
    new_rhs <- lab
  }
  #new_nms <- paste0(lhs, "=", new_rhs)
  #names(survfit_object$strata) <- new_nms
  names(survfit_object$strata) <- new_rhs
  return(survfit_object)
}

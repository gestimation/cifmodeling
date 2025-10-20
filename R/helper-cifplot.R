style_theme_classic <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 4, family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 4, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

style_theme_bold <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

style_theme_framed <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(fill = "grey90", color = "black"),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 0.8)
    )
}

style_theme_monochrome <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 1, face = "bold"),
      axis.text         = ggplot2::element_text(size = font.size, color = "grey20"),
      legend.text       = ggplot2::element_text(size = font.size - 1, color = "grey20"),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

scale_monochrome <- function(n_strata = 6) {
  ltys_all   <- c("dashed","solid","dotted","longdash","dotdash","twodash","dashed","solid","dotted","longdash","dotdash","twodash")
  shapes_all <- c(16, 1, 3, 4, 15, 17, 16, 1, 3, 4, 15, 17)
  n_use <- min(n_strata, length(ltys_all), length(shapes_all))
  list(
    ggplot2::scale_color_manual(values = rep("black", n_use)),
    ggplot2::scale_fill_manual(values = gray(seq(0.85, 0.30, length.out = n_use))),
    ggplot2::scale_linetype_manual(values = ltys_all[seq_len(n_use)]),
    ggplot2::scale_shape_manual(values = shapes_all[seq_len(n_use)])
  )
}

applyStyle <- function(
    p,
    style = c("CLASSIC", "BOLD", "FRAMED", "MONOCHROME"),
    font.family = "sans",
    font.size = 14,
    legend.position = "top",
    n_strata = 6
) {
  style <- match.arg(style)
  style_theme <- switch(
    style,
    CLASSIC    = style_theme_classic(font.family, font.size, legend.position),
    BOLD       = style_theme_bold(font.family, font.size, legend.position),
    FRAMED     = style_theme_framed(font.family, font.size, legend.position),
    MONOCHROME = style_theme_monochrome(font.family, font.size, legend.position)
  )
  p <- p + style_theme
  if (identical(style, "MONOCHROME")) {
    p <- p + scale_monochrome(n_strata = n_strata)
  }
  return(p)
}

drawMarks <- function(p, survfit_object, marks, type.y, shape, size) {
  time <- y <- strata <- NULL
  if (is.null(marks) || !length(marks)) return(p)
  mark_df <- makeMarkDataFrame(survfit_object, marks, type.y, extend = TRUE)
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

alignMarkKeys <- function(fit, marks) {
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

makeMarkDataFrame <- function(fit, marks, type.y, extend = TRUE) {
  if (is.null(marks) || length(marks) == 0) return(NULL)
  out_alignMarkKeys <- alignMarkKeys(fit, marks)
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

theme_risktable_font <- function(
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

.plot_make_x_max <- function(sf) {
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

.plot_make_label.strata.map <- function(fit, label.strata) {
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




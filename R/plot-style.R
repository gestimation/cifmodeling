# --- helpers (internal) -------------------------------------------------------

# 既定色
.default_fallback_colors <- function(k) {
  hue_fn <- scales::hue_pal(h = c(15, 375), c = 100, l = 65, h.start = 0)
  hue_fn(k)
}

.validate_or_fix_color <- function(x) {
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

.resolve_palette_colors <- function(levels_final, palette, n, fallback_colors = NULL) {
  if (is.null(fallback_colors)) fallback_colors <- .default_fallback_colors(max(1, n))
  if (is.null(palette)) {
    cols <- rep_len(fallback_colors, max(1, n))
    return(list(values = unname(cols), supplied = FALSE))
  }

  stopifnot(is.character(palette))
  pal <- .validate_or_fix_color(palette)

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

#' Apply presentation style theme to CIF plots
#'
#' @param p A ggplot object produced by \code{call_ggsurvfit()}.
#' @param style Character scalar matching one of the supported styles.
#' @param font.family Base font family used in theme definitions.
#' @param font.size Base font size used in theme definitions.
#' @param legend.position Legend placement.
#'
#' @return The plot with the requested theme applied.
#' @keywords internal
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

#' Apply stratified scales for CIF plots in a single pass
#'
#' @param p A ggplot object.
#' @param style Plot style requested by the user.
#' @param palette Optional color palette supplied to `cifplot()`.
#' @param n_strata Number of strata present in the plotted object.
#' @param strata_levels_final Character vector of final strata levels after label/order adjustments.
#' @param strata_labels_final Character vector of final strata labels.
#'
#' @return A ggplot object with color, fill, linetype, and shape scales applied.
#' @keywords internal
apply_all_scales_once <- function(
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
  palette_info <- .resolve_palette_colors(lvls, palette, n_effective)
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

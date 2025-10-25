# --- helpers (internal) -------------------------------------------------------

# 既定色
.default_fallback_colors <- function(k) {
  # ggplot2 標準の hue パレットを使用して k 色生成
  hue_fn <- scales::hue_pal(h = c(15, 375), c = 100, l = 65, h.start = 0)
  hue_fn(k)
}

# 色名の検証（未知は警告して black にフォールバック）
.validate_or_fix_color <- function(x) {
  v <- tolower(grDevices::colors())
  out <- tolower(x)
  bad <- which(!out %in% v)
  if (length(bad)) {
    warning("Unknown color name: ", paste(unique(out[bad]), collapse = ", "),
            " (defaulting to 'black')")
    out[bad] <- "black"
  }
  out
}

# palette 入力（色名のみ）を最終凡例ラベル順の色ベクトルに解決
# levels_final: label.strata 適用後の最終レベル
# palette: NULL | 無名/名前付き 文字ベクトル（色名のみ）
.resolve_colors_from_palette <- function(levels_final, palette, fallback_colors = NULL) {
  k <- length(levels_final)
  if (k == 0) return(character(0))
  if (is.null(fallback_colors)) fallback_colors <- .default_fallback_colors(k)

  if (is.null(palette)) {
    cols <- rep_len(fallback_colors, k)
    cols <- unname(cols)
    names(cols) <- NULL
  }

  stopifnot(is.character(palette))

  if (is.null(names(palette))) {
    pal <- .validate_or_fix_color(palette)
    cols <- rep_len(pal, k)
    names(cols) <- levels_final
    return(cols)
  } else {
    pal <- .validate_or_fix_color(palette)
    cols <- vapply(levels_final, function(s) pal[[s]] %||% NA_character_, character(1))
    if (anyNA(cols)) {
      miss <- which(is.na(cols))
      cols[miss] <- rep_len(fallback_colors, k)[miss]
    }
    cols <- unname(cols)
    names(cols) <- NULL
  }
}


# 線種の自動判定：
# - 全色が同一（unique(colors)==1）なら ltys_all を順に割当
# - それ以外は "solid" に統一
.resolve_linetypes_auto <- function(colors, ltys_all) {
  k <- length(colors)
  if (k == 0) return(character(0))
  stopifnot(is.character(colors), length(colors) == k)

  if (length(unique(colors)) == 1L) {
    lt <- rep_len(ltys_all, k)
  } else {
    lt <- rep_len("solid", k)
  }
  names(lt) <- names(colors)
  lt
}

plot_apply_style <- function(
    p,
    style = c("CLASSIC", "BOLD", "FRAMED", "MONOCHROME"),
    font.family = "sans",
    font.size = 14,
    legend.position = "top",
    n_strata = 6,
    palette_colors = NULL,
    strata_levels_final = NULL,
    strata_labels_final = NULL
) {
  print(strata_levels_final)
  print(strata_labels_final)
  style <- match.arg(style)
  style_theme <- switch(
    style,
    CLASSIC    = plot_style_classic(font.family, font.size, legend.position),
    BOLD       = plot_style_bold(font.family, font.size, legend.position),
    FRAMED     = plot_style_framed(font.family, font.size, legend.position),
    MONOCHROME = plot_style_monochrome(font.family, font.size, legend.position)
  )
  p <- p + style_theme

  # MONOCHROME は既存のスケール関数を使用（ここでは変更しない）
  if (identical(style, "MONOCHROME")) {
    p <- p + plot_scale_monochrome(n_strata = n_strata)
    return(p)
  }

  # デフォルト維持：palette_colors が NULL のときはスケールを追加しない
  if (is.null(palette_colors) || is.null(strata_levels_final)) {
    return(p)
  }

# plot_apply_style() 内のスケール付与部分を差し替え
cols <- .resolve_colors_from_palette(strata_levels_final, palette_colors)

ltys_all <- c("dashed","solid","dotted","longdash","dotdash","twodash",
              "dashed","solid","dotted","longdash","dotdash","twodash",
              "solid","dotted","longdash","dotdash","twodash")
lts <- .resolve_linetypes_auto(cols, ltys_all)

breaks <- strata_levels_final
labels <- strata_labels_final  # NULL でもOK（その場合は既定表示）

p +
  ggplot2::scale_color_manual(
    values = unname(cols),
    drop = FALSE, guide = "legend"
  ) +
  ggplot2::scale_linetype_manual(
    values = unname(lts),  breaks = breaks, labels = labels,
    drop = FALSE, guide = "legend"
  ) +
  ggplot2::scale_fill_manual(
    values = unname(cols), breaks = breaks, labels = labels,
    drop = FALSE, guide = "legend"
  )
}

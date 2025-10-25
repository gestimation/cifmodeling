#' Internal helpers for color/linetype resolution
#' @keywords internal

# 既定パレット（Okabe-Ito が使えればそれを優先）
.default_fallback_colors <- function(k) {
  if ("Okabe-Ito" %in% grDevices::palette.pals()) {
    pal <- grDevices::palette.colors(palette = "Okabe-Ito")
    return(rep_len(unname(pal), k))
  } else {
    return(scales::hue_pal()(k))
  }
}

# 有効な色名チェック（未知の色名は警告して black にフォールバック）
.validate_or_fix_color <- function(x) {
  v <- tolower(grDevices::colors())
  out <- tolower(x)
  bad <- which(!out %in% v)
  if (length(bad)) {
    warning(
      "Unknown color name: ", paste(unique(out[bad]), collapse = ", "),
      " (defaulting to 'black')"
    )
    out[bad] <- "black"
  }
  out
}

# "red", "blue-dot", "green-broken", "dot" を色＋線種に展開
# palette: 文字ベクトル（名前付き/無名どちらも可）
# levels_final: label.strata 適用後の最終凡例ラベル（順序確定後）
# 返り値: list(color=..., linetype=...) with names(levels_final)
parse_colorlinetype_patterns <- function(palette, levels_final, fallback_colors = NULL) {
  if (is.null(palette)) return(NULL)
  stopifnot(is.character(palette))
  if (is.null(levels_final)) return(NULL)

  k <- length(levels_final)
  if (k == 0L) return(NULL)

  if (is.null(fallback_colors)) fallback_colors <- .default_fallback_colors(k)

  # 名前付き → 最終ラベル順に並べ替え、無名 → そのまま並び使用
  target <- if (is.null(names(palette))) {
    stats::setNames(rep_len(palette, k), levels_final)
  } else {
    stats::setNames(vapply(levels_final, function(s) palette[[s]] %||% NA_character_, character(1)), levels_final)
  }

  col <- character(k)
  lt  <- character(k)

  i <- 0L
  for (lab in names(target)) {
    i <- i + 1L
    spec <- tolower(target[[lab]])

    if (is.na(spec) || !nzchar(spec)) {
      # 未指定 → フォールバック
      col[i] <- fallback_colors[i]
      lt[i]  <- "solid"
      next
    }

    if (identical(spec, "dot")) {
      col[i] <- "black"
      lt[i]  <- "dotted"
      next
    }

    parts <- strsplit(spec, "-", fixed = TRUE)[[1]]
    base_col <- parts[1]
    pattern  <- if (length(parts) >= 2) parts[2] else "solid"

    base_col <- .validate_or_fix_color(base_col)

    lt[i] <- switch(
      pattern,
      dot    = "dotted",
      broken = "dashed",
      dash   = "longdash",
      solid  = "solid",
      fill   = "solid",
      # 未知パターンは solid にフォールバック
      "solid"
    )

    col[i] <- base_col
  }

  list(
    color = stats::setNames(col, names(target)),
    linetype = stats::setNames(lt, names(target))
  )
}

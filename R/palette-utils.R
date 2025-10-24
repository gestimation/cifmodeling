#' Internal palettes and helpers for color resolution
#'
#' @keywords internal
#' @importFrom grDevices colors palette.colors palette.pals
#' @importFrom scales hue_pal
.palette.cif <- list(
  basic  = c("red", "blue", "green", "orange", "purple"),
  pastel = c("lightblue", "pink", "palegreen", "gold", "lightgray"),
  bw     = c("black", "gray50")
)

# Return default fallback palette of length k
.default_fallback <- function(k) {
  pals_fun    <- get0("palette.pals",    envir = asNamespace("grDevices"), mode = "function")
  colors_fun  <- get0("palette.colors",  envir = asNamespace("grDevices"), mode = "function")
  okabe_ito   <- NULL
  if (!is.null(pals_fun) && !is.null(colors_fun)) {
    pals <- tryCatch(pals_fun(), error = function(...) NULL)
    if (!is.null(pals) && "Okabe-Ito" %in% pals) {
      okabe_ito <- tryCatch(colors_fun(palette = "Okabe-Ito"), error = function(...) NULL)
    }
  }
  if (!is.null(okabe_ito)) {
    return(rep_len(unname(okabe_ito), k))
  }
  scales::hue_pal()(k)
}

# Validate color names (character vector), warn unknown names
.validate_color_names <- function(cols) {
  if (is.null(cols)) return(invisible(NULL))
  stopifnot(is.character(cols))
  valid <- tolower(grDevices::colors())
  bad <- setdiff(tolower(unique(cols)), valid)
  if (length(bad)) {
    warning("Unknown color names: ", paste(bad, collapse = ", "))
  }
  invisible(NULL)
}

# Resolve user 'palette' argument to a named vector aligned to final strata levels.
# - levels_final: levels (after label.strata applied)
# - palette: NULL | scalar preset name | unnamed vector | named vector (by final labels)
# - fallback: character vector used for missing colors
resolve_palette <- function(levels_final, palette, fallback = NULL) {
  if (is.null(levels_final)) stop("levels_final must not be NULL")
  k <- length(levels_final)
  if (is.null(fallback)) fallback <- get_default_palette(k)

  # 1) preset string (length=1, matches .palette.cif)
  if (is.character(palette) && length(palette) == 1L && palette %in% names(.palette.cif)) {
    palette <- .palette.cif[[palette]]
  }

  # 2) NULL -> fallback
  if (is.null(palette)) {
    pal <- rep_len(fallback, k)
    names(pal) <- levels_final
    return(pal)
  }

  stopifnot(is.character(palette))

  # 3) named vector -> map by final labels; missing filled by fallback
  if (!is.null(names(palette)) && any(nzchar(names(palette)))) {
    .validate_color_names(palette)
    pal <- vapply(levels_final, function(s) {
      if (s %in% names(palette)) palette[[s]] else NA_character_
    }, character(1))
    if (anyNA(pal)) {
      miss <- which(is.na(pal))
      fb   <- rep_len(fallback, k)[miss]
      pal[miss] <- fb
    }
    names(pal) <- levels_final
    return(pal)
  }

  # 4) unnamed vector -> assign in order
  .validate_color_names(palette)
  pal <- rep_len(palette, k)
  names(pal) <- levels_final
  pal
}

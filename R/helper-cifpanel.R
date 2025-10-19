.panel_as_formula_global <- function(f) {
  if (is.character(f))   return(stats::as.formula(f, env = .GlobalEnv))
  if (inherits(f, "formula")) return(stats::as.formula(f, env = .GlobalEnv))
  .err("formula_must_be")
}

.panel_recycle_to <- function(x, n) {
  if (length(x) == n) return(x)
  if (length(x) == 0L) .err("recycle_empty")
  rep_len(x, n)
}
.panel_to_list <- function(x) if (is.null(x)) NULL else if (is.list(x)) x else as.list(x)
.panel_drop_nulls <- function(lst) lst[!vapply(lst, is.null, logical(1))]
.panel_strip_overrides_from_dots <- function(dots, override_names) {
  if (length(override_names) == 0L) return(dots)
  dots[setdiff(names(dots), override_names)]
}
.panel_is_surv <- function(x) {
  x <- toupper(as.character(x %||% ""))
  x %in% c("S", "SURVIVAL")
}
.panel_is_comp <- function(x) {
  x <- toupper(as.character(x %||% ""))
  x %in% c("C", "COMPETING-RISK", "COMPETING_RISK", "COMPETINGRISK")
}
.panel_norm_outcome <- function(x) {
  if (.panel_is_surv(x)) return("S")
  if (.panel_is_comp(x)) return("C")
  stop("Unknown outcome.type: ", x, " (use 'S'/'SURVIVAL' or 'C'/'COMPETING-RISK').")
}
.panel_validate_code_events <- function(code_events_list, outcome_flags) {
  stopifnot(length(code_events_list) == length(outcome_flags))
  for (i in seq_along(code_events_list)) {
    pair <- code_events_list[[i]]
    if (outcome_flags[i] == "S") {
      if (!(is.numeric(pair) && length(pair) == 2L))
        .err("code_events_len_surv", i = i)
    } else {
      if (!(is.numeric(pair) && length(pair) == 3L))
        .err("code_events_len_cr", i = i)
    }
  }
  invisible(TRUE)
}

.panel_extract_fonts <- function(dots) {
  list(
    family = dots$font.family %||% "sans",
    size   = dots$font.size   %||% 12
  )
}
.panel_build_theme <- function(font.family = "sans", font.size = 12) {
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

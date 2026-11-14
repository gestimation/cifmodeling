# Display metadata are keyed by the original values and never used to recode data.
.cif_variable_label <- function(x, variable, explicit = NULL) {
  if (.cif_has_label(explicit)) return(explicit)
  label <- attr(x, "label", exact = TRUE)
  if (.cif_has_label(label)) label else variable
}

.cif_has_label <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.cif_value_labels <- function(x) {
  labels <- attr(x, "labels", exact = TRUE)
  if (!length(labels) || is.null(names(labels))) return(character())
  keep <- !is.na(labels) & !is.na(names(labels)) & nzchar(names(labels))
  stats::setNames(names(labels)[keep], as.character(labels[keep]))
}

.cif_display_levels <- function(keys, values = NULL, explicit = NULL, fallback = keys) {
  out <- stats::setNames(fallback, keys)
  for (map in list(values, explicit)) {
    if (is.null(names(map)) && length(map) == length(keys)) names(map) <- keys
    matched <- match(keys, names(map))
    labels <- if (length(map)) as.character(map)[matched] else rep(NA_character_, length(keys))
    use <- !is.na(labels) & nzchar(labels)
    out[use] <- labels[use]
  }
  out
}

.cif_curve_display <- function(formula, data, terms, model, strata) {
  lhs <- if (length(formula) == 3L && is.call(formula[[2L]])) formula[[2L]]
  time <- if (length(lhs) >= 2L && is.symbol(lhs[[2L]])) as.character(lhs[[2L]])
  event <- if (length(lhs) >= 3L && is.symbol(lhs[[3L]])) as.character(lhs[[3L]])
  variables <- attr(terms, "term.labels")
  if (!length(variables)) variables <- character()
  codes <- levels(strata)
  labels <- NULL
  if (length(codes) && length(variables) && all(variables %in% names(model))) {
    maps <- stats::setNames(lapply(variables, function(variable) {
      .cif_value_labels(data[[variable]] %||% model[[variable]])
    }), variables)
    if (any(lengths(maps))) {
      rows <- match(codes, as.character(strata))
      parts <- lapply(variables, function(variable) {
        key <- as.character(model[[variable]][rows])
        unname(.cif_display_levels(key, maps[[variable]]))
      })
      labels <- stats::setNames(do.call(paste, c(parts, sep = " / ")), codes)
    }
  }
  list(time = if (!is.null(time)) .cif_variable_label(data[[time]], time),
       event = if (!is.null(event)) .cif_variable_label(data[[event]], event),
       event.levels = if (!is.null(event)) .cif_value_labels(data[[event]]),
       variable = if (length(variables)) paste(vapply(variables, function(variable) {
         .cif_variable_label(data[[variable]], variable)
       }, character(1L)), collapse = " / "),
       levels = labels, codes = codes)
}

.cif_plot_labels <- function(fit, explicit, expected = NULL) {
  keys <- names(fit$strata)
  values <- fit$reporting$display$levels
  if (!length(keys) || (is.null(explicit) && is.null(values))) return(NULL)
  short <- unname(fit$reporting$display$codes[keys] %||% sub("^[^=]*=", "", keys))
  short[is.na(short)] <- keys[is.na(short)]
  .cif_display_levels(keys, values, .cif_strata_overrides(keys, short, explicit, expected),
                      fallback = short)
}

.cif_strata_overrides <- function(keys, short, explicit, expected = NULL) {
  if (!length(explicit)) return(NULL)
  if (is.null(names(explicit))) {
    expected <- expected %||% keys
    if (length(expected) != length(explicit))
      stop("`label.strata` length must match the number of expected strata.")
    names(explicit) <- as.character(expected)
  }
  index <- match(names(explicit), keys)
  missing <- is.na(index)
  index[missing] <- match(names(explicit)[missing], short)
  stats::setNames(explicit[!is.na(index)], keys[index[!is.na(index)]])
}

.cif_scale_labeller <- function(mapping) {
  force(mapping)
  function(breaks) unname(.cif_display_levels(breaks, mapping))
}

.cif_test_display <- function(input) {
  data <- input$data.original
  variable <- input$exposure
  codes <- levels(factor(input$data[[variable]]))
  list(variable = .cif_variable_label(data[[variable]], variable),
       levels = .cif_display_levels(codes, .cif_value_labels(data[[variable]])))
}

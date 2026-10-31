# Reporting information describes exposure contrasts, not nuisance regressions.
.polyreg_reporting_metadata <- function(x) {
  levels <- as.character(x$estimand$exposure.labels)
  if (length(levels) < 2L) return(NULL)
  k <- length(levels) - 1L
  events <- if (x$outcome.type %in% c("competing-risk", "proportional-competing-risk"))
    c("event1", "event2") else "event1"
  value <- x$diagnostics[[x$exposure]]
  label <- attr(value, "label")
  parts <- lapply(events, function(event) {
    part <- if (length(events) == 2L) x$summary[[event]] else x$summary[[1L]]
    rows <- utils::tail(seq_len(nrow(part$tidy)), k)
    if (length(rows) != k) stop("Exposure reporting metadata do not match the fitted contrasts")
    list(rows = rows, term = as.character(part$tidy$term[rows]),
         coefficient.index = attr(part$tidy, "coefficient.index")[rows],
         effect.measure = x$estimand[[if (event == "event1") "effect.measure1" else "effect.measure2"]])
  })
  names(parts) <- events
  list(exposure = x$exposure, levels = levels,
       label = if (is.null(label)) x$exposure else as.character(label), events = parts)
}

.polyreg_reporting_info <- function(model) {
  meta <- model$reporting %||% .polyreg_reporting_metadata(model)
  if (is.null(meta)) stop("Exposure reporting metadata are unavailable for this polyreg object")
  event <- attr(model, "cifmodeling.reporting.event") %||% "event1"
  part <- meta$events[[event]]
  if (is.null(part)) stop("No competing event in this polyreg model.", call. = FALSE)
  c(meta[c("exposure", "levels", "label")], part, list(event = event))
}

#' @rdname polyreg-methods
#' @exportS3Method gtsummary::tbl_regression
tbl_regression.polyreg <- function(x, ..., event = c("event1", "event2")) {
  if (utils::packageVersion("gtsummary") < "2.6.1")
    stop("Reporting polyreg objects requires gtsummary >= 2.6.1.", call. = FALSE)
  if (identical(event, "both"))
    stop("Report one event at a time, then combine tables with gtsummary::tbl_merge().", call. = FALSE)
  event <- match.arg(event)
  x$reporting <- x$reporting %||% .polyreg_reporting_metadata(x)
  attr(x, "cifmodeling.reporting.event") <- event
  info <- .polyreg_reporting_info(x)
  # Filter a local reporting copy. The original fit and default tidy output
  # retain all the coefficients requested when the model was fitted.
  key <- if (length(x$reporting$events) == 1L) 1L else event
  part <- x$summary[[key]]
  part$tidy <- part$tidy[info$rows, , drop = FALSE]
  attr(part$tidy, "coefficient.index") <- info$coefficient.index
  x$summary[[key]] <- part
  x$reporting$events[[event]]$rows <- seq_len(nrow(part$tidy))
  args <- list(...)
  if (!"conf.level" %in% names(args)) args$conf.level <- x$conf.int
  default <- utils::getS3method("tbl_regression", "default", envir = asNamespace("gtsummary"))
  out <- do.call(default, c(list(x = x), args, list(event = event)))
  exponentiate <- isTRUE(out$inputs$exponentiate)
  label <- if (exponentiate) info$effect.measure else paste0("log(", info$effect.measure, ")")
  out$table_body$coefficients_label <- label
  out <- gtsummary::modify_header(out, estimate = paste0("**", label, "**"))
  horizon <- if (x$outcome.type %in% c("proportional-survival", "proportional-competing-risk"))
    "over time" else paste("at time", paste(x$estimand$time.point, collapse = ", "))
  gtsummary::modify_footnote_header(out, columns = "estimate",
    footnote = paste(label, "for", sub("event", "event ", event), horizon))
}

# These broom.helpers methods supply explicit contrast mappings. No stats
# terms(), model.frame(), or model.matrix() method is manufactured for the
# nuisance model or the full estimating-equation coefficient vector.

#' @exportS3Method broom.helpers::model_identify_variables
model_identify_variables.polyreg <- function(model) {
  info <- .polyreg_reporting_info(model)
  data.frame(term = info$term, variable = info$exposure, var_class = "factor",
    var_type = if (length(info$levels) == 2L) "dichotomous" else "categorical",
    var_nlevels = length(info$levels), stringsAsFactors = FALSE)
}

#' @exportS3Method broom.helpers::model_list_variables
model_list_variables.polyreg <- function(model, labels = NULL, only_variable = FALSE,
                                         add_var_type = FALSE, instrumental_suffix = " (instrumental)") {
  info <- .polyreg_reporting_info(model)
  if (only_variable) return(info$exposure)
  label <- if (info$exposure %in% names(labels)) as.character(labels[[info$exposure]]) else info$label
  out <- data.frame(variable = info$exposure, var_class = "factor",
                    label_attr = info$label, var_label = label, stringsAsFactors = FALSE)
  if (add_var_type) {
    out$var_type <- if (length(info$levels) == 2L) "dichotomous" else "categorical"
    out$var_nlevels <- length(info$levels)
  }
  out
}

#' @exportS3Method broom.helpers::model_list_contrasts
model_list_contrasts.polyreg <- function(model) {
  info <- .polyreg_reporting_info(model)
  data.frame(variable = info$exposure, contrasts = "contr.treatment",
             contrasts_type = "treatment", reference = 1L, stringsAsFactors = FALSE)
}

#' @exportS3Method broom.helpers::model_list_terms_levels
model_list_terms_levels.polyreg <- function(model, label_pattern = "{level}",
                                            variable_labels = NULL, sdif_term_level = c("diff", "ratio")) {
  info <- .polyreg_reporting_info(model)
  label <- model_list_variables.polyreg(model, labels = variable_labels)$var_label
  out <- data.frame(variable = info$exposure, contrasts_type = "treatment",
    term = c(paste0(info$exposure, ", ", info$levels[1L], " vs ", info$levels[1L]), info$term),
    level = info$levels, level_rank = seq_along(info$levels),
    reference = seq_along(info$levels) == 1L, reference_level = info$levels[1L],
    var_label = label, var_nlevels = length(info$levels),
    dichotomous = length(info$levels) == 2L, stringsAsFactors = FALSE)
  out$label <- as.character(glue::glue_data(out, label_pattern))
  out
}

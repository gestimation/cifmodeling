# Shared, model-independent output contracts. No fitting or optional table
# rendering package is invoked by these adapters.
output_check_dots <- function(...) {
  if (length(list(...))) stop("Unsupported additional arguments; use the explicit arguments")
}

output_model_tidy <- function(x, ...) generics::tidy(x, ...)
output_model_glance <- function(x, ...) generics::glance(x, ...)

output_check_level <- function(level, argument = "conf.int") {
  if (!is.numeric(level) || length(level) != 1L || !is.null(dim(level)) ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop(argument, " must be one number strictly between zero and one")
  invisible(NULL)
}

#' Extract a ggplot for customization or export
#'
#' @param x A ggplot or a single-plot cifplot object.
#' @param ... No additional arguments are supported.
#' @return A ggplot. This function does not estimate or predict anything.
#' @export
as_ggplot <- function(x, ...) UseMethod("as_ggplot")

#' @export
as_ggplot.default <- function(x, ...) {
  output_check_dots(...)
  if (!inherits(x, "ggplot")) stop("Supply a ggplot or a single-plot cifplot object")
  x
}

#' @export
as_ggplot.cifplot <- function(x, ...) {
  output_check_dots(...)
  if (is.null(x$plot)) stop("as_ggplot() requires a single plot")
  as_ggplot(x$plot)
}

#' Extract prediction data in the common long format
#'
#' @param x Prediction data or a supported fitted model.
#' @param ... Model-specific arguments. GAM/CIF methods take newdata and
#'   inference settings. polyreg extracts stored training-cohort potential CIFs.
#' @return A data frame with .row, estimate, optional category, and optional
#'   conf.low/conf.high columns. No newdata prediction is added to polyreg.
#' @export
as_prediction_data <- function(x, ...) UseMethod("as_prediction_data")

#' @export
as_prediction_data.default <- function(x, ...) {
  output_check_dots(...)
  if (!is.data.frame(x) || !all(c(".row", "estimate") %in% names(x)))
    stop("Supply prediction data with .row and estimate, or a supported fitted model")
  x
}

#' @rdname as_prediction_data
#' @param exposure.level For polyreg, one fitted exposure label; defaults to
#'   the reference. Values are counterfactual predictions for this level.
#' @param event For polyreg, event1, event2, or both (the default).
#' @export
as_prediction_data.polyreg <- function(x, exposure.level = x$estimand$code.exposure.ref,
                                       event = c("both", "event1", "event2"), ...) {
  output_check_dots(...)
  event <- match.arg(event)
  labels <- x$estimand$exposure.labels
  if (is.null(labels) || is.null(x$predictors))
    stop("Prediction labels are unavailable; refit this legacy polyreg model")
  if (length(exposure.level) != 1L || is.na(exposure.level) || !exposure.level %in% labels)
    stop("exposure.level must be one fitted exposure label")
  values <- x$diagnostics$potential.CIFs
  if (!is.matrix(values) && !is.data.frame(values))
    stop("Stored potential CIFs are unavailable for this model")
  competing <- identical(x$outcome.type, "competing-risk")
  events <- if (event == "both") if (competing) c("event1", "event2") else "event1" else event
  if (!competing && "event2" %in% events) stop("No competing event in this polyreg model")
  k <- length(labels)
  if (ncol(values) != k * (if (competing) 2L else 1L) || nrow(values) != nrow(x$diagnostics))
    stop("Stored potential CIFs do not match their labels or analysis rows")
  d <- x$diagnostics
  d[c("potential.CIFs", "influence.function", "ip.weight")] <- NULL
  if (any(names(d) %in% c(".row", "estimate", "category", "time.point", "conf.low", "conf.high")))
    stop("Training data contains reserved prediction column names")
  label <- as.character(exposure.level)
  original <- d[[x$exposure]]
  d[[x$exposure]] <- if (is.factor(original)) factor(rep(label, nrow(d)), levels = levels(original)) else
    if (is.numeric(original)) rep(as.numeric(label), nrow(d)) else rep(label, nrow(d))
  j <- match(label, labels)
  out <- do.call(rbind, lapply(events, function(ev) {
    column <- j + if (ev == "event2") k else 0L
    cbind(d, .row = seq_len(nrow(d)), category = ev,
      time.point = x$estimand$time.point, estimate = as.numeric(values[, column]))
  }))
  rownames(out) <- NULL
  attr(out, "predictors") <- x$predictors
  attr(out, "interval") <- "none"
  attr(out, "exposure.level") <- label
  attr(out, "prediction.target") <- "training-cohort potential CIF"
  out
}

#' Convert regression results to modelsummary input
#'
#' @param x A supported fitted model.
#' @param ... Model-specific inference and coefficient-selection arguments.
#' @return A modelsummary_list, or a named list of these for multiple outcomes.
#'   Each contains tidy and glance data frames. Extraction does not require
#'   modelsummary. Inference settings are recorded in glance and metadata.
#' @export
as_modelsummary <- function(x, ...) UseMethod("as_modelsummary")

#' @export
as_modelsummary.default <- function(x, ...) {
  stop("No as_modelsummary() method for this model")
}

output_modelsummary <- function(tidy, glance, metadata = list()) {
  if (!is.data.frame(tidy) || !all(c("term", "estimate") %in% names(tidy)) ||
      !is.data.frame(glance) || nrow(glance) != 1L)
    stop("Model summaries require a coefficient table and one model-information row")
  for (nm in names(metadata)) {
    value <- metadata[[nm]]
    if (length(value) == 1L && is.atomic(value)) glance[[nm]] <- value
  }
  note <- if (is.null(metadata$conf.level)) character() else if
    (identical(metadata$interval.type, "none")) "Intervals unavailable." else paste0(
    format(100 * metadata$conf.level, trim = TRUE), "% ", metadata$interval.type,
    " intervals; variance: ", metadata$variance, "; scale: ", metadata$coefficient.scale, ".")
  structure(list(tidy = tidy, glance = glance), class = "modelsummary_list",
            inference = metadata, notes = note)
}

output_inference_metadata <- function(table) {
  info <- attr(table, "variance.info")
  list(interval.type = attr(table, "interval"), conf.level = attr(table, "conf.int"),
    variance = info$type, variance.source = info$source,
    interval.approximation = info$interval.approximation,
    smoothing.uncertainty = info$smoothing.uncertainty,
    censoring.uncertainty.included = info$censoring.uncertainty.included,
    coefficient.scale = "link")
}

output_coefficient_table <- function(beta, covariance, limits) {
  se <- sqrt(pmax(0, diag(covariance)))
  statistic <- unname(beta / se)
  posterior <- identical(attr(limits, "interval"), "credible")
  out <- data.frame(term = names(beta), estimate = unname(beta),
    std.error = unname(se), statistic = statistic,
    p.value = if (posterior) rep(NA_real_, length(beta)) else
      2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
    conf.low = as.numeric(limits[, 1L]), conf.high = as.numeric(limits[, 2L]),
    stringsAsFactors = FALSE)
  for (nm in c("interval", "conf.int", "variance.info")) attr(out, nm) <- attr(limits, nm)
  out
}

#' @rdname as_modelsummary
#' @param event For polyreg, event1, event2, or both (the default).
#' @param conf.int Interval level. For polyreg, defaults to its fitted level.
#' @param variance Covariance estimator: if or bootstrap for polyreg.
#'   Omitting it preserves the stored inference.
#' @export
as_modelsummary.polyreg <- function(x, event = c("both", "event1", "event2"),
                                    conf.int = x$conf.int, variance = NULL, ...) {
  output_check_dots(...)
  event <- match.arg(event)
  events <- if (event == "both") {
    if (x$outcome.type %in% c("competing-risk", "proportional-competing-risk"))
      c("event1", "event2") else "event1"
  } else event
  models <- stats::setNames(lapply(events, function(label) {
    table <- output_model_tidy(x, event = label, conf.int = conf.int, variance = variance)
    metadata <- output_inference_metadata(table)
    metadata$coefficient.scale <- "log effect"
    output_modelsummary(table, output_model_glance(x, event = label), metadata)
  }), events)
  if (length(models) == 1L) models[[1L]] else models
}

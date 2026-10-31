# Keep the native inference grid independent of survfit0()'s inserted rows.
# The component is copied by survival; it contains no original data or calls.
.cifcurve_reporting <- function(out, input, error, conf.type, level, conf.dist, logse) {
  counts <- if (length(out$strata)) as.integer(out$strata) else length(out$time)
  labels <- if (length(out$strata)) names(out$strata) else "All"
  if (is.null(labels)) labels <- paste0("strata=", seq_along(counts))
  groups <- sort(unique(as.integer(input$strata)))
  out$error <- error; out$conf.int <- level; out$conf.dist <- conf.dist
  out$logse <- isTRUE(logse)
  n <- length(out$time)
  out$df <- out$df %||% rep(Inf, n)
  out$df.reason <- out$df.reason %||% rep("normal", n)
  out$critical.value <- out$critical.value %||% rep(stats::qnorm((1 + level) / 2), n)
  se <- out$std.err
  curves <- vector("list", length(counts)); names(curves) <- labels
  offset <- 0L
  for (g in seq_along(counts)) {
    j <- offset + seq_len(counts[g]); offset <- offset + counts[g]
    rows <- which(as.integer(input$strata) == groups[g] & input$w > 0)
    support <- if (length(rows)) max(input$t[rows]) else -Inf
    status <- if (is.null(out$conf.status)) .cifcurve_ci_validity(
      out$surv[j], se[j], out$critical.value[j], out$df.reason[j], out$time[j] <= support)$status else out$conf.status[j]
    if (is.null(conf.type) || tolower(conf.type) %in% c("none", "n")) status[] <- "not_requested"
    curves[[g]] <- list(time = out$time[j], se = se[j], support = support,
      lower = if (length(out$lower)) out$lower[j] else rep(NA_real_, length(j)),
      upper = if (length(out$upper)) out$upper[j] else rep(NA_real_, length(j)),
      df = out$df[j], df.reason = out$df.reason[j],
      critical.value = out$critical.value[j], conf.status = status)
  }
  out$conf.status <- unlist(lapply(curves, `[[`, "conf.status"), use.names = FALSE)
  out$reporting <- list(curves = curves, conf.dist = conf.dist,
    has.ci = length(out$lower) > 0L, outcome.type = input$outcome.type,
    event.code = input$code.event1)
  if (!length(out$strata)) {
    # survfit0 drops NULL components, so no other field may partially match strata.
    out$strata <- out$strata.levels <- NULL
  }
  if (out$logse) out$std.err <- se / out$surv
  # These are the same derived quantities supplied by survfit0 for ordinary
  # survfit objects; the subclass must make their interpretation explicit.
  out$cumhaz <- -log(out$surv)
  out$std.chaz <- se / out$surv
  class(out) <- c("cifcurve", "survfit")
  out
}

.cifcurve_reporting_at <- function(info, time, strata = NULL) {
  labels <- if (is.null(strata)) rep(names(info$curves)[1L], length(time)) else as.character(strata)
  out <- list(se = numeric(length(time)), lower = rep(1, length(time)),
    upper = rep(1, length(time)), df = rep(NA_real_, length(time)),
    df.reason = rep("zero_contrast_information", length(time)),
    critical.value = rep(NA_real_, length(time)), conf.status = rep("endpoint_degenerate", length(time)))
  for (label in unique(labels)) {
    rows <- which(labels == label); curve <- info$curves[[label]]
    if (is.null(curve)) stop("Curve reporting metadata do not match the summary strata")
    index <- findInterval(time[rows], curve$time)
    has <- index > 0L
    for (nm in names(out)) out[[nm]][rows[has]] <- curve[[nm]][index[has]]
    if (info$conf.dist == "normal") {
      out$df[rows[!has]] <- Inf
      out$df.reason[rows[!has]] <- "normal"
      out$critical.value[rows[!has]] <- curve$critical.value[1L]
    }
    outside <- rows[time[rows] > curve$support]
    out$lower[outside] <- out$upper[outside] <- NA_real_
    out$conf.status[outside] <- "outside_observed_followup"
    if (!info$has.ci) out$conf.status[rows] <- "not_requested"
  }
  out
}

#' @rdname cifcurve
#' @name cifcurve
#' @export
summary.cifcurve <- function(object, ...) {
  info <- object$reporting
  args <- list(...)
  # Match the standard method's named and positional scale argument.
  matched <- match.call(utils::getS3method("summary", "survfit"),
    as.call(c(list(quote(summary), object = quote(object)), args)))
  scale <- matched$scale %||% 1
  out <- NextMethod("summary")
  if (!length(out$time)) return(out)
  # Reversing the division by scale can move an event time to its left.
  # Let survival select the same rows in native units instead; retain its
  # scaled output (including the restricted-mean table) unchanged.
  time <- if (scale == 1) out$time else NextMethod("summary", scale = 1)$time
  values <- .cifcurve_reporting_at(info, time, out$strata)
  out$std.err <- values$se
  if (info$has.ci) {
    out$lower <- values$lower; out$upper <- values$upper
  } else {
    # cardx expects these standard columns even for a point-estimate table.
    out$lower <- out$upper <- rep(NA_real_, length(out$time))
  }
  for (nm in c("df", "df.reason", "critical.value", "conf.status")) out[[nm]] <- values[[nm]]
  out
}

#' @export
`[.cifcurve` <- function(x, i, ...) {
  out <- NextMethod("[")
  selected <- stats::setNames(seq_along(x$reporting$curves), names(x$reporting$curves))
  if (!missing(i)) selected <- selected[i]
  info <- x$reporting
  info$curves <- info$curves[selected]
  out$reporting <- info
  for (nm in c("error", "conf.int", "conf.dist", "df.method", "fay.bound", "logse"))
    if (!is.null(x[[nm]])) out[[nm]] <- x[[nm]]
  if (!length(out$strata)) out$strata <- out$strata.levels <- NULL
  labels <- if (length(out$strata)) rep(names(out$strata), as.integer(out$strata)) else NULL
  values <- .cifcurve_reporting_at(info, out$time, labels)
  for (nm in c("df", "df.reason", "critical.value", "conf.status")) out[[nm]] <- values[[nm]]
  out
}

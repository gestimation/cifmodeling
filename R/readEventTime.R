
#' Extract per-stratum event times from (formula, data)
#'
#' @description
#' Build per-stratum event-time lists that you can pass to `competing.risk.time` or
#' `intercurrent.event.time` in `cifplot()` / `cifpanel()`.
#'
#' @param formula Event(time, status) ~ strata(...)
#' @param data Data frame.
#' @param which One of "event2", "event1", "censor", "code".
#' @param event.code Used when which="code": status value to pick (e.g., 3 for an intercurrent event).
#' @param code.event1,code.event2,code.censoring Integers describing your coding (used to build d1/d2/d0 inside readSurv).
#' @param subset.condition Optional character expression to subset `data`.
#' @param na.action Function to handle NA (default na.omit).
#' @param unique_times,drop_empty See `cif_event_times_from_frame()`.
#' @return Named list of numeric vectors (times per stratum).
#' @export
readEventTime <- function(
    formula, data,
    which_event = c("event2", "event1", "censor", "censoring", "user_specified"),
    code.event1 = 1, code.event2 = 2, code.censoring = 0, user_specified_code = NULL,
    subset.condition = NULL, na.action = stats::na.omit,
    unique_times = TRUE, drop_empty = TRUE
){
  which_event <- match.arg(which_event)
  out_readSurv <- readSurv(
    formula = formula, data = data, weights = NULL,
    code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
    subset.condition = subset.condition, na.action = na.action  # ← 修正：stats::na.action ではなく引数をそのまま
  )
  getEventTime(
    out_readSurv = out_readSurv,
    which_event = which_event, user_specified_code = user_specified_code,
    unique_times = unique_times, drop_empty = drop_empty
  )
}

getEventTime <- function(
    out_readSurv,
    which_event = c("event2", "event1", "censor", "censoring", "user_specified"),
    user_specified_code = NULL,
    unique_times = TRUE,
    drop_empty = TRUE
){
  which_event <- match.arg(which_event)

  if (is.null(out_readSurv) || !is.list(out_readSurv)) {
    .err("req", arg = "out_readSurv (list)")
  }

  strata <- out_readSurv$strata
  if (is.null(strata)) strata <- factor(rep("all", length(out_readSurv$t)))
  if (is.factor(strata)) strata <- as.character(strata)

  tvec    <- suppressWarnings(as.numeric(out_readSurv$t))
  epsilon <- suppressWarnings(as.numeric(out_readSurv$epsilon))
  if (anyNA(tvec))       .err("na", arg = "out_readSurv$t")
  if (any(tvec < 0))     .err("time_nonneg", arg = "out_readSurv$t")

  pick <- switch(
    which_event,
    event1 = as.integer(out_readSurv$d1),
    event2 = as.integer(out_readSurv$d2),
    censor = as.integer(out_readSurv$d0),
    censoring = as.integer(out_readSurv$d0),
    user_specified = {
      if (is.null(user_specified_code)) .err("req", arg = "user_specified_code")
      as.integer(epsilon == as.numeric(user_specified_code))
    }
  )

  labs <- unique(strata)
  out  <- stats::setNames(vector("list", length(labs)), labs)
  for (s in labs) {
    idx <- (strata == s) & (pick > 0L) & is.finite(tvec)
    tt  <- tvec[idx]
    if (isTRUE(unique_times)) tt <- sort(unique(tt)) else tt <- sort(tt)
    if (length(tt) > 0L || !isTRUE(drop_empty)) out[[s]] <- tt
  }
  if (isTRUE(drop_empty)) out <- out[vapply(out, length, integer(1)) > 0L]
  return(out)
}

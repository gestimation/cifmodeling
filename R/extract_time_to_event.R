#' Extract per-stratum event times from a formula and data
#'
#' @description
#' Creates a list of event times that can be passed to downstream
#' visualization or analysis functions such as `competing.risk.time` or
#' `intercurrent.event.time` in [cifplot()] and [cifpanel()].
#' Event types are specified by event 1, event 2, censoring, or user-specified codes.
#'
#' @param formula A model formula specifying the outcome and (optionally) \code{strata()}.
#' @param data A data frame containing variables in \code{formula}.
#' @param subset.condition Optional expression (as a character string) defining a
#'   subset of \code{data} to analyse. Defaults to \code{NULL}.
#' @param na.action Function to handle missing values (default: \code{na.omit} in \pkg{stats}).
#' @param which_event One of \code{"event1"}, \code{"event2"}, \code{"censor"},
#'   \code{"censoring"}, or \code{"user_specified"}, indicating which event type
#'   to extract times for.
#' @param code.event1,code.event2,code.censoring Integer codes representing the
#'   event and censoring categories. Defaults are \code{1}, \code{2}, and
#'   \code{0}, respectively.
#' @param user_specified_code When \code{which_event = "user_specified"},
#'   the integer event code to extract (e.g., 3 for an intercurrent event).
#' @param readUniqueTime Logical if \code{TRUE}, only unique and sorted time points
#'   are returned for each stratum.
#' @param dropEmpty Logical if \code{TRUE} (default), strata with no events are
#'   dropped from the returned list. Set to \code{FALSE} to retain empty strata
#'   as \code{numeric(0)} vectors (useful for diagnostics or consistent list length).
#'
#' @return
#' A named list of numeric vectors, where each element corresponds to a stratum
#' and contains the event times of the selected type.
#'
#' @details
#' This function is typically used internally by plotting and model functions,
#' but can also be called directly to inspect the per-stratum event-time
#' structure of a dataset.
#'
#' @examples
#' data(diabetes.complications)
#' output <- extract_time_to_event(Event(t,epsilon) ~ fruitq,
#'                                 data = diabetes.complications,
#'                                 which_event = "event2")
#' cifplot(Event(t,epsilon) ~ fruitq,
#'         data = diabetes.complications,
#'         outcome.type="COMPETING-RISK",
#'         addConfidenceInterval=FALSE,
#'         addRiskTable=FALSE,
#'         addCensorMark=FALSE,
#'         addCompetingRiskMark=TRUE,
#'         competing.risk.time=output,
#'         label.y='CIF of diabetic retinopathy',
#'         label.x='Years from registration')
#' @name extract_time_to_event
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifcurve()] for KM/AJ estimators; [cifplot()] for display of a CIF; [cifpanel()] for display of multiple CIFs; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
extract_time_to_event <- function(
    formula, data,
    subset.condition = NULL, na.action = na.omit,
    which_event = c("event2", "event1", "censor", "censoring", "user_specified"),
    code.event1 = 1, code.event2 = 2, code.censoring = 0, user_specified_code = NULL,
    readUniqueTime = TRUE, dropEmpty = TRUE
){
  which_event <- match.arg(which_event)
  out_readSurv <- curve_read_surv(
    formula = formula, data = data, weights = NULL,
    code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
    subset.condition = subset.condition, na.action = na.action
  )
  getTimeToEvent(
    out_readSurv = out_readSurv,
    which_event = which_event, user_specified_code = user_specified_code,
    readUniqueTime = readUniqueTime, dropEmpty = dropEmpty
  )
}

getTimeToEvent <- function(
    out_readSurv,
    which_event = c("event2", "event1", "censor", "censoring", "user_specified"),
    user_specified_code = NULL,
    readUniqueTime = TRUE,
    dropEmpty = TRUE
){
  which_event <- match.arg(which_event)

  if (is.null(out_readSurv) || !is.list(out_readSurv)) .err("req", arg = "out_readSurv (list)")

  strata <- out_readSurv$strata
  if (is.null(strata)) strata <- factor(rep("all", length(out_readSurv$t)))
  if (is.factor(strata)) strata <- as.character(strata)

  tvec    <- suppressWarnings(as.numeric(out_readSurv$t))
  epsilon <- suppressWarnings(as.numeric(out_readSurv$epsilon))
  if (anyNA(tvec))   .err("na",          arg = "out_readSurv$t")
  if (any(tvec < 0)) .err("nonneg",      arg = "out_readSurv$t")  # 既存の .msg$nonneg を使用

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

  # --- 単一ストラタなら numeric を返す ---
  if (length(labs) == 1L) {
    idx <- (strata == labs[[1L]]) & (pick > 0L) & is.finite(tvec)
    tt  <- tvec[idx]
    tt  <- if (isTRUE(readUniqueTime)) sort(unique(tt)) else sort(tt)
    if (length(tt) == 0L && isTRUE(dropEmpty)) return(numeric(0))
    return(tt)
  }

  # --- 複数ストラタなら list を返す（従来挙動） ---
  out  <- stats::setNames(vector("list", length(labs)), labs)
  for (s in labs) {
    idx <- (strata == s) & (pick > 0L) & is.finite(tvec)
    tt  <- tvec[idx]
    tt  <- if (isTRUE(readUniqueTime)) sort(unique(tt)) else sort(tt)
    if (length(tt) > 0 || !isTRUE(dropEmpty)) out[[s]] <- tt
  }
  if (isTRUE(dropEmpty)) out <- out[vapply(out, length, integer(1)) > 0]
  out
}

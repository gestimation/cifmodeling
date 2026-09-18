#' Internal Kaplan-Meier interface
#'
#' The numerical engine always retains fractional weighted counts. The public
#' package interface defaults to integer reporting for plotting compatibility;
#' internal estimation code can request the unrounded weighted totals.
#'
#' @param t Follow-up times.
#' @param d Event indicator.
#' @param w Optional non-negative case weights.
#' @param strata Optional integer stratum codes.
#' @param error Variance calculation method.
#' @param count.type Return weighted `n.risk`, `n.event`, and `n.censor` as
#'   upward-rounded integers (`"integer"`) or exact numeric totals
#'   (`"numeric"`).
#' @keywords internal
calculateKM <- function(
    t,
    d,
    w = numeric(),
    strata = integer(),
    error = "greenwood",
    count.type = c("integer", "numeric")
) {
  count.type <- match.arg(count.type)
  out <- calculateKM_engine(t, d, w, strata, error)
  if (identical(count.type, "integer")) {
    count_fields <- c("n.risk", "n.event", "n.censor")
    for (field in count_fields) {
      out[[field]] <- as.integer(ceiling(out[[field]]))
    }
  }
  out
}

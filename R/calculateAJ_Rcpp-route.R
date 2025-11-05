#' Route to Rcpp AJ engine with normalized vectors
#' Internal helper: expects normalized vectors already prepared upstream.
#' @keywords internal
calculateAJ_Rcpp_route <- function(t, epsilon, w = NULL, strata = NULL,
                                   error = "greenwood",
                                   conf.type = "arcsin",
                                   return_if = FALSE,
                                   conf.int = 0.95) {
  calculateAJ_Rcpp(
    t = t,
    epsilon = epsilon,
    w = w %||% numeric(),
    strata = strata %||% integer(),
    error = error,
    conf_type = conf.type,
    return_if = isTRUE(return_if),
    conf_int = conf.int
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

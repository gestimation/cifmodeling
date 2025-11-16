#' These arguments are shared by \code{cifplot()}, \code{cifpanel()}, and
#' \code{cifcurve()}.
#'
#' @name cif-stat-arguments
#' @keywords internal
#'
#' @param data A data frame containing variables in the formula.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative.
#' @param subset.condition Optional character string giving a logical condition to subset
#' \code{data} (default \code{NULL}).
#' @param na.action A function specifying the action to take on missing values (default \code{na.omit}).
#' @param outcome.type Character string specifying the type of time-to-event outcome.
#' One of \code{"survival"} (Kaplan–Meier) or \code{"competing-risk"} (Aalen–Johansen).
#' If \code{NULL} (default), the function automatically infers the outcome type from the data:
#' if the event variable has more than two unique levels, \code{"competing-risk"} is assumed;
#' otherwise, \code{"survival"} is used. You can also use abbreviations such as \code{"S"} or \code{"C"}.
#' Mixed or ambiguous inputs (e.g., \code{"c("S", "C")"} trigger automatic detection based on the event coding.
#' @param code.event1 Integer code of the event of interest (default \code{1}).
#' @param code.event2 Integer code of the competing event (default \code{2}).
#' @param code.censoring Integer code of censoring (default \code{0}).
#' @param error Character string specifying the method for SEs and CIs used internally.
#' For \code{"survival"}, choose one of \code{"greenwood"} (default), \code{"tsiatis"}, or \code{"if"}.
#' For \code{"competing-risk"}, choose one of \code{"delta"} (default), \code{"aalen"}, or \code{"if"}.
#' @param conf.type Character specifying the method of transformation for confidence intervals
#' used internally (default \code{arcsine-square root}).
#' @param conf.int Numeric two-sided confidence level (default \code{0.95}).
#'
NULL

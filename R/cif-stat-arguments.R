#' Common data / outcome arguments for CIF functions
#'
#' These arguments are shared by \code{cifplot()}, \code{cifpanel()}, and
#' \code{cifcurve()}.
#'
#' @param data A data frame containing variables in the formula.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative.
#' @param subset.condition Optional expression (as a character string) defining a
#'   subset of \code{data} to analyze (default \code{NULL}).
#' @param na.action A function specifying the action to take on missing values (default \code{na.omit}).
#' @param outcome.type Character string specifying the type of time-to-event outcome.
#'   One of \code{"SURVIVAL"} (Kaplan–Meier) or \code{"COMPETING-RISK"} (Aalen–Johansen).
#'   Abbreviations such as \code{"S"} or \code{"C"} are also accepted.
#'
#' @param outcome.type Character string specifying the type of time-to-event outcome.
#' One of \code{"SURVIVAL"} (Kaplan–Meier type) or \code{"COMPETING-RISK"} (Aalen–Johansen type).
#' If \code{NULL} (default), the function automatically infers the outcome type from the data:
#' if the event variable has more than two unique levels, \code{"COMPETING-RISK"} is assumed;
#' otherwise, \code{"SURVIVAL"} is used. You can also use abbreviations such as \code{"S"} or \code{"C"}.
#' Mixed or ambiguous inputs (e.g., \code{"c("S", "C")"} trigger automatic detection based on the event coding.
#' @param code.event1 Integer code of the event of interest (default \code{1}).
#' @param code.event2 Integer code of the competing event (default \code{2}).
#' @param code.censoring Integer code of censoring (default \code{0}).
#' @param error Character specifying the method for variance and standard error used internally.
#' @param conf.type Character specifying the method of transformation for confidence intervals used internally.
#' @param conf.int Numeric two-sided confidence level (default \code{0.95}).
#'
#' @name cif-stat-arguments
#' @keywords internal
NULL

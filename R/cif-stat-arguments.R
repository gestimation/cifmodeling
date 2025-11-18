#' These arguments are shared by [cifplot()], [cifpanel()], and
#' [cifcurve()].
#'
#' @name cif-stat-arguments
#' @keywords internal
#'
#' @param data A data frame containing variables in the formula.
#' @param weights Optional name of the weight variable in `data`. Weights must be nonnegative.
#' @param subset.condition Optional character string giving a logical condition to subset
#' `data` (default `NULL`).
#' @param na.action A function specifying the action to take on missing values (default `na.omit`).
#' @param outcome.type Character string specifying the type of time-to-event outcome.
#' One of `"survival"` (Kaplan-Meier) or `"competing-risk"` (Aalen-Johansen).
#' If `NULL` (default), the function automatically infers the outcome type from the data:
#' if the event variable has more than two unique levels, `"competing-risk"` is assumed;
#' otherwise, `"survival"` is used. You can also use abbreviations such as `"S"` or `"C"`.
#' Mixed or ambiguous inputs (e.g., `c("S", "C")`) trigger automatic detection based on the event coding.
#' @param code.event1 Integer code of the event of interest (default `1`).
#' @param code.event2 Integer code of the competing risk (default `2`).
#' @param code.censoring Integer code of censoring (default `0`).
#' @param error Character string specifying the method for SEs and CIs used internally.
#' For `"survival"`, choose one of `"greenwood"` (default), `"tsiatis"`, or `"if"`.
#' For `"competing-risk"`, choose one of `"delta"` (default), `"aalen"`, or `"if"`.
#' @param conf.type Character specifying the method of transformation for CIs
#' used internally (default `arcsine-square root`).
#' @param conf.int Numeric two-sided level of CIs (default `0.95`).
#'
NULL

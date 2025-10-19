#' @title Plot survival/CIF curves (ggsurvfit)
#' @description
#' Draw a publication-ready plot. Accepts a \code{survfit} object or a \code{formula+data}
#' (in which case it computes a \code{survfit} via \code{cifcurve()} first).
#' @return A ggplot object.
#' @export
cifplot <- function(
    x,
    data = NULL,
    # ----- display options (moved from cifcurve) -----
    type.y = NULL,
    label.x = "Time",
    label.y = NULL,
    label.strata = NULL,
    limits.x = NULL,
    limits.y = NULL,
    breaks.x = NULL,
    breaks.y = NULL,
    use_coord_cartesian = FALSE,
    addConfidenceInterval = TRUE,
    addRiskTable = TRUE,
    addEstimateTable = FALSE,
    addCensorMark = TRUE,
    shape.censor.mark = 3,
    size.censor.mark = 2,
    addCompetingRiskMark = FALSE,
    competing.risk.time = list(),
    shape.competing.risk.mark = 16,
    size.competing.risk.mark = 2,
    addIntercurrentEventMark = FALSE,
    intercurrent.event.time = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark = 2,
    addQuantileLine = FALSE,
    quantile = 0.5,
    style = "CLASSIC",
    font.family = "sans",
    font.size = 12,
    legend.position = "top",
    conf.type = NULL,
    filename = NULL,
    width = 6, height = 6, dpi = 300
) {
  # compute survfit if formula
  if (!inherits(x, "survfit")) {
    if (is.null(data)) stop("When `x` is a formula, `data` must be provided.")
    x <- cifcurve(x, data = data)  # estimation only
  }

  if (!requireNamespace("ggsurvfit", quietly = TRUE))
    stop("Package 'ggsurvfit' is required. Please install.packages('ggsurvfit').")

  p <- call_ggsurvfit(
    survfit_object                = x,
    out_readSurv                  = NULL,    # limitsは下でcoord等で与える
    conf.type                     = conf.type,
    addConfidenceInterval         = addConfidenceInterval,
    addRiskTable                  = addRiskTable,
    addEstimateTable              = addEstimateTable,
    addCensorMark                 = addCensorMark,
    shape.censor.mark             = shape.censor.mark,
    size.censor.mark              = size.censor.mark,
    addCompetingRiskMark          = addCompetingRiskMark,
    competing.risk.time           = competing.risk.time,
    shape.competing.risk.mark     = shape.competing.risk.mark,
    size.competing.risk.mark      = size.competing.risk.mark,
    addIntercurrentEventMark      = addIntercurrentEventMark,
    intercurrent.event.time       = intercurrent.event.time,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark,
    size.intercurrent.event.mark  = size.intercurrent.event.mark,
    addQuantileLine               = addQuantileLine,
    quantile                      = quantile,
    type.y                        = type.y,
    label.x                       = label.x,
    label.y                       = label.y,
    label.strata                  = label.strata,
    limits.x                      = limits.x,
    limits.y                      = limits.y,
    breaks.x                      = breaks.x,
    breaks.y                      = breaks.y,
    use_coord_cartesian           = use_coord_cartesian,
    style                         = style,
    font.family                   = font.family,
    font.size                     = font.size,
    legend.position               = legend.position
  )

  if (!is.null(filename)) ggplot2::ggsave(filename, plot = p, width = width, height = height, dpi = dpi)
  return(p)
}

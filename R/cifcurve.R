#' @title Calculate the Kaplan–Meier estimator and the Aalen–Johansen estimator with various SE and CI methods.
#' @description
#' Estimate and plot survival curves using the Kaplan–Meier estimator or
#' cumulative incidence curves under competing risks using the Aalen–Johansen estimator.
#' Returns a \code{survfit}-compatible object and, by default, draws a publication-ready plot via \pkg{ggsurvfit}.
#'
#' @details
#' \strong{Estimation:}
#' \itemize{
#'   \item \code{outcome.type = "SURVIVAL"}: Kaplan–Meier estimator with Greenwood-type variance.
#'   \item \code{outcome.type = "COMPETING-RISK"}: Aalen–Johansen estimator for CIF of \code{code.event1}
#'         using IPCW for the censoring distribution. The returned \code{surv} corresponds to \code{1 - CIF}.
#' }
#' \strong{Confidence intervals:}
#' Constructed on the probability scale with the specified \code{conf.type}.
#' If \code{conf.type \%in\% c("none","n")}, the plot suppresses CI bands.
#'
#' @param formula A model formula specifying the outcome and (optionally) \code{strata()}.
#' @param data A data frame containing variables in \code{formula}.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative; strictly positive is recommended.
#' @param subset.condition Optional character expression to subset \code{data} before analysis.
#' @param na.action Function to handle missing values (default: \code{\link[stats]{na.omit}}).
#' @param outcome.type \code{"SURVIVAL"} (KM) or \code{"COMPETING-RISK"} (AJ).
#' @param code.event1 Integer code of the event of interest (default \code{1}).
#' @param code.event2 Integer code of the competing event (default \code{2}).
#' @param code.censoring Integer code of censoring (default \code{0}).
#' @param error Character specifying variance type used internally. For \code{"SURVIVAL"} typically \code{"greenwood"}.
#'   For \code{"COMPETING-RISK"} pass options supported by \code{calculateAalenDeltaSE()} (\code{"aalen"}, \code{"delta"}, \code{"none"}).
#' @param conf.type Character transformation for CI on the probability scale (default \code{"arcsine-square root"}).
#' @param conf.int numeric two-sided confidence level (default \code{0.95}).
#' @param report.survfit.std.err If \code{TRUE}, report SE on the log-survival scale (survfit's convention). Otherwise SE is on the probability scale.

#' @returns A \code{survfit} object. For \code{outcome.type="SURVIVAL"}, \code{$surv} is the survival function.
#' For \code{outcome.type="COMPETING-RISK"}, \code{$surv} equals \code{1 - CIF} for \code{code.event1}.
#' Standard error and CIs are provided per \code{conf.type}. Note that some methods for \code{survfit} (e.g., \code{residuals.survfit}) may not be supported.
#'
#' @seealso \code{\link{cifplot}} for display of a plot; \code{\link{cifpanel}} for display of multiple plots; \pkg{ggsurvfit} for plotting helpers; \code{\link{polyreg}} for log-odds product models of CIFs.
#'
#' @examples
#' data(diabetes.complications)
#' survfit_by_group <- cifcurve(Event(t,epsilon) ~ fruitq, data = diabetes.complications,
#'                     outcome.type='COMPETING-RISK', error='delta')
#' cifplot(survfit_by_group, type.y = 'risk', label.y = 'CIF of diabetic retinopathy', label.x = 'Years from registration')
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib cifmodeling, .registration = TRUE
#' @export
cifcurve <- function(
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = c("SURVIVAL","COMPETING-RISK"),
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    error = NULL,
    conf.type = "arcsine-square root",
    conf.int = 0.95,
    report.survfit.std.err = FALSE
) {
  outcome.type  <- check_outcome.type(match.arg(outcome.type))
  out_readSurv  <- readSurv(formula, data, weights, code.event1, code.event2, code.censoring, subset.condition, na.action)
  error <- check_error(error, outcome.type)

  if (identical(outcome.type, "SURVIVAL")) {
    out_km <- calculateKM(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), error)
    out_km$std.err <- out_km$surv * out_km$std.err
    out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower = NULL)
    if (isTRUE(report.survfit.std.err)) out_km$std.err <- out_km$std.err / out_km$surv

    survfit_object <- list(
      time      = out_km$time,
      surv      = out_km$surv,
      n         = out_km$n,
      n.risk    = out_km$n.risk,
      n.event   = out_km$n.event,
      n.censor  = out_km$n.censor,
      std.err   = out_km$std.err,
      upper     = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower     = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type = conf.type,
      call      = match.call(),
      type      = "kaplan-meier",
      method    = "Kaplan-Meier"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) {
      names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
      survfit_object$strata <- out_km$strata
    }
    class(survfit_object) <- "survfit"

  } else {
    out_aj <- calculateAJ(out_readSurv)
    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))

    if (any(as.integer(out_readSurv$strata) != 1)) {
      n <- table(as.integer(out_readSurv$strata))
      rep_list <- mapply(rep, n, out_aj$strata1, SIMPLIFY = FALSE)
      n.risk <- do.call(c, rep_list) - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    } else {
      n <- length(out_readSurv$strata)
      n.risk <- n - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    }
    out_aj$std.err <- calculateAalenDeltaSE(out_aj$time1, out_aj$aj1, out_aj$n.event1, out_aj$n.event2,
                                            n.risk, out_aj$time0, out_aj$km0, out_aj$strata1, error)
    out_aj$surv <- 1 - out_aj$aj1
    out_ci <- calculateCI(out_aj, conf.int, conf.type, conf.lower = NULL)
    if (isTRUE(report.survfit.std.err)) out_aj$std.err <- out_aj$std.err / out_aj$surv

    survfit_object <- list(
      time        = out_aj$time1,
      surv        = out_aj$surv,
      n           = n,
      n.risk      = n.risk,
      n.event     = out_aj$n.event1,
      n.censor    = out_aj$n.censor,
      std.err     = out_aj$std.err,
      upper       = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower       = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type   = conf.type,
      call        = match.call(),
      type        = "aalen-johansen",
      method      = "aalen-johansen"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) survfit_object$strata <- out_aj$strata1
    class(survfit_object) <- "survfit"
  }
  return(survfit_object)
}

calculateAJ <- function(data) {
  out_km0 <- calculateKM(data$t, data$d0, data$w, as.integer(data$strata), "none")
  km0 <- get_surv(data$t, out_km0$surv, out_km0$time, as.integer(data$strata), out_km0$strata, out_km0$strata.levels)
  ip.weight <- (data$d0 == 0) * ifelse(km0 > 0, 1 / km0, 0)
  d1_ipw <- as.matrix(data$w * data$d1 * ip.weight)

  aj1 <- time1 <- integer(0)
  n.cum.event1 <- n.cum.event2 <- n.cum.censor <- numeric(0)
  n.event1 <- n.event2 <- n.censor <- numeric(0)
  strata1 <- integer(0)
  strata_vec <- as.integer(data$strata)

  for (level in sort(unique(strata_vec))) {
    idx <- which(strata_vec == level)

    sub_t  <- data$t[idx]
    sub_d0 <- data$d0[idx]
    sub_d1 <- data$d1[idx]
    sub_d2 <- data$d2[idx]
    sub_d1_ipw <- d1_ipw[idx, , drop = FALSE]

    o <- order(sub_t)
    sub_t  <- sub_t[o]
    sub_d0 <- sub_d0[o]
    sub_d1 <- sub_d1[o]
    sub_d2 <- sub_d2[o]
    sub_d1_ipw <- sub_d1_ipw[o, , drop = FALSE]

    not_atrisk <- outer(sub_t, sub_t, FUN = ">=")
    sub_aj1 <- as.vector(not_atrisk %*% sub_d1_ipw) / length(sub_t)
    sub_n.censor <- as.vector(not_atrisk %*% as.matrix(sub_d0))
    sub_n.event1 <- as.vector(not_atrisk %*% as.matrix(sub_d1))
    sub_n.event2 <- as.vector(not_atrisk %*% as.matrix(sub_d2))

    keep <- !duplicated(rev(sub_t))
    keep <- rev(keep)

    u_t  <- sub_t[keep]
    u_aj1 <- sub_aj1[keep]
    u_nc  <- sub_n.censor[keep]
    u_ne1 <- sub_n.event1[keep]
    u_ne2 <- sub_n.event2[keep]

    oo <- order(u_t)
    u_t  <- u_t[oo]
    u_aj1 <- u_aj1[oo]
    u_nc  <- u_nc[oo]
    u_ne1 <- u_ne1[oo]
    u_ne2 <- u_ne2[oo]

    inc_nc  <- c(u_nc[1],  diff(u_nc))
    inc_ne1 <- c(u_ne1[1], diff(u_ne1))
    inc_ne2 <- c(u_ne2[1], diff(u_ne2))

    time1 <- c(time1, u_t)
    aj1   <- c(aj1, u_aj1)

    n.cum.censor <- c(n.cum.censor, u_nc)
    n.cum.event1 <- c(n.cum.event1, u_ne1)
    n.cum.event2 <- c(n.cum.event2, u_ne2)

    n.censor <- c(n.censor, inc_nc)
    n.event1 <- c(n.event1, inc_ne1)
    n.event2 <- c(n.event2, inc_ne2)

    strata1 <- c(strata1, length(u_t))
  }

  list(
    time1 = time1,
    aj1 = aj1,
    n.event1 = n.event1,
    n.event2 = n.event2,
    n.censor = n.censor,
    n.cum.event1 = n.cum.event1,
    n.cum.event2 = n.cum.event2,
    n.cum.censor = n.cum.censor,
    strata1 = strata1,
    time0 = out_km0$time,
    km0 = out_km0$surv
  )
}

get_surv <- function(
    predicted.time,
    estimated.surv,
    estimated.time,
    predicted.strata = NULL,
    estimated.strata = NULL,
    strata.levels = NULL
){
  if (anyNA(predicted.time)) stop("Invalid predicted.time: contains NA.")
  if (length(estimated.surv) != length(estimated.time))
    stop("estimated.surv and estimated.time must have the same length.")

  prepareSeries <- function(time_vec, surv_vec) {
    ok <- !(is.na(time_vec) | is.na(surv_vec))
    time_vec <- time_vec[ok]; surv_vec <- surv_vec[ok]
    if (!length(time_vec)) return(list(t = numeric(0), s = numeric(0)))
    o <- order(time_vec)
    t2 <- time_vec[o]; s2 <- surv_vec[o]
    keep <- !duplicated(t2, fromLast = TRUE)
    list(t = t2[keep], s = s2[keep])
  }

  n_pred <- length(predicted.time)
  predicted.surv <- numeric(n_pred)

  strata_mode <- !(
    is.null(predicted.strata) || is.null(estimated.strata) || is.null(strata.levels) ||
      length(estimated.strata) == 0L || length(strata.levels) == 0L
  )

  if (!strata_mode) {
    ser <- prepareSeries(estimated.time, estimated.surv)
    if (!length(ser$t)) return(rep(1.0, n_pred))
    for (i in seq_len(n_pred)) {
      idx <- findInterval(predicted.time[i], ser$t, left.open = TRUE)
      predicted.surv[i] <- if (idx > 0L) ser$s[idx] else 1.0
    }
    return(predicted.surv)
  }

  if (!is.numeric(estimated.strata) || any(estimated.strata < 0))
    stop("'estimated.strata' must be a non-negative integer vector.")
  if (sum(estimated.strata) != length(estimated.time))
    stop("sum(estimated.strata) must equal length(estimated.time).")
  K <- length(estimated.strata)
  if (length(strata.levels) != K)
    stop("'strata.levels' must have length K = length(estimated.strata).")

  if (length(predicted.strata) == 1L) {
    predicted.strata <- rep(predicted.strata, n_pred)
  } else if (length(predicted.strata) != n_pred) {
    stop("Length of predicted.strata must be 1 or match length(predicted.time).")
  }

  mapped <- if (is.factor(predicted.strata)) {
    match(as.character(predicted.strata), as.character(strata.levels))
  } else {
    match(predicted.strata, strata.levels)
  }
  if (any(is.na(mapped))) {
    bad <- unique(predicted.strata[is.na(mapped)])
    stop("Some values in predicted.strata are not found in 'strata.levels': ",
         paste(bad, collapse = ", "))
  }

  strata_start <- c(1L, head(cumsum(estimated.strata), -1L) + 1L)
  strata_end   <- cumsum(estimated.strata)

  for (i in seq_len(n_pred)) {
    s <- mapped[i]
    if (estimated.strata[s] == 0L) { out[i] <- NA_real_; next }
    idx <- strata_start[s]:strata_end[s]
    ser <- prepareSeries(estimated.time[idx], estimated.surv[idx])
    if (!length(ser$t)) { out[i] <- 1.0; next }
    j <- findInterval(predicted.time[i], ser$t, left.open = TRUE)
    predicted.surv[i] <- if (j > 0L) ser$s[j] else 1.0
  }
  return(predicted.surv)
}

calculateIPCW <- function(formula, data, code.censoring, strata_name, specific.time,
                          return_censoring = FALSE) {
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  t <- Y[, 1L]
  epsilon <- Y[, 2L]
  strata <- if (is.null(strata_name)) rep.int(1L, nrow(data)) else
    as.integer(factor(data[[strata_name]]))

  km <- estimate_censoring_km(t, epsilon, code.censoring, strata,
    tie.order = "censoring_first")
  weights <- vapply(specific.time, function(horizon) {
    calculate_ipcw_weights(km, pmin(t, horizon),
      observed = t > horizon | epsilon != code.censoring, strata = strata)$weights
  }, numeric(length(t)))
  weights <- matrix(weights, nrow = length(t), ncol = length(specific.time))
  if (return_censoring) list(weights = weights, km = km, epsilon = epsilon, strata = strata) else weights
}

#calculateIPCW_phreg <- function(formula, data, code.censoring, strata_name, specific.time) {
#  cl <- match.call()
#  mf <- match.call(expand.dots = TRUE)[1:3]
#  special <- c("strata", "cluster", "offset")
#  Terms <- terms(formula, special, data = data)
#  mf$formula <- Terms
#  mf[[1]] <- as.name("model.frame")
#  mf <- eval(mf, parent.frame())
#  Y <- model.extract(mf, "response")
#  if (!inherits(Y, c("Event", "Surv")))
#    stop("Expected a 'Surv' or 'Event'-object")
#  if (ncol(Y) == 2) {
#    t <- Y[, 1]
#    epsilon <- Y[, 2]
#  }

#  censoring.model <- createCensoringFormula(formula=formula, code.censoring.updated=code.censoring, strata_name=strata_name)
#  resC <- mets::phreg(censoring.model, data)
#  if (resC$p > 0) kmt <- FALSE
#  kmt <- TRUE
#  s <- rep(specific.time, length(t))
#  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = t, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
#  out_predict2 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = s, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
#  km1 <- out_predict1$surv
#  km2 <- out_predict2$surv

#  tmp1 <- ifelse(km1 > 0, 1 / km1, 0)
#  tmp2 <- ifelse(km2 > 0, 1 / km2, 0)
#  censoring_status <- as.numeric(epsilon==code.censoring)
#  tmp3 <- (t <= specific.time) * (censoring_status==0) * tmp1
#  tmp4 <- (t > specific.time) * tmp2
#  ip.weight <- tmp3 + tmp4
#  if (any(is.na(ip.weight)))
#    stop("Inverse probability weights contain NA values")
#  return(ip.weight)
#}

createCensoringFormula <- function(formula, code.censoring.updated, strata_name = NULL) {
  lhs <- formula[[2]]
  if (is.call(lhs) && length(lhs) >= 3) {
    lhs[[3]] <- call("==", lhs[[3]], code.censoring.updated)
  }
  if (is.null(strata_name)) {
    rhs <- quote(+1)  # Intercept-only model
  } else {
    if (!is.character(strata_name)) {
      strata_name <- deparse(substitute(strata_name))
    }
    rhs <- as.call(list(as.symbol("strata"), as.symbol(strata_name)))
  }
  censoring.model <- as.formula(call("~", lhs, rhs))
  return(censoring.model)
}

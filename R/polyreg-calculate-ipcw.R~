calculateIPCW <- function(formula, data, code.censoring, strata_name, specific.time) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (ncol(Y) == 2) {
    t <- Y[, 1]
    epsilon <- Y[, 2]
  }

  d <- as.integer(epsilon==0)
  if (is.null(strata_name)) {
    strata <- rep(1,nrow(data))
  } else {
    strata <- data[[strata_name]]
    strata <- as.integer(strata)
  }

  out_km <- calculateKM(t=t, d=d, strata=strata, error="none")
  s <- rep(specific.time, length(t))
  km1 <- get_surv(t, out_km$surv, out_km$time, strata, out_km$strata, out_km$strata.levels)
  km2 <- get_surv(s, out_km$surv, out_km$time, strata, out_km$strata, out_km$strata.levels)

  tmp1 <- ifelse(km1 > 0, 1 / km1, 0)
  tmp2 <- ifelse(km2 > 0, 1 / km2, 0)
  censoring_status <- as.numeric(epsilon==code.censoring)
  tmp3 <- (t <= specific.time) * (censoring_status==0) * tmp1
  tmp4 <- (t > specific.time) * tmp2
  ip.weight <- tmp3 + tmp4

  if (any(is.na(ip.weight)))
    stop("Inverse probability weights contain NA values")
  return(as.matrix(ip.weight))
}

calculateIPCWMatrix <- function(formula, data, code.censoring, strata_name, estimand, out_normalizeCovariate) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv")))
    stop("Expected a 'Surv' or 'Event'-object")
  if (ncol(Y) == 2) {
    t <- Y[, 1]
    epsilon <- Y[, 2]
  }
  i_time <- 0
  ip.weight.matrix <- matrix(NA, nrow=out_normalizeCovariate$n, ncol=length(estimand$time.point))
  for (specific.time in estimand$time.point) {
    i_time <- i_time + 1
    ip.weight.matrix[,i_time] <- calculateIPCW(formula, data, code.censoring, strata_name, specific.time)
  }
  return(ip.weight.matrix)
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

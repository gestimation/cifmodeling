calculateCI <- function(survfit_object, conf.int, conf.type, conf.lower=NULL,
                        critical_value=NULL) {
  if (conf.int <= 0 | conf.int >= 1)
    stop("Confidence level must be between 0 and 1")
  alpha <- 1 - conf.int
  if (is.null(critical_value)) critical_value <- qnorm(1 - alpha / 2)
  if (is.null(conf.type) | conf.type == "none" | conf.type == "n") {
    lower <- NULL
    upper <- NULL
  } else if (conf.type == "arcsine-square root" | conf.type == "arcsin" | conf.type == "a") {
    #    se <- survfit_object$surv*survfit_object$std.err/2/sqrt(survfit_object$surv * (1 - survfit_object$surv))
    se <- survfit_object$std.err/2/sqrt(survfit_object$surv * (1 - survfit_object$surv))
    lower <- sin(pmax(asin(sqrt(survfit_object$surv)) - critical_value*se, 0))^2
    upper <- sin(pmin(asin(sqrt(survfit_object$surv)) + critical_value*se, pi/2))^2
  } else if (conf.type == "plain" | conf.type == "p" | conf.type == "linear") {
    #    se <- survfit_object$surv*survfit_object$std.err
    se <- survfit_object$std.err
    lower <- pmax(survfit_object$surv - critical_value*se, 0)
    upper <- pmin(survfit_object$surv + critical_value*se, 1)
  } else if (conf.type == "log") {
    #    se <- survfit_object$std.err
    se <- survfit_object$std.err/survfit_object$surv
    lower <- survfit_object$surv * exp(-critical_value*se)
    upper <- pmin(survfit_object$surv * exp(critical_value*se), 1)
  } else if (conf.type == "log-log") {
    #    se <- survfit_object$std.err / log(survfit_object$surv)
    se <- survfit_object$std.err/survfit_object$surv/log(survfit_object$surv)
    lower <- survfit_object$surv^exp(-critical_value*se)
    upper <- survfit_object$surv^exp(critical_value*se)
  } else if (conf.type == "logit") {
    #    se <- survfit_object$std.err/(1 - survfit_object$surv)
    se <- survfit_object$std.err/survfit_object$surv/(1 - survfit_object$surv)
    lower <- survfit_object$surv / (survfit_object$surv + (1 - survfit_object$surv)*exp(critical_value*se))
    upper <- survfit_object$surv / (survfit_object$surv + (1 - survfit_object$surv)*exp(-critical_value*se))
  }
  lower <- sapply(lower, function(x) ifelse(is.nan(x), NA, x))
  upper <- sapply(upper, function(x) ifelse(is.nan(x), NA, x))
  lower <- sapply(lower, function(x) ifelse(x>=1, 1, x))
  upper <- sapply(upper, function(x) ifelse(x>=1, 1, x))
  lower <- sapply(lower, function(x) ifelse(x<=0, 0, x))
  upper <- sapply(upper, function(x) ifelse(x<=0, 0, x))
  return(list(upper=upper, lower=lower))
}

calculateJackKnife <- function(data, estimator) {
  theta <- estimator(data)
  n <- nrow(data)
  theta_i <- rep(0, n)
  for (i in 1:n){
    theta_i[i] <- estimator(data[-i,])
  }
  pseudo_observations <- (n*theta) - (n-1)*theta_i
  std.err <- sqrt(var(pseudo_observations)/n)
  list(std.err = std.err, pseudo_observations = pseudo_observations)
}

calculateJackKnifeSE <- function(data, estimator, outputPseudo=FALSE) {
  estimate <- estimator(data)
  std.err <- NULL
  n <- table(data$strata)
  n_strata <- length(n)
  if (n_strata == 1) estimate$strata <- length(estimate$time)
  strata_start <- 1
  for (i in 1:n_strata) {
    if (i>1) strata_start <- strata_start + estimate$strata[i-1]
    strata_end <- strata_start + estimate$strata[i] -1
    index_strata <- strata_start:strata_end
    surv <- estimate$surv[index_strata]
    time <- estimate$time[index_strata]
    m <- length(time)
    pseudo_observations <- matrix(0, nrow = n[i], ncol = m)
    v <- rep(0, m)
    if (n_strata > 1) {
      t <- data$t[data$strata == levels(data$strata)[i]]
      d <- data$d[data$strata == levels(data$strata)[i]]
      w <- data$w[data$strata == levels(data$strata)[i]]
    } else {
      t <- data$t
      d <- data$d
      w <- data$w
    }
    for (j in 1:m) {
      for (k in 1:n[i]) {
        data_k <- list(t = t[-k], d = d[-k], w = w[-k], strata = rep(1, (n[i]-1)))
        estimate_k <- estimator(data_k)
        surv_k <- util_get_surv(time[j], estimate_k$surv, estimate_k$time)
        pseudo_observations[k, j] <- (n[i]*surv[j]) - (n[i]-1)*surv_k
      }
      v[j] <- var(pseudo_observations[, j], na.rm = TRUE)
    }
    std.err <- c(std.err, sqrt(v/n[i])/surv)
  }
  if (outputPseudo) {
    results <- cbind(time, std.err, t(pseudo_observations))
    results_sorted <- results[order(results[,1]), ]
    return(list(std.err = results_sorted[,2], pseudo.obs = results_sorted[,-c(1, 2)]))
  }
  return(std.err)
}

getInitialValues <- function(formula, data, outcome.type, exposure, estimand,
                             specific.time, data.initial.values = NULL, prob.bound = 1e-8) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")

  if (!inherits(Y, c("Event", "Surv"))) {
    t <- rep(0, length(Y))
    epsilon <- Y
    if (!all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
      stop("BINOMIAL requires event codes {censoring,event1}.")
    }
  } else {
    t <- as.numeric(Y[, 1])
    epsilon <- as.numeric(Y[, 2])
  }

  if (!is.null(attributes(Terms)$specials$offset)) {
    ts <- untangle.specials(Terms, "offset")
    if (length(ts$vars) > 0) {
      Terms2 <- Terms[-ts$terms]
      offset <- mf[[ts$vars]]
      x_l <- model.matrix(Terms2, mf)
    } else {
      offset <- rep(0, nrow(mf))
      x_l <- model.matrix(Terms, mf)
    }
  } else {
    offset <- rep(0, nrow(mf))
    x_l <- model.matrix(Terms, mf)
  }
  if (any(is.na(x_l))) stop("Covariate design (x_l) contains NA.")

  out_defineExposureDesign <- defineExposureDesign(data, exposure, estimand$code.exposure.ref)
  x_a <- out_defineExposureDesign$x_a
  p_l <- ncol(x_l)
  p_a <- ncol(x_a)

  if (!is.null(data.initial.values)) {
    expected_len <- if (outcome.type == "SURVIVAL" ||
                        all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
      p_l + p_a
    } else {
      2 * (p_l + p_a)
    }
    if (length(data.initial.values) != expected_len) {
      stop(sprintf("Invalid initial values length: expected %d, got %d.",
                   expected_len, length(data.initial.values)))
    }
    return(data.initial.values)
  }

  link_from_measure <- function(measure) {
    switch(toupper(measure),
           "RR"  = "log",
           "OR"  = "logit",
           "SHR" = "cloglog",
           stop("effect.measure must be RR, OR, or SHR"))
  }

  fit_init <- function(y01, x_l, x_a, link, offset = NULL, pb = 1e-8) {
    y_tilde <- as.numeric(y01)
    y_tilde <- y_tilde * (1 - 2 * pb) + pb

    X <- cbind(x_l, x_a)
    df <- data.frame(y = y_tilde, X)
    fam <- binomial(link = link)

    fit <- suppressWarnings(
      try(glm(y ~ . - 1, data = df,
              family = fam,
              weights = rep(1, nrow(df)),
              offset = offset),
          silent = TRUE)
    )

    if (inherits(fit, "try-error")) {
      co <- rep(0, ncol(X))
      names(co) <- colnames(X)
    } else {
      co <- coef(fit)
      co[is.na(co)] <- 0
    }
    list(alpha = unname(co[seq_len(p_l)]),
         beta  = unname(co[p_l + seq_len(p_a)]))
  }

  y1 <- as.integer(epsilon == estimand$code.event1 & t <= specific.time)

  if (outcome.type == "SURVIVAL" ||
      all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
    link1 <- link_from_measure(estimand$effect.measure1)
    est1  <- fit_init(y1, x_l, x_a, link1, offset = offset, pb = prob.bound)
    return(c(est1$alpha, est1$beta))
  } else {
    y2    <- as.integer(epsilon == estimand$code.event2 & t <= specific.time)
    link1 <- link_from_measure(estimand$effect.measure1)
    link2 <- link_from_measure(estimand$effect.measure2)

    est1  <- fit_init(y1, x_l, x_a, link1, offset = offset, pb = prob.bound)
    est2  <- fit_init(y2, x_l, x_a, link2, offset = offset, pb = prob.bound)

    return(c(est1$alpha, est1$beta, est2$alpha, est2$beta))
  }
}

getInitialValuesProportional <- function(formula, data, outcome.type, exposure, estimand, data.initial.values, prob.bound, out_normalizeCovariate) {
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

  n_para_1 <- out_normalizeCovariate$n_covariate + 1
  n_para_2 <- out_normalizeCovariate$n_covariate + 2
  n_para_3 <- out_normalizeCovariate$n_covariate + 3
  n_para_4 <- 2*out_normalizeCovariate$n_covariate + 3
  n_para_5 <- 2*out_normalizeCovariate$n_covariate + 4
  n_para_6 <- length(estimand$time.point) * (n_para_5 - 2) + 2

  if (outcome.type == "PROPORTIONAL") {
    alpha_beta_0 <- rep(NA_real_, n_para_6/2)
  } else {
    alpha_beta_0 <- rep(NA_real_, n_para_6)
  }

  sum1 <- 0
  sum2 <- 0
  i_time <- 0

  for (specific.time in estimand$time.point) {
    specific.time <- as.numeric(specific.time)
    out_getInitialValues <- getInitialValues(
      formula = formula,
      data = data,
      outcome.type = outcome.type,
      exposure = exposure,
      estimand = estimand,
      specific.time = specific.time,
      data.initial.values = data.initial.values,
      prob.bound = prob.bound
    )

    i_time <- i_time + 1
    i_para <- n_para_1*(i_time - 1) + 1

    tmp1 <- out_getInitialValues[seq.int(1, length.out = n_para_1)]
    if (outcome.type == "PROPORTIONAL") {
      idx1 <- seq.int(i_para, length.out = n_para_1)
      alpha_beta_0[idx1] <- tmp1
    } else {
      idx1 <- seq.int(i_para, length.out = n_para_1)
      alpha_beta_0[idx1] <- tmp1
    }
    sum1 <- sum1 + out_getInitialValues[n_para_2]

    if (outcome.type == "POLY-PROPORTIONAL") {
      tmp2 <- out_getInitialValues[seq.int(n_para_3, n_para_4)]
      idx2_start <- (n_para_6 %/% 2) + i_para
      idx2 <- seq.int(idx2_start, length.out = n_para_1)
      alpha_beta_0[idx2] <- tmp2
      sum2 <- sum2 + out_getInitialValues[n_para_5]
    }
  }
  if (outcome.type == "PROPORTIONAL") {
    alpha_beta_0[n_para_6 %/% 2] <- sum1 / length(estimand$time.point)
  } else {
    alpha_beta_0[n_para_6 %/% 2] <- sum1 / length(estimand$time.point)
    alpha_beta_0[n_para_6]       <- sum2 / length(estimand$time.point)
  }
  return(alpha_beta_0)
}

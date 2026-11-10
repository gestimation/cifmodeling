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

  out_readExposureDesign <- reg_read_exposure_design(data, exposure, estimand$code.exposure.ref)
  x_a <- out_readExposureDesign$x_a
  p_l <- ncol(x_l)
  p_a <- ncol(x_a)

  competing <- outcome.type %in% c("competing-risk", "proportional-competing-risk")
  if (!is.null(data.initial.values)) {
    expected_len <- (p_l + p_a) * if (competing) 2L else 1L
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

  if (!competing) {
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

getInitialValuesProportional <- function(formula, data, outcome.type, exposure, estimand,
                                         data.initial.values, prob.bound) {
  layout <- reg_proportional_index(estimand)
  iv <- estimand$index.vector
  n_causes <- if (outcome.type == "proportional-competing-risk") 2L else 1L
  values <- vapply(estimand$time.point, function(time) {
    getInitialValues(formula, data, outcome.type, exposure, estimand, time,
      data.initial.values, prob.bound)
  }, numeric(n_causes * iv[3L]))
  initial <- numeric(n_causes * layout$block)
  for (cause in seq_len(n_causes)) {
    local_offset <- (cause - 1L) * iv[3L]
    offset <- (cause - 1L) * layout$block
    initial[layout$nuisance + offset] <- values[seq_len(iv[1L]) + local_offset, , drop = FALSE]
    initial[layout$exposure + offset] <- rowMeans(
      values[seq.int(iv[2L], iv[3L]) + local_offset, , drop = FALSE]
    )
  }
  initial
}

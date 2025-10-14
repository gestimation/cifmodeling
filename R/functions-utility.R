clampP <- function(p, eps = 1e-5) pmin(pmax(p, eps), 1 - eps)

clampLogP <- function(x, eps = 1e-5) {
  if (!is.numeric(x)) stop("`x` must be numeric")
  low <- log(eps)
  high <- log1p(-eps)
  pmin(pmax(x, low), high)
}

calculateIndexForParameter <- function(i_parameter,x_l,x_a,length.time.point) {
  i_parameter[1] <- ncol(x_l)
  i_parameter[2] <- i_parameter[1] + 1
  i_parameter[3] <- i_parameter[1] + ncol(x_a)
  i_parameter[4] <- i_parameter[1] + ncol(x_a) + 1
  i_parameter[5] <- 2 * i_parameter[1] + ncol(x_a)
  i_parameter[6] <- 2 * i_parameter[1] + ncol(x_a) + 1
  i_parameter[7] <- 2 * i_parameter[1] + 2 * ncol(x_a)
  i_parameter[8] <- length.time.point*(2 * i_parameter[1]) + 2 * ncol(x_a)
  return(i_parameter)
}

normalizeCovariate <- function(formula, data, should.normalize.covariate,
                                   outcome.type, exposure.levels) {
  mf <- model.frame(formula, data)
  Y  <- model.extract(mf, "response")
  if (inherits(Y, c("Surv","Event"))) {
    response_vars  <- all.vars(formula[[2]])
    covariate_cols <- setdiff(all.vars(formula), response_vars)
  } else {
    covariate_cols <- all.vars(formula)[-1]
  }

  normalized_data <- data

  robust_scale <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(1)
    s <- IQR(x, na.rm=TRUE)
    if (!is.finite(s) || s==0) {
      s <- mad(x, center = median(x, na.rm=TRUE), constant = 1.4826, na.rm=TRUE)
    }
    if (!is.finite(s) || s==0) s <- sd(x, na.rm=TRUE)
    if (!is.finite(s) || s==0) s <- 1
    s
  }

  if (length(covariate_cols) > 0) {
    mm <- model.matrix(reformulate(covariate_cols), data = data)
    mm_cols <- colnames(mm)
  } else {
    mm <- cbind(`(Intercept)` = 1)
    mm_cols <- "(Intercept)"
  }

  is_num <- if (length(covariate_cols) > 0)
    vapply(data[covariate_cols], is.numeric, logical(1)) else logical(0)
  num_covars <- covariate_cols[is_num]

  scales_num <- numeric(0)
  if (isTRUE(should.normalize.covariate) && length(num_covars) > 0) {
    for (col in num_covars) {
      s <- robust_scale(normalized_data[[col]])
      normalized_data[[col]] <- normalized_data[[col]] / s
      scales_num[col] <- s
    }
  } else if (length(num_covars) > 0) {
    scales_num <- setNames(rep(1, length(num_covars)), num_covars)
  }

  scales_in_mm_order <- c()
  if (length(mm_cols) > 1) {
    for (cn in mm_cols[-1]) {
      scales_in_mm_order <- c(scales_in_mm_order, scales_num[cn] %||% 1)
    }
  }

  k_ex <- max(0L, as.integer(exposure.levels) - 1L)

  block <- c(1, as.numeric(scales_in_mm_order), rep(1, k_ex))

  if (outcome.type %in% c("PROPORTIONAL","POLY-PROPORTIONAL")) {
    range_for_params <- NULL
  } else if (outcome.type %in% c("SURVIVAL","BINOMIAL")) {
    range_for_params <- block
  } else if (outcome.type == "COMPETING-RISK") {
    range_for_params <- c(block, block)
  } else {
    stop("Unknown outcome.type: ", outcome.type)
  }

  list(
    normalized_data = normalized_data,
    range           = range_for_params,
    n_covariate_all = length(covariate_cols),
    n               = nrow(data),
    k_ex            = k_ex,
    mm_cols         = mm_cols,
    scaled_vars     = names(scales_num)
  )
}

normalizeEstimate <- function(
    outcome.type,
    report.sandwich.conf,
    should.normalize.covariate,
    current_params,
    out_getResults,
    estimand,
    prob.bound,
    out_normalizeCovariate,
    out_calculateCov = NULL
) {
  if (isFALSE(report.sandwich.conf)) {
    if (isTRUE(should.normalize.covariate) && !is.null(out_normalizeCovariate$range)) {
      adj <- 1 / as.vector(out_normalizeCovariate$range)
      if (length(adj) != length(current_params)) {
        stop(sprintf(
          "Length mismatch: adj=%d vs params=%d. outcome.type=%s, p_num=%d, k_ex=%d",
          length(adj), length(current_params), outcome.type,
          out_normalizeCovariate$p_num %||% NA_integer_,
          out_normalizeCovariate$k_ex  %||% NA_integer_
        ))
      }
      alpha_beta_estimated <- adj * current_params
    } else {
      alpha_beta_estimated <- current_params
    }
    return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = NULL))
  }

  if (isTRUE(should.normalize.covariate) && !is.null(out_normalizeCovariate$range)) {
    adj <- 1 / as.vector(out_normalizeCovariate$range)

    if (length(adj) != length(current_params)) {
      stop(sprintf(
        "Length mismatch: adj=%d vs params=%d. outcome.type=%s, p_num=%d, k_ex=%d",
        length(adj), length(current_params), outcome.type,
        out_normalizeCovariate$p_num %||% NA_integer_,
        out_normalizeCovariate$k_ex  %||% NA_integer_
      ))
    }

    alpha_beta_estimated <- adj * current_params
    if (is.null(out_calculateCov) || is.null(out_calculateCov$cov_estimated)) {
      stop("out_calculateCov$cov_estimated is required when report.sandwich.conf=TRUE.")
    }
    A <- diag(adj, length(adj))
    cov_estimated <- A %*% out_calculateCov$cov_estimated %*% A
  } else {
    alpha_beta_estimated <- current_params
    if (is.null(out_calculateCov) || is.null(out_calculateCov$cov_estimated)) {
      stop("out_calculateCov$cov_estimated is required when report.sandwich.conf=TRUE.")
    }
    cov_estimated <- out_calculateCov$cov_estimated
  }

  list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = cov_estimated)
}


normalizeCovariate_old <- function(formula, data, should.normalize.covariate, outcome.type, exposure.levels) {
  mf <- model.frame(formula, data)
  Y <- model.extract(mf, "response")
  response_term <- formula[[2]]
  if (inherits(mf[[1]], "Surv") || inherits(mf[[1]], "Event")) {
    response_vars <- all.vars(response_term)
    covariate_cols <- setdiff(all.vars(formula), response_vars)
  } else {
    covariate_cols <- all.vars(formula)[-1]
  }
  normalized_data <- data
  range_vector <- 1
  exposure.range <- matrix(1, 1, exposure.levels-1)
  if (should.normalize.covariate == TRUE & length(covariate_cols)>0) {
    for (col in covariate_cols) {
      x <- normalized_data[[col]]
      range <- max(x)-min(x)
      normalized_data[[col]] <- x/range
      range_vector <- cbind(range_vector, range)
    }
    if (outcome.type == "PROPORTIONAL" || outcome.type == "POLY-PROPORTIONAL") {
      range_vector <- NULL
    } else if (outcome.type == "SURVIVAL" || outcome.type == "BINOMIAL") {
      range_vector <- cbind(range_vector,exposure.range)
    } else {
      range_vector <- cbind(range_vector,exposure.range,range_vector,exposure.range)
    }
  } else {
    if (outcome.type == "PROPORTIONAL" || outcome.type == "POLY-PROPORTIONAL") {
      range_vector <- NULL
    } else if (outcome.type == "SURVIVAL" || outcome.type == "BINOMIAL") {
      range_vector <- rep(1, (length(covariate_cols)+exposure.levels))
    } else {
      range_vector <- rep(1, (2*length(covariate_cols)+2*exposure.levels))
    }
  }
  n_covariate <- length(covariate_cols)
  out <- list(normalized_data=normalized_data, range=range_vector, n_covariate=n_covariate, n=nrow(data))
  return(out)
}


normalizeEstimate_old <- function(
    outcome.type,
    report.sandwich.conf,
    should.normalize.covariate,
    current_params,
    out_getResults,
    estimand,
    prob.bound,
    out_normalizeCovariate
) {
  if (report.sandwich.conf == FALSE) {
    alpha_beta_estimated <- if (should.normalize.covariate) {
      adj <- 1 / as.vector(out_normalizeCovariate$range)
      if (length(adj) != length(current_params)) stop("Length of adj (range) must match length of current_params.")
      adj * current_params
    } else {
      current_params
    }
    return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = NULL))
  }
  if (should.normalize.covariate) { #本来はoutcome.typeで分岐が必要
    adj <- 1 / as.vector(out_normalizeCovariate$range)
    if (length(adj) != length(current_params)) stop("Length of adj (range) must match length of current_params.")
    alpha_beta_estimated <- adj * current_params
    adj_matrix <- diag(adj, length(adj))
    cov_estimated <- adj_matrix %*% out_calculateCov$cov_estimated %*% adj_matrix
  } else {
    alpha_beta_estimated <- current_params
    cov_estimated <- out_calculateCov$cov_estimated
  }
  return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = cov_estimated))
}

clampP <- function(p, eps = 1e-5) pmin(pmax(p, eps), 1 - eps)

clampLogP <- function(x, eps = 1e-5) {
  ex <- exp(x)
  x[ex < eps]         <- log(eps)
  x[(1 - ex) < eps]   <- log(1 - eps)
  x
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

normalizeCovariate <- function(formula, data, should.normalize.covariate, outcome.type, exposure.levels) {
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

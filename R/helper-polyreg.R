#readStrata <- function(out_readSurv, out_aj, label.strata=NULL) {
#  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
#    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))
#  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
#    names(out_aj$strata1) <- label.strata
#  }
#  return(out_aj)
#}

#check_label.strata <- function(out_readSurv, label.strata) {
#  strata_levels <- levels(as.factor(out_readSurv$strata))
#  n_strata <- length(strata_levels)
#  if (!is.null(label.strata)) {
#    if (!is.character(label.strata)) {
#      stop("`label.strata` must be a character vector.", call. = FALSE)
#    }
#    if (length(label.strata) != n_strata) {
#      stop(sprintf("`label.strata` must have length %d (number of strata), but got %d.", n_strata, length(label.strata)), call. = FALSE)
#    }
#  }
#}


clampP <- function(p, eps = 1e-5) pmin(pmax(p, eps), 1 - eps)

clampLogP <- function(x, eps = 1e-5) {
  if (!is.numeric(x)) stop("`x` must be numeric")
  low <- log(eps)
  high <- log1p(-eps)
  pmin(pmax(x, low), high)
}

reg_index_for_parameter <- function(i_parameter,x_l,x_a,length.time.point) {
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

reg_normalize_covariate <- function(formula, data, should.normalize.covariate,
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

reg_normalize_estimate <- function(
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

reg_check_input <- function(data, formula, exposure, code.event1, code.event2, code.censoring,
                       code.exposure.ref, outcome.type, conf.level, report.sandwich.conf,
                       report.boot.conf, nleqslv.method, should.normalize.covariate,
                       strata = NULL, subset.condition = NULL, na.action = na.omit) {

  other_vars <- c(exposure, strata)
  other_vars <- other_vars[!is.null(other_vars) & nzchar(other_vars)]
  data <- createAnalysisDataset(formula, data,
                                other.variables.analyzed = other_vars,
                                subset.condition = subset.condition,
                                na.action = na.action)

  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata)) .err("no_strata_in_formula")
  if (!is.null(attr(out_terms, "specials")$offset))  .err("no_offset_in_formula")
  if (!is.null(attr(out_terms, "specials")$cluster)) .err("no_cluster_in_formula")

  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- na.action
  mf <- eval(mf, parent.frame())

  Y <- model.extract(mf, "response")
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL","PROPORTIONAL","POLY-PROPORTIONAL")) {
    if (!inherits(Y, c("Event","Surv"))) .err("surv_expected")
    t <- as.numeric(Y[,1]); epsilon <- as.numeric(Y[,2])
    if (any(t < 0, na.rm = TRUE)) .err("time_nonneg", arg = "time")
    if (outcome.type == "SURVIVAL") {
      ok <- all(epsilon %in% c(code.event1, code.censoring), na.rm = TRUE)
      if (!ok) .err("codes_required_surv",
                    censoring = code.censoring, event1 = code.event1,
                    found = paste(sort(unique(na.omit(epsilon))), collapse = ", "))
    } else if (outcome.type == "COMPETING-RISK") {
      ok <- all(epsilon %in% c(code.event1, code.event2, code.censoring), na.rm = TRUE)
      if (!ok) .err("codes_required_cr",
                    censoring = code.censoring, event1 = code.event1, event2 = code.event2,
                    found = paste(sort(unique(na.omit(epsilon))), collapse = ", "))
    }
  }

  out_readExposureDesign <- reg_read_exposure_design(data, exposure, code.exposure.ref)
  x_a <- out_readExposureDesign$x_a
  x_l <- model.matrix(out_terms, mf)

  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1) .err("conf_level")

  if (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL") {
    should.normalize.covariate.corrected <- FALSE
    report.sandwich.conf.corrected <- FALSE
    if (is.null(report.boot.conf)) {
      report.boot.conf.corrected <- TRUE
    } else {
      report.boot.conf.corrected <- report.boot.conf
    }
  } else {
    should.normalize.covariate.corrected <- should.normalize.covariate
    if (is.null(report.boot.conf)) {
      report.boot.conf.corrected <- FALSE
    } else {
      report.boot.conf.corrected <- report.boot.conf
    }
    if (report.boot.conf == FALSE || is.null(report.boot.conf)) {
      report.sandwich.conf.corrected <- report.sandwich.conf
    } else {
      report.sandwich.conf.corrected <- FALSE
    }
  }

  outer_choices <- c("nleqslv","Newton","Broyden")
  nleqslv.method <- match.arg(nleqslv.method, choices = outer_choices)
  return(list(should.normalize.covariate = should.normalize.covariate.corrected, report.sandwich.conf = report.sandwich.conf.corrected, report.boot.conf = report.boot.conf.corrected, out_readExposureDesign=out_readExposureDesign, x_a=x_a, x_l=x_l))
}

reg_read_exposure_design <- function(data, exposure, code.exposure.ref = NULL, prefix = "a") {
  stopifnot(is.data.frame(data))
  if (!exposure %in% names(data)) {
    stop("exposure = '", exposure, "' is not found in data.")
  }

  a_ <- data[[exposure]]
  a_ <- factor(a_)
  a_ <- base::droplevels(a_)
  lev <- levels(a_)
  K <- length(lev)

  ref_lab <- NULL
  if (!is.null(code.exposure.ref)) {
    ref_lab <- if (is.numeric(code.exposure.ref)) as.character(code.exposure.ref) else code.exposure.ref
    if (length(lev) > 0 && ref_lab %in% lev) {
      a_ <- stats::relevel(a_, ref = ref_lab)
    } else {
      warning("code.exposure.ref = ", ref_lab," is not found among factor levels. The first level is used as reference.")
      ref_lab <- NULL
    }
  }
  if (is.null(ref_lab)) ref_lab <- lev[1L]
  if (K < 1 || K == 1) stop("Exposure has only one level (", lev, ") or no valid levels. Effect estimation is not possible.")
  X <- stats::model.matrix(~ a_)[, -1, drop = FALSE]
  cn <- colnames(X)
  cn <- sub("^a_", paste0(prefix, "_"), cn, perl = TRUE)
  colnames(X) <- cn
  return(list(
    x_a = as.matrix(X),
    exposure.levels = K,
    exposure.labels = lev,
    ref = ref_lab
  ))
}

reg_read_time.point <- function(formula, data, x_a, outcome.type, code.censoring, should.terminate.time.point, time.point) {
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL")) {
    if (is.null(time.point) || !length(time.point)) .err("timepoint_required")
    tp <- suppressWarnings(max(time.point, na.rm = TRUE))
    if (!is.finite(tp) || tp < 0) .err("timepoint_nonneg_finite")
    return(tp)
  } else if (outcome.type == "BINOMIAL") {
    tp <- Inf
    return(tp)
  } else if (outcome.type %in% c("PROPORTIONAL","POLY-PROPORTIONAL") & is.null(time.point)) {
    cl <- match.call()
    mf <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster", "offset")
    Terms <- terms(formula, special, data = data)
    mf$formula <- Terms
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.extract(mf, "response")
    t <- Y[, 1]
    epsilon <- Y[, 2]
    tp <- t[epsilon != code.censoring]
    tp <- sort(unique(tp[is.finite(tp)]))
    if (should.terminate.time.point) {
      valid <- is.finite(t) & !is.na(epsilon) & (epsilon != code.censoring)
      if (ncol(x_a) <= 1L) {
        idx0 <- valid & (x_a[, 1L] == 0)
        idx1 <- valid & (x_a[, 1L] != 0)
        maxs <- c(if (any(idx0)) max(t[idx0]) else Inf,
                  if (any(idx1)) max(t[idx1]) else Inf)
        cutoff <- min(maxs)
      } else {
        rs   <- rowSums(x_a != 0, na.rm = TRUE)
        m0   <- if (any(valid & rs == 0L)) max(t[valid & rs == 0L]) else Inf
        mj   <- vapply(seq_len(ncol(x_a)), function(j) {
          idx <- valid & (x_a[, j] != 0)
          if (any(idx)) max(t[idx]) else Inf
        }, numeric(1))
        cutoff <- min(c(m0, mj))
      }
      tp     <- tp[tp <= cutoff]
    }
    return(tp)
  } else {
    return(time.point)
  }
}

reg_choose_nleqslv_method <- function(nleqslv.method) {
  if (nleqslv.method == "nleqslv" || nleqslv.method == "Broyden") {
    "Broyden"
  } else if (nleqslv.method == "Newton") {
    "Newton"
  } else {
    stop("Unsupported nleqslv.method: ", nleqslv.method)
  }
}

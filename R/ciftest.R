#data(diabetes.complications)
#ciftest(Event(t, epsilon)~fruitq1, data=diabetes.complications, test.rho=1, test.gamma=1)

ciftest <- function(
    formula,
    data,
    subset.condition = NULL,
    na.action = stats::na.omit,
    outcome.type = "COMPETING-RISK",
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.strata.ref = 0,
    code.weight.ref = 0,
    test.type = c("cumulative incidence function", "cause-specific hazard"),
    effect.measure = c("DIFF","RR","OR"),
    test.contrast = NULL,
    test.weight   = NULL,
    test.rho      = NULL,
    test.gamma    = NULL,
    conf.int = 0.95,
    n_simulation = 10000
){
  outcome.type <- util_check_outcome_type(outcome.type, formula=formula, data=data)
  type <- test_check_type(match.arg(test.type))
  effect.measure <- match.arg(effect.measure)

  out <- util_read_surv(
    formula = formula, data = data, weights = NULL,
    code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
    subset.condition = subset.condition, na.action = na.action
  )

  strata_fac <- if (is.factor(out$strata)) out$strata else factor(out$strata)
  lev_codes  <- levels(strata_fac)
  K          <- length(lev_codes)
  if (K < 2L) stop("the number of strata should be 2 or more")
  group_idx  <- as.integer(strata_fac)
  code_map   <- data.frame(group = 1:K, code = lev_codes, stringsAsFactors = FALSE)

  if (identical(code.strata.ref, 0)) {
    pos0 <- match("0", lev_codes)
    ref_group <- if (is.na(pos0)) 1L else as.integer(pos0)
  } else {
    m <- match(as.character(code.strata.ref), lev_codes)
    if (is.na(m)) stop("code.strata.ref is not founc in strata levels")
    ref_group <- as.integer(m)
  }

  weight_base <- if (identical(code.weight.ref, 0)) 0L else {
    m <- match(as.character(code.weight.ref), lev_codes)
    if (is.na(m)) stop("code.weight.ref is not founc in strata levels")
    as.integer(m)
  }

  Cmat <- build_contrasts(K, test.contrast, ref_group)

  weights_list <- test.weight
  if (is.null(weights_list)) weights_list <- list()
  need_beta <- (length(weights_list) == 0L)

  if (need_beta) {
    if (is.null(test.rho) || is.null(test.gamma)) stop("test.rho and test.gamma are required when test.weight=NULL")
    if (length(test.rho)!=length(test.gamma) || length(test.rho)==0L) stop("test.rho and test.gamma should have the same length>0")
  }

  if (type == "cumulative incidence function") {
    res_if <- calculateIFofAJ(out$t, out$epsilon, as.integer(group_idx))
    time   <- res_if$time
    aj1    <- res_if$aj1
    if_aj1 <- res_if$if_aj1

    out_wa <- calculateWeightedAverage(
      time        = time,
      aj1         = aj1,
      if_aj1      = if_aj1,
      contrasts   = Cmat,
      weights     = weights_list,
      rho         = if (need_beta) test.rho else NULL,
      gamma       = if (need_beta) test.gamma else NULL,
      weight_base = weight_base,
      conf_int    = conf.int,
      ref         = ref_group
    )

    if (effect.measure == "DIFF") {
      df  <- out_wa$diff_aj1
      cov <- out_wa$cov_diff
    } else if (effect.measure == "RR") {
      df  <- out_wa$rr_aj1
      cov <- out_wa$cov_rr
    } else { # "OR"
      df  <- out_wa$or_aj1
      cov <- out_wa$cov_or
    }
    z <- df$z
    mc <- calculateMaxCombo(z, cov, n_simulation)

    return(structure(list(
      call      = match.call(),
      type      = type,
      effect    = effect.measure,
      mapping   = list(
        strata_levels = lev_codes,
        code_map      = code_map,
        ref_group     = ref_group,
        weight_base   = weight_base
      ),
      grid      = list(time = time),
      weights   = out_wa$weight,
      estimates = df,
      cov       = cov,
      maxcombo  = mc
    ), class="ciftest"))

  } else {
    if (effect.measure != "DIFF") {
      stop("Only effect.measure=DIFF is supported for type='cause-specific hazard'")
    }

    res_na <- calculateIFofNA(out$t, out$epsilon, as.integer(group_idx))
    time   <- res_na$time
    H1     <- res_na$cumhaz1
    IFH    <- res_na$if_cumhaz1

    if (need_beta) {
      # Weights are calculated using rho, gamma and "base function", KM pooled on time grid calculated using any event=(epsilon>0)
      base <- km_pooled_on_grid(t, as.integer(out$epsilon>0), time)
      out_lr <- calculateWeightedAverageLinear(
        time      = time,
        curve     = H1,
        if_list   = IFH,
        contrasts = Cmat,
        weights   = list(),
        rho       = test.rho,
        gamma     = test.gamma,
        base      = base,
        conf_int  = conf.int
      )
    } else {
      out_lr <- calculateWeightedAverageLinear(
        time      = time,
        curve     = H1,
        if_list   = IFH,
        contrasts = Cmat,
        weights   = weights_list,
        conf_int  = conf.int
      )
    }

    df  <- out_lr$table
    cov <- out_lr$cov
    z   <- df$z
    mc  <- calculateMaxCombo(z, cov, n_simulation)

    return(structure(list(
      call      = match.call(),
      type      = type,
      effect    = effect.measure,
      mapping   = list(
        strata_levels = lev_codes,
        code_map      = code_map,
        ref_group     = ref_group,
        weight_base   = weight_base
      ),
      grid      = list(time = time),
      weights   = out_lr$weight,
      estimates = df,
      cov       = cov,
      maxcombo  = mc
    ), class="ciftest"))
  }
}

km_pooled_on_grid <- function(t, event_any, time_grid) {
  m <- length(time_grid)-1L
  Y <- d <- numeric(m+1L)
  out <- numeric(m+1L); out[1] <- 1
  for (j in 2:(m+1L)) {
    tj <- time_grid[j]
    Y[j] <- sum(t + 1e-12 >= tj)
    d[j] <- sum(abs(t - tj) <= 1e-12 & event_any==1L)
    if (Y[j] > 0) {
      out[j] <- out[j-1] * (1 - d[j]/Y[j])
    } else {
      out[j] <- out[j-1]
    }
  }
  return(out)
}

build_contrasts <- function(K, C, ref_group){
  if (is.null(C)) {
    M <- matrix(0, nrow=K-1, ncol=K)
    r <- 0L
    for (g in seq_len(K)) if (g!=ref_group){ r <- r+1L; M[r,g] <- 1; M[r,ref_group] <- -1 }
    return(M)
  }
  if (is.matrix(C)) {
    M <- C
  } else {
    v <- as.numeric(C)
    if (length(v)!=K) stop("test.contrast should have length same as the number of strata")
    M <- matrix(v, nrow=1)
  }
  if (ncol(M)!=K) stop("test.contrast should have length same as the number of stratum")
  if (any(abs(rowSums(M))>1e-8)) stop("Each row of test.contrast should sum to zero")
  M
}

test_check_type <- function(type) {
  x <- tolower(trimws(type))
  if (x %in% c("cumulative incidence function","cumulative incidence","cif")) return("cumulative incidence function")
  if (x %in% c("cause-specific hazard","cause specific hazard","csh"))     return("cause-specific hazard")
  stop("type must be one of {'cumulative incidence function','cause-specific hazard'} (case-insensitive).")
}

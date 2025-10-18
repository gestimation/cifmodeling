reportEffects <- function(outcome.type,
                          report.nuisance.parameter,
                          report.optim.convergence,
                          report.sandwich.conf,
                          report.boot.conf,
                          nuisance.model,
                          exposure,
                          estimand,
                          alpha_beta_estimated,
                          cov_estimated,
                          out_bootstrap,
                          out_getResults,
                          iteration,
                          converged.by,
                          objective.function,
                          max.absolute.difference,
                          relative.difference,
                          sol,
                          conf.level,
                          nleqslv.method) {

  out_getCoefTerm <- getCoefTerm(report.nuisance.parameter, report.sandwich.conf, report.boot.conf, nuisance.model, outcome.type, exposure, estimand, alpha_beta_estimated, cov_estimated, out_bootstrap, conf.level, out_getResults)
  coef1 <- out_getCoefTerm$coef1
  coef2 <- out_getCoefTerm$coef2
  terms_text <- out_getCoefTerm$terms_text

  median_followup <- paste(round(median(out_getResults$t), 2))
  range_followup <- paste("[",round(min(out_getResults$t), 2), ",",round(max(out_getResults$t), 2), "]")

  tg_event1_tidy <- tidy_df(coef1, terms_text)
  if (outcome.type == "COMPETING-RISK" | outcome.type == "SURVIVAL" | outcome.type == "BINOMIAL") {
    effect1_text <- paste(estimand$effect.measure1, "at", estimand$time.point)
  } else if (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL") {
    effect1_text <- paste(estimand$effect.measure1, "of", exposure, "over time")
  }

  tg_event1_glance <- glance_df(
    report.optim.convergence,
    effect_text                  = effect1_text,
    event_text                   = paste(sum(out_getResults$y_1), "in N =", length(out_getResults$y_1)),
    median_followup              = median_followup,
    range_followup               = range_followup,
    param_len                    = length(alpha_beta_estimated),
    iteration                    = iteration,
    converged.by                 = converged.by,
    objective.function           = objective.function,
    max.absolute.difference      = max.absolute.difference,
    relative.difference          = relative.difference,
    nleqslv.value.at.solution    = max(sol$fvec),
    nleqslv.iteration            = sol$iter,
    nleqslv.message              = sol$message
  )

  if (outcome.type == "COMPETING-RISK" | outcome.type == "POLY-PROPORTIONAL") {
    if (outcome.type == "COMPETING-RISK") {
      effect2_text <- paste(estimand$effect.measure2, "at", estimand$time.point)
    } else if (outcome.type == "POLY-PROPORTIONAL") {
      effect2_text <- paste(estimand$effect.measure2, "of", exposure, "over time")
    }
    tg_event2_tidy <- tidy_df(coef2, terms_text)
    tg_event2_glance <- glance_df(
      report.optim.convergence,
      effect_text                  = effect2_text,
      event_text                   = paste(sum(out_getResults$y_2), "in N =", length(out_getResults$y_2)),
      median_followup              = "-",
      range_followup               = "-",
      param_len                    = "-",
      iteration                    = "-",
      converged.by                 = "-",
      objective.function           = "-",
      max.absolute.difference      = "-",
      relative.difference          = "-",
      nleqslv.value.at.solution    = "-",
      nleqslv.iteration            = "-",
      nleqslv.message              = "-"
    )
    tg <- list(
      event1 = list(tidy = tg_event1_tidy, glance = tg_event1_glance),
      event2 = list(tidy = tg_event2_tidy, glance = tg_event2_glance)
    )
    class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
    return(tg)
  }

  tg1 <- list(
    tidy   = tg_event1_tidy,
    glance = tg_event1_glance
  )
  class(tg1) <- "modelsummary_list"
  tg <- list("event 1 (no competing risk)" = tg1)
  return(tg)
}

tidy_df <- function(coef, text) {
  data.frame(
    term = text,
    estimate = coef$coef,
    std.error = coef$coef_se,
    p.value = coef$p_value,
    conf.low = coef$conf_low,
    conf.high = coef$conf_high
  )
}

glance_df <- function(
    report.optim.convergence, effect_text, event_text, median_followup, range_followup,
    param_len, iteration, converged.by, objective.function, max.absolute.difference,relative.difference,
    nleqslv.value.at.solution,nleqslv.iteration,nleqslv.message) {
  if (report.optim.convergence == FALSE) {
    data.frame(
      effect.measure = effect_text,
      n.events = event_text,
      median.follow.up = median_followup,
      range.follow.up = range_followup,
      n.parameters = param_len,
      converged.by = converged.by,
      nleqslv.message = nleqslv.message
    )
  } else {
    data.frame(
      effect.measure = effect_text,
      n.events  = event_text,
      median.follow.up = median_followup,
      range.follow.up = range_followup,
      n.parameters = param_len,
      n.iteration = iteration,
      converged.by = converged.by,
      objective.function = objective.function,
      max.absolute.difference = max.absolute.difference,
      relative.difference = relative.difference,
      nleqslv.value.at.solution = nleqslv.value.at.solution,
      nleqslv.iteration = nleqslv.iteration,
      nleqslv.message = nleqslv.message
    )
  }
}


getCoefTerm <- function(report.nuisance.parameter, report.sandwich.conf, report.boot.conf, nuisance.model, outcome.type, exposure, estimand, alpha_beta_estimated, cov_estimated, out_bootstrap, conf.level, out_getResults) {
  iv <- estimand$index.vector
  code.exposure.ref <- as.character(estimand$code.exposure.ref)

  if (outcome.type == "COMPETING-RISK") {
    if (report.nuisance.parameter == TRUE) {
      index1 <- seq_len(iv[3])
      index2 <- seq.int(iv[4], iv[7])
      nuisance_terms <- c("Intercept", attr(terms(nuisance.model), "term.labels"))
    } else {
      index1 <- seq.int(iv[2], iv[3])
      index2 <- seq.int(iv[6], iv[7])
      nuisance_terms <- character(0)
    }
  }
  if (outcome.type == "SURVIVAL" || outcome.type == "BINOMIAL") {
    if (report.nuisance.parameter == TRUE) {
      index1 <- seq_len(iv[3])
      index2 <- NULL
      nuisance_terms <- c("Intercept", attr(terms(nuisance.model), "term.labels"))
    } else {
      index1 <- seq.int(iv[2], iv[3])
      index2 <- NULL
      nuisance_terms <- character(0)
    }
  }
  if (outcome.type == "POLY-PROPORTIONAL") {
    if (report.nuisance.parameter == TRUE) {
      index1 <- seq_len(iv[8]/2)
      index2 <- seq.int(iv[8]/2+1, iv[8])
      nuisance_terms <- c("Intercept", attr(terms(nuisance.model), "term.labels"))
    } else {
      index1 <- seq.int(iv[8]/2-estimand$exposure.levels+2, iv[8]/2)
      index2 <- seq.int(iv[8]-estimand$exposure.levels+2, iv[8])
      nuisance_terms <- character(0)
    }
  }
  if (outcome.type == "PROPORTIONAL") {
    if (report.nuisance.parameter == TRUE) {
      index1 <- seq_len(iv[8]/2)
      index2 <- NULL
      nuisance_terms <- c("Intercept", attr(terms(nuisance.model), "term.labels"))
    } else {
      index1 <- seq.int(iv[8]/2-estimand$exposure.levels+2, iv[8]/2)
      index2 <- NULL
      nuisance_terms <- character(0)
    }
  }
  coef1 <- getCoef(index1, alpha_beta_estimated, cov_estimated, report.sandwich.conf, report.boot.conf, out_bootstrap, conf.level)
  coef2 <- getCoef(index2, alpha_beta_estimated, cov_estimated, report.sandwich.conf, report.boot.conf, out_bootstrap, conf.level)
  code.exposure.nonref <- colnames(out_getResults$x_a)
  terms_text_exposure <- if (length(code.exposure.nonref)) {
    paste0(exposure, ", ", sub("^.*_", "", code.exposure.nonref), " vs ", code.exposure.ref)
  } else {
    character(0)
  }
  terms_text <- c(nuisance_terms, terms_text_exposure)
  return(list(coef1=coef1, coef2=coef2, terms_text=terms_text))
}

getCoef <- function(
    index,
    alpha_beta_estimated,
    cov_estimated,
    report.sandwich.conf,
    report.boot.conf,
    out_bootstrap,
    conf.level) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  if (report.sandwich.conf == TRUE) {
    coef <- alpha_beta_estimated[index]
    coef_se <- sqrt(diag(cov_estimated)[index])
    conf_low <- coef - critical_value * coef_se
    conf_high <- coef + critical_value * coef_se
    p_value <- 2 * (1 - pnorm(abs(coef) / coef_se))
  } else if (report.boot.conf == TRUE) {
    coef <- out_bootstrap$boot.coef[index]
    coef_se <- out_bootstrap$boot.coef_se[index]
    conf_low <- out_bootstrap$boot.conf_low[index]
    conf_high <- out_bootstrap$boot.conf_high[index]
    p_value <- out_bootstrap$boot.p_value[index]
  } else {
    coef <- alpha_beta_estimated[index]
    coef_se <- "N.A."
    conf_low <- "N.A."
    conf_high <- "N.A."
    p_value <- "N.A."
  }
  list(coef = coef, coef_se = coef_se, conf_low = conf_low, conf_high = conf_high, p_value = p_value)
}

create_rr_text <- function(coefficient, cov, index, omit.conf.int=TRUE, conf.int=0.95) {
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  coef <- coefficient[index]
  coef_se <- sqrt(diag(cov)[index])
  conf_low <- coef - critical_value * coef_se
  conf_high <- coef + critical_value * coef_se
  p_value <- floor(2 * (1 - pnorm(abs(coef) / coef_se)))
  if (omit.conf.int==TRUE) {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), ", p<0.01")
    else text <- paste0("RR=", round(exp(coef), digit=2), ", p=", p_value)
  } else {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p<0.01", ")")
    else text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p=", p_value, ")")
  }
  return(text)
}

check_ggsurvfit <- function(
    survfit_object,
    type.y,
    conf.type,
    label.y,
    limits.x, limits.y,
    breaks.x,  breaks.y,
    addConfidenceInterval,
    addCensorMark,
    addCompetingRiskMark,
    addIntercurrentEventMark,
    shape.censor.mark,
    shape.competing.risk.mark,
    shape.intercurrent.event.mark,
    out_readSurv,
    use_coord_cartesian,
    style
){
  if (isTRUE(addCensorMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.censor.mark, shape.intercurrent.event.mark)) {
    warning("shape.censor.mark and shape.intercurrent.event.mark specify an idencical type of symbol", call. = FALSE)
  }
  if (isTRUE(addCensorMark) && isTRUE(addCompetingRiskMark) &&
      identical(shape.censor.mark, shape.competing.risk.mark)) {
    warning("shape.censor.mark and shape.competing.risk.mark specify an idencical type of symbol", call. = FALSE)
  }
  if (isTRUE(addCompetingRiskMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.competing.risk.mark, shape.intercurrent.event.mark)) {
    warning("shape.competing.risk.mark and shape.intercurrent.event.mark specify an idencical type of symbol", call. = FALSE)
  }

  is_len2_num <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  is_increasing <- function(x) all(diff(x) > 0, na.rm = TRUE)
  is_nondec <- function(x) all(diff(x) >= 0, na.rm = TRUE)

  if (!is.null(limits.x)) {
    if (!is_len2_num(limits.x)) {
      warning("`limits.x` must be a numeric length-2 vector.", call. = FALSE)
    } else if (!is_increasing(limits.x)) {
      warning("`limits.x` must be strictly increasing, e.g., c(xmin, xmax).", call. = FALSE)
    }
    tmax <- suppressWarnings(max(survfit_object$time, na.rm = TRUE))
    if (is.finite(tmax) && is_len2_num(limits.x)) {
      if (tmax < limits.x[1] || tmax > limits.x[2]) {
        warning(sprintf("Max of `survfit_object$time` (%.4g) lies outside `limits.x` = [%.4g, %.4g].",
                        tmax, limits.x[1], limits.x[2]), call. = FALSE)
      }
    }
  } else if (!is.null(out_readSurv) && !is.null(out_readSurv$t)) {
    tmax <- suppressWarnings(max(out_readSurv$t, na.rm = TRUE))
    if (!is.finite(tmax) || tmax <= 0) {
      warning("`out_readSurv$t` provided but has no finite positive max; x-limits fallback may fail.", call. = FALSE)
    }
  }

  if (!is.null(limits.y)) {
    if (!is_len2_num(limits.y)) {
      warning("`limits.y` must be a numeric length-2 vector.", call. = FALSE)
    } else if (!is_increasing(limits.y)) {
      warning("`limits.y` must be strictly increasing, e.g., c(ymin, ymax).", call. = FALSE)
    }
    surv  <- survfit_object$surv
    upper <- survfit_object$upper
    lower <- survfit_object$lower
    if (identical(type.y, "risk")) {
      surv  <- 1 - surv
      if (!is.null(upper)) upper <- 1 - upper
      if (!is.null(lower)) lower <- 1 - lower
    }
    if (any(surv < limits.y[1] | surv > limits.y[2], na.rm = TRUE)) {
      warning(sprintf("Some point estimates fall outside `limits.y` = [%.4g, %.4g].",
                      limits.y[1], limits.y[2]), call. = FALSE)
    }
    if (isTRUE(addConfidenceInterval)) {
      if (!is.null(upper) && any(upper < limits.y[1] | upper > limits.y[2], na.rm = TRUE)) {
        warning(sprintf("Some upper CI values fall outside `limits.y` = [%.4g, %.4g].",
                        limits.y[1], limits.y[2]), call. = FALSE)
      }
      if (!is.null(lower) && any(lower < limits.y[1] | lower > limits.y[2], na.rm = TRUE)) {
        warning(sprintf("Some lower CI values fall outside `limits.y` = [%.4g, %.4g].",
                        limits.y[1], limits.y[2]), call. = FALSE)
      }
    }
  }

  check_breaks <- function(bk, nm, lims) {
    if (is.null(bk) || is.function(bk)) return(invisible())
    if (!is.numeric(bk)) {
      warning(sprintf("`%s` should be numeric (or a function).", nm), call. = FALSE); return(invisible())
    }
    if (!is_nondec(bk)) {
      warning(sprintf("`%s` must be non-decreasing.", nm), call. = FALSE)
    }
    if (!is.null(lims) && is_len2_num(lims)) {
      if (any(bk < lims[1] | bk > lims[2], na.rm = TRUE)) {
        warning(sprintf("Some `%s` are outside plotting range [%g, %g].", nm, lims[1], lims[2]), call. = FALSE)
      }
    }
    invisible()
  }
  check_breaks(breaks.x, "breaks.x", limits.x)
  check_breaks(breaks.y, "breaks.y", limits.y)

  type.y <- if (identical(tolower(type.y), "risk")) "risk"
  else if (identical(tolower(type.y), "r")) "risk"
  else if (identical(tolower(type.y), "survival")) "survival"
  else if (identical(tolower(type.y), "s")) "survival"
  else type.y

  label.y <- if (is.null(label.y) && !identical(type.y, "risk") && identical(survfit_object$type, "kaplan-meier")) "Survival"
  else if (is.null(label.y) && identical(type.y, "risk") && identical(survfit_object$type, "kaplan-meier")) "Risk"
  else if (is.null(label.y) && !identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) "Cumulative incidence"
  else if (is.null(label.y) && identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) "1 - cumulative incidence"
  else label.y

  coerce_conf <- function(survfit_object, conf.type) {
    if (!is.null(survfit_object$lower) && !is.null(survfit_object$upper)) return(survfit_object)
    if (conf.type %in% c("none","n") || length(survfit_object$strata) > 2) {
      x <- survfit_object
      x$lower <- x$surv
      x$upper <- x$surv
      return(x)
    }
    return(survfit_object)
  }
  survfit_object <- coerce_conf(survfit_object, conf.type)

  if (identical(style, "MONOCHROME")) options("ggsurvfit.switch-color-linetype" = TRUE)
  if (       !identical(type.y,     "risk") && identical(survfit_object$type, "kaplan-meier")) {
    survfit_object <- ggsurvfit(survfit_object, type = "survival")
    type.y <- "survival"
  } else if ( identical(type.y,     "risk") && identical(survfit_object$type, "kaplan-meier")) {
    survfit_object <- ggsurvfit(survfit_object, type = "risk")
    type.y <- "risk"
  } else if (!identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) {
    survfit_object <- ggsurvfit(survfit_object, type = "risk")
    type.y <- "risk"
  } else if ( identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) {
    survfit_object <- ggsurvfit(survfit_object, type = "survival")
    type.y <- "survival"
  }
  if (identical(style, "MONOCHROME")) options("ggsurvfit.switch-color-linetype" = FALSE)
  return(list(out_survfit_object=survfit_object, label.y=label.y, type.y=type.y))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
extractOptimizationInfo <- function(sol, method) {
  out <- list(solver = method)
  if (method %in% c("nleqslv","Broyden","Newton")) {
    out$code    <- sol$termcd %||% NA_integer_
    out$message <- sol$message %||% NA_character_
    out$fn.evals <- sol$feval %||% NA_integer_
    out$jac.evals <- sol$jeval %||% NA_integer_
  } else if (method == "multiroot") {
    out$message <- sol$estim.precis %||% NA_character_
    out$iter    <- sol$iter %||% NA_integer_
    out$estim.precis <- sol$estim.precis %||% NA_character_
  } else if (method %in% c("optim","SANN","BFGS")) {
    out$code       <- sol$convergence %||% NA_integer_
    out$message    <- sol$message %||% NA_character_
    if (!is.null(sol$counts)) {
      out$fn.evals <- sol$counts["function"] %||% NA_integer_
      out$gr.evals <- sol$counts["gradient"] %||% NA_integer_
    }
  }
  return(out)
}

append_trace <- function(trace_df, iteration, computation.time.second = NA_real_, nleqslv.method, objective.function, relative.difference,
                         max.absolute.difference, converged.by, nleqslv.info, store_params = FALSE,
                         coefficient = NULL) {
  row <- data.frame(
    iteration = iteration,
    computation.time.second = computation.time.second,
    nleqslv.method = nleqslv.method,
    objective.function = objective.function,
    relative.difference = relative.difference,
    max.absolute.difference = max.absolute.difference,
    converged.by = if (isTRUE(converged.by)) NA_character_ else as.character(converged.by),
    code = nleqslv.info$code %||% NA_integer_,
    msg  = nleqslv.info$message %||% NA_character_,
    fn_evals = nleqslv.info$fn.evals %||% NA_integer_,
    gr_evals = nleqslv.info$gr.evals %||% NA_integer_,
    jac_evals = nleqslv.info$jac.evals %||% NA_integer_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(store_params)) {
    row$coefficient <- list(coefficient)
  }
  if (is.null(trace_df)) return(row)
  rbind(trace_df, row)
}



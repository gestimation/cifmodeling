Surv <- function(time, event) {
  if (missing(time))
    stop("Must have a time argument")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  #  if (any(is.na(time)))
  #    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("Must have an event argument")
  if (is.numeric(event)) {
  #  if (any(is.na(event)))
  #    warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
  #  if (any(is.na(event)))
  #    warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
  #  warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
  #  warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  #if (nlevels(as.factor(event)) > 3)
  #  warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Surv"
  return(ss)
}

#' Create a survival or competing risks response
#'
#' A lightweight response constructor used in \code{survival.curve()} and \code{polyreg()}
#' to pass survival and competing risks data via a model formula.
#'
#' @param time Numeric vector of follow-up times (non-negative).
#' @param event Integer or factor vector of event codes
#'   (e.g., 0 = censored, 1 = event of interest, 2 = competing event).
#'
#' @return An object of class \code{"Event"} with fields \code{time}, \code{event}
#' @examples
#' data(diabetes.complications)
#' polyreg(
#'   nuisance.model = Event(t, epsilon) ~ +1,
#'   exposure = "fruitq1",
#'   data = diabetes.complications,
#'   effect.measure1 = "RR",
#'   effect.measure2 = "RR",
#'   time.point = 8,
#'   outcome.type = "COMPETING-RISK"
#' )
#' @export
Event <- function(time, event) {
  if (missing(time))
    stop("Must have a time argument")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  #  if (any(is.na(time)))
  #    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("Must have an event argument")
  if (is.numeric(event)) {
    #  if (any(is.na(event)))
    #    warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    #  if (any(is.na(event)))
    #    warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
    #  warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
    #  warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  #if (nlevels(as.factor(event)) > 3)
  #  warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Event"
  return(ss)
}

createAnalysisDataset <- function(formula, data, other.variables.analyzed=NULL, subset.condition=NULL, na.action=na.pass) {
  if (!is.null(subset.condition)) {
    analysis_dataset <- subset(data, eval(parse(text = subset.condition)))
  } else {
    analysis_dataset <- data
  }
  all_vars <- c(all.vars(formula), other.variables.analyzed)
  analysis_dataset <- analysis_dataset[, all_vars, drop = FALSE]
  return(na.action(analysis_dataset))
}

readStrata <- function(out_readSurv, out_aj, label.strata=NULL) {
  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))
  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
    names(out_aj$strata1) <- label.strata
  }
  return(out_aj)
}

readSurv <- function(formula, data, weights = NULL, code.event1 = 1, code.event2 = 2, code.censoring = 0, subset.condition = NULL, na.action = na.omit) {
  data <- createAnalysisDataset(formula, data, weights, subset.condition, na.action)
  cl <- match.call()
  if (missing(formula))
    stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv"))) {
    stop("A 'Surv' or 'Event' object is expected")
  } else {
    t <- Y[, 1]
    if (any(t<0)) {
      stop("Invalid time variable. Expected non-negative values. ")
    }
    if (!all(Y[, 2] %in% c(code.event1, code.event2, code.censoring))) {
      stop("Invalid event codes. Must be 0 or 1 for survival and 0, 1 or 2 for competing risks, with 0 representing censoring, if event codes are not specified. ")
    } else {
      epsilon <- Y[, 2]
      d <- ifelse(Y[, 2] == code.censoring, 0, 1)
      d0 <- ifelse(Y[, 2] == code.censoring, 1, 0)
      d1 <- ifelse(Y[, 2] == code.event1, 1, 0)
      d2 <- ifelse(Y[, 2] == code.event2, 1, 0)
    }
  }
  if (is.na(all.vars(out_terms)[3])) {
    strata <- rep(1, nrow(data))
    strata_name <- NULL
  } else {
    strata_name <- all.vars(out_terms)[3]
    strata <- as.factor(data[[strata_name]])
  }
  if (is.null(weights)) {
    w <- rep(1, nrow(data))
  } else {
    w <- data[[weights]]
    if (!is.numeric(w))
      stop("Weights must be numeric")
    if (any(!is.finite(w)))
      stop("Weights must be finite")
    if (any(w < 0))
      stop("Weights must be non-negative")
    if (any(is.na(w)))
      stop("Weights contain NA values")
  }
  return(list(t = t, epsilon = epsilon, d = d, d0 = d0, d1 = d1, d2 = d2, strata = strata, strata_name = strata_name, w=w))
}

defineExposureDesign <- function(data, exposure, code.exposure.ref = NULL, prefix = "a") {
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


read_time.point <- function(formula, data, x_a, outcome.type, code.censoring, should.terminate.time.point, time.point) {
  #  read_time.point <- function(formula, data, outcome.type, exposure, code.censoring, code.exposure.ref, time.point) {
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL")) {
    if (is.null(time.point) || !length(time.point)) stop("time.point is required when outcome.type is COMPETING-RISK or SURVIVAL.")
    tp <- suppressWarnings(max(time.point, na.rm = TRUE))
    if (!is.finite(tp) || tp < 0) stop("time.point must be non-negative and finite when outcome.type is COMPETING-RISK or SURVIVAL.")
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

checkInput <- function(data, formula, exposure, code.event1, code.event2, code.censoring, code.exposure.ref, outcome.type, conf.level, report.sandwich.conf, report.boot.conf, nleqslv.method, should.normalize.covariate) {
  cl <- match.call()
  if (missing(formula)) stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL","POLY-PROPORTIONAL","POLY-PROPORTIONAL")) {
    if (!inherits(Y, c("Event", "Surv"))) {
      stop("Surv- or Event-object is expected")
    } else {
      t <- as.numeric(Y[, 1])
      if (any(t<0)) stop("Invalid time variable. Expected non-negative values. ")
      if (any(is.na(t))) stop("Time variable contains NA.")

      epsilon <- as.numeric(Y[, 2])
      if (any(is.na(epsilon))) stop("Event variable contains NA.")
      if (outcome.type == "SURVIVAL") {
        if (!all(epsilon %in% c(code.event1, code.censoring)))
          stop("SURVIVAL requires event codes {censoring,event1}.")
      } else {
        if (!all(epsilon %in% c(code.event1, code.event2, code.censoring)))
          stop("COMPETING-RISK requires event codes {censoring,event1,event2}.")
      }
    }
  }

  out_defineExposureDesign <- defineExposureDesign(data, exposure, code.exposure.ref)
  x_a <- out_defineExposureDesign$x_a
  x_l <- model.matrix(out_terms, mf)

  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be a single number between 0 and 1")

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
  return(list(should.normalize.covariate = should.normalize.covariate.corrected, report.sandwich.conf = report.sandwich.conf.corrected, report.boot.conf = report.boot.conf.corrected, out_defineExposureDesign=out_defineExposureDesign, x_a=x_a, x_l=x_l))
}

check_effect.measure <- function(effect.measure1, effect.measure2) {
  if (effect.measure1 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure1.corrected <- "RR"
  } else if (effect.measure1 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure1.corrected <- "OR"
  } else if (effect.measure1 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure1.corrected <- "SHR"
  } else {
    stop("Invalid input for effect.measure1, Choose 'RR', 'OR', or 'SHR'.")
  }
  if (effect.measure2 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure2.corrected <- "RR"
  } else if (effect.measure2 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure2.corrected <- "OR"
  } else if (effect.measure2 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure2.corrected <- "SHR"
  } else {
    stop("Invalid input for effect.measure2, Choose 'RR', 'OR', or 'SHR'.")
  }
  return(list(effect.measure1 = effect.measure1.corrected, effect.measure2 = effect.measure2.corrected))
}

check_outcome.type <- function(outcome.type) {
  if (outcome.type %in% c("COMPETING-RISK", "COMPETINGRISK", "C", "CR", "COMPETING RISK", "COMPETING-RISKS", "COMPETINGRISKS", "COMPETING RISKS", "Competingrisk", "Competing-risk", "Competing risk", "Competingrisks", "Competing-risks", "Competing risks", "competing-risk", "competingrisk", "competing risk", "competing-risks", "competingrisks", "competing risks")) {
    outcome.type.corrected <- "COMPETING-RISK"
  } else if (outcome.type %in% c("SURVIVAL", "S", "Survival", "Survival")) {
    outcome.type.corrected <- "SURVIVAL"
  } else if (outcome.type %in% c("POLY-PROPORTIONAL", "PP", "Poly-proportional", "poly-proportional")) {
    outcome.type.corrected <- "POLY-PROPORTIONAL"
  } else if (outcome.type %in% c("PROPORTIONAL", "P", "Proportional", "proportional")) {
    outcome.type.corrected <- "PROPORTIONAL"
  } else if (outcome.type %in% c("BINOMIAL", "B", "Binomial", "binomial")) {
    outcome.type.corrected <- "BINOMIAL"
  } else {
    stop("Invalid input for outcome.type, Choose 'COMPETING-RISK', 'SURVIVAL', 'BINOMIAL', 'PROPORTIONAL', or 'POLY-PROPORTIONAL'.")
  }
  return(outcome.type.corrected)
}

check_error <- function(error, outcome.type) {
  if (outcome.type == "SURVIVAL") {
    if (is.null(error)) error <- "greenwood"
  } else {
    if (is.null(error)) error <- "delta"
  }
  if (error %in% c("g", "G", "greenwood", "Greenwood", "GREENWOOD")) error <- "greenwood"
  if (error %in% c("t", "T", "tsiatis", "Tsiatis", "TSIATIS")) error <- "tsiatis"
  if (error %in% c("j", "J", "jackknife", "Jackknife", "JACKKNIFE")) error <- "jackknife"
  if (error %in% c("a", "A", "aalen", "Aalen", "AALEN")) error <- "aalen"
  if (error %in% c("d", "D", "delta", "Delta", "DELTA")) error <- "delta"
  if (outcome.type == "SURVIVAL") {
    if (!error %in% c("greenwood", "tsiatis", "jackknife")) {
      error <- "greenwood"
      warning("Invalid input for standard error for SURVIVAL outcome. Supported options are 'greenwood', 'tsiatis', and 'jackknife'. 'greenwood' is selected. ")
    }
  } else {
    if (!error %in% c("delta", "delta", "jackknife")) {
      error <- "greenwood"
      warning("Invalid input for standard error for COMPETING-RISK outcome. Supported options are 'aalen', 'delta', and 'jackknife'. 'delta' is selected. ")
    }
  }
  return(error)
}

check_label.strata <- function(out_readSurv, label.strata) {
  strata_levels <- levels(as.factor(out_readSurv$strata))
  n_strata <- length(strata_levels)
  if (!is.null(label.strata)) {
    if (!is.character(label.strata)) {
      stop("`label.strata` must be a character vector.", call. = FALSE)
    }
    if (length(label.strata) != n_strata) {
      stop(sprintf("`label.strata` must have length %d (number of strata), but got %d.", n_strata, length(label.strata)), call. = FALSE)
    }
  }
}

untangle.specials <- function(tt, special, order = 1) {
  spc <- attr(tt, "specials")[[special]]
  if (length(spc) == 0)
    return(list(vars = character(0), terms = numeric(0)))
  facs <- attr(tt, "factors")
  fname <- dimnames(facs)
  ff <- apply(facs[spc, , drop = FALSE], 2, sum)
  list(vars = (fname[[1]])[spc], tvar = spc - attr(tt, "response"),
       terms = seq(ff)[ff & match(attr(tt, "order"), order,
                                  nomatch = 0)])
}

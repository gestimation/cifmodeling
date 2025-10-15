#' Create a survival or competing-risks response
#'
#' A lightweight response constructor used in \code{cifcurve()} and \code{polyreg()}
#' to pass survival and competing-risks data via a model formula.
#'
#' @param time Numeric vector of follow-up times (non-negative).
#' @param event Integer (0=censor, 1,2,...) or a character/factor vector whose levels
#'   are numeric codes "0","1","2",... for competing events.
#'
#' @return An object of class \code{"Event"} (a 2-column matrix) with columns \code{time}, \code{event}.
#' @examples
#' ## event: 0=censor, 1=primary, 2=competing
#' # polyreg(nuisance.model = Event(t, epsilon) ~ 1, data = df, outcome.type="COMPETING-RISK")
#'
#' @export
Event <- function(time, event) {
  te <- normalize_time_event(time, event)
  ss <- cbind(time = te$time, event = te$event)
  dimnames(ss) <- list(NULL, c("time","event"))
  attr(ss, "type") <- "right"
  class(ss) <- c("Event", class(ss))
  ss
}

untangle.specials <- function(tt, special, order = 1) {
  spc <- attr(tt, "specials")[[special]]
  if (length(spc) == 0)
    return(list(vars = character(0), terms = numeric(0)))
  facs <- attr(tt, "factors")
  fname <- dimnames(facs)
  ff <- apply(facs[spc, , drop = FALSE], 2, sum)
  list(vars = (fname[[1]])[spc], tvar = spc - attr(tt, "response"),
       terms = seq(ff)[ff & match(attr(tt, "order"), order, nomatch = 0)])
}

createAnalysisDataset <- function(formula,
                                  data,
                                  other.variables.analyzed = NULL,
                                  subset.condition = NULL,
                                  na.action = na.pass,
                                  fill_missing = FALSE) {
  stopifnot(is.data.frame(data))
  if (!is.null(subset.condition)) {
    cond <- if (is.character(subset.condition)) parse(text = subset.condition)[[1]] else subset.condition
    analysis_dataset <- subset(data, eval(cond, envir = data, enclos = parent.frame()))
  } else {
    analysis_dataset <- data
  }

  all_vars <- unique(c(all.vars(formula), other.variables.analyzed))
  missing_cols <- setdiff(all_vars, names(analysis_dataset))

  if (length(missing_cols)) {
    if (isTRUE(fill_missing)) {
      warning(sprintf("The following columns are not in `data` and will be filled with NA: %s",
                      paste(missing_cols, collapse = ", ")))
      for (v in missing_cols) analysis_dataset[[v]] <- NA
    } else {
      stop(sprintf("Undefined columns selected: %s",
                   paste(missing_cols, collapse = ", ")))
    }
  }

  analysis_dataset <- analysis_dataset[, all_vars, drop = FALSE]
  na.action(analysis_dataset)
}

readStrata <- function(out_readSurv, out_aj, label.strata=NULL) {
  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))
  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
    names(out_aj$strata1) <- label.strata
  }
  return(out_aj)
}

readSurv <- function(formula, data, weights = NULL,
                        code.event1 = 1, code.event2 = 2, code.censoring = 0,
                        subset.condition = NULL, na.action = na.omit) {
  stopifnot(is.data.frame(data))
  allowed <- c(code.censoring, code.event1, code.event2)
  data <- createAnalysisDataset(formula, data, weights, subset.condition, na.action)

  Terms <- terms(formula, specials = c("strata","offset","cluster"), data = data)
  mf    <- model.frame(Terms, data = data, na.action = na.action)

  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event","Surv"))) stop("A 'Surv' or 'Event' object is expected")

  te <- normalize_time_event(Y[,1], Y[,2], allowed = allowed)
  t <- te$time
  epsilon <- te$event
  d  <- as.integer(epsilon != code.censoring)
  d0 <- as.integer(epsilon == code.censoring)
  d1 <- as.integer(epsilon == code.event1)
  d2 <- as.integer(epsilon == code.event2)

  mf_rows <- rownames(mf)
  idx <- suppressWarnings(as.integer(mf_rows))
  if (any(is.na(idx))) {
    rn <- rownames(data)
    if (!is.null(rn)) idx <- match(mf_rows, rn)
  }
  if (any(is.na(idx))) stop("Failed to align analysis rows to original data (cannot compute numeric index).")

  data_sync <- data[idx, , drop = FALSE]
  term_labels <- attr(Terms, "term.labels")
  if (length(term_labels) == 0L) {
    strata_name <- NULL
    strata <- factor(rep(1, nrow(mf)))
  } else if (length(term_labels) == 1) {
    strata_name <- term_labels[1]
    strata <- factor(mf[[strata_name]])
  } else {
    strata_name <- paste(term_labels, collapse = ":")
    strata <- interaction(mf[term_labels], drop = TRUE)
  }

  if (is.null(weights)) {
    w <- rep(1, nrow(mf))
  } else if (is.character(weights) && length(weights) == 1) {
    if (!weights %in% names(data)) {
      w <- rep(1, nrow(data))
    }
    else {
      w <- data[[weights]]
      if (!is.numeric(w)) stop("Weights must be numeric.")
      if (any(!is.finite(w))) stop("Weights must be finite.")
      if (any(w < 0)) stop("Weights must be non-negative.")
      if (any(is.na(w))) stop("Weights contain NA values.")
    }
  }
  list(t=t, epsilon=epsilon, d=d, d0=d0, d1=d1, d2=d2, strata=strata, strata_name=strata_name, w=w, data_sync=data_sync)
}

make_competingrisk_marks <- function(
    out_readSurv,
    event_code = 2,
    strata_label_style = c("auto","full","plain"),
    na_rm = TRUE,
    unique = TRUE,
    sort = TRUE
){
  stopifnot(is.list(out_readSurv))
  strata_label_style <- match.arg(strata_label_style)

  if (!is.numeric(out_readSurv$t)) stop("Time vector must be numeric.")
  if (length(out_readSurv$t) != length(out_readSurv$epsilon)) stop("Length of time and event vectors must match.")

  sfac <- out_readSurv$strata
  has_strata <- !is.null(sfac) && is.factor(sfac) && nlevels(sfac) > 1L
  if (has_strata) sfac <- droplevels(sfac)

  clean_vec <- function(v){
    if (na_rm) v <- v[is.finite(v)]
    if (unique) v <- base::unique(v)
    if (sort)   v <- base::sort(v)
    unname(v)
  }

  mk_key <- function(level = NULL) {
    if (!has_strata) return("(all)")
    sname <- out_readSurv$strata_name
    if (identical(strata_label_style, "plain") || is.null(sname)) as.character(level)
    else paste0(sname, "=", level)
  }

  pick <- (out_readSurv$epsilon == event_code)

  if (!has_strata) {
    return( setNames(list(clean_vec(out_readSurv$t[pick])), mk_key()) )
  }

  levs <- levels(sfac)
  lst  <- lapply(levs, function(lv) out_readSurv$t[pick & (sfac == lv)])
  names(lst) <- vapply(levs, mk_key, FUN.VALUE = character(1))
  lapply(lst, clean_vec)
}

readExposureDesign <- function(data, exposure, code.exposure.ref = NULL, prefix = "a") {
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

checkInput <- function(data, formula, exposure, code.event1, code.event2, code.censoring,
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
  if (!is.null(attr(out_terms, "specials")$strata)) stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))  stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster)) stop("cluster() cannot appear in formula")

  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- na.action
  mf <- eval(mf, parent.frame())

  Y <- model.extract(mf, "response")
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL","PROPORTIONAL","POLY-PROPORTIONAL")) {
    if (!inherits(Y, c("Event","Surv"))) stop("A 'Surv' or 'Event' object is expected")
    t <- as.numeric(Y[,1]); epsilon <- as.numeric(Y[,2])
    if (any(t < 0, na.rm = TRUE)) stop("Invalid time variable. Expected non-negative values.")
    if (outcome.type == "SURVIVAL") {
      if (!all(epsilon %in% c(code.event1, code.censoring), na.rm = TRUE))
        stop("SURVIVAL requires event codes {censoring,event1}.")
    } else {
      if (!all(epsilon %in% c(code.event1, code.event2, code.censoring), na.rm = TRUE))
        stop("COMPETING-RISK requires event codes {censoring,event1,event2}.")
    }
  }

  out_readExposureDesign <- readExposureDesign(data, exposure, code.exposure.ref)
  x_a <- out_readExposureDesign$x_a
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
  return(list(should.normalize.covariate = should.normalize.covariate.corrected, report.sandwich.conf = report.sandwich.conf.corrected, report.boot.conf = report.boot.conf.corrected, out_readExposureDesign=out_readExposureDesign, x_a=x_a, x_l=x_l))
}

checkInput_old <- function(data, formula, exposure, code.event1, code.event2, code.censoring, code.exposure.ref, outcome.type, conf.level, report.sandwich.conf, report.boot.conf, nleqslv.method, should.normalize.covariate) {
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

  out_readExposureDesign <- readExposureDesign(data, exposure, code.exposure.ref)
  x_a <- out_readExposureDesign$x_a
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
  return(list(should.normalize.covariate = should.normalize.covariate.corrected, report.sandwich.conf = report.sandwich.conf.corrected, report.boot.conf = report.boot.conf.corrected, out_readExposureDesign=out_readExposureDesign, x_a=x_a, x_l=x_l))
}


#check_outcome.type <- normalize_outcome_type
#check_error <- normalize_error_method
check_effect.measure <- function(effect.measure1, effect.measure2) {
  list(
    effect.measure1 = normalize_effect_measure(effect.measure1, "effect.measure1"),
    effect.measure2 = normalize_effect_measure(effect.measure2, "effect.measure2")
  )
}


check_effect.measure_old <- function(effect.measure1, effect.measure2) {
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

check_outcome.type_old <- function(outcome.type) {
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

check_error_old <- function(error, outcome.type) {
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

.msg <- list(
  req           = "`{arg}` is required.",
  numeric       = "`{arg}` must be numeric.",
  nonneg        = "`{arg}` must be non-negative.",
  na            = "`{arg}` contains NA.",
  len_mismatch  = "`{x}` and `{y}` must have the same length (got {nx} vs {ny}).",
  ev_type       = "`event` must be numeric, logical, factor, or character.",
  ev_codes      = "`event` codes must be non-negative integers: {allowed}. Found: {found}.",
  weights_num   = "`weights` must be numeric.",
  weights_fin   = "`weights` must be finite.",
  weights_pos   = "`weights` must be non-negative.",
  outcome_type  = "Invalid input for outcome.type. Choose one of: {choices}.",
  effect_meas   = "Invalid input for {which}. Choose 'RR', 'OR', or 'SHR'.",
  error_surv    = "Invalid SE method for SURVIVAL. Use 'greenwood','tsiatis','jackknife'. Defaulting to 'greenwood'.",
  error_cr      = "Invalid SE method for COMPETING-RISK. Use 'aalen','delta','jackknife'. Defaulting to 'delta'."
)

.err <- function(key, ..., .class = NULL, .messages = .msg) {
  args <- list(...)

  tmpl <- .messages[[key]]
  if (is.null(tmpl)) {
    stop(sprintf("Unknown error key: %s", key), call. = FALSE)
  }

  msg <- tmpl
  if (length(args)) {
    for (nm in names(args)) {
      pat <- paste0("\\{", nm, "\\}")
      msg <- gsub(pat, as.character(args[[nm]]), msg, perl = TRUE)
    }
  }

  cls <- unique(c(.class, paste0("cifmodeling_", key), "cifmodeling_error"))
  cond <- structure(
    list(message = msg, call = sys.call(-1L)),
    class = c(cls, "simpleError", "error", "condition")
  )
  stop(cond)
}

.chk_numeric_nonneg <- function(x, name) {
  if (!is.numeric(x)) .err("numeric", arg = name)
  if (anyNA(x))       .err("na",      arg = name)
  if (any(x < 0))     .err("nonneg",  arg = name)
  invisible(TRUE)
}

normalize_time_event <- function(time, event, allowed = NULL) {
  if (missing(time))  .err("req", arg = "time")
  if (missing(event)) .err("req", arg = "event")
  if (!is.numeric(time)) .err("numeric", arg = "time")
  if (any(time < 0, na.rm = TRUE)) .err("nonneg", arg = "time")

  if (is.numeric(event)) {
    if (any(event < 0, na.rm = TRUE) || any(event != floor(event), na.rm = TRUE)) {
      .err("ev_codes", allowed = "{0,1,2,...}",
           found = paste(unique(event[!is.na(event)]), collapse = ", "))
    }
    status <- suppressWarnings(as.integer(event))
  } else if (is.logical(event)) {
    status <- ifelse(is.na(event), NA_integer_, as.integer(event))  # FALSE=0, TRUE=1
  } else if (is.factor(event) || is.character(event)) {
    ev_chr <- as.character(event)  # NA は NA のまま
    ok <- !is.na(ev_chr)
    if (!all(grepl("^[0-9]+$", ev_chr[ok]))) {
      .err("ev_codes", allowed = "'0','1','2',...",
           found = paste(unique(ev_chr[ok]), collapse = ", "))
    }
    status <- suppressWarnings(as.integer(ev_chr))  # "NA" 等は NA になる
  } else {
    .err("ev_type")
  }

  if (length(status) != length(time)) {
    .err("len_mismatch", x = "time", y = "event",
         nx = length(time), ny = length(status))
  }

  if (!is.null(allowed)) {
    ok <- is.na(status) | status %in% allowed
    if (!all(ok)) {
      .err("ev_codes",
           allowed = paste0("{", paste(allowed, collapse = ","), "}"),
           found   = paste(sort(unique(status[!ok])), collapse = ", "))
    }
  }
  list(time = as.numeric(time), event = status)
}


normalize_time_event_old <- function(time, event, allowed = NULL) {
  if (missing(time))  .err("req", arg = "time")
  if (missing(event)) .err("req", arg = "event")
  .chk_numeric_nonneg(time, "time")
  status <- NULL
  if (is.numeric(event)) {
    if (anyNA(event)) .err("na", arg = "event")
    if (any(event < 0) || any(event != floor(event))) {
      .err("ev_codes", allowed = "{0,1,2,...}", found = paste(unique(event), collapse=", "))
    }
    status <- as.integer(event)
  } else if (is.logical(event)) {
    if (anyNA(event)) .err("na", arg = "event")
    status <- as.integer(event)  # FALSE=0/TRUE=1
  } else if (is.factor(event) || is.character(event)) {
    ev_chr <- as.character(event)
    if (anyNA(ev_chr)) .err("na", arg = "event")
    if (!all(grepl("^[0-9]+$", ev_chr))) {
      .err("ev_codes", allowed = "'0','1','2',...", found = paste(unique(ev_chr), collapse=", "))
    }
    status <- as.integer(ev_chr)
  } else {
    .err("ev_type")
  }

  if (length(status) != length(time)) {
    .err("len_mismatch", x = "time", y = "event", nx = length(time), ny = length(status))
  }
  if (!is.null(allowed) && !all(status %in% allowed)) {
    .err("ev_codes", allowed = paste0("{", paste(allowed, collapse=","), "}"),
         found = paste(sort(unique(status)), collapse=", "))
  }
  list(time = as.numeric(time), event = status)
}

check_weights <- function(w) {
  if (!is.numeric(w))     .err("weights_num")
  if (any(!is.finite(w))) .err("weights_fin")
  if (any(w < 0))         .err("weights_pos")
  if (anyNA(w))           .err("na", arg = "weights")
  invisible(TRUE)
}

normalize_effect_measure <- function(x, which = "effect.measure") {
  ux <- toupper(x)
  if (ux %in% c("RR","RISK RATIO")) return("RR")
  if (ux %in% c("OR","ODDS RATIO")) return("OR")
  if (ux %in% c("SHR","HR","SUBDISTRIBUTION HAZARD RATIO")) return("SHR")
  .err("effect_meas", which = which)
}

check_outcome.type <- function(x) {
  map <- list(
    "COMPETING-RISK"   = c("competing-risk","competing risk","competingrisks","competing-risks","cr","c"),
    "SURVIVAL"         = c("survival","s"),
    "POLY-PROPORTIONAL"= c("poly-proportional","pp"),
    "PROPORTIONAL"     = c("proportional","p"),
    "BINOMIAL"         = c("binomial","b")
  )
  ux <- toupper(gsub("[[:space:]]+", " ", x))
  for (k in names(map)) {
    if (ux == k || tolower(ux) %in% tolower(map[[k]])) return(k)
  }
  .err("outcome_type", choices = paste(names(map), collapse = ", "))
}

check_error <- function(x, outcome.type) {
  ot <- toupper(as.character(outcome.type))
  out <- if (is.null(x)) if (ot == "SURVIVAL") "greenwood" else "delta" else tolower(x)

  if (ot == "SURVIVAL") {
    if (!out %in% c("greenwood", "tsiatis", "jackknife")) {
      warning(.msg$error_surv, call. = FALSE); out <- "greenwood"
    }
  } else if (ot == "COMPETING-RISK") {
    if (!out %in% c("aalen", "delta", "jackknife")) {
      warning(.msg$error_cr, call. = FALSE); out <- "delta"
    }
  } else {
    stop(sprintf("Invalid outcome.type: %s", outcome.type), call. = FALSE)
  }
  out
}


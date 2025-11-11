.msg <- list(
  align_rows_fail = "Failed to align analysis rows to original data (cannot compute numeric index).",
  codes_required_cr = "COMPETING-RISK requires event codes {censoring,event1,event2} (got: {found}).",
  codes_required_surv = "SURVIVAL requires event codes {censoring,event1} (got: {found}).",
  code_events_len_vec = "code.events must be a numeric length-3 vector c(event1, event2, censoring).",
  code_events_integer = "code.events entries must be integer-like (whole numbers).",
  code_events_distinct = "code.events[1]` and `code.events[2]` must be different event codes.",
  code_events_len_surv = "code.events[[{i}]] must be c(event.code1, censoring) for SURVIVAL.",
  code_events_len_cr = "code.events[[{i}]] must be c(event.code1, event.code2, censoring) for COMPETING-RISK.",
  conf_level = "conf.level must be a single number in (0, 1).",
  effect_meas   = "Invalid input for {which}. Choose 'RR', 'OR', or 'SHR'.",
  error_cr      = "Invalid SE method for COMPETING-RISK. Use aalen, delta, jackknife. Defaulting to delta.",
  error_surv    = "Invalid SE method for SURVIVAL. Use greenwood, tsiatis, jackknife. Defaulting to greenwood.",
  est_outside_limits_y = "Some point estimates fall outside `{arg}` = [{a}, {b}].",
  ev_codes      = "event codes must be non-negative integers: {allowed}. Found: {found}.",
  ev_type       = "event must be numeric, logical, factor, or character.",
  formula_must_be = "Each formula must be a character string or a formula object.",
  finite = "{arg} must be finite.",
  incompatible_flags = "{which} cannot be used together.",
  inset_need_two = "inset.panel=TRUE requires at least two plots.",
  inset_extra_drop = "inset.panel=TRUE: only the first two plots will be used.",
  infer_outcome_fail = "Failed to infer outcome.type from code.events; each must be length 2 (S) or 3 (C).",
  na            = "{arg} contains NA.",
  nonneg        = "{arg} must be non-negative.",
  no_cluster_in_formula = "cluster() cannot appear in formula.",
  no_offset_in_formula = "offset() cannot appear in formula.",
  no_strata_in_formula = "strata() cannot appear in formula.",
  need_code_events = "Provide non-empty `code.events` as a list of c(event1, event2[, censor]) per panel.",
  need_data = "`data` must be provided.",
  need_formula_or_formulas = "Provide either formula (single) or formulas (list).",
  numeric       = "`{arg}` must be numeric.",
  outcome_type  = "Invalid input for outcome.type. Choose one of: {choices}.",
  panel_disables_tables = "add.risktable/add.estimate.table are ignored in panel mode (panel.per.variable/panel.per.event/cifpanel); set to FALSE internally.",
  panel_disables_labelstrata = "label.strata is ignored in panel mode (panel.per.variable/panel.per.event/cifpanel).",
  recycle_empty = "Cannot recycle an empty object to nonzero length.",
  req           = "`{arg}` is required.",
  style         = "Invalid input for style. Choose one of: {choices}.",
  surv_expected = "A 'Surv' or 'Event' object is expected.",
  timepoint_required = "`time.point` is required when outcome.type is COMPETING-RISK or SURVIVAL.",
  timepoint_nonneg_finite = "`time.point` must be a finite non-negative number.",
  time_nonneg = "`{arg}` must be non-negative.",
  limits_len2 = "`{arg}` must be a numeric length-2 vector.",
  limits_increasing = "`{arg}` must be strictly increasing, e.g., c(min, max).",
  limits_x_outside = "Max of `survfit_object$time` ({tmax}) lies outside `{arg}` = [{a}, {b}].",
  len_mismatch  = "`{x}` and `{y}` must have the same length (got {nx} vs {ny}).",
  upper_outside_limits_y = "Some upper CI values fall outside `{arg}` = [{a}, {b}].",
  lower_outside_limits_y = "Some lower CI values fall outside `{arg}` = [{a}, {b}].",
  ors_tmax_bad = "`out_readSurv$t` provided but has no finite positive max; x-limits fallback may fail.",
  both_formula_forms = "Both `formula` and `formulas` provided. Using `formulas` only.",
  plots_extra_dropped = "There are {n_plots} plots but grid holds {n_slots}. Extra plots are dropped.",
  shape_identical = "{a} and {b} specify an identical type of symbol.",
  order_strata_type = "`order.strata` must be a character vector or a named list.",
  need_formula_for_panel.per.variable = "panel.per.variable=TRUE requires a formula interface (not a survfit object).",
  need_rhs_vars_for_panel.per.variable = "panel.per.variable=TRUE requires one or more variables on the RHS.",
  no_transform_for_panel.per.variable = "panel.per.variable=TRUE does not support transformations/interactions. Remove: {which}.",
  need_formula_or_formulas = "Provide a formula or a fitted survfit-like object.",
  weights_num   = "weights must be numeric.",
  weights_fin   = "weights must be finite.",
  weights_pos   = "weights must be non-negative."
)

.warn <- function(key, ..., .messages = NULL) {
  if (is.null(.messages)) {
    .messages <- tryCatch(get(".msg", inherits = TRUE), error = function(e) NULL)
  }
  tmpl <- if (!is.null(.messages)) .messages[[key]] else NULL
  if (is.null(tmpl)) {
    warning(sprintf("Unknown warning key: %s", key), call. = FALSE)
    return(invisible(FALSE))
  }
  args <- list(...)
  msg <- tmpl
  if (length(args)) {
    for (nm in names(args)) {
      pat <- paste0("\\{", nm, "\\}")
      msg <- gsub(pat, as.character(args[[nm]]), msg, perl = TRUE)
    }
  }
  warning(msg, call. = FALSE)
  invisible(TRUE)
}


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

.assert <- function(cond, key, ...) {
  if (!isTRUE(cond)) .err(key, ...)
  invisible(TRUE)
}

.assert_choice <- function(x, choices, key, ...) {
  .assert(length(x) == 1L && is.character(x), key, ...)
  .assert(x %in% choices, key, choices = paste(choices, collapse = ", "), ...)
  x
}

.assert_limits <- function(x, arg) {
  .assert(is.numeric(x) && length(x) == 2L && all(is.finite(x)), "limits_len2", arg = arg)
  .assert(x[1L] < x[2L], "limits_increasing", arg = arg)
  invisible(TRUE)
}

util_check_outcome_type <- function(
    x = NULL,
    formula = NULL,
    data = NULL,
    na.action = stats::na.omit,
    auto_message = TRUE,
    allow_multiple = FALSE
) {
  normalize_na_action <- function(na) {
    if (is.function(na)) return(na)
    if (is.character(na) && length(na) == 1L) {
      if (na %in% c("na.omit", "na.exclude", "na.pass", "na.fail")) {
        return(get(na, envir = asNamespace("stats"), inherits = FALSE))
      }
    }
    stats::na.omit
  }
  na.action <- normalize_na_action(na.action)

  ok <- c("competing-risk", "survival", "binomial",
          "proportional-competing-risk", "proportional-survival")

  normalize <- function(z) {
    if (is.null(z)) return(character())
    z <- as.character(z)
    z <- tolower(trimws(z))
    z <- gsub("[[:space:]_.]+", "-", z)
    z <- gsub("-{2,}", "-", z)
    z <- gsub("^-|-$", "", z)
    map <- c(
      "c"   = "competing-risk", "cr"  = "competing-risk", "competing" = "competing-risk",
      "s"   = "survival",       "surv" = "survival",
      "b"   = "binomial",       "bin" = "binomial",
      "pc"  = "proportional-competing-risk", "pcr" = "proportional-competing-risk",
      "ps"  = "proportional-survival"
    )
    ifelse(z %in% names(map), unname(map[z]), z)
  }

  if (!is.null(x)) {
    canon <- normalize(x)
    invalid <- is.na(canon) | !nzchar(canon) | !(canon %in% ok)
    if (any(invalid)) {
      bad <- unique(canon[invalid])
      bad <- bad[!is.na(bad)]
      stop(sprintf("Invalid outcome.type: '%s'. Allowed: %s",
                   paste(bad, collapse = ", "), paste(ok, collapse = ", ")),
           call. = FALSE)
    }
    u <- unique(canon)
    if (length(u) == 1L) return(u)
    if (isTRUE(allow_multiple)) return(canon)
    stop("`outcome.type` is ambiguous and matched multiple types: ",
         paste(u, collapse = ", "), call. = FALSE)
  }

  if (is.null(formula) || is.null(data)) {
    stop("Provide either `outcome.type` or both `formula` and `data` for automatic detection.", call. = FALSE)
  }

  Terms <- stats::terms(formula, specials = c("strata", "offset", "cluster"), data = data)
  mf <- tryCatch(
    stats::model.frame(Terms, data = data, na.action = na.action),
    error = function(e) stats::model.frame(Terms, data = data, na.action = stats::na.omit)
  )

  Y <- stats::model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv")))
    stop("Response must be Event() or Surv().", call. = FALSE)

  status <- suppressWarnings(as.numeric(Y[, 2]))
  n_levels <- length(unique(stats::na.omit(status)))
  if (n_levels > 2L) {
    if (isTRUE(auto_message)) message("Detected >2 status levels; outcome.type set to 'competing-risk'.")
    "competing-risk"
  } else {
    if (isTRUE(auto_message)) message("Detected <= 2 status levels; outcome.type set to 'survival'.")
    "survival"
  }
}

util_check_type_y <- function(x = NULL) {
  .norm <- function(s) gsub("[^A-Z0-9]+", "", toupper(trimws(as.character(s))))

  if (!is.null(x)) {
    buckets <- list(
      surv = c("surv", "s", "survival", "km", "kaplan-meier", "kaplanmeier",
               "survival-probability", "survivalprobability"),
      risk = c("risk", "r", "cif", "ci", "cumulative-incidence", "cumulativeincidence",
               "cuminc", "failure", "incidence")
    )

    ali_norm <- lapply(buckets, .norm)
    names(ali_norm) <- names(buckets)
    alias_rev <- stats::setNames(
      rep(names(ali_norm), lengths(ali_norm)),
      unlist(ali_norm, use.names = FALSE)
    )

    ux <- .norm(x)
    ux <- ux[!is.na(ux) & nzchar(ux)]

    canon <- unique(stats::na.omit(vapply(
      ux,
      function(u) {
        if (u %in% .norm(names(buckets))) {
          idx <- match(u, .norm(names(buckets)))
          names(buckets)[idx]
        } else {
          cn <- unname(alias_rev[u])
          if (length(cn) == 0L || is.na(cn)) NA_character_ else cn
        }
      },
      FUN.VALUE = character(1)
    )))

    if (length(canon) == 1L) {
      return(canon)
    } else if (length(canon) > 1L) {
      stop("`type.y` is ambiguous and matched multiple types: ",
           paste(canon, collapse = ", "), call. = FALSE)
    } else {
      allowed <- sort(unique(c(
        names(buckets),
        unlist(buckets, use.names = FALSE)
      )))
      stop("Invalid `type.y`. Allowed values are: ",
           paste(allowed, collapse = ", "), call. = FALSE)
    }
    return(canon)
  }
  return(NULL)
}

check_weights <- function(w) {
  if (!is.numeric(w))     .err("weights_num")
  if (any(!is.finite(w))) .err("weights_fin")
  if (any(w < 0))         .err("weights_pos")
  if (anyNA(w))           .err("na", arg = "weights")
  invisible(TRUE)
}

reg_check_effect.measure <- function(effect.measure1, effect.measure2) {
  normalize_effect_measure <- function(x, which = "effect.measure") {
    ux <- toupper(x)
    if (ux %in% c("RR","RISK RATIO")) return("RR")
    if (ux %in% c("OR","ODDS RATIO")) return("OR")
    if (ux %in% c("SHR","HR","SUBDISTRIBUTION HAZARD RATIO")) return("SHR")
    .err("effect_meas", which = which)
  }
  list(
    effect.measure1 = normalize_effect_measure(effect.measure1, "effect.measure1"),
    effect.measure2 = normalize_effect_measure(effect.measure2, "effect.measure2")
  )
}

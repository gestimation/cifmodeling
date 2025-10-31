.msg <- list(
  panel_disables_tables = "addRiskTable/addEstimateTable are ignored in panel mode (printEachVar/printEachEvent); set to FALSE internally.",
  panel_disables_labelstrata = "`label.strata` is ignored in panel mode (printEachVar/printEachEvent).",
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
  style         = "Invalid input for style. Choose one of: {choices}.",
  effect_meas   = "Invalid input for {which}. Choose 'RR', 'OR', or 'SHR'.",
  error_surv    = "Invalid SE method for SURVIVAL. Use 'greenwood','tsiatis','jackknife'. Defaulting to 'greenwood'.",
  error_cr      = "Invalid SE method for COMPETING-RISK. Use 'aalen','delta','jackknife'. Defaulting to 'delta'.",
  surv_expected = "A 'Surv' or 'Event' object is expected.",
  time_nonneg = "`{arg}` must be non-negative.",
  timepoint_required = "`time.point` is required when outcome.type is COMPETING-RISK or SURVIVAL.",
  timepoint_nonneg_finite = "`time.point` must be a finite non-negative number.",
  conf_level = "`conf.level` must be a single number in (0, 1).",
  no_strata_in_formula = "`strata()` cannot appear in formula.",
  no_offset_in_formula = "`offset()` cannot appear in formula.",
  no_cluster_in_formula = "`cluster()` cannot appear in formula.",
  align_rows_fail = "Failed to align analysis rows to original data (cannot compute numeric index).",
  codes_required_surv = "SURVIVAL requires event codes {censoring,event1} (got: {found}).",
  codes_required_cr = "COMPETING-RISK requires event codes {censoring,event1,event2} (got: {found}).",
  need_data = "`data` must be provided.",
  need_code_events = "Provide non-empty `code.events` as a list of c(event1, event2[, censor]) per panel.",
  need_formula_or_formulas = "Provide either `formula` (single) or `formulas` (list).",
  formula_must_be = "Each formula must be a character string or a formula object.",
  recycle_empty = "Cannot recycle an empty object to nonzero length.",
  limits_len2 = "`{arg}` must be a numeric length-2 vector.",
  limits_increasing = "`{arg}` must be strictly increasing, e.g., c(min, max).",
  limits_x_outside = "Max of `survfit_object$time` ({tmax}) lies outside `{arg}` = [{a}, {b}].",
  est_outside_limits_y = "Some point estimates fall outside `{arg}` = [{a}, {b}].",
  upper_outside_limits_y = "Some upper CI values fall outside `{arg}` = [{a}, {b}].",
  lower_outside_limits_y = "Some lower CI values fall outside `{arg}` = [{a}, {b}].",
  ors_tmax_bad = "`out_readSurv$t` provided but has no finite positive max; x-limits fallback may fail.",
  both_formula_forms = "Both `formula` and `formulas` provided. Using `formulas` only.",
  inset_need_two = "`use_inset_element=TRUE` requires at least two plots.",
  inset_extra_drop = "`use_inset_element=TRUE`: only the first two plots will be used.",
  plots_extra_dropped = "There are {n_plots} plots but grid holds {n_slots}. Extra plots are dropped.",
  code_events_len_vec = "`code.events` must be a numeric length-3 vector c(event1, event2, censoring).",
  code_events_integer = "`code.events` entries must be integer-like (whole numbers).",
  code_events_distinct = "`code.events[1]` and `code.events[2]` must be different event codes.",
  code_events_len_surv = "`code.events[[{i}]]` must be c(event.code1, censoring) for SURVIVAL.",
  code_events_len_cr = "`code.events[[{i}]]` must be c(event.code1, event.code2, censoring) for COMPETING-RISK.",
  infer_outcome_fail = "Failed to infer outcome.type from code.events; each must be length 2 (S) or 3 (C).",
  shape_identical = "`{a}` and `{b}` specify an identical type of symbol.",
  finite = "`{arg}` must be finite.",
  incompatible_flags = "`{which}` cannot be used together.",
  need_formula_for_printEachVar = "printEachVar=TRUE requires a formula interface (not a survfit object).",
  need_rhs_vars_for_printEachVar = "printEachVar=TRUE requires one or more variables on the RHS.",
  no_transform_for_printEachVar = "printEachVar=TRUE does not support transformations/interactions. Remove: {which}.",
  order_strata_type = "`order.strata` must be a character vector or a named list.",
  need_formula_or_formulas = "Provide a formula or a fitted survfit-like object."
)

.warn <- function(key, ..., .messages = .msg) {
  tmpl <- .messages[[key]]
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
    auto_message = TRUE
) {
  if (!is.null(x)) {
    map <- list(
      "COMPETING-RISK"    = c("competing-risk","competing risk","competingrisks","competing-risks","cr","c"),
      "SURVIVAL"          = c("survival","s"),
      "PROPORTIONAL-COMPETING-RISK" = c("proportional-competing-risk","pc","pcr"),
      "PROPORTIONAL-SURVIVAL"       = c("proportional-survival","ps"),
      "BINOMIAL"          = c("binomial","b")
    )
    alias_rev <- {
      ali <- lapply(map, toupper)
      stats::setNames(rep(names(ali), lengths(ali)), unlist(ali, use.names = FALSE))
    }

    ux <- toupper(trimws(as.character(x)))
    ux <- ux[!is.na(ux) & nzchar(ux)]

    canon <- unique(stats::na.omit(vapply(
      ux,
      function(u) {
        if (u %in% names(map)) return(u)
        cn <- alias_rev[[u]]
        if (is.null(cn)) NA_character_ else cn
      },
      FUN.VALUE = character(1)
    )))

    if (length(canon) == 1L) {
      return(canon)
    } else {
      x <- NULL
    }
  }

  if (is.null(formula) || is.null(data))
    .err("req", arg = "formula and data (for automatic detection)")

  Terms <- stats::terms(formula, specials = c("strata","offset","cluster"), data = data)
  mf    <- stats::model.frame(Terms, data = data, na.action = na.action)

  Y <- stats::model.extract(mf, "response")
  if (!inherits(Y, c("Event","Surv"))) .err("surv_expected")

  status <- suppressWarnings(as.numeric(Y[, 2]))
  n_levels <- length(unique(stats::na.omit(status)))

  if (n_levels > 2L) {
    if (isTRUE(auto_message)) message("Detected >2 status levels; outcome.type set to 'COMPETING-RISK'.")
    return("COMPETING-RISK")
  } else {
    if (isTRUE(auto_message)) message("Detected <= 2 status levels; outcome.type set to 'SURVIVAL'.")
    return("SURVIVAL")
  }
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

#plot_check_style <- function(x) {
#  map <- list(
#    "CLASSIC"     = c("classic","c"),
#    "BOLD"        = c("bold","b"),
#    "FRAMED"      = c("framed","f"),
#    "GRID"        = c("grid","g"),
#    "GRAY"        = c("gray"),
#    "GGSURVFIT"   = c("ggsurvfit")
#  )
#  ux <- toupper(gsub("[[:space:]]+", " ", x))
#  for (k in names(map)) {
#    if (ux == k || tolower(ux) %in% tolower(map[[k]])) return(k)
#  }
#  .err("style", choices = paste(names(map), collapse = ", "))
#}

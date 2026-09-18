`%||%` <- function(x, y) if (is.null(x)) y else x

util_resolve_weights <- function(weights.expr, data, envir = parent.frame(), missing = FALSE) {
  if (isTRUE(missing) || identical(weights.expr, quote(NULL))) return(NULL)
  if (is.name(weights.expr)) {
    nm <- as.character(weights.expr)
    if (nm %in% names(data)) return(nm)
  }
  eval(weights.expr, envir = envir)
}

createAnalysisDataset <- function(formula,
                                  data,
                                  other.variables.analyzed = NULL,
                                  subset.condition = NULL,
                                  na.action = na.pass,
                                  fill_missing = FALSE,
                                  return.info = FALSE) {
  stopifnot(is.data.frame(data))
  stopifnot(is.function(na.action))

  index <- rep_len(TRUE, nrow(data))
  if (!is.null(subset.condition)) {
    if (is.logical(subset.condition)) {
      if (length(subset.condition) != nrow(data))
        stop("`subset.condition` logical length must equal nrow(data).")
      index <- subset.condition & !is.na(subset.condition)
    } else {
      cond_expr <- if (inherits(subset.condition, "formula")) {
        if (length(subset.condition) != 2)
          stop("Use a one-sided formula like `~ condition` for `subset.condition`.")
        subset.condition[[2]]
      } else if (is.character(subset.condition)) {
        parse(text = subset.condition)[[1]]
      } else if (is.expression(subset.condition)) {
        subset.condition[[1]]
      } else if (is.language(subset.condition)) {
        subset.condition
      } else stop("Unsupported `subset.condition` type.")
      val <- eval(cond_expr, envir = data, enclos = parent.frame())
      if (!is.logical(val)) stop("Evaluated `subset.condition` is not logical.")
      val[is.na(val)] <- FALSE
      index <- val
    }
  }
  row_id_name <- ".cifmodeling_row_id"
  while (row_id_name %in% names(data)) row_id_name <- paste0(row_id_name, "_")
  data[[row_id_name]] <- seq_len(nrow(data))
  analysis_dataset <- data[index, , drop = FALSE]

  all_vars <- unique(c(all.vars(formula), other.variables.analyzed))
  missing_cols <- setdiff(all_vars, names(analysis_dataset))

  if (length(missing_cols)) {
    if (isTRUE(fill_missing)) {
      warning(sprintf("The following columns are not in `data` and will be filled with NA: %s", paste(missing_cols, collapse = ", ")))
      for (v in missing_cols) analysis_dataset[[v]] <- NA
    } else {
      stop(sprintf("Undefined columns selected: %s", paste(missing_cols, collapse = ", ")))
    }
  }
  analysis_dataset <- analysis_dataset[, c(all_vars, row_id_name), drop = FALSE]
  analysis_dataset <- na.action(analysis_dataset)

  row.index <- as.integer(analysis_dataset[[row_id_name]])
  analysis_dataset[[row_id_name]] <- NULL
  if (!isTRUE(return.info)) return(analysis_dataset)

  included <- rep_len(FALSE, nrow(data))
  included[row.index] <- TRUE
  exclusion.reason <- rep_len(NA_character_, nrow(data))
  exclusion.reason[!index] <- "subset"
  exclusion.reason[index & !included] <- "missing"

  list(
    data = analysis_dataset,
    row.index = row.index,
    included = included,
    exclusion.reason = exclusion.reason,
    subset.index = index
  )
}

util_get_surv <- function(
    predicted.time,
    estimated.surv,
    estimated.time,
    predicted.strata = NULL,
    estimated.strata = NULL,
    strata.levels = NULL
){
  if (anyNA(predicted.time)) stop("Invalid predicted.time: contains NA.")
  if (length(estimated.surv) != length(estimated.time))
    stop("estimated.surv and estimated.time must have the same length.")

  prepareSeries <- function(time_vec, surv_vec) {
    ok <- !(is.na(time_vec) | is.na(surv_vec))
    time_vec <- time_vec[ok]; surv_vec <- surv_vec[ok]
    if (!length(time_vec)) return(list(t = numeric(0), s = numeric(0)))
    o <- order(time_vec)
    t2 <- time_vec[o]; s2 <- surv_vec[o]
    keep <- !duplicated(t2, fromLast = TRUE)
    list(t = t2[keep], s = s2[keep])
  }
  n_pred <- length(predicted.time)
  predicted.surv <- numeric(n_pred)

  strata_mode <- !(
    is.null(predicted.strata) || is.null(estimated.strata) || is.null(strata.levels) ||
      length(estimated.strata) == 0L || length(strata.levels) == 0L
  )
  if (!strata_mode) {
    ser <- prepareSeries(estimated.time, estimated.surv)
    if (!length(ser$t)) return(rep(1.0, n_pred))
    for (i in seq_len(n_pred)) {
      idx <- findInterval(predicted.time[i], ser$t, left.open = TRUE)
      predicted.surv[i] <- if (idx > 0L) ser$s[idx] else 1.0
    }
    return(predicted.surv)
  }

  if (!is.numeric(estimated.strata) || any(estimated.strata < 0))
    stop("'estimated.strata' must be a non-negative numeric vector of counts.")
  if (sum(estimated.strata) != length(estimated.time))
    stop("sum(estimated.strata) must equal length(estimated.time).")

  K <- length(estimated.strata)
  if (length(strata.levels) != K)
    stop("'strata.levels' must have length K = length(estimated.strata).")

  if (length(predicted.strata) == 1L) {
    predicted.strata <- rep(predicted.strata, n_pred)
  } else if (length(predicted.strata) != n_pred) {
    stop("Length of predicted.strata must be 1 or match length(predicted.time).")
  }

  mapped <- if (is.factor(predicted.strata)) {
    match(as.character(predicted.strata), as.character(strata.levels))
  } else {
    match(predicted.strata, strata.levels)
  }
  if (any(is.na(mapped))) {
    bad <- unique(predicted.strata[is.na(mapped)])
    stop("Some values in predicted.strata are not found in 'strata.levels': ",
         paste(bad, collapse = ", "))
  }

  cs <- cumsum(estimated.strata)
  strata_start <- c(1L, cs[-K] + 1L)
  strata_end   <- cs

  series_per_stratum <- vector("list", K)
  for (s in seq_len(K)) {
    if (estimated.strata[s] == 0L) {
      series_per_stratum[[s]] <- list(t = numeric(0), s = numeric(0))
    } else {
      idx <- strata_start[s]:strata_end[s]
      series_per_stratum[[s]] <- prepareSeries(estimated.time[idx], estimated.surv[idx])
    }
  }

  for (i in seq_len(n_pred)) {
    s <- mapped[i]
    ser <- series_per_stratum[[s]]
    if (!length(ser$t)) {
      predicted.surv[i] <- 1.0
    } else {
      j <- findInterval(predicted.time[i], ser$t, left.open = TRUE)
      predicted.surv[i] <- if (j > 0L) ser$s[j] else 1.0
    }
  }

  predicted.surv
}

util_validate_event_codes <- function(
    code.event1,
    code.event2,
    code.censoring,
    outcome.type = NULL
) {
  codes <- c(code.event1, code.event2, code.censoring)
  if (!is.numeric(codes) || length(codes) != 3L || anyNA(codes) ||
      any(!is.finite(codes)) || any(codes < 0) || any(codes != floor(codes)) ||
      anyDuplicated(codes)) {
    stop(
      "`code.event1`, `code.event2`, and `code.censoring` must be distinct non-negative integers.",
      call. = FALSE
    )
  }
  as.integer(codes)
}

cif_prepare_input <- function(
    formula,
    data,
    weights = NULL,
    other.variables.analyzed = NULL,
    subset.condition = NULL,
    na.action = stats::na.omit,
    outcome.type = NULL,
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    auto_message = TRUE,
    validate.observed.codes = TRUE
) {
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.", call. = FALSE)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)

  codes <- util_validate_event_codes(
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring,
    outcome.type = outcome.type
  )
  code.event1 <- codes[1L]
  code.event2 <- codes[2L]
  code.censoring <- codes[3L]
  na.action <- util_normalize_na_action(na.action)

  data_work <- data
  weight_name <- NULL
  if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data_work)) {
      stop("weights = '", weights, "' is not found in data.", call. = FALSE)
    }
    weight_name <- weights
  } else if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(data_work)) {
      stop("Numeric `weights` must have length nrow(data).", call. = FALSE)
    }
    weight_name <- ".cifmodeling_weights"
    while (weight_name %in% names(data_work)) weight_name <- paste0(weight_name, "_")
    data_work[[weight_name]] <- weights
  }

  other_vars <- unique(c(other.variables.analyzed, weight_name))
  prep <- createAnalysisDataset(
    formula = formula,
    data = data_work,
    other.variables.analyzed = other_vars,
    subset.condition = subset.condition,
    na.action = na.action,
    return.info = TRUE
  )
  analysis_data <- prep$data
  if (nrow(analysis_data) == 0L) {
    stop("No observations remain after subsetting and missing-data handling.", call. = FALSE)
  }

  allowed <- unique(c(code.censoring, code.event1, code.event2))
  old_opt <- getOption("cifmodeling.allowed", NULL)
  on.exit(options(cifmodeling.allowed = old_opt), add = TRUE)
  options(cifmodeling.allowed = allowed)

  outcome.type <- util_check_outcome_type(
    x = outcome.type,
    formula = formula,
    data = analysis_data,
    na.action = stats::na.pass,
    auto_message = auto_message
  )
  if (!outcome.type %in% c("survival", "competing-risk")) {
    stop("This interface supports only survival and competing-risk outcomes.", call. = FALSE)
  }

  Terms <- stats::terms(
    formula,
    specials = c("strata", "offset", "cluster"),
    data = analysis_data
  )
  mf <- stats::model.frame(Terms, data = analysis_data, na.action = stats::na.pass)
  Y <- stats::model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv"))) .err("surv_expected")

  te <- util_normalize_time_event(Y[, 1L], Y[, 2L], allowed = allowed)
  t <- te$time
  epsilon <- te$event
  if (anyNA(t) || anyNA(epsilon)) {
    stop("Missing time or event values remain after `na.action`.", call. = FALSE)
  }

  allowed_observed <- if (identical(outcome.type, "survival")) {
    c(code.censoring, code.event1)
  } else {
    c(code.censoring, code.event1, code.event2)
  }
  unexpected <- setdiff(unique(epsilon), allowed_observed)
  if (isTRUE(validate.observed.codes) && length(unexpected)) {
    stop(
      "Observed status codes are incompatible with `outcome.type`: ",
      paste(sort(unexpected), collapse = ", "),
      call. = FALSE
    )
  }

  term_labels <- attr(Terms, "term.labels")
  if (length(term_labels) == 0L) {
    strata_name <- NULL
    formula_strata <- factor(rep.int(1L, nrow(mf)))
  } else if (length(term_labels) == 1L) {
    strata_name <- term_labels[1L]
    formula_strata <- factor(mf[[strata_name]])
  } else {
    strata_name <- paste(term_labels, collapse = ":")
    formula_strata <- interaction(mf[term_labels], drop = TRUE)
  }

  if (is.null(weight_name)) {
    w <- rep.int(1, nrow(analysis_data))
  } else {
    w <- analysis_data[[weight_name]]
    check_weights(w)
  }

  if (!is.null(weight_name) && startsWith(weight_name, ".cifmodeling_weights")) {
    analysis_data[[weight_name]] <- NULL
  }

  structure(
    list(
      formula = formula,
      terms = Terms,
      data = analysis_data,
      data.original = data,
      row.index = prep$row.index,
      included = prep$included,
      exclusion.reason = prep$exclusion.reason,
      outcome.type = outcome.type,
      code.event1 = code.event1,
      code.event2 = code.event2,
      code.censoring = code.censoring,
      t = t,
      epsilon = epsilon,
      d = as.integer(epsilon != code.censoring),
      d0 = as.integer(epsilon == code.censoring),
      d1 = as.integer(epsilon == code.event1),
      d2 = as.integer(epsilon == code.event2),
      strata = formula_strata,
      strata_name = strata_name,
      w = as.numeric(w)
    ),
    class = "cif_prepared_input"
  )
}

util_read_surv <- function(formula, data, weights = NULL,
                           code.event1 = 1, code.event2 = 2, code.censoring = 0,
                           subset.condition = NULL, na.action = stats::na.omit,
                           outcome.type = NULL, auto_message = TRUE,
                           other.variables.analyzed = NULL,
                           validate.observed.codes = TRUE) {

  # --- resolve weights without forcing evaluation ---
  weights_expr <- substitute(weights)
  weights_resolved <- NULL

  if (missing(weights) || identical(weights_expr, quote(NULL))) {
    weights_resolved <- NULL

  } else if (is.name(weights_expr)) {
    nm <- as.character(weights_expr)
    if (nm %in% names(data)) {
      weights_resolved <- nm            # data列名として扱う（"w" と同等）
    } else {
      # 親フレームのオブジェクトとして評価（ip.weight など）
      weights_resolved <- eval(weights_expr, parent.frame(), parent.frame(2))
    }

  } else if (is.character(weights_expr) && length(weights_expr) == 1) {
    # "w" のような文字列リテラル
    weights_resolved <- as.character(weights_expr)

  } else {
    # df$w や get("w") 等の式は親フレームで評価
    weights_resolved <- eval(weights_expr, parent.frame(), parent.frame(2))
  }
  # -------------------------------------------------

  prepared <- cif_prepare_input(
    formula = formula,
    data = data,
    weights = weights_resolved,
    other.variables.analyzed = other.variables.analyzed,
    subset.condition = subset.condition,
    na.action = na.action,
    outcome.type = outcome.type,
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring,
    auto_message = auto_message,
    validate.observed.codes = validate.observed.codes
  )

  list(
    t = prepared$t,
    epsilon = prepared$epsilon,
    d = prepared$d,
    d0 = prepared$d0,
    d1 = prepared$d1,
    d2 = prepared$d2,
    strata = prepared$strata,
    strata_name = prepared$strata_name,
    w = prepared$w,
    data_sync = prepared$data,
    outcome.type = prepared$outcome.type,
    row.index = prepared$row.index,
    included = prepared$included,
    exclusion.reason = prepared$exclusion.reason,
    prepared = prepared
  )
}

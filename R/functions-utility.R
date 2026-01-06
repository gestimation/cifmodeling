`%||%` <- function(x, y) if (is.null(x)) y else x

createAnalysisDataset <- function(formula,
                                  data,
                                  other.variables.analyzed = NULL,
                                  subset.condition = NULL,
                                  na.action = na.pass,
                                  fill_missing = FALSE) {
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
  analysis_dataset <- analysis_dataset[, all_vars, drop = FALSE]
  return(na.action(analysis_dataset))
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

util_read_surv <- function(formula, data, weights = NULL,
                           code.event1 = 1, code.event2 = 2, code.censoring = 0,
                           subset.condition = NULL, na.action = stats::na.omit) {

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

  data <- createAnalysisDataset(formula, data, weights_resolved, subset.condition, na.action)
  weights <- weights_resolved

  # 以下、あなたの既存コードをそのまま
  allowed <- c(code.censoring, code.event1, code.event2)
  allowed <- unique(stats::na.omit(allowed))
  old_opt <- getOption("cifmodeling.allowed", NULL)
  on.exit(options(cifmodeling.allowed = old_opt), add = TRUE)
  options(cifmodeling.allowed = allowed)

  Terms <- terms(formula, specials = c("strata","offset","cluster"), data = data)
  mf    <- model.frame(Terms, data = data, na.action = na.action)

  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event","Surv"))) .err("surv_expected")

  te <- util_normalize_time_event(Y[,1], Y[,2], allowed = allowed)
  t <- te$time
  epsilon <- te$event
  if (any(t < 0, na.rm = TRUE)) .err("time_nonneg", arg = "time")

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
  if (any(is.na(idx))) .err("align_rows_fail")

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
    } else {
      w <- data[[weights]]
      check_weights(w)
    }
  } else {
    check_weights(weights)
    w <- weights
  }

  list(t=t, epsilon=epsilon, d=d, d0=d0, d1=d1, d2=d2,
       strata=strata, strata_name=strata_name, w=w, data_sync=data_sync)
}

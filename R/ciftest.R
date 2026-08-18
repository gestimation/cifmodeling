#' Tests for survival and cumulative incidence curves
#'
#' @description
#' Provides a formula-based interface for log-rank, Fleming-Harrington, Gray,
#' and augmented score tests. `test = "auto"` selects log-rank for survival
#' outcomes and the augmented score test for competing-risk outcomes. The
#' `"early"` and `"late"` presets select that same outcome-specific test with
#' Fleming--Harrington parameters `(rho, gamma) = (1, 0)` and `(0, 1)`,
#' respectively. Analysis strata can be specified independently with `strata`.
#' For augmented tests, working Aalen-Johansen distributions are estimated
#' within exposure-by-`strata.competing.risk` cells, censoring distributions
#' are estimated within `strata.censor`, and the closed-form augmentation is
#' used. Setting `iteration` to a positive integer applies that many Scheike
#' fixed-point refinements. Finite iterates are anchored as the closed-form
#' score plus the change in the AIPWCC image from the seed to the refined
#' direction, preserving the zeroth-step identity in the observed sample.
#'
#' @param formula A two-sided formula with an `Event()` or `Surv()` response
#'   and one grouping variable on the right-hand side.
#' @param data A data frame containing variables in the formula.
#' @param weights Optional numeric case weights or the name of a weight column.
#'   The standard Gray branch currently accepts integer frequency weights.


#' @param subset.condition Optional character string giving a logical condition to subset
#' `data` (default `NULL`).
#' @param na.action A function specifying the action to take on missing values (default `na.omit`).
#' @param outcome.type One of `"auto"`, `"competing-risk"`, or `"survival"`.
#'   The aliases `"C"` and `"S"` are accepted.
#' @param code.event1 Integer code of the event of interest (default `1`).
#' @param code.event2 Integer code of the competing risk (default `2`).
#' @param code.censoring Integer code of censoring (default `0`).
#' @param test One of `"auto"`, `"early"`, `"late"`, `"logrank"`, `"gray"`,
#'   `"score"`, `"augmented"`, or `"multiple"`. The aliases `"L"`, `"LR"`,
#'   and `"log-rank"` select log-rank; `"G"` selects Gray; and `"A"`, `"aug"`,
#'   and `"augmentation"` select the augmented score test. The early and late
#'   choices are outcome-specific single-direction presets. `"multiple"`
#'   (aliases `"multi"` and `"m"`) returns [ciftest_mdir()] using the default
#'   directions `(rho, gamma) = (2, 0), (0, 2), (0, 0)`. `"score"` uses the
#'   null score-IID variance; for competing-risk outcomes it may use censoring
#'   Kaplan--Meier strata.
#' @param rho,gamma Optional non-negative Fleming-Harrington weight
#'   parameters. Omitted values default to zero. They cannot be combined with
#'   the fixed `"early"` and `"late"` presets.
#' @param iteration Non-negative integer giving the number of Scheike
#'   fixed-point refinements. `0` returns the closed-form augmented score;
#'   positive values use the closed-form-anchored AIPWCC difference.
#' @param tolerance Optional positive convergence tolerance. `NULL` performs
#'   exactly the requested number of finite iterations.
#' @param strata Optional character vector of one or more column names
#'   defining analysis strata. Multiple columns define their observed
#'   interaction. The grouping variable cannot be included.
#'   Event risk sets, null event distributions, and Fleming--Harrington weight
#'   processes are constructed separately within these strata, and their score
#'   vectors are summed. This is distinct from the two nuisance-model strata.
#' @param strata.censor Optional character vector of one or more column names
#'   defining censoring Kaplan-Meier strata. The grouping variable may be
#'   included explicitly when group-specific censoring distributions are
#'   required.
#' @param strata.competing.risk Optional character vector of one or more column
#'   names defining the exposure-by-stratum working Aalen-Johansen models used
#'   for the competing-risk nuisance distribution. The grouping variable must
#'   not be included because it is crossed with these strata automatically.
#' @param tau Optional finite non-negative analysis horizon.
#' @param prob.bound Strictly positive numerical/positivity bound. Required
#'   nuisance probabilities at or below this value produce an error rather
#'   than being silently replaced.
#' @param prob.truncation Optional probability lower truncation strictly above
#'   `prob.bound`. It is applied only to positive censoring and working-survival
#'   denominators used by score and augmented tests.
#' @param ... Reserved for future extensions.
#'
#' @return An object inheriting from `"ciftest"` and `"htest"`, with a
#'   test-specific subclass. An ordinary unweighted log-rank result also
#'   inherits from `"survdiff"` and contains a fully compatible object in
#'   `survdiff`. Every score-based branch returns an analysis-row by contrast
#'   matrix of individual null-score contributions in `score.iid`.
#'   Augmented results include `score.iid.base`,
#'   `score.iid.censor`, and `score.iid.augment` matrices and use their summed
#'   empirical cross-product as `vcov.score`. Standard Gray results retain the
#'   classical Gray covariance. If the optional standard-Gray score-IID
#'   diagnostic cannot be constructed, the Gray test is still returned; its
#'   score-IID matrices contain `NA` and the reason is recorded in
#'   `diagnostics$score.iid.error`. Positive `iteration` values return the
#'   requested finite-iteration score, its subject-level score matrix,
#'   anchoring diagnostics, and fixed-point diagnostics. The event-time FH
#'   process actually used by the fit is returned in
#'   `diagnostics$fh.weight.process`, including whether it came from the pooled
#'   or event-stratified AJ left limit or the native Gray subdistribution
#'   construction.
#' @export
ciftest <- function(
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = "auto",
    test = "auto",
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    rho = NULL,
    gamma = NULL,
    iteration = 0L,
    tolerance = NULL,
    strata = NULL,
    strata.censor = NULL,
    strata.competing.risk = NULL,
    tau = NULL,
    prob.bound = 1e-7,
    prob.truncation = NULL,
    ...
) {
  call <- match.call()
  dots <- list(...)
  # Temporary development bridge. These names are no longer part of the
  # public signature and will be removed after package tests and simulation
  # profiles have been migrated to `test` and `strata`.
  legacy_augmentation <- dots$augmentation
  legacy_strata <- dots$strata.event
  dots$augmentation <- NULL
  dots$strata.event <- NULL
  if (!is.null(legacy_strata)) {
    if (!is.null(strata)) {
      stop("Use only `strata`; the removed `strata.event` alias was also supplied.",
           call. = FALSE)
    }
    strata <- legacy_strata
  }
  if (length(dots)) {
    dot_names <- names(dots)
    displayed <- if (is.null(dot_names)) {
      rep.int("<unnamed>", length(dots))
    } else {
      ifelse(nzchar(dot_names), dot_names, "<unnamed>")
    }
    stop(
      "Unused argument", if (length(dots) == 1L) "" else "s", ": ",
      paste(displayed, collapse = ", "),
      call. = FALSE
    )
  }
  weights.expr <- substitute(weights)
  weights.resolved <- util_resolve_weights(
    weights.expr = weights.expr,
    data = data,
    envir = parent.frame(),
    missing = missing(weights)
  )

  test_key <- if (is.character(test) && length(test) == 1L && !is.na(test)) {
    gsub("[-_. ]", "", tolower(trimws(test)))
  } else {
    ""
  }
  if (test_key %in% c("multiple", "multi", "m")) {
    if (!is.null(rho) || !is.null(gamma)) {
      stop(
        "`test = \"multiple\"` uses the fixed directions (2, 0), (0, 2), ",
        "and (0, 0); use `ciftest_mdir()` for custom directions.",
        call. = FALSE
      )
    }
    if (!is.null(legacy_augmentation)) {
      stop(
        "Use `test = \"multiple\"` without the removed `augmentation` argument.",
        call. = FALSE
      )
    }
    out <- ciftest_mdir(
      formula = formula,
      data = data,
      directions = c("early", "late", "unweighted"),
      weights = weights.resolved,
      subset.condition = subset.condition,
      outcome.type = outcome.type,
      test = "auto",
      code.event1 = code.event1,
      code.event2 = code.event2,
      code.censoring = code.censoring,
      iteration = iteration,
      tolerance = tolerance,
      strata = strata,
      strata.censor = strata.censor,
      strata.competing.risk = strata.competing.risk,
      tau = tau,
      prob.bound = prob.bound,
      prob.truncation = prob.truncation,
      na.action = na.action
    )
    out$call <- call
    out$test.requested <- "multiple"
    return(out)
  }

  input <- ciftest_prepare(
    formula = formula,
    data = data,
    weights = weights.resolved,
    subset.condition = subset.condition,
    outcome.type = outcome.type,
    test = test,
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring,
    rho = rho,
    gamma = gamma,
    iteration = iteration,
    tolerance = tolerance,
    strata = strata,
    strata.censor = strata.censor,
    strata.competing.risk = strata.competing.risk,
    tau = tau,
    prob.bound = prob.bound,
    prob.truncation = prob.truncation,
    legacy.augmentation = legacy_augmentation,
    na.action = na.action
  )
  ciftest_fit_prepared(
    input = input,
    formula = formula,
    call = call,
    nuisance.cache = NULL
  )
}

ciftest_resolve_test_spec <- function(
    outcome.type,
    test = "auto",
    rho = NULL,
    gamma = NULL,
    iteration = 0L,
    legacy.augmentation = NULL
) {
  if (!is.character(test) || length(test) != 1L || is.na(test) ||
      !nzchar(trimws(test))) {
    stop("`test` must be one non-missing character value.", call. = FALSE)
  }
  requested <- tolower(trimws(test))
  requested <- gsub("[-_. ]", "", requested)
  aliases <- c(
    auto = "auto", early = "early", late = "late",
    logrank = "logrank", l = "logrank", lr = "logrank",
    gray = "gray", g = "gray", score = "score",
    augmented = "augmented", a = "augmented", aug = "augmented",
    augmentation = "augmented"
  )
  if (!requested %in% names(aliases)) {
    stop(
      "`test` must be one of 'auto', 'early', 'late', 'logrank', ",
      "'gray', 'score', or 'augmented' (documented aliases are accepted).",
      call. = FALSE
    )
  }
  requested <- unname(aliases[requested])

  if (!is.null(legacy.augmentation)) {
    if (!is.logical(legacy.augmentation) ||
        length(legacy.augmentation) != 1L || is.na(legacy.augmentation)) {
      stop("The removed `augmentation` argument must be TRUE or FALSE.",
           call. = FALSE)
    }
    if (!identical(requested, "auto")) {
      stop("Use `test` instead of combining it with the removed `augmentation` argument.",
           call. = FALSE)
    }
    requested <- if (isTRUE(legacy.augmentation)) {
      "augmented"
    } else if (identical(outcome.type, "survival")) {
      "logrank"
    } else {
      "gray"
    }
  }

  if (requested %in% c("early", "late")) {
    if (!is.null(rho) || !is.null(gamma)) {
      stop(
        "`test = \"", requested,
        "\"` is a fixed preset; specify the underlying test with `rho` and `gamma` for a custom weight.",
        call. = FALSE
      )
    }
    resolved <- if (identical(outcome.type, "survival")) {
      "logrank"
    } else {
      "augmented"
    }
    rho <- if (identical(requested, "early")) 1 else 0
    gamma <- if (identical(requested, "late")) 1 else 0
    weight.label <- requested
  } else {
    resolved <- if (identical(requested, "auto")) {
      if (identical(outcome.type, "survival")) "logrank" else "augmented"
    } else {
      requested
    }
    rho <- if (is.null(rho)) 0 else rho
    gamma <- if (is.null(gamma)) 0 else gamma
    weight.label <- if (rho == 0 && gamma == 0) "unweighted" else "custom"
  }

  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0) {
    stop("`rho` must be NULL or one finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma) != 1L ||
      !is.finite(gamma) || gamma < 0) {
    stop("`gamma` must be NULL or one finite non-negative number.", call. = FALSE)
  }

  allowed <- if (identical(outcome.type, "survival")) {
    c("logrank", "score")
  } else {
    c("gray", "score", "augmented")
  }
  if (!resolved %in% allowed) {
    stop(
      "`test = \"", resolved, "\"` is not available for outcome.type = \"",
      outcome.type, "\".",
      call. = FALSE
    )
  }
  if (iteration > 0L && !identical(resolved, "augmented")) {
    stop("Positive `iteration` requires `test = \"augmented\"`.",
         call. = FALSE)
  }

  list(
    requested = if (is.null(legacy.augmentation)) requested else "legacy",
    test = resolved,
    rho = as.numeric(rho),
    gamma = as.numeric(gamma),
    weight.label = weight.label,
    augmentation = identical(resolved, "augmented"),
    score.construction = if (identical(resolved, "score") &&
      identical(outcome.type, "competing-risk")) "fine-gray" else "standard"
  )
}

# Normalize the public analysis/nuisance strata specification.  The returned
# value is always NULL or a non-empty character vector of distinct data-column
# names.  The exposure restrictions are role-specific: it is a valid censoring
# stratum but cannot define an analysis stratum or be repeated in the working
# AJ strata because working AJ cells are already exposure-specific.
ciftest_normalize_strata_columns <- function(
    value,
    argument,
    data,
    exposure = NULL,
    role = c("event", "censor", "competing-risk")) {
  role <- match.arg(role)
  if (is.null(value)) return(NULL)
  if (!is.character(value) || !length(value) || anyNA(value) ||
      any(!nzchar(value))) {
    stop(
      "`", argument,
      "` must be NULL or a non-empty character vector of column names.",
      call. = FALSE
    )
  }
  if (anyDuplicated(value)) {
    duplicated_names <- unique(value[duplicated(value)])
    stop(
      "`", argument, "` contains duplicated column names: ",
      paste(duplicated_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  missing_names <- setdiff(value, names(data))
  if (length(missing_names)) {
    stop(
      "`", argument, "` column", if (length(missing_names) == 1L) "" else "s",
      " not found in data: ", paste(missing_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.null(exposure) && exposure %in% value) {
    if (identical(role, "event")) {
      stop(
        "`", argument, "` must not include the grouping variable `", exposure,
        "`; within-stratum exposure contrasts would be lost.",
        call. = FALSE
      )
    }
    if (identical(role, "competing-risk")) {
      stop(
        "`strata.competing.risk` must not include the grouping variable `",
        exposure, "`; working Aalen-Johansen models are already fitted ",
        "within exposure-by-stratum cells.",
        call. = FALSE
      )
    }
  }
  value
}

ciftest_make_strata_info <- function(
    data,
    columns = NULL,
    role = c("event", "censor", "competing-risk"),
    weights = NULL) {
  role <- match.arg(role)
  n <- nrow(data)
  if (is.null(weights)) weights <- rep.int(1, n)
  if (length(weights) != n || anyNA(weights) || any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("Internal ciftest strata weights are invalid.", call. = FALSE)
  }
  if (is.null(columns)) {
    values <- factor(rep.int("pooled", n), levels = "pooled")
    mapping <- data.frame(stratum = "pooled", stringsAsFactors = FALSE)
  } else {
    parts <- lapply(data[columns], function(value) {
      droplevels(factor(value))
    })
    if (any(vapply(parts, anyNA, logical(1L)))) {
      stop("Missing strata values remain after preprocessing.", call. = FALSE)
    }
    if (length(columns) == 1L) {
      values <- parts[[1L]]
      mapping <- data.frame(
        stratum = levels(values), stringsAsFactors = FALSE,
        check.names = FALSE
      )
      mapping[[columns[[1L]]]] <- levels(values)
    } else {
      raw <- do.call(
        interaction,
        c(parts, list(drop = TRUE, lex.order = TRUE, sep = "\r"))
      )
      first <- match(levels(raw), as.character(raw))
      source_values <- data[first, columns, drop = FALSE]
      labels <- vapply(seq_len(nrow(source_values)), function(index) {
        paste(
          paste0(columns, "=", vapply(
            source_values[index, , drop = FALSE], as.character,
            character(1L)
          )),
          collapse = " | "
        )
      }, character(1L))
      labels <- make.unique(labels, sep = " #")
      values <- factor(
        as.integer(raw), levels = seq_along(levels(raw)), labels = labels
      )
      mapping <- data.frame(
        stratum = labels, source_values,
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }
  level_names <- levels(values)
  count_n <- tabulate(as.integer(values), nbins = length(level_names))
  weight_sum <- vapply(seq_along(level_names), function(index) {
    sum(weights[as.integer(values) == index])
  }, numeric(1L))
  counts <- data.frame(
    stratum = level_names,
    n = as.integer(count_n),
    weight = as.numeric(weight_sum),
    stringsAsFactors = FALSE
  )
  list(
    name = if (is.null(columns)) NULL else paste(columns, collapse = ":"),
    columns = if (is.null(columns)) character() else columns,
    key = if (is.null(columns)) NULL else paste(columns, collapse = "\r"),
    values = values,
    mapping = mapping,
    counts = counts,
    role = role
  )
}

ciftest_score_fh_process <- function(score.parts) {
  if (isTRUE(score.parts$diagnostics$event.stratified)) {
    out <- do.call(rbind, lapply(seq_along(score.parts$event.strata), function(i) {
      data.frame(
        time = score.parts$event.time,
        weight = score.parts$fh.weight[, i],
        stratum = score.parts$event.strata[i],
        source = "event-stratified-aj-left",
        stringsAsFactors = FALSE
      )
    }))
    rownames(out) <- NULL
    return(out)
  }
  data.frame(
    time = score.parts$event.time,
    weight = score.parts$fh.weight,
    stratum = "pooled",
    source = "pooled-aj-left",
    stringsAsFactors = FALSE
  )
}

# Internal fit engine shared by the scalar UI and simulation batch path.
ciftest_fit_prepared <- function(
    input,
    formula,
    call,
    nuisance.cache = NULL,
    precomputed = NULL
) {
  score_parts <- NULL
  augmentation_parts <- NULL
  score_iid_error <- NULL
  fh_weight_process <- NULL

  iteration_parts <- NULL
  iteration_anchor <- NULL
  fixed_point_solver <- if (is.null(input$fixed.point.solver)) {
    "finite"
  } else {
    input$fixed.point.solver
  }
  score_construction <- if (is.null(input$score.construction)) {
    "standard"
  } else {
    input$score.construction
  }
  if (identical(input$outcome.type, "competing-risk")) {
    exposure_design <- reg_read_exposure_design(
      input$data,
      exposure = input$exposure
    )

    if (identical(score_construction, "fine-gray")) {
      censoring <- ciftest_cache_censoring(
        nuisance.cache,
        input = input,
        strata = input$strata.censor.info$values,
        key = input$strata.censor.info$key
      )
      score_parts <- build_fg_score_iid(
        t = input$t,
        epsilon = input$epsilon,
        x = exposure_design$x_a,
        code.event1 = input$code.event1,
        code.event2 = input$code.event2,
        code.censoring = input$code.censoring,
        strata.event = input$strata.event.info$values,
        strata = input$strata.censor.info$values,
        weights = input$weights,
        rho = input$rho,
        gamma = input$gamma,
        censoring = censoring,
        prob.bound = input$prob.bound,
        prob.truncation = input$prob.truncation
      )
      fh_weight_process <- ciftest_score_fh_process(score_parts)
      comp <- list(
        score = score_parts$score,
        var = crossprod(score_parts$score.iid),
        exposure.labels = exposure_design$exposure.labels
      )
      method <- "Fine-Gray score test"
      variance_method <- "score-iid"
    } else if (isTRUE(input$augmentation)) {
      if (!is.null(precomputed)) {
        score_parts <- precomputed$score.parts
        augmentation_parts <- precomputed$augmentation.parts
      } else {
        censoring <- ciftest_cache_censoring(
          nuisance.cache,
          input = input,
          strata = input$strata.censor.info$values,
          key = input$strata.censor.info$key
        )
        score_parts <- build_fg_score_iid(
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          strata.event = input$strata.event.info$values,
          strata = input$strata.censor.info$values,
          weights = input$weights,
          rho = input$rho,
          gamma = input$gamma,
          censoring = censoring,
          prob.bound = input$prob.bound,
          prob.truncation = input$prob.truncation
        )
        working_aj <- ciftest_cache_working_aj(
          nuisance.cache,
          input = input,
          strata = input$strata.competing.risk.info$values,
          key = input$strata.competing.risk.info$key
        )
        augmentation_parts <- build_closed_form_augmentation(
          base = score_parts,
          working = working_aj,
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          exposure = input$data[[input$exposure]],
          strata.event = input$strata.event.info$values,
          strata.censor = input$strata.censor.info$values,
          strata.competing.risk = input$strata.competing.risk.info$values,
          weights = input$weights,
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          prob.bound = input$prob.bound
        )
      }
      fh_weight_process <- ciftest_score_fh_process(score_parts)
      closed_form_iid <- score_parts$score.iid +
        augmentation_parts$score.iid.augment
      if (identical(fixed_point_solver, "seed-map")) {
        iteration_parts <- build_seed_aipwcc_score_iid(
          base = score_parts,
          working = augmentation_parts$working.aj,
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          exposure = input$data[[input$exposure]],
          weights = input$weights,
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          prob.bound = input$prob.bound
        )
      } else if (identical(fixed_point_solver, "direct")) {
        iteration_setup <- ciftest_iteration_setup(
          base = score_parts,
          working = augmentation_parts$working.aj,
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          exposure = input$data[[input$exposure]],
          weights = input$weights,
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          prob.bound = input$prob.bound
        )
        fixed_point_operator <- ciftest_cache_fixed_point_operator(
          nuisance.cache,
          key = "pooled",
          compute = function() ciftest_appendix_e_operator(
            iteration_setup, score_parts, augmentation_parts$working.aj
          )
        )
        iteration_parts <- build_direct_fixed_point_score_iid(
          base = score_parts,
          working = augmentation_parts$working.aj,
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          exposure = input$data[[input$exposure]],
          weights = input$weights,
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          prob.bound = input$prob.bound,
          tolerance = if (is.null(input$tolerance)) 1e-8 else input$tolerance,
          setup = iteration_setup,
          operator = fixed_point_operator
        )
      } else if (input$iteration > 0L) {
        iteration_parts <- build_iterated_score_iid(
          base = score_parts,
          working = augmentation_parts$working.aj,
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          exposure = input$data[[input$exposure]],
          weights = input$weights,
          iterations = input$iteration,
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          prob.bound = input$prob.bound,
          tolerance = input$tolerance
        )
      }
      total_iid <- if (is.null(iteration_parts)) {
        closed_form_iid
      } else if (identical(fixed_point_solver, "seed-map")) {
        iteration_anchor <- list(
          applied = FALSE,
          construction = "raw AIPWCC(seed) diagnostic",
          closed.form.score = colSums(closed_form_iid),
          seed.aipwcc.score = colSums(iteration_parts$score.iid.seed),
          raw.aipwcc.score = colSums(iteration_parts$score.iid),
          score.adjustment = rep.int(0, ncol(closed_form_iid)),
          seed.score.difference =
            colSums(iteration_parts$score.iid.seed) -
              colSums(closed_form_iid),
          seed.iid.rms.difference = sqrt(mean(
            (iteration_parts$score.iid.seed - closed_form_iid)^2
          )),
          identity.error = 0
        )
        iteration_parts$score.iid
      } else {
        if (is.null(iteration_parts$score.iid.seed)) {
          stop("The iterated score is missing its AIPWCC seed image.",
               call. = FALSE)
        }
        anchored <- ciftest_anchor_aipwcc_iid(
          closed.form.iid = closed_form_iid,
          value.aipwcc.iid = iteration_parts$score.iid,
          seed.aipwcc.iid = iteration_parts$score.iid.seed
        )
        anchored_iid <- anchored$score.iid
        iteration_anchor <- list(
          applied = TRUE,
          construction =
            "closed-form + AIPWCC(value) - AIPWCC(seed)",
          closed.form.score = colSums(closed_form_iid),
          seed.aipwcc.score = colSums(iteration_parts$score.iid.seed),
          raw.aipwcc.score = colSums(iteration_parts$score.iid),
          score.adjustment = colSums(anchored$score.adjustment),
          seed.score.difference = anchored$seed.score.difference,
          seed.iid.rms.difference = anchored$seed.iid.rms.difference,
          identity.error = anchored$identity.error
        )
        anchored_iid
      }
      total_score <- colSums(total_iid)
      comp <- list(
        score = stats::setNames(total_score, colnames(exposure_design$x_a)),
        var = crossprod(total_iid),
        exposure.labels = exposure_design$exposure.labels
      )
      method <- if (is.null(iteration_parts)) {
        "Closed-form augmented Fine-Gray score test"
      } else if (identical(fixed_point_solver, "seed-map")) {
        "AIPWCC seed-map diagnostic score test"
      } else if (identical(fixed_point_solver, "direct")) {
        "Direct fixed-point time-weighted Fine-Gray score test"
      } else {
        "Iterated time-weighted Fine-Gray score test"
      }
      variance_method <- "score-iid"
    } else {
      if (!is.null(input$strata.censor.info$name) ||
          !is.null(input$strata.competing.risk.info$name)) {
        stop(
          "Nuisance-model strata are not available for the standard Gray test; ",
          "use `augmentation = TRUE`.",
          call. = FALSE
        )
      }
      comp <- calculate_gray(
        t = input$t,
        epsilon = as.integer(input$epsilon == input$code.event1) +
          2L * as.integer(input$epsilon == input$code.event2),
        exposure = input$exposure,
        weights = input$weights,
        strata = input$strata.event.info$values,
        data = input$data,
        rho = input$rho,
        gamma = input$gamma,
        prob.bound = input$prob.bound
      )
      fh_weight_process <- do.call(
        rbind,
        lapply(names(comp$fh.weight.process), function(level) {
          process <- comp$fh.weight.process[[level]]
          if (is.null(process) || !nrow(process)) return(NULL)
          data.frame(
            time = process$time,
            weight = process$weight,
            stratum = level,
            source = if (is.null(input$strata.event.info$name)) {
              "gray-pooled-subdistribution-left"
            } else {
              "gray-event-stratified-subdistribution-left"
            },
            stringsAsFactors = FALSE
          )
        })
      )
      if (!is.null(fh_weight_process)) rownames(fh_weight_process) <- NULL
      score_attempt <- if (!is.null(input$strata.event.info$name)) {
        simpleError(
          "The stratified standard-Gray score-IID diagnostic is not yet available."
        )
      } else {
        gray_censoring_strata <- factor(
          input$data[[input$exposure]], levels = comp$exposure.labels
        )
        gray_censoring <- ciftest_cache_censoring(
          nuisance.cache,
          input = input,
          strata = gray_censoring_strata,
          key = paste0(".gray-exposure:", input$exposure)
        )
        tryCatch(
          build_fg_score_iid(
            t = input$t,
            epsilon = input$epsilon,
            x = exposure_design$x_a,
            code.event1 = input$code.event1,
            code.event2 = input$code.event2,
            code.censoring = input$code.censoring,
            strata = gray_censoring_strata,
            weights = input$weights,
            rho = input$rho,
            gamma = input$gamma,
            fh.weight = comp$fh.weight.process[[1L]]$weight,
            censoring = gray_censoring,
            prob.bound = input$prob.bound
          ),
          error = identity
        )
      }
      if (inherits(score_attempt, "error")) {
        score_iid_error <- conditionMessage(score_attempt)
      } else {
        decomposition_tolerance <- 1e-8 * max(1, max(abs(comp$score)))
        decomposition_error <- max(abs(score_attempt$score - comp$score))
        if (decomposition_error > decomposition_tolerance) {
          score_iid_error <- paste0(
            "Internal Gray score decomposition check failed: error = ",
            format(decomposition_error, digits = 6),
            ", tolerance = ", format(decomposition_tolerance, digits = 6),
            "."
          )
        } else {
          score_parts <- score_attempt
        }
      }
      method <- if (input$rho == 0 && input$gamma == 0) {
        "Gray's test"
      } else {
        "Fleming-Harrington-type weighted Gray test"
      }
      variance_method <- "gray"
    }
  } else {
    comp <- calculate_log_rank(
      t = input$t,
      epsilon = as.integer(input$epsilon == input$code.event1),
      exposure = input$exposure,
      weights = input$weights,
      strata = input$strata.event.info$values,
      data = input$data,
      rho = input$rho,
      gamma = input$gamma,
      prob.bound = input$prob.bound
    )
    zero_iid <- matrix(
      0, nrow = nrow(comp$score.iid), ncol = ncol(comp$score.iid),
      dimnames = dimnames(comp$score.iid)
    )
    score_parts <- list(
      score = comp$score,
      score.iid = comp$score.iid,
      score.iid.base = comp$score.iid,
      score.iid.censor = zero_iid,
      censoring = NULL,
      diagnostics = list(
        score.decomposition.error = max(abs(
          colSums(comp$score.iid) - comp$score
        )),
        engine = "R"
      )
    )
    if (identical(input$test, "score")) {
      comp$var <- crossprod(comp$score.iid)
      method <- "Null robust log-rank score test"
      variance_method <- "score-iid"
    } else if (input$rho == 0 && input$gamma == 0) {
      method <- "Log-rank test"
      variance_method <- "hypergeometric"
    } else {
      method <- "Fleming-Harrington weighted log-rank test"
      variance_method <- "hypergeometric"
    }
  }
  test_stat <- ciftest_quadratic_form(comp$score, comp$var)

  if (is.null(score_parts)) {
    score_iid <- matrix(
      NA_real_,
      nrow = length(input$t),
      ncol = length(comp$score),
      dimnames = list(NULL, names(comp$score))
    )
    score_iid_base <- score_iid_censor <- score_iid
    score_iid_augment <- score_iid_iterated <- score_iid
  } else {
    score_iid_augment <- if (is.null(augmentation_parts)) {
      matrix(
        0,
        nrow = nrow(score_parts$score.iid),
        ncol = ncol(score_parts$score.iid),
        dimnames = dimnames(score_parts$score.iid)
      )
    } else {
      augmentation_parts$score.iid.augment
    }
    closed_form_iid <- score_parts$score.iid + score_iid_augment
    score_iid <- if (exists("total_iid", inherits = FALSE)) {
      total_iid
    } else if (is.null(iteration_parts)) {
      closed_form_iid
    } else {
      iteration_parts$score.iid
    }
    score_iid_base <- score_parts$score.iid.base
    score_iid_censor <- score_parts$score.iid.censor
    score_iid_iterated <- if (is.null(iteration_parts)) {
      matrix(
        NA_real_,
        nrow = nrow(score_iid),
        ncol = ncol(score_iid),
        dimnames = dimnames(score_iid)
      )
    } else {
      score_iid
    }
  }

  censoring_truncation <- if (is.null(score_parts)) NULL else
    score_parts$censoring
  working_truncation <- if (is.null(augmentation_parts)) NULL else
    augmentation_parts$working.aj
  censoring_truncation_count <- if (is.null(censoring_truncation)) 0L else
    as.integer(censoring_truncation$truncation.count)
  working_truncation_count <- if (is.null(working_truncation)) 0L else
    as.integer(working_truncation$truncation.count)
  truncation_diagnostics <- list(
    requested = !is.null(input$prob.truncation),
    applied = censoring_truncation_count > 0L ||
      working_truncation_count > 0L,
    probability = input$prob.truncation,
    censoring.count = censoring_truncation_count,
    working.survival.count = working_truncation_count,
    censoring.minimum.raw = if (is.null(censoring_truncation)) NA_real_ else
      censoring_truncation$minimum.survival,
    censoring.minimum.used = if (is.null(censoring_truncation)) NA_real_ else
      censoring_truncation$minimum.survival.used,
    working.survival.minimum.raw = if (is.null(working_truncation)) {
      NA_real_
    } else {
      working_truncation$minimum.survival
    },
    working.survival.minimum.used = if (is.null(working_truncation)) {
      NA_real_
    } else {
      working_truncation$minimum.survival.used
    }
  )

  out <- list(
    statistic = stats::setNames(test_stat$statistic, "X-squared"),
    parameter = stats::setNames(test_stat$rank, "df"),
    p.value = stats::pchisq(test_stat$statistic, df = test_stat$rank, lower.tail = FALSE),
    method = method,
    data.name = paste(deparse(formula), collapse = " "),
    call = call,
    outcome.type = input$outcome.type,
    test.requested = input$test.requested,
    test = input$test,
    weight.label = input$weight.label,
    code.event1 = input$code.event1,
    code.event2 = input$code.event2,
    code.censoring = input$code.censoring,
    rho = input$rho,
    gamma = input$gamma,
    augmentation = input$augmentation,
    iteration = input$iteration,
    tolerance = input$tolerance,
    tau = input$tau,
    prob.bound = input$prob.bound,
    prob.truncation = input$prob.truncation,
    variance.method = variance_method,
    score = comp$score,
    vcov.score = comp$var,
    observed = if (is.null(comp$observed)) NULL else comp$observed,
    expected = if (is.null(comp$expected)) NULL else comp$expected,
    z = if (length(comp$score) == 1L && comp$var[1L, 1L] > 0) {
      as.numeric(comp$score / sqrt(comp$var[1L, 1L]))
    } else {
      NA_real_
    },
    score.iid = score_iid,
    score.iid.base = score_iid_base,
    score.iid.censor = score_iid_censor,
    score.iid.augment = score_iid_augment,
    score.iid.iterated = score_iid_iterated,
    iterations = if (is.null(iteration_parts)) input$iteration else
      iteration_parts$iterations,
    converged = if (is.null(iteration_parts)) NA else
      iteration_parts$converged,
    fixed.point.residual = if (is.null(iteration_parts)) NA_real_ else
      iteration_parts$fixed.point.residual,
    last.increment = if (is.null(iteration_parts)) NA_real_ else
      iteration_parts$last.increment,
    contraction.ratio = if (is.null(iteration_parts)) NA_real_ else
      iteration_parts$contraction.ratio,
    n = length(input$t),
    n.events = table(factor(input$epsilon,
                            levels = c(input$code.censoring,
                                       input$code.event1,
                                       input$code.event2))),
    strata.info = input$strata.event.info,
    strata.event.info = input$strata.event.info,
    strata.censor.info = input$strata.censor.info,
    strata.competing.risk.info = input$strata.competing.risk.info,
    diagnostics = list(
      exposure = input$exposure,
      iteration.history = if (is.null(iteration_parts)) NULL else
        iteration_parts$history,
      iteration.components = if (is.null(iteration_parts) ||
        is.null(iteration_parts$components)) NULL else
        iteration_parts$components,
      fixed.point = if (is.null(iteration_parts)) NULL else
        iteration_parts$fixed.point,
      iteration.support = if (is.null(iteration_parts)) NULL else
        iteration_parts$diagnostics,
      iteration.anchor = iteration_anchor,
      fixed.point.solver = fixed_point_solver,
      score.construction = score_construction,
      score.iid.available = !is.null(score_parts),
      score.iid.error = score_iid_error,
      score.iid.variance.role = if (is.null(score_parts)) {
        "not available"
      } else if (identical(input$test, "score") ||
                 !is.null(augmentation_parts) ||
                 identical(score_construction, "fine-gray")) {
        "empirical score-iid cross-product used as the test variance"
      } else if (identical(input$test, "logrank")) {
        "diagnostic only; hypergeometric covariance remains the test variance"
      } else {
        "diagnostic only; Gray covariance remains the test variance"
      },
      fh.weight.process = fh_weight_process,
      censoring.km = if (is.null(score_parts)) NULL else score_parts$censoring,
      working.aj = if (is.null(augmentation_parts)) NULL else
        augmentation_parts$working.aj,
      augmentation.cells = if (is.null(augmentation_parts)) NULL else
        augmentation_parts$augmentation.cells,
      h.process = if (is.null(augmentation_parts)) NULL else
        augmentation_parts$h.process,
      augmentation.centering.error = if (is.null(augmentation_parts)) {
        NA_real_
      } else {
        augmentation_parts$diagnostics$augmentation.centering.error
      },
      score.decomposition.error = if (is.null(score_parts)) {
        NA_real_
      } else {
        score_parts$diagnostics$score.decomposition.error
      },
      score.engine = if (is.null(score_parts)) {
        NA_character_
      } else {
        score_parts$diagnostics$engine
      },
      augmentation.engine = if (is.null(augmentation_parts)) {
        NA_character_
      } else {
        augmentation_parts$diagnostics$engine
      },
      truncation = truncation_diagnostics,
      analysis.row.index = input$row.index,
      analysis.included = input$included,
      exclusion.reason = input$exclusion.reason
    ),
    data = input$data.original
  )
  survdiff_compatible <- identical(input$outcome.type, "survival") &&
    identical(input$test, "logrank") && input$rho == 0 && input$gamma == 0 &&
    all(input$weights == 1) && !is.null(comp$n.group)
  if (survdiff_compatible) {
    out$survdiff <- structure(
      list(
        n = comp$n.group,
        obs = comp$observed,
        exp = comp$expected,
        var = comp$var.full,
        chisq = unname(test_stat$statistic),
        pvalue = out$p.value,
        call = call
      ),
      class = "survdiff"
    )
  } else {
    out$survdiff <- NULL
  }
  subclass <- switch(
    input$test,
    logrank = "ciftest_logrank",
    gray = "ciftest_gray",
    score = "ciftest_score",
    augmented = "ciftest_augmented",
    "ciftest_result"
  )
  class(out) <- c(
    subclass, "ciftest", "htest",
    if (survdiff_compatible) "survdiff" else NULL
  )
  out
}

new_ciftest_nuisance_cache <- function() {
  cache <- new.env(parent = emptyenv())
  cache$censoring <- new.env(parent = emptyenv())
  cache$working.aj <- new.env(parent = emptyenv())
  cache$fixed.point.operator <- new.env(parent = emptyenv())
  cache$prepare.seconds <- 0
  cache$hits <- 0L
  cache$misses <- 0L
  class(cache) <- "ciftest_nuisance_cache"
  cache
}

ciftest_cache_fixed_point_operator <- function(cache, key, compute) {
  if (is.null(cache)) return(compute())
  if (!inherits(cache, "ciftest_nuisance_cache")) {
    stop("Internal ciftest nuisance cache has an invalid class.",
         call. = FALSE)
  }
  cache_key <- ciftest_cache_key("fixed-point", key)
  if (exists(cache_key, envir = cache$fixed.point.operator,
             inherits = FALSE)) {
    cache$hits <- cache$hits + 1L
    return(get(cache_key, envir = cache$fixed.point.operator,
               inherits = FALSE))
  }
  started <- proc.time()[["elapsed"]]
  value <- compute()
  cache$prepare.seconds <- cache$prepare.seconds +
    proc.time()[["elapsed"]] - started
  cache$misses <- cache$misses + 1L
  assign(cache_key, value, envir = cache$fixed.point.operator)
  value
}

ciftest_cache_key <- function(prefix, key) {
  label <- if (is.null(key) || !length(key) || is.na(key) || !nzchar(key)) {
    "pooled"
  } else {
    as.character(key)
  }
  paste(prefix, label, sep = "::")
}

ciftest_cache_censoring <- function(cache, input, strata, key = NULL) {
  compute <- function() {
    estimate_censoring_km(
      t = input$t,
      epsilon = input$epsilon,
      code.censoring = input$code.censoring,
      strata = strata,
      weights = input$weights,
      censoring.event = input$censoring.event,
      prob.bound = input$prob.bound,
      prob.truncation = input$prob.truncation
    )
  }
  if (is.null(cache)) return(compute())
  if (!inherits(cache, "ciftest_nuisance_cache")) {
    stop("Internal ciftest nuisance cache has an invalid class.",
         call. = FALSE)
  }
  probability_key <- paste(
    format(input$prob.bound, digits = 17L),
    if (is.null(input$prob.truncation)) "none" else
      format(input$prob.truncation, digits = 17L),
    sep = ":"
  )
  cache_key <- ciftest_cache_key(
    paste0("censoring[", probability_key, "]"), key
  )
  if (exists(cache_key, envir = cache$censoring, inherits = FALSE)) {
    cache$hits <- cache$hits + 1L
    return(get(cache_key, envir = cache$censoring, inherits = FALSE))
  }
  started <- proc.time()[["elapsed"]]
  value <- compute()
  cache$prepare.seconds <- cache$prepare.seconds +
    proc.time()[["elapsed"]] - started
  cache$misses <- cache$misses + 1L
  assign(cache_key, value, envir = cache$censoring)
  value
}

ciftest_cache_working_aj <- function(cache, input, strata, key = NULL) {
  compute <- function() {
    estimate_working_aj(
      t = input$t,
      epsilon = input$epsilon,
      exposure = input$data[[input$exposure]],
      strata = strata,
      weights = input$weights,
      code.event1 = input$code.event1,
      code.event2 = input$code.event2,
      code.censoring = input$code.censoring,
      prob.bound = input$prob.bound,
      prob.truncation = input$prob.truncation
    )
  }
  if (is.null(cache)) return(compute())
  if (!inherits(cache, "ciftest_nuisance_cache")) {
    stop("Internal ciftest nuisance cache has an invalid class.",
         call. = FALSE)
  }
  probability_key <- paste(
    format(input$prob.bound, digits = 17L),
    if (is.null(input$prob.truncation)) "none" else
      format(input$prob.truncation, digits = 17L),
    sep = ":"
  )
  cache_key <- ciftest_cache_key(
    paste0("working-aj[", probability_key, "]"), key
  )
  if (exists(cache_key, envir = cache$working.aj, inherits = FALSE)) {
    cache$hits <- cache$hits + 1L
    return(get(cache_key, envir = cache$working.aj, inherits = FALSE))
  }
  started <- proc.time()[["elapsed"]]
  value <- compute()
  cache$prepare.seconds <- cache$prepare.seconds +
    proc.time()[["elapsed"]] - started
  cache$misses <- cache$misses + 1L
  assign(cache_key, value, envir = cache$working.aj)
  value
}

ciftest_batch_strata_value <- function(method, name) {
  if (!name %in% names(method)) return(NULL)
  column <- method[[name]]
  value <- if (is.list(column)) column[[1L]] else column[1L]
  if (is.null(value) ||
      (length(value) == 1L && is.na(value))) return(NULL)
  value
}

ciftest_batch_has_strata <- function(methods, name) {
  vapply(seq_len(nrow(methods)), function(index) {
    !is.null(ciftest_batch_strata_value(
      methods[index, , drop = FALSE], name
    ))
  }, logical(1L))
}

ciftest_batch_strata_key <- function(method, name) {
  value <- ciftest_batch_strata_value(method, name)
  if (is.null(value)) "pooled" else paste(value, collapse = "\r")
}

ciftest_batch_method_input <- function(input, method) {
  scalar <- function(name, default = NULL) {
    if (!name %in% names(method)) return(default)
    value <- method[[name]][1L]
    if (!length(value) || is.na(value)) default else value
  }
  rho <- scalar("rho", 0)
  gamma <- scalar("gamma", 0)
  augmentation <- scalar("augmentation", FALSE)
  iteration <- scalar("iteration", 0L)
  fixed_point_solver <- scalar("fixed.point.solver", "finite")
  score_construction <- scalar("score.construction", "standard")
  strata.event <- ciftest_batch_strata_value(method, "strata.event")
  strata.censor <- ciftest_batch_strata_value(method, "strata.censor")
  strata.competing.risk <- ciftest_batch_strata_value(
    method, "strata.competing.risk"
  )
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0 ||
      !is.numeric(gamma) || length(gamma) != 1L ||
      !is.finite(gamma) || gamma < 0) {
    stop("Batch method rho and gamma must be finite non-negative numbers.",
         call. = FALSE)
  }
  if (!is.logical(augmentation) || length(augmentation) != 1L ||
      is.na(augmentation)) {
    stop("Batch method augmentation must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(iteration) || is.logical(iteration) ||
      length(iteration) != 1L || is.na(iteration) ||
      !is.finite(iteration) || iteration < 0 ||
      iteration != floor(iteration)) {
    stop("Batch method iteration must be one non-negative integer.",
         call. = FALSE)
  }
  if (iteration > 0L && !augmentation) {
    stop("Positive batch iteration requires augmentation = TRUE.",
         call. = FALSE)
  }
  if (!is.character(fixed_point_solver) ||
      length(fixed_point_solver) != 1L ||
      !fixed_point_solver %in% c("finite", "seed-map", "direct")) {
    stop("Batch fixed.point.solver must be 'finite', 'seed-map', or 'direct'.",
         call. = FALSE)
  }
  if (fixed_point_solver %in% c("seed-map", "direct") && !augmentation) {
    stop("Seed-map and direct fixed-point methods require augmentation = TRUE.",
         call. = FALSE)
  }
  if (fixed_point_solver %in% c("seed-map", "direct") && iteration != 0L) {
    stop("Seed-map and direct fixed-point methods require iteration = 0.",
         call. = FALSE)
  }
  if (!is.character(score_construction) ||
      length(score_construction) != 1L ||
      !score_construction %in% c("standard", "fine-gray")) {
    stop("Batch score.construction must be 'standard' or 'fine-gray'.",
         call. = FALSE)
  }
  if (identical(score_construction, "fine-gray") && augmentation) {
    stop("Fine-Gray score construction requires augmentation = FALSE.",
         call. = FALSE)
  }
  if (identical(score_construction, "fine-gray") &&
      !is.null(strata.competing.risk)) {
    stop("The Fine-Gray score does not use `strata.competing.risk`.",
         call. = FALSE)
  }
  if (!augmentation &&
      (!is.null(strata.competing.risk) ||
       (!identical(score_construction, "fine-gray") &&
        !is.null(strata.censor)))) {
    stop("Batch nuisance strata require augmentation = TRUE.", call. = FALSE)
  }
  if ((iteration > 0L || fixed_point_solver %in% c("seed-map", "direct")) &&
      (!is.null(strata.event) || !is.null(strata.censor) ||
       !is.null(strata.competing.risk))) {
    stop("Batch iteration currently requires pooled analysis and nuisance models.",
         call. = FALSE)
  }
  strata.event <- ciftest_normalize_strata_columns(
    strata.event, "strata.event", input$data, input$exposure, role = "event"
  )
  strata.censor <- ciftest_normalize_strata_columns(
    strata.censor, "strata.censor", input$data, input$exposure,
    role = "censor"
  )
  strata.competing.risk <- ciftest_normalize_strata_columns(
    strata.competing.risk, "strata.competing.risk", input$data,
    input$exposure, role = "competing-risk"
  )
  input$rho <- as.numeric(rho)
  input$gamma <- as.numeric(gamma)
  input$augmentation <- augmentation
  input$iteration <- as.integer(iteration)
  input$fixed.point.solver <- fixed_point_solver
  input$score.construction <- score_construction
  input$test <- if (identical(score_construction, "fine-gray")) {
    "score"
  } else if (isTRUE(augmentation)) {
    "augmented"
  } else if (identical(input$outcome.type, "survival")) {
    "logrank"
  } else {
    "gray"
  }
  input$test.requested <- "batch"
  input$weight.label <- if (rho == 0 && gamma == 0) {
    "unweighted"
  } else if (rho == 1 && gamma == 0) {
    "early"
  } else if (rho == 0 && gamma == 1) {
    "late"
  } else {
    "custom"
  }
  input$strata.event.info <- ciftest_make_strata_info(
    input$data, strata.event, role = "event", weights = input$weights
  )
  input$strata.censor.info <- ciftest_make_strata_info(
    input$data, strata.censor, role = "censor", weights = input$weights
  )
  input$strata.competing.risk.info <- ciftest_make_strata_info(
    input$data, strata.competing.risk, role = "competing-risk",
    weights = input$weights
  )
  input
}

ciftest_batch_precompute_multi <- function(base_input, methods, cache) {
  values <- vector("list", nrow(methods))
  available <- identical(base_input$outcome.type, "competing-risk") &&
    any(methods$augmentation) &&
    !any(ciftest_batch_has_strata(methods, "strata.event")) &&
    ciftest_use_cpp_kernel(
      "_cifmodeling_ciftest_fg_iid_multi_kernel_cpp"
    ) &&
    ciftest_use_cpp_kernel(
      "_cifmodeling_ciftest_augmentation_iid_multi_kernel_cpp"
    )
  if (!available) {
    return(list(values = values, elapsed = 0, errors = character()))
  }
  started <- proc.time()[["elapsed"]]
  augmented <- which(methods$augmentation)
  group_key <- vapply(augmented, function(index) {
    paste(
      ciftest_batch_strata_key(
        methods[index, , drop = FALSE], "strata.censor"
      ),
      ciftest_batch_strata_key(
        methods[index, , drop = FALSE], "strata.competing.risk"
      ),
      sep = "\034"
    )
  }, character(1L))
  groups <- split(augmented, group_key)
  score_cache <- new.env(parent = emptyenv())
  errors <- character()
  for (indices in groups) {
    if (length(indices) < 2L) next
    attempt <- tryCatch({
      input <- ciftest_batch_method_input(
        base_input, methods[indices[1L], , drop = FALSE]
      )
      exposure_design <- reg_read_exposure_design(
        input$data, exposure = input$exposure
      )
      censoring <- ciftest_cache_censoring(
        cache, input, input$strata.censor.info$values,
        input$strata.censor.info$key
      )
      working <- ciftest_cache_working_aj(
        cache, input, input$strata.competing.risk.info$values,
        input$strata.competing.risk.info$key
      )
      score_key <- paste(
        if (is.null(input$strata.censor.info$key)) {
          "pooled"
        } else input$strata.censor.info$key,
        paste(methods$rho[indices], methods$gamma[indices], sep = ":",
              collapse = ","),
        sep = "::"
      )
      if (exists(score_key, envir = score_cache, inherits = FALSE)) {
        score_parts <- get(score_key, envir = score_cache, inherits = FALSE)
      } else {
        score_parts <- build_fg_score_iid_multi(
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          rho = methods$rho[indices],
          gamma = methods$gamma[indices],
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          strata = input$strata.censor.info$values,
          weights = input$weights,
          censoring = censoring
        )
        assign(score_key, score_parts, envir = score_cache)
      }
      augmentation_parts <- build_closed_form_augmentation_multi(
        bases = score_parts,
        working = working,
        t = input$t,
        epsilon = input$epsilon,
        x = exposure_design$x_a,
        exposure = input$data[[input$exposure]],
        strata.censor = input$strata.censor.info$values,
        strata.competing.risk = input$strata.competing.risk.info$values,
        weights = input$weights,
        code.event1 = input$code.event1,
        code.event2 = input$code.event2,
        code.censoring = input$code.censoring
      )
      list(score = score_parts, augmentation = augmentation_parts)
    }, error = identity)
    if (inherits(attempt, "error")) {
      errors <- c(errors, conditionMessage(attempt))
      next
    }
    for (local_index in seq_along(indices)) {
      values[[indices[local_index]]] <- list(
        score.parts = attempt$score[[local_index]],
        augmentation.parts = attempt$augmentation[[local_index]]
      )
    }
  }
  list(
    values = values,
    elapsed = proc.time()[["elapsed"]] - started,
    errors = unique(errors)
  )
}

# Internal batch path. It intentionally is not exported: public ciftest()
# remains scalar while simulations can share nuisance fits within a replicate.
ciftest_batch_internal <- function(
    formula,
    data,
    methods,
    weights = NULL,
    subset.condition = NULL,
    outcome.type = "competing-risk",
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    tau = NULL,
    na.action = stats::na.omit
) {
  if (!is.data.frame(methods) || !nrow(methods)) {
    stop("`methods` must be a non-empty data frame.", call. = FALSE)
  }
  required <- c("augmentation", "strata.censor",
                "strata.competing.risk", "rho", "gamma")
  if (!all(required %in% names(methods))) {
    stop("`methods` is missing required batch columns: ",
         paste(setdiff(required, names(methods)), collapse = ", "),
         call. = FALSE)
  }
  if (!is.logical(methods$augmentation) || anyNA(methods$augmentation)) {
    stop("Batch augmentation values must be non-missing logical values.",
         call. = FALSE)
  }
  if (!"iteration" %in% names(methods)) {
    methods$iteration <- 0L
  }
  if (!"fixed.point.solver" %in% names(methods)) {
    methods$fixed.point.solver <- "finite"
  }
  if (!"score.construction" %in% names(methods)) {
    methods$score.construction <- "standard"
  }
  if (!"strata.event" %in% names(methods)) {
    methods$strata.event <- NA_character_
  }
  has_event_strata <- ciftest_batch_has_strata(methods, "strata.event")
  has_censor_strata <- ciftest_batch_has_strata(methods, "strata.censor")
  has_competing_strata <- ciftest_batch_has_strata(
    methods, "strata.competing.risk"
  )
  if (!is.numeric(methods$iteration) || is.logical(methods$iteration) ||
      anyNA(methods$iteration) || any(!is.finite(methods$iteration)) ||
      any(methods$iteration < 0) ||
      any(methods$iteration != floor(methods$iteration))) {
    stop("Batch iteration values must be non-negative integers.",
         call. = FALSE)
  }
  if (!is.numeric(methods$rho) || anyNA(methods$rho) ||
      any(!is.finite(methods$rho)) || any(methods$rho < 0) ||
      !is.numeric(methods$gamma) || anyNA(methods$gamma) ||
      any(!is.finite(methods$gamma)) || any(methods$gamma < 0)) {
    stop("Batch rho and gamma values must be finite and non-negative.",
         call. = FALSE)
  }
  if (!is.character(methods$fixed.point.solver) ||
      anyNA(methods$fixed.point.solver) ||
      any(!methods$fixed.point.solver %in% c("finite", "seed-map", "direct"))) {
    stop(
      "Batch fixed.point.solver values must be 'finite', 'seed-map', or 'direct'.",
         call. = FALSE)
  }
  if (!is.character(methods$score.construction) ||
      anyNA(methods$score.construction) ||
      any(!methods$score.construction %in% c("standard", "fine-gray"))) {
    stop("Batch score.construction values must be 'standard' or 'fine-gray'.",
         call. = FALSE)
  }
  fine_gray <- methods$score.construction == "fine-gray"
  if (any(fine_gray & methods$augmentation)) {
    stop("Fine-Gray score methods require augmentation = FALSE.",
         call. = FALSE)
  }
  if (any(fine_gray & (has_censor_strata | has_competing_strata))) {
    stop("Fine-Gray score methods currently require pooled nuisance models.",
         call. = FALSE)
  }
  fixed_point <- methods$fixed.point.solver %in% c("seed-map", "direct")
  if (any(fixed_point & !methods$augmentation)) {
    stop("Seed-map and direct fixed-point methods require augmentation = TRUE.",
         call. = FALSE)
  }
  if (any(fixed_point & methods$iteration != 0L)) {
    stop("Seed-map and direct fixed-point methods require iteration = 0.",
         call. = FALSE)
  }
  if (any(fixed_point & (has_event_strata | has_censor_strata |
                         has_competing_strata))) {
    stop("Seed-map and direct methods currently require pooled analysis and nuisance models.",
         call. = FALSE)
  }
  method_ids <- if ("method_id" %in% names(methods)) {
    as.character(methods$method_id)
  } else {
    paste0("method", seq_len(nrow(methods)))
  }
  if (anyNA(method_ids) || any(!nzchar(method_ids)) || anyDuplicated(method_ids)) {
    stop("Batch method IDs must be unique and non-missing.", call. = FALSE)
  }
  collect_strata_names <- function(name) {
    unique(unlist(lapply(seq_len(nrow(methods)), function(index) {
      ciftest_batch_strata_value(methods[index, , drop = FALSE], name)
    }), use.names = FALSE))
  }
  event_names <- collect_strata_names("strata.event")
  censor_names <- collect_strata_names("strata.censor")
  competing_names <- collect_strata_names("strata.competing.risk")
  strata_names <- unique(c(event_names, censor_names, competing_names))
  if (length(strata_names) &&
      any(!strata_names %in% names(data))) {
    stop("A batch analysis or nuisance stratum is not present in `data`.",
         call. = FALSE)
  }
  weights.resolved <- if (is.null(weights)) {
    NULL
  } else if (is.character(weights) && length(weights) == 1L &&
             weights %in% names(data)) {
    data[[weights]]
  } else {
    weights
  }
  input_started <- proc.time()[["elapsed"]]
  base_input <- ciftest_prepare(
    formula = formula,
    data = data,
    weights = weights.resolved,
    subset.condition = subset.condition,
    outcome.type = outcome.type,
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring,
    rho = 0,
    gamma = 0,
    iteration = 0L,
    strata = if (length(event_names)) event_names else NULL,
    strata.censor = if (length(censor_names)) censor_names else NULL,
    strata.competing.risk = if (length(competing_names)) {
      competing_names
    } else NULL,
    tau = tau,
    legacy.augmentation = any(methods$augmentation),
    na.action = na.action
  )
  input_prepare_seconds <- proc.time()[["elapsed"]] - input_started
  cache <- new_ciftest_nuisance_cache()
  batch_started <- proc.time()[["elapsed"]]
  precomputed_bundle <- ciftest_batch_precompute_multi(
    base_input, methods, cache
  )
  precomputed <- precomputed_bundle$values
  results <- vector("list", nrow(methods))
  names(results) <- method_ids
  timing <- vector("list", nrow(methods))
  batch_call <- match.call()
  for (i in seq_len(nrow(methods))) {
    method_started <- proc.time()[["elapsed"]]
    prepare_before <- cache$prepare.seconds
    input <- ciftest_batch_method_input(
      base_input,
      methods[i, , drop = FALSE]
    )
    results[[i]] <- tryCatch(
      ciftest_fit_prepared(
        input = input,
        formula = formula,
        call = batch_call,
        nuisance.cache = cache,
        precomputed = precomputed[[i]]
      ),
      error = identity
    )
    total <- proc.time()[["elapsed"]] - method_started
    prepare_delta <- cache$prepare.seconds - prepare_before
    timing[[i]] <- data.frame(
      method_id = method_ids[[i]],
      elapsed_prepare_seconds = prepare_delta,
      elapsed_method_seconds = max(0, total - prepare_delta),
      elapsed_seconds = total,
      stringsAsFactors = FALSE
    )
  }
  attr(results, "timing") <- do.call(rbind, timing)
  attr(results, "batch.timing") <- list(
    elapsed_input_prepare_seconds = input_prepare_seconds,
    elapsed_nuisance_prepare_seconds = cache$prepare.seconds,
    elapsed_multiweight_precompute_seconds = precomputed_bundle$elapsed,
    elapsed_batch_seconds = proc.time()[["elapsed"]] - batch_started,
    cache_hits = cache$hits,
    cache_misses = cache$misses,
    multiweight_errors = precomputed_bundle$errors
  )
  class(results) <- c("ciftest_batch", "list")
  results
}

ciftest_prepare <- function(
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    outcome.type = "auto",
    test = "auto",
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    rho = NULL,
    gamma = NULL,
    iteration = 0L,
    tolerance = NULL,
    strata = NULL,
    strata.censor = NULL,
    strata.competing.risk = NULL,
    tau = NULL,
    prob.bound = 1e-7,
    prob.truncation = NULL,
    legacy.augmentation = NULL,
    augmentation = NULL,
    na.action = stats::na.omit
) {
  if (!is.null(augmentation)) {
    if (!is.null(legacy.augmentation)) {
      stop("Supply only one internal legacy augmentation argument.",
           call. = FALSE)
    }
    legacy.augmentation <- augmentation
  }
  if (!is.numeric(iteration) || is.logical(iteration) ||
      length(iteration) != 1L || is.na(iteration) ||
      !is.finite(iteration) || iteration < 0 ||
      iteration != floor(iteration)) {
    stop("`iteration` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.null(tolerance) &&
      (!is.numeric(tolerance) || length(tolerance) != 1L ||
       !is.finite(tolerance) || tolerance <= 0)) {
    stop("`tolerance` must be NULL or one positive finite number.",
         call. = FALSE)
  }
  if (!is.null(tolerance) && iteration < 1L) {
    stop("`tolerance` requires a positive `iteration`.", call. = FALSE)
  }
  if (!is.numeric(prob.bound) || length(prob.bound) != 1L ||
      !is.finite(prob.bound) || prob.bound <= 0 || prob.bound >= 1) {
    stop("`prob.bound` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.null(prob.truncation) &&
      (!is.numeric(prob.truncation) || length(prob.truncation) != 1L ||
       !is.finite(prob.truncation) || prob.truncation <= prob.bound ||
       prob.truncation >= 1)) {
    stop("`prob.truncation` must be NULL or one number strictly between `prob.bound` and one.",
         call. = FALSE)
  }
  if (!is.null(tau) &&
      (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0)) {
    stop("`tau` must be NULL or one finite non-negative number.", call. = FALSE)
  }
  Terms0 <- stats::terms(formula, specials = c("strata", "offset", "cluster"), data = data)
  term_labels <- attr(Terms0, "term.labels")
  if (length(term_labels) != 1L || grepl("[:*(]", term_labels, fixed = FALSE)) {
    stop("The right-hand side of `formula` must be one untransformed grouping variable.", call. = FALSE)
  }
  exposure_vars <- all.vars(formula[[3L]])
  if (length(exposure_vars) != 1L || !exposure_vars %in% names(data)) {
    stop("The grouping variable must be one column in `data`.", call. = FALSE)
  }
  exposure <- exposure_vars[[1L]]

  strata <- ciftest_normalize_strata_columns(
    strata, "strata", data, exposure, role = "event"
  )
  strata.censor <- ciftest_normalize_strata_columns(
    strata.censor, "strata.censor", data, exposure, role = "censor"
  )
  strata.competing.risk <- ciftest_normalize_strata_columns(
    strata.competing.risk, "strata.competing.risk", data, exposure,
    role = "competing-risk"
  )

  prepared <- cif_prepare_input(
    formula = formula,
    data = data,
    weights = weights,
    other.variables.analyzed = unique(c(
      strata, strata.censor, strata.competing.risk
    )),
    subset.condition = subset.condition,
    na.action = na.action,
    outcome.type = if (is.character(outcome.type) && length(outcome.type) == 1L &&
      identical(tolower(trimws(outcome.type)), "auto")) NULL else outcome.type,
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring
  )
  reg_read_exposure_design(prepared$data, exposure = exposure)

  spec <- ciftest_resolve_test_spec(
    outcome.type = prepared$outcome.type,
    test = test,
    rho = rho,
    gamma = gamma,
    iteration = as.integer(iteration),
    legacy.augmentation = legacy.augmentation
  )
  if (!is.null(strata.censor) &&
      !(identical(spec$test, "augmented") ||
        (identical(spec$test, "score") &&
         identical(prepared$outcome.type, "competing-risk")))) {
    stop("`strata.censor` is available only for competing-risk score and augmented tests.",
         call. = FALSE)
  }
  if (!is.null(strata.competing.risk) && !identical(spec$test, "augmented")) {
    stop("`strata.competing.risk` requires `test = \"augmented\"`.",
         call. = FALSE)
  }
  if (!is.null(prob.truncation) &&
      !(identical(prepared$outcome.type, "competing-risk") &&
        spec$test %in% c("score", "augmented"))) {
    stop("`prob.truncation` is available only for competing-risk score and augmented tests.",
         call. = FALSE)
  }
  if (iteration > 0L &&
      (!is.null(strata) || !is.null(strata.censor) ||
       !is.null(strata.competing.risk))) {
    stop(
      "Analysis and nuisance-model strata are not available with iteration in the initial release.",
      call. = FALSE
    )
  }
  t <- prepared$t
  epsilon <- prepared$epsilon
  tau.used <- if (is.null(tau)) max(t) else tau
  horizon.complete <- rep.int(FALSE, length(t))
  if (!is.null(tau)) {
    horizon.complete <- t >= tau
    after_tau <- t > tau
    t[after_tau] <- tau
    epsilon[after_tau] <- prepared$code.censoring
  }
  censoring.event <- epsilon == prepared$code.censoring & !horizon.complete

  strata_event_info <- ciftest_make_strata_info(
    prepared$data, strata, role = "event", weights = prepared$w
  )
  strata_censor_info <- ciftest_make_strata_info(
    prepared$data, strata.censor, role = "censor", weights = prepared$w
  )
  strata_competing_risk_info <- ciftest_make_strata_info(
    prepared$data, strata.competing.risk, role = "competing-risk",
    weights = prepared$w
  )

  list(
    formula = formula,
    data = prepared$data,
    data.original = data,
    row.index = prepared$row.index,
    included = prepared$included,
    exclusion.reason = prepared$exclusion.reason,
    exposure = exposure,
    outcome.type = prepared$outcome.type,
    code.event1 = prepared$code.event1,
    code.event2 = prepared$code.event2,
    code.censoring = prepared$code.censoring,
    t = t,
    epsilon = epsilon,
    horizon.complete = horizon.complete,
    censoring.event = censoring.event,
    weights = prepared$w,
    test.requested = spec$requested,
    test = spec$test,
    weight.label = spec$weight.label,
    rho = spec$rho,
    gamma = spec$gamma,
    augmentation = spec$augmentation,
    score.construction = spec$score.construction,
    iteration = as.integer(iteration),
    tolerance = tolerance,
    tau = tau.used,
    prob.bound = as.numeric(prob.bound),
    prob.truncation = prob.truncation,
    strata.event.info = strata_event_info,
    strata.censor.info = strata_censor_info,
    strata.competing.risk.info = strata_competing_risk_info
  )
}

ciftest_quadratic_form <- function(score, variance) {
  variance <- (variance + t(variance)) / 2
  eig <- eigen(variance, symmetric = TRUE)
  tol <- max(1, max(abs(eig$values))) * sqrt(.Machine$double.eps)
  keep <- eig$values > tol
  rank <- sum(keep)
  if (rank < 1L) stop("The score covariance matrix has rank zero.", call. = FALSE)
  z <- crossprod(eig$vectors[, keep, drop = FALSE], as.numeric(score))
  statistic <- sum((z^2) / eig$values[keep])
  list(statistic = as.numeric(statistic), rank = as.integer(rank))
}

#' Compute Gray's K-sample test
#'
#' Implements the predictable-weight estimating equations and covariance from
#' Gray (1988). The pooled weight is
#' `(1 - F_1(t-))^rho * F_1(t-)^gamma`; `gamma = 0` recovers the family used
#' by `cmprsk::cuminc()`.
#'
#' `weights` are frequency weights. `strata` defines independent analysis
#' strata; stratum-specific scores and covariance matrices are summed.
#'
#' @keywords internal
calculate_gray <- function(
    t,
    epsilon,
    exposure,
    code.exposure.ref = NULL,
    prefix = "a",
    weights = rep.int(1, length(t)),
    strata = rep.int(1L, length(t)),
    data,
    rho = 0,
    gamma = 0,
    prob.bound = 1e-10
) {
  n <- length(t)
  stopifnot(
    length(epsilon) == n,
    length(weights) == n,
    length(strata) == n
  )
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!is.character(exposure) || length(exposure) != 1L ||
      !exposure %in% names(data)) {
    stop("`exposure` must name one column in `data`.", call. = FALSE)
  }
  if (any(!is.finite(t)) || any(t < 0)) {
    stop("Gray's test requires finite non-negative follow-up times.", call. = FALSE)
  }
  if (anyNA(epsilon) || any(!epsilon %in% 0:2)) {
    stop("`epsilon` must use 0 = censoring, 1 = event of interest, 2 = competing event.",
         call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0) ||
      any(abs(weights - round(weights)) > sqrt(.Machine$double.eps))) {
    stop("Standard Gray currently supports non-negative integer frequency weights only.",
         call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0 ||
      !is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("`rho` and `gamma` must be finite non-negative numbers.", call. = FALSE)
  }

  exp_info <- reg_read_exposure_design(
    data = data,
    exposure = exposure,
    code.exposure.ref = code.exposure.ref,
    prefix = prefix
  )
  K <- exp_info$exposure.levels
  if (K < 2L) stop("Exposure must have at least two levels.", call. = FALSE)

  group <- factor(data[[exposure]], levels = exp_info$exposure.labels)
  gid <- as.integer(group)
  if (anyNA(gid)) {
    stop("Missing or unknown exposure levels remain after preprocessing.", call. = FALSE)
  }

  q <- K - 1L
  score_base <- numeric(q)
  variance_base <- matrix(0, q, q)
  stratum <- droplevels(factor(strata))
  weight_process <- vector("list", nlevels(stratum))
  names(weight_process) <- levels(stratum)

  for (lev in levels(stratum)) {
    use <- which(stratum == lev & weights > 0)
    if (!length(use)) next
    part <- gray_stratum_components(
      t = t[use],
      epsilon = epsilon[use],
      gid = gid[use],
      weights = as.numeric(weights[use]),
      K = K,
      rho = rho,
      gamma = gamma,
      prob.bound = prob.bound
    )
    score_base <- score_base + part$score
    variance_base <- variance_base + part$var
    weight_process[[lev]] <- data.frame(
      time = part$event.time,
      weight = part$fh.weight
    )
  }

  # Gray's covariance is naturally parameterized by the first K-1 groups.
  # Transform it to the package convention: level 1 is the reference and the
  # returned coordinates correspond to levels 2,...,K.
  transform <- matrix(0, q, q)
  if (q > 1L) {
    transform[cbind(seq_len(q - 1L), 2:q)] <- 1
  }
  transform[q, ] <- -1
  score <- as.vector(transform %*% score_base)
  variance <- transform %*% variance_base %*% t(transform)

  score_names <- colnames(exp_info$x_a)
  names(score) <- score_names
  dimnames(variance) <- list(score_names, score_names)

  list(
    score = score,
    var = (variance + t(variance)) / 2,
    df = q,
    exposure.levels = K,
    exposure.labels = exp_info$exposure.labels,
    ref = exp_info$ref,
    fh.weight.process = weight_process
  )
}

# One-stratum Gray score and covariance in the first K-1 group coordinates.
gray_stratum_components <- function(
    t,
    epsilon,
    gid,
    weights,
    K,
    rho,
    gamma,
    prob.bound
) {
  q <- K - 1L
  times <- sort(unique(t))
  time_id <- match(t, times)
  m <- length(times)

  all_count <- event1_count <- event2_count <- matrix(0, m, K)
  for (i in seq_along(t)) {
    j <- time_id[i]
    k <- gid[i]
    all_count[j, k] <- all_count[j, k] + weights[i]
    if (epsilon[i] == 1L) {
      event1_count[j, k] <- event1_count[j, k] + weights[i]
    } else if (epsilon[i] == 2L) {
      event2_count[j, k] <- event2_count[j, k] + weights[i]
    }
  }

  risk <- colSums(all_count)
  survival_left <- rep.int(1, K)
  cif1_left <- numeric(K)
  pooled_cif_left <- 0
  score <- numeric(q)
  variance <- matrix(0, q, q)
  censor_link <- matrix(0, q, K)
  cross_term <- matrix(0, q, K)
  censor_variance <- numeric(K)
  event_time <- numeric()
  event_weight <- numeric()

  for (j in seq_len(m)) {
    d1 <- event1_count[j, ]
    d2 <- event2_count[j, ]
    d_all <- all_count[j, ]
    total_d1 <- sum(d1)
    total_d2 <- sum(d2)

    if (total_d1 + total_d2 > 0) {
      invalid <- risk > 0 & survival_left <= prob.bound
      if (any(invalid)) {
        stop("Gray covariance is undefined after a group survival estimate reaches zero.",
             call. = FALSE)
      }

      risk_over_survival <- numeric(K)
      active <- risk > 0
      risk_over_survival[active] <- risk[active] / survival_left[active]
      pooled_risk <- sum(risk_over_survival)
      if (!is.finite(pooled_risk) || pooled_risk <= 0) {
        stop("Gray's test encountered an empty pooled risk set.", call. = FALSE)
      }

      pooled_cif_right <- pooled_cif_left + total_d1 / pooled_risk
      gray_weight <- (1 - pooled_cif_left)^rho * pooled_cif_left^gamma

      influence_map <- -gray_weight *
        outer(risk_over_survival[seq_len(q)], risk_over_survival) / pooled_risk
      influence_map[cbind(seq_len(q), seq_len(q))] <-
        influence_map[cbind(seq_len(q), seq_len(q))] +
        gray_weight * risk_over_survival[seq_len(q)]

      if (total_d1 > 0) {
        event_time <- c(event_time, times[j])
        event_weight <- c(event_weight, gray_weight)
        cif_risk <- numeric(K)
        cif_risk[active] <- risk[active] * (1 - cif1_left[active]) /
          survival_left[active]
        cif_risk_total <- sum(cif_risk)
        if (!is.finite(cif_risk_total) || cif_risk_total <= 0) {
          stop("Gray's test encountered an empty subdistribution risk set.",
               call. = FALSE)
        }
        score <- score + gray_weight * (
          d1[seq_len(q)] - total_d1 * cif_risk[seq_len(q)] / cif_risk_total
        )

        remaining_cif <- 1 - pooled_cif_left
        if (remaining_cif <= prob.bound) {
          stop("Gray covariance is undefined when the pooled CIF reaches one.",
               call. = FALSE)
        }
        censor_link <- censor_link +
          influence_map * total_d1 / (pooled_risk * remaining_cif)
      }

      survival_right <- survival_left
      cif1_right <- cif1_left
      survival_right[active] <- survival_left[active] *
        (1 - (d1[active] + d2[active]) / risk[active])
      cif1_right[active] <- cif1_left[active] +
        survival_left[active] * d1[active] / risk[active]

      if (total_d1 > 0) {
        for (k in which(active)) {
          censor_ratio <- 1
          if (survival_right[k] > 0) {
            censor_ratio <- 1 - (1 - pooled_cif_right) / survival_right[k]
          }
          tie_adjust <- 1
          if (total_d1 > 1) {
            denominator <- pooled_risk * survival_left[k] - 1
            if (denominator <= 0) {
              stop("Gray target-event tie correction is undefined.", call. = FALSE)
            }
            tie_adjust <- 1 - (total_d1 - 1) / denominator
          }
          increment <- tie_adjust * survival_left[k] * total_d1 /
            (pooled_risk * risk[k])
          residual <- influence_map[, k] - censor_ratio * censor_link[, k]
          variance <- variance + tcrossprod(residual) * increment
          cross_term[, k] <- cross_term[, k] +
            residual * censor_ratio * increment
          censor_variance[k] <- censor_variance[k] +
            censor_ratio^2 * increment
        }
      }

      if (total_d2 > 0) {
        for (k in which(d2 > 0 & survival_right > 0)) {
          censor_ratio <- (1 - pooled_cif_right) / survival_right[k]
          tie_adjust <- 1
          if (d2[k] > 1) {
            if (risk[k] <= 1) {
              stop("Gray competing-event tie correction is undefined.", call. = FALSE)
            }
            tie_adjust <- 1 - (d2[k] - 1) / (risk[k] - 1)
          }
          increment <- tie_adjust * survival_left[k]^2 * d2[k] / risk[k]^2
          linked <- censor_ratio * censor_link[, k]
          variance <- variance + tcrossprod(linked) * increment
          cross_term[, k] <- cross_term[, k] -
            linked * censor_ratio * increment
          censor_variance[k] <- censor_variance[k] +
            censor_ratio^2 * increment
        }
      }

      pooled_cif_left <- pooled_cif_right
      survival_left <- survival_right
      cif1_left <- cif1_right
    }

    risk <- risk - d_all
    risk[abs(risk) < prob.bound] <- 0
  }

  variance <- variance +
    (censor_link * rep(censor_variance, each = q)) %*% t(censor_link) +
    censor_link %*% t(cross_term) + cross_term %*% t(censor_link)

  list(
    score = score,
    var = (variance + t(variance)) / 2,
    event.time = event_time,
    fh.weight = event_weight
  )
}

#' @export
print.ciftest <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("Outcome: ", x$outcome.type, "\n", sep = "")
  cat("Test: ", x$test, "\n", sep = "")
  cat("FH weights: rho = ", x$rho, ", gamma = ", x$gamma, "\n", sep = "")
  strata_text <- function(info) {
    if (is.null(info) || is.null(info$name)) return("pooled")
    columns <- info$columns
    if (is.null(columns) || !length(columns)) columns <- info$name
    paste0(
      paste(columns, collapse = " x "), " (",
      nlevels(info$values),
      if (nlevels(info$values) == 1L) {
        " observed stratum)"
      } else {
        " observed strata)"
      }
    )
  }
  cat("Analysis strata: ", strata_text(x$strata.info), "\n", sep = "")
  if (identical(x$test, "augmented")) {
    cat("Censoring strata: ", strata_text(x$strata.censor.info), "\n",
        sep = "")
    working_columns <- x$strata.competing.risk.info$columns
    if (is.null(working_columns) || !length(working_columns)) {
      working_text <- paste0(x$diagnostics$exposure, " only")
    } else {
      working_text <- paste(
        c(x$diagnostics$exposure, working_columns), collapse = " x "
      )
    }
    working_cells <- x$diagnostics$working.aj$cells
    cat(
      "Working AJ strata: ", working_text,
      if (!is.null(working_cells)) paste0(" (", nrow(working_cells),
                                          " fitted cells)") else "",
      "\n", sep = ""
    )
  } else if (identical(x$test, "score") &&
             identical(x$outcome.type, "competing-risk")) {
    cat("Censoring strata: ", strata_text(x$strata.censor.info), "\n",
        sep = "")
    cat("Working AJ strata: not used by this test\n")
  } else {
    cat("Censoring strata: not used by this test\n")
    cat("Working AJ strata: not used by this test\n")
  }
  fixed_point_solver <- x$diagnostics$fixed.point.solver
  if (identical(fixed_point_solver, "seed-map")) {
    cat("Diagnostic state: Appendix-E seed mapped to observed data\n")
  } else if (identical(fixed_point_solver, "direct")) {
    cat("Fixed-point solver: direct\n")
    cat("Iteration-converged: ", if (isTRUE(x$converged)) "yes" else "no",
        "\n", sep = "")
    cat("Fixed-point residual: ",
        format(x$fixed.point.residual, digits = 5L), "\n", sep = "")
  } else if (isTRUE(x$iteration > 0L)) {
    cat("Iteration: ", x$iterations, " of ", x$iteration,
        " requested refinement", if (x$iteration == 1L) "" else "s",
        "\n", sep = "")
    if (!is.null(x$tolerance)) {
      cat("Iteration tolerance: ", format(x$tolerance, digits = 5L),
          " (converged: ", if (isTRUE(x$converged)) "yes" else "no",
          ")\n", sep = "")
    }
    cat("Fixed-point residual: ",
        format(x$fixed.point.residual, digits = 5L), "\n", sep = "")
  } else if (identical(x$test, "augmented")) {
    cat("Iteration: 0 (closed-form augmentation)\n")
  }
  if (isTRUE(x$diagnostics$truncation$requested)) {
    truncation <- x$diagnostics$truncation
    cat(
      "Probability truncation: ", truncation$probability,
      " (censoring ", truncation$censoring.count,
      ", working survival ", truncation$working.survival.count,
      " replacements)\n", sep = ""
    )
  }
  cat("Variance: ", x$variance.method, "\n\n", sep = "")
  cat(
    "Chi-squared = ", format(unname(x$statistic), digits = 5L),
    ", df = ", unname(x$parameter),
    ", p-value = ", format.pval(x$p.value, digits = 4L), "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
summary.ciftest <- function(object, ...) object

#' @export
nobs.ciftest <- function(object, ...) object$n

#' @export
tidy.ciftest <- function(x, ...) {
  data.frame(
    statistic = unname(x$statistic),
    p.value = x$p.value,
    parameter = unname(x$parameter),
    method = x$method,
    outcome.type = x$outcome.type,
    test = x$test,
    weight = x$weight.label,
    iteration = x$iteration,
    fixed.point.residual = x$fixed.point.residual,
    n = x$n,
    stringsAsFactors = FALSE
  )
}

#' @export
glance.ciftest <- function(x, ...) tidy.ciftest(x, ...)

#' @export
augment.ciftest <- function(x, ...) {
  n0 <- nrow(x$data)
  p <- ncol(x$score.iid)
  restore <- function(value) {
    out <- matrix(NA_real_, nrow = n0, ncol = p, dimnames = list(NULL, colnames(value)))
    out[x$diagnostics$analysis.row.index, ] <- value
    I(out)
  }
  out <- x$data
  out$.score_iid <- restore(x$score.iid)
  out$.score_base <- restore(x$score.iid.base)
  out$.score_censoring <- restore(x$score.iid.censor)
  out$.score_augmentation <- restore(x$score.iid.augment)
  out$.score_iterated <- restore(x$score.iid.iterated)
  out$.censoring_weight <- NA_real_
  out$.analysis_included <- x$diagnostics$analysis.included
  out
}

#' @keywords internal
calculate_log_rank <- function(
    t,
    epsilon,
    exposure,
    code.exposure.ref = NULL,
    prefix = "a",
    weights,
    strata,
    data,
    rho = 0,
    gamma = 0,
    prob.bound = 1e-7
) {
  # --- basic checks ---
  n <- length(t)
  stopifnot(length(epsilon) == n, length(weights) == n, length(strata) == n)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(exposure) || length(exposure) != 1L) stop("`exposure` must be a single character string.")
  if (!exposure %in% names(data)) stop("exposure = '", exposure, "' is not found in data.")

  # --- exposure design (K levels -> K-1 dummies) ---
  exp_info <- reg_read_exposure_design(
    data = data,
    exposure = exposure,
    code.exposure.ref = code.exposure.ref,
    prefix = prefix
  )
  K <- exp_info$exposure.levels
  if (K < 2L) stop("Exposure must have >= 2 levels.")

  # Reconstruct the same factor coding/order as exp_info to get group IDs (1..K)
  a_ <- factor(data[[exposure]])
  a_ <- base::droplevels(a_)
  a_ <- factor(a_, levels = exp_info$exposure.labels) # ensures ref is level 1
  gid <- as.integer(a_)
  if (anyNA(gid)) stop("Missing values in exposure after preprocessing; handle via na.action before calling.")

  # event of interest indicator (epsilon==1)
  is_event1 <- (as.integer(epsilon) == 1L)

  # unique strata
  strata_fac <- factor(strata)
  L <- nlevels(strata_fac)

  p <- K - 1L
  score_total <- rep.int(0, p)
  var_total <- matrix(0, nrow = p, ncol = p)
  score_iid <- matrix(0, nrow = n, ncol = p)
  observed_total <- expected_total <- numeric(K)
  var_full_total <- matrix(0, nrow = K, ncol = K)

  # helper: compute FH weights from pooled survival in a stratum
  fh_weights <- function(Yw, dNw, rho, gamma, prob.bound) {
    m <- length(Yw)
    Kt <- numeric(m)
    S <- 1
    for (j in seq_len(m)) {
      S_clip <- min(max(S, prob.bound), 1 - prob.bound)
      Kt[j] <- (S_clip^rho) * ((1 - S_clip)^gamma)
      if (is.finite(Yw[j]) && Yw[j] > 0) {
        S <- S * (1 - dNw[j] / Yw[j])
      }
    }
    Kt
  }

  # --- loop strata ---
  for (ll in seq_len(L)) {
    idx <- which(strata_fac == levels(strata_fac)[ll])
    if (length(idx) == 0L) next

    tt <- t[idx]
    ww <- weights[idx]
    ee <- is_event1[idx]
    gg <- gid[idx]

    # event times in this stratum
    times <- sort(unique(tt[ee]))
    M <- length(times)
    if (M == 0L) next

    # order by time for risk set suffix sums
    o <- order(tt)
    tt_o <- tt[o]
    ww_o <- ww[o]
    gg_o <- gg[o]

    # position of each event time within ordered times (first occurrence)
    pos <- match(times, tt_o)
    if (anyNA(pos)) stop("Internal error: could not match event times in ordered times.")

    # --- risk set sums by group: Yw_mat (M x K) ---
    Yw_mat <- matrix(0, nrow = M, ncol = K)
    for (k in seq_len(K)) {
      suf <- rev(cumsum(rev(ww_o * (gg_o == k))))
      Yw_mat[, k] <- suf[pos]
    }
    Yw <- rowSums(Yw_mat)

    # --- event sums by group at each event time: dNw_mat (M x K) ---
    dNw_mat <- matrix(0, nrow = M, ncol = K)
    ev_idx <- which(ee)
    if (length(ev_idx) > 0L) {
      m_idx <- match(tt[ev_idx], times)
      k_idx <- gg[ev_idx]
      for (j in seq_along(ev_idx)) {
        dNw_mat[m_idx[j], k_idx[j]] <- dNw_mat[m_idx[j], k_idx[j]] + ww[ev_idx[j]]
      }
    }
    dNw <- rowSums(dNw_mat)

    # --- FH weights K(t-) computed from pooled survival in this stratum ---
    Kt <- fh_weights(Yw = Yw, dNw = dNw, rho = rho, gamma = gamma, prob.bound = prob.bound)

    # --- score: U = sum_t K(t) (dN_g - Y_g/Y * dN) for non-ref groups ---
    # Avoid division by zero (shouldn't happen at event times, but guard anyway)
    P <- Yw_mat
    for (j in seq_len(M)) {
      if (Yw[j] > 0) P[j, ] <- Yw_mat[j, ] / Yw[j] else P[j, ] <- 0
    }
    U_full <- dNw_mat - P * dNw
    U_red  <- U_full[, -1, drop = FALSE]             # drop reference group
    score_l <- as.vector(crossprod(Kt, U_red))       # (K-1)-vector
    observed_total <- observed_total +
      as.vector(crossprod(Kt, dNw_mat))
    expected_total <- expected_total +
      as.vector(crossprod(Kt, P * dNw))

    # Subject-level null score contributions. The compensator term sums to
    # zero at every event time, so the column sums reproduce the score while
    # retaining the individual martingale variation needed by robust score
    # covariance calculations.
    x_red <- matrix(0, nrow = length(idx), ncol = p)
    for (column in seq_len(p)) {
      x_red[, column] <- as.numeric(gg == column + 1L)
    }
    for (j in seq_len(M)) {
      if (!is.finite(Yw[j]) || Yw[j] <= 0 || dNw[j] <= 0) next
      centered <- sweep(x_red, 2L, P[j, -1L], "-")
      event_here <- ee & tt == times[j]
      at_risk <- tt >= times[j]
      martingale <- ww * (
        as.numeric(event_here) - as.numeric(at_risk) * dNw[j] / Yw[j]
      )
      score_iid[idx, ] <- score_iid[idx, , drop = FALSE] +
        centered * (Kt[j] * martingale)
    }

    # --- variance-covariance: sum_t K(t)^2 * Cov( O-E at t ) ---
    # Multigroup log-rank covariance increment (survdiff-style):
    # c = d*(Y-d)/(Y*(Y-1));  Cov = c*( diag(Yg) - (Yg Yh)/Y )
    var_l <- matrix(0, nrow = p, ncol = p)
    for (j in seq_len(M)) {
      Y <- Yw[j]
      d <- dNw[j]
      if (!is.finite(Y) || !is.finite(d) || Y <= 1 || d <= 0) next
      if (Y - d < 0) next  # numerical guard for extreme weights
      cfac <- d * (Y - d) / (Y * (Y - 1))

      Yvec <- Yw_mat[j, ]
      C_full <- cfac * (diag(Yvec, nrow = K, ncol = K) - tcrossprod(Yvec) / Y)
      C_red  <- C_full[-1, -1, drop = FALSE]
      var_l  <- var_l + (Kt[j]^2) * C_red
      var_full_total <- var_full_total + (Kt[j]^2) * C_full
    }

    score_total <- score_total + score_l
    var_total   <- var_total + var_l
  }

  colnames(var_total) <- rownames(var_total) <- colnames(exp_info$x_a)
  names(score_total) <- colnames(exp_info$x_a)
  colnames(score_iid) <- colnames(exp_info$x_a)
  names(observed_total) <- names(expected_total) <- exp_info$exposure.labels
  dimnames(var_full_total) <- list(
    exp_info$exposure.labels, exp_info$exposure.labels
  )

  list(
    score = score_total,
    var = var_total,
    score.iid = score_iid,
    observed = observed_total,
    expected = expected_total,
    var.full = var_full_total,
    n.group = stats::setNames(
      tabulate(gid, nbins = K), exp_info$exposure.labels
    ),
    df = length(score_total),
    exposure.levels = exp_info$exposure.levels,
    exposure.labels = exp_info$exposure.labels,
    ref = exp_info$ref
  )
}

#' Extract Mparts from a weightit object (robustly)
#' @keywords internal
get_weightit_mparts <- function(weightit) {
  if (!is.null(weightit$Mparts)) return(weightit$Mparts)
  mp <- attr(weightit, "Mparts", exact = TRUE)
  if (!is.null(mp)) return(mp)
  # NOTE: weightitMSM uses "Mparts.list" attribute (not for standard weightit)
  mp_list <- attr(weightit, "Mparts.list", exact = TRUE)
  if (!is.null(mp_list)) return(mp_list)
  stop("No Mparts found in `weightit`. Make sure the fit stored M-estimation parts.",
       call. = FALSE)
}

#' Choose a safe finite-difference step so w +/- h*delta stays >= 0
#' @keywords internal
fd_step_safe <- function(w, delta, rel = 1e-6, max_tries = 20L) {
  stopifnot(length(w) == length(delta))
  # target: max|h*delta| ≈ rel * max(1, median(w))
  scale_w <- max(1, stats::median(abs(w), na.rm = TRUE))
  max_abs_delta <- max(abs(delta), na.rm = TRUE)
  if (!is.finite(max_abs_delta) || max_abs_delta == 0) return(0)

  h <- rel * scale_w / max_abs_delta

  # ensure nonnegativity for both plus and minus
  idx <- which(abs(delta) > 0)
  if (length(idx)) {
    h_cap <- min(w[idx] / abs(delta[idx]), na.rm = TRUE)
    if (is.finite(h_cap)) h <- min(h, 0.49 * h_cap)
  }

  # backoff if still violates (due to numerical issues)
  for (k in seq_len(max_tries)) {
    if (all(w + h * delta >= 0) && all(w - h * delta >= 0)) return(h)
    h <- h * 0.5
  }
  stop("Could not find a safe finite-difference step (weights hit negative).",
       call. = FALSE)
}

#' A12 for log-rank-type score using directional finite differences with dw/dB from WeightIt
#'
#' Returns A12 on the "mean score" scale: (1/n) * dU_total/dB^T
#' @keywords internal
calculate_A12_logrank_weightit <- function(
    t,
    epsilon,
    strata,
    data,              # data_sync used by reg_read_exposure_design inside calculate_log_rank_components()
    exposure,
    weightit,
    code.exposure.ref = NULL,
    prefix = "a",
    rho = 0,
    gamma = 0,
    prob.bound = 1e-7,
    fd_rel_step = 1e-6
) {
  stopifnot(length(t) == length(epsilon),
            length(t) == length(strata))

  w0 <- weightit$weights
  if (is.null(w0)) stop("weightit$weights is NULL.", call. = FALSE)
  if (length(w0) != length(t)) {
    stop("Length mismatch: length(weightit$weights) != length(t).",
         call. = FALSE)
  }

  mp <- get_weightit_mparts(weightit)
  dw <- mp$dw_dBtreat
  if (is.null(dw)) stop("Mparts$dw_dBtreat is missing.", call. = FALSE)

  dw <- as.matrix(dw)
  if (nrow(dw) != length(t)) {
    stop("nrow(dw_dBtreat) != length(t). Check row alignment after subsetting/na.omit.",
         call. = FALSE)
  }

  # baseline score (total score, not averaged)
  base <- calculate_log_rank(
    t = t, epsilon = epsilon, weights = w0, strata = strata,
    data = data, exposure = exposure,
    code.exposure.ref = code.exposure.ref, prefix = prefix,
    rho = rho, gamma = gamma, prob.bound = prob.bound
  )
  if (is.null(base$score)) stop("calculate_log_rank() must return $score.", call. = FALSE)

  score0 <- as.numeric(base$score)
  p <- length(score0)
  n <- length(t)
  q <- ncol(dw)

  A12 <- matrix(NA_real_, nrow = p, ncol = q)
  rownames(A12) <- names(base$score) %||% paste0("score", seq_len(p))
  colnames(A12) <- colnames(dw) %||% paste0("Btreat", seq_len(q))

  steps <- numeric(q)

  for (j in seq_len(q)) {
    delta <- dw[, j]
    h <- fd_step_safe(w0, delta, rel = fd_rel_step)
    steps[j] <- h

    if (h == 0) {
      A12[, j] <- 0
      next
    }

    w_plus  <- w0 + h * delta
    w_minus <- w0 - h * delta

    # Compute scores at perturbed weights
    s_plus <- calculate_log_rank(
      t = t, epsilon = epsilon, weights = w_plus, strata = strata,
      data = data, exposure = exposure,
      code.exposure.ref = code.exposure.ref, prefix = prefix,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    s_minus <- calculate_log_rank(
      t = t, epsilon = epsilon, weights = w_minus, strata = strata,
      data = data, exposure = exposure,
      code.exposure.ref = code.exposure.ref, prefix = prefix,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    dU_dB_j <- (as.numeric(s_plus) - as.numeric(s_minus)) / (2 * h)

    # Put on mean-score scale: A12 = (1/n) * dU_total/dB^T
    A12[, j] <- dU_dB_j / n
  }

  list(A12 = A12, fd_steps = steps, score0 = base$score)
}

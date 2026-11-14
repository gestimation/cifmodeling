.cifreport_run <- function(data, plan, characteristics, overall, conf.int, table1.args,
                           subset.env = parent.frame()) {
  log <- list()
  record <- function(message, status, analysis, timepoint = NA_real_, outcome_id = NA_character_,
                     group_id = NA_character_, group_level = NA_character_) {
    log[[length(log) + 1L]] <<- data.frame(
      population_id, group_id, outcome_id, timepoint, group_level, analysis, status, message)
  }
  attempt <- function(expr, analysis, failure = "estimation_failed", ...) {
    context <- list(...)
    note <- function(message, status) {
      do.call(record, c(list(message = message, status = status, analysis = analysis), context))
    }
    tryCatch(withCallingHandlers(force(expr), warning = function(w) {
      note(conditionMessage(w), "warning")
      invokeRestart("muffleWarning")
    }), error = function(e) {
      note(conditionMessage(e), failure)
      NULL
    })
  }

  units <- list()
  exposure_titles <- list()
  population_id <- "All"
  cohort <- attempt(.cifreport_population(plan$subset.condition, data, subset.env),
                    "Population", "invalid_input")
  populations <- list()
  if (!is.null(cohort)) {
    if (!nrow(cohort)) {
      record("The population contains no subjects.", "skipped", "Population")
    } else {
      populations$All <- cohort
      if (!is.null(plan$subgroup)) {
        subgroups <- attempt(.cifreport_subgroups(cohort, plan$subgroup), "Subgroup", "invalid_input")
        if (!is.null(subgroups)) {
          if (anyNA(cohort[[plan$subgroup]]))
            record("Subjects with missing subgroup values appear in the overall population only.",
                   "skipped", "Subgroup")
          populations <- c(populations, subgroups)
        }
      }
    }
  }
  strata <- list()
  for (group_id in names(plan$exposure)) {
    if (!length(populations)) break
    definition <- plan$exposure[[group_id]]
    if (!definition$variable %in% names(data)) next
    strata[group_id] <- list(attempt(.cifreport_strata(definition, data, populations$All),
                                    "Grouping", "invalid_input", group_id = group_id))
  }
  for (population_id in names(populations)) {
    cohort <- populations[[population_id]]
    units[[population_id]] <- list()
    for (group_id in names(plan$exposure)) {
      definition <- plan$exposure[[group_id]]
      group <- attempt(.cifreport_group(cohort, definition$variable), "Grouping",
                       "invalid_input", group_id = group_id)
      if (is.null(group)) next
      scope <- function(expr, analysis, failure = "estimation_failed", ...) {
        attempt(expr, analysis, failure, group_id = group_id, ...)
      }
      note <- function(message, status, analysis, ...) {
        record(message, status, analysis, group_id = group_id, ...)
      }
      level <- levels(group)
      ref <- definition$code.exposure.ref %||% level[1L]
      labels <- strata[[group_id]]
      if (is.null(labels)) next
      labels <- labels[names(labels) %in% level]
      group <- factor(group, levels = names(labels))
      exposure_title <- .cif_variable_label(data[[definition$variable]], definition$variable)
      if (length(.cif_value_labels(data[[definition$variable]]))) {
        exposure_title <- paste0(exposure_title, " (", paste(unname(labels), collapse = ", "), ")")
      }
      exposure_titles[[population_id]][[group_id]] <- exposure_title
      counts <- stats::setNames(as.integer(table(group)), unname(labels))
      table1 <- scope(.cifreport_baseline(cohort, definition$variable,
                       factor(group, levels = names(labels), labels = unname(labels)),
                       characteristics, overall, table1.args, source = data), "Table 1", "invalid_input")

      valid_ref <- (is.character(ref) || is.numeric(ref) || is.logical(ref)) &&
        length(ref) == 1L && !is.na(ref) && as.character(ref) %in% level
      if (valid_ref) ref <- as.character(ref)
      if (length(level) == 1L)
        note("Only one group level is observed.", "skipped", "Tables 4 and 5")
      if (!valid_ref)
        note("The specified reference is not an observed group level.", "invalid_input", "Table 4")

      results <- list()
      followup <- list()
      for (outcome_id in names(plan$outcomes)) {
        specification <- plan$outcomes[[outcome_id]]
        time <- if (.cifreport_string(specification$time)) cohort[[specification$time]]
        if (is.numeric(time) && !is.complex(time) && is.null(dim(time)) && any(is.finite(time)))
          followup[[outcome_id]] <- range(time[is.finite(time)])
        outcome <- scope(.cifreport_outcome(specification, cohort, group),
                         "Outcome", "invalid_input", outcome_id = outcome_id)
        if (is.null(outcome)) next
        attr(outcome$data$followup, "label") <-
          .cif_variable_label(data[[outcome$time.variable]], outcome$time.variable)
        attr(outcome$data$group, "label") <-
          .cif_variable_label(data[[definition$variable]], definition$variable)
        results[[outcome_id]] <- .cifreport_analyse_outcome(outcome, conf.int, ref, valid_ref, labels,
          subgroup = population_id != "All",
          attempt = function(expr, analysis, ...) scope(expr, analysis, outcome_id = outcome_id, ...),
          note = function(message, status, analysis, ...) {
            note(message, status, analysis, outcome_id = outcome_id, ...)
          })
      }

      combine <- function(name, columns) {
        rows <- lapply(names(results), function(id) {
          x <- results[[id]][[name]]
          if (is.null(x) || !nrow(x)) return(NULL)
          data.frame(population_id = population_id, group_id = group_id, outcome_id = id, x)
        })
        rows <- Filter(Negate(is.null), rows)
        if (length(rows)) return(do.call(rbind, c(rows, list(make.row.names = FALSE))))
        data.frame(population_id = character(), group_id = character(), outcome_id = character(),
                   columns)
      }
      interval <- data.frame(estimate = numeric(), conf.low = numeric(), conf.high = numeric())
      units[[population_id]][[group_id]] <- list(
        narrative = list(exposure = .cif_variable_label(data[[definition$variable]], definition$variable),
                         counts = counts, followup = followup,
                         reference = if (valid_ref && length(level) > 1L) unname(labels[ref])),
        tables = list(table1 = table1,
          table2 = combine("counts", data.frame(time = character(), event = character(),
            group_level = character(), group_label = character(), primary = integer(), code.event1 = integer(),
            competing = integer(), code.event2 = integer(), censored = integer(),
            code.censoring = integer(), total = integer())),
          table3 = combine("probabilities", data.frame(timepoint = numeric(), group_level = character(),
            group_label = character(),
            primary = integer(), competing = integer(), interval, estimator = character(),
            status = character())),
          table4 = combine("ratios", data.frame(timepoint = numeric(), event_type = character(),
            comparison = character(), reference = character(), comparison_label = character(),
            reference_label = character(),
            interval, p.value = numeric(), status = character())),
          table5 = combine("tests", data.frame(method = character(), statistic = numeric(),
            parameter = numeric(), p.value = numeric()))),
        figures = lapply(results, `[[`, "plot"),
        fits = lapply(results, `[[`, "fits"))
    }
  }
  collect <- function(name) lapply(units, function(x) lapply(x, `[[`, name))
  subgroup <- plan$subgroup
  subgroup_titles <- if (!is.null(subgroup) && subgroup %in% names(data)) {
    variable <- data[[subgroup]]
    ids <- setdiff(names(populations), "All")
    values <- substring(ids, nchar(subgroup) + 4L)
    stats::setNames(paste(.cif_variable_label(variable, subgroup),
      .cifreport_labels_at(variable, values)), ids)
  }
  list(tables = collect("tables"), figures = collect("figures"), fits = collect("fits"),
       narrative = collect("narrative"),
       exposure_titles = exposure_titles, subgroup_titles = subgroup_titles,
       status = if (length(log)) do.call(rbind, log) else data.frame(
         population_id = character(), group_id = character(), outcome_id = character(),
         timepoint = numeric(), group_level = character(), analysis = character(),
         status = character(), message = character()))
}

.cifreport_labels_at <- function(x, codes) {
  unname(.cif_display_levels(codes, .cif_value_labels(x)))
}

.cifreport_strata <- function(definition, data, population) {
  variable <- data[[definition$variable]]
  original <- levels(factor(variable))
  observed <- unique(as.character(population[[definition$variable]]))
  observed <- observed[!is.na(observed)]
  codes <- function(x, arg) {
    if (!is.atomic(x) || !is.null(dim(x)) || !length(x) || anyNA(x))
      stop("`", arg, "` must contain non-missing exposure codes.")
    x <- as.character(x)
    if (anyDuplicated(x)) stop("`", arg, "` must contain distinct exposure codes.")
    x
  }
  level <- if (is.null(definition$level.strata)) original else
    codes(definition$level.strata, "level.strata")
  if (!all(observed %in% level)) stop("`level.strata` must include every observed exposure code.")
  order <- if (is.null(definition$order.strata)) level else codes(definition$order.strata, "order.strata")
  if (!all(order %in% level)) stop("`order.strata` includes an unknown exposure code.")
  labels <- definition$label.strata
  if (!is.null(labels)) {
    if (!is.character(labels) || !is.null(dim(labels)) || !length(labels) ||
        anyNA(labels) || any(!nzchar(labels)))
      stop("`label.strata` must contain nonempty character labels.")
    if (is.null(names(labels))) {
      if (length(labels) != length(level))
        stop("Unnamed `label.strata` must have one label per exposure code.")
      names(labels) <- level
    } else if (anyNA(names(labels)) || anyDuplicated(names(labels)) || !all(names(labels) %in% level)) {
      stop("Names of `label.strata` must be distinct exposure codes.")
    }
  }
  mapping <- .cif_display_levels(level, .cif_value_labels(variable), labels)
  if (anyDuplicated(mapping[observed])) stop("Exposure display labels must be distinct.")
  mapping[c(order, setdiff(level, order))]
}

.cifreport_subgroups <- function(data, variable) {
  if (!variable %in% names(data)) stop("Subgroup column is absent: ", variable, ".")
  value <- data[[variable]]
  if (!is.null(dim(value)) || !(is.factor(value) || is.character(value) || is.logical(value) ||
      (is.numeric(value) && !is.complex(value) && all(is.finite(value[!is.na(value)])))))
    stop("The subgroup column must contain categorical values or finite numeric category codes.")
  groups <- split(data, value, drop = TRUE)
  if (length(groups)) names(groups) <- paste0(variable, " = ", names(groups))
  groups
}

.cifreport_baseline <- function(data, variable, group, characteristics, overall, args, source = data) {
  missing <- setdiff(characteristics, names(data))
  if (length(missing)) stop("Baseline columns are absent: ", paste(missing, collapse = ", "), ".")
  source <- source[, unique(c(variable, characteristics)), drop = FALSE]
  data <- data[, unique(c(variable, characteristics)), drop = FALSE]
  label <- attr(source[[variable]], "label", exact = TRUE)
  attr(group, "label") <- label
  data[[variable]] <- group
  include <- setdiff(characteristics, variable)
  if (!length(include)) {
    # Include a single constant characteristic to show group sample sizes.
    sample <- make.unique(c(names(data), "Participants"))[ncol(data) + 1L]
    data[[sample]] <- factor(rep("Included", nrow(data)))
    attr(data[[sample]], "label") <- "Participants"
    include <- sample
  }
  args <- utils::modifyList(list(missing = "ifany"), args)
  labelled <- names(Filter(function(x) length(.cif_value_labels(x)) > 0L,
                           source[intersect(include, names(source))]))
  for (name in labelled) {
    x <- data[[name]]
    values <- as.character(x)
    values[is.na(x)] <- NA_character_
    codes <- unique(c(names(.cif_value_labels(source[[name]])), levels(factor(values))))
    data[[name]] <- factor(values, levels = codes, labels = .cifreport_labels_at(source[[name]], codes))
    attr(data[[name]], "label") <- attr(source[[name]], "label", exact = TRUE)
  }
  if (length(labelled)) {
    type <- args$type
    if (inherits(type, "formula")) type <- list(type)
    args$type <- c(type, stats::setNames(as.list(rep("categorical", length(labelled))), labelled))
  }
  table <- do.call(gtsummary::tbl_summary, c(list(data = data, by = variable, include = include), args))
  if (nlevels(group) > 1L) {
    comparable <- vapply(data[include], function(x) {
      present <- !is.na(x)
      length(unique(group[present])) > 1L && length(unique(x[present])) > 1L
    }, logical(1L))
    if (any(comparable)) {
      table <- gtsummary::add_p(table, include = names(comparable)[comparable],
        pvalue_fun = function(p) ifelse(is.na(p), NA_character_, sprintf("%.3f", p)))
    }
  }
  if (overall && nlevels(group) > 1L) table <- gtsummary::add_overall(table)
  table
}

.cifreport_analyse_outcome <- function(outcome, conf.int, ref, valid_ref, labels, subgroup, attempt, note) {
  d <- outcome$data
  args <- outcome$arguments
  competing <- args$outcome.type == "competing-risk"
  level <- levels(d$group)
  times <- outcome$times
  estimator <- outcome$estimator
  risk_scale <- estimator %in% c("AJ", "1-KM")
  scale <- switch(estimator, KM = "Survival probability", "1-KM" = "Event probability (1 - survival)",
    AJ = "Primary-event cumulative incidence", "1-AJ" = "1 - Primary-event cumulative incidence")
  groups <- split(d, d$group)
  counts <- do.call(rbind, lapply(groups, function(x) {
    data.frame(time = outcome$time.variable, event = outcome$event.variable,
      group_level = as.character(x$group[1L]),
      group_label = unname(labels[as.character(x$group[1L])]),
      primary = sum(x$status == args$code.event1),
      code.event1 = args$code.event1,
      competing = if (competing) sum(x$status == args$code.event2) else NA_integer_,
      code.event2 = if (competing) args$code.event2 else NA_integer_,
      censored = sum(x$status == args$code.censoring), code.censoring = args$code.censoring,
      total = nrow(x))
  }))
  if (subgroup && sum(counts$primary, counts$competing, na.rm = TRUE) == 0L) {
    note("No events were observed in this subgroup over follow-up.", "skipped",
         "Tables 3-5 and Figure 1")
    return(list(counts = counts))
  }

  fit_args <- c(args, list(data = quote(report_data)))
  # The shared curve/test input reader validates all three codes, even for survival.
  if (!competing) fit_args$code.event2 <- setdiff(0:2, c(args$code.event1, args$code.censoring))[1L]
  f <- if (length(level) > 1L) Event(followup, status) ~ group else Event(followup, status) ~ 1
  # cardx evaluates the data in the formula environment to recover group labels.
  report_data <- d
  environment(f) <- list2env(list(report_data = report_data), parent = environment(cifreport))
  curve <- attempt(do.call(cifcurve, c(list(formula = f, conf.int = conf.int), fit_args)), "Curve")
  support <- vapply(groups, function(x) max(x$followup), numeric(1L))
  probabilities <- do.call(rbind, lapply(groups, function(x) {
    cumulative <- function(code) {
      findInterval(times, sort(x$followup[x$status == code]))
    }
    data.frame(timepoint = times, group_level = as.character(x$group[1L]),
      group_label = unname(labels[as.character(x$group[1L])]),
      primary = cumulative(args$code.event1),
      competing = if (competing) cumulative(args$code.event2) else NA_integer_,
      estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, estimator = estimator,
      status = "estimate_unavailable")
  }))
  outside <- probabilities$timepoint > support[probabilities$group_level]
  for (i in which(outside)) {
    note("Time is beyond the last observed follow-up in this group.", "unsupported", "Table 3",
         timepoint = probabilities$timepoint[i], group_level = probabilities$group_level[i])
  }
  if (!is.null(curve)) {
    estimates <- attempt(.cifreport_probabilities(curve, times, level, risk_scale),
                         "Table 3 (tbl_survfit)")
    if (is.null(estimates)) {
      estimates <- attempt(.cifreport_curve_summary(curve, times, level, risk_scale), "Table 3")
    }
    if (!is.null(estimates)) {
      probabilities[names(estimates)] <- estimates
      probabilities[outside, c("estimate", "conf.low", "conf.high")] <- NA_real_
      probabilities$status <- ifelse(is.finite(probabilities$estimate), "ok", "estimate_unavailable")
      for (i in which(!outside & probabilities$status == "estimate_unavailable")) {
        note("The fitted curve did not provide a finite probability.",
             "estimation_failed", "Table 3", timepoint = probabilities$timepoint[i],
             group_level = probabilities$group_level[i])
      }
      no_ci <- !outside & probabilities$status == "ok" &
        (!is.finite(probabilities$conf.low) | !is.finite(probabilities$conf.high))
      for (i in which(no_ci)) {
        probabilities$status[i] <- "interval_unavailable"
        note("The fitted curve did not provide a finite confidence interval.",
             "unsupported", "Table 3", timepoint = probabilities$timepoint[i],
             group_level = probabilities$group_level[i])
      }
    }
  }
  probabilities$status[outside] <- "outside_observed_followup"

  plot <- if (!is.null(curve)) attempt(.cifreport_plot(curve, labels, scale, args$outcome.type,
                                                       conf.int, risk_scale), "Figure 1")

  test <- tests <- NULL
  ratios <- models <- vector("list", length(times))
  names(models) <- format(times, digits = 17L, trim = TRUE)
  observed_primary <- any(d$status == args$code.event1)
  if (length(level) > 1L) {
    if (!observed_primary) {
      note("No primary events were observed.", "unsupported", "Table 5")
    } else {
      test <- attempt(do.call(ciftest, c(list(formula = f,
        test = if (competing) "gray" else "logrank"), fit_args)), "Table 5")
      if (!is.null(test)) {
        tests <- attempt({
          x <- generics::tidy(test)
          if (!is.finite(x$statistic) || !is.finite(x$parameter) || x$parameter < 1 ||
              !is.finite(x$p.value)) stop("The curve test did not provide valid inference.")
          x[, c("method", "statistic", "parameter", "p.value"), drop = FALSE]
        }, "Table 5")
      }
    }
    if (valid_ref) {
      # Give polyreg a factor with the reference first so its exposure metadata
      # and coefficient ordering agree, including for multi-level exposures.
      model_data <- d
      model_data$group <- stats::relevel(d$group, ref)
      model_args <- c(args, list(data = quote(model_data)))
      for (i in seq_along(times)) {
        t <- times[i]
        observed <- d$status[d$followup <= t]
        reason <- if (t > min(support)) {
          "Time is beyond observed follow-up in at least one group."
        } else if (!args$code.event1 %in% observed) {
          "No primary events were observed by this time."
        } else if (competing && !args$code.event2 %in% observed) {
          "Competing-risk polyreg requires both event types by this time."
        }
        if (!is.null(reason)) {
          note(reason, "unsupported", "Table 4", timepoint = t)
          next
        }
        model <- attempt(do.call(polyreg, c(list(nuisance.model = Event(followup, status) ~ 1,
          exposure = "group", code.exposure.ref = ref, effect.measure1 = "RR", effect.measure2 = "RR",
          time.point = t, conf.int = conf.int, report.optim.convergence = TRUE), model_args)),
          "Table 4", timepoint = t)
        models[i] <- list(model)
        if (is.null(model)) next
        rows <- attempt(.cifreport_ratios(model, t), "Table 4", timepoint = t)
        if (!is.null(rows)) {
          rows$comparison_label <- unname(labels[rows$comparison])
          rows$reference_label <- unname(labels[rows$reference])
          for (j in which(rows$status != "ok")) {
            note(paste(rows$event_type[j], rows$status[j], sep = ": "),
                 "estimation_failed", "Table 4", timepoint = t, group_level = rows$comparison[j])
          }
        }
        ratios[i] <- list(rows)
      }
    }
  }
  list(counts = counts, probabilities = probabilities, ratios = do.call(rbind, ratios),
       tests = tests, plot = plot, fits = list(curve = curve, risk.ratios = models, test = test))
}

.cifreport_probabilities <- function(fit, times, levels, risk_scale) {
  table <- gtsummary::tbl_survfit(fit, times = times,
                                 type = if (risk_scale) "risk" else "survival")
  card <- table$cards$tbl_survfit[[1L]]
  get_stat <- function(name) {
    rows <- card$variable == "time" & card$stat_name == name
    time <- vapply(card$variable_level[rows], as.numeric, numeric(1L))
    group <- if (length(levels) == 1L) rep(levels, length(time)) else
      vapply(card$group1_level[rows], as.character, character(1L))
    value <- vapply(card$stat[rows], as.numeric, numeric(1L))
    index <- (match(group, levels) - 1L) * length(times) + match(time, times)
    if (length(index) != length(levels) * length(times) || anyNA(index) || anyDuplicated(index))
      stop("The numeric survival table did not match the requested groups and times.")
    value[order(index)]
  }
  data.frame(estimate = get_stat("estimate"), conf.low = get_stat("conf.low"),
             conf.high = get_stat("conf.high"))
}

.cifreport_curve_summary <- function(fit, times, levels, risk_scale) {
  x <- summary(fit, times = times, extend = TRUE)
  group <- if (length(levels) == 1L) rep(levels, length(x$time)) else
    sub("^group=", "", as.character(x$strata))
  index <- (match(group, levels) - 1L) * length(times) + match(x$time, times)
  if (length(index) != length(levels) * length(times) || anyNA(index) || anyDuplicated(index))
    stop("The curve summary did not match the requested groups and times.")
  result <- if (risk_scale) data.frame(estimate = 1 - x$surv,
    conf.low = 1 - x$upper, conf.high = 1 - x$lower) else
    data.frame(estimate = x$surv, conf.low = x$lower, conf.high = x$upper)
  result[order(index), , drop = FALSE]
}

.cifreport_ratios <- function(model, time) {
  convergence <- utils::tail(model$optim.inofo$converged.by, 1L)
  if (length(convergence) && (is.na(convergence) || convergence == "FALSE"))
    stop("The risk-ratio fit did not converge.")
  tidy <- generics::tidy(model, event = "both", exponentiate = TRUE)
  info <- model$reporting
  do.call(rbind, lapply(names(info$events), function(event) {
    terms <- info$events[[event]]$term
    x <- tidy[tidy$event == event, , drop = FALSE]
    index <- match(terms, x$term)
    if (anyNA(index) || length(index) != length(info$levels) - 1L)
      stop("The risk-ratio coefficients did not match the exposure levels.")
    x <- x[index, , drop = FALSE]
    x <- x[, c("estimate", "conf.low", "conf.high", "p.value")]
    x[] <- lapply(x, function(value) replace(value, !is.finite(value), NA_real_))
    status <- ifelse(is.na(x$estimate), "Risk ratio unavailable",
      ifelse(is.na(x$conf.low) | is.na(x$conf.high) | is.na(x$p.value),
             "Confidence interval or p-value unavailable", "ok"))
    data.frame(timepoint = time, event_type = if (event == "event2") "Competing event" else "Primary event",
      comparison = info$levels[-1L], reference = info$levels[1L], x, status)
  }))
}

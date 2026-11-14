.cifreport_english_list <- function(x) {
  n <- length(x)
  if (n <= 1L) return(paste(x, collapse = ""))
  if (n == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-n], collapse = ", "), ", and ", x[n])
}

.cifreport_baseline_text <- function(table, number, dataset, info) {
  counts <- .cifreport_english_list(sprintf("%s (n=%d)", names(info$counts), info$counts))
  leading <- c(sprintf("Table %d summarizes baseline characteristics of participants in %s.",
                       number, dataset),
               sprintf("The distribution of %s was %s.", info$exposure, counts))
  if (is.null(table)) return(leading[2L])
  body <- table$table_body
  p <- body$p.value
  conclusion <- if (any(is.finite(p))) {
    significant <- unique(body$var_label[which(is.finite(p) & p < 0.05)])
    if (length(significant)) {
      sprintf("Statistically significant differences across %s groups (p < 0.05) were observed for %s.",
              info$exposure, .cifreport_english_list(significant))
    } else {
      sprintf("No tested baseline characteristics differed significantly across %s groups (p < 0.05).",
              info$exposure)
    }
  } else {
    sprintf("No p-values were available for comparisons across %s groups.", info$exposure)
  }
  c(leading, conclusion)
}

.cifreport_outcome_text <- function(table, number, exposure) {
  if (!nrow(table)) return(character())
  text <- sprintf("Table %d shows the definitions of the outcomes used in the analysis.", number)
  sparse <- unique(table$outcome_id[which(table$primary < 5L)])
  if (length(sparse)) {
    text <- c(text, sprintf(
      "Interpret results for %s with caution because at least one %s group had fewer than five primary events.",
      .cifreport_english_list(sparse), exposure))
  }
  text
}

.cifreport_time_text <- function(table, number, level, followup, figures) {
  ranges <- vapply(names(followup), function(outcome) {
    times <- format(followup[[outcome]], trim = TRUE)
    sprintf("Observed follow-up for %s ranged from %s to %s.", outcome, times[1L], times[2L])
  }, character(1L))
  nonpositive <- names(followup)[vapply(followup, function(x) x[1L] <= 0, logical(1L))]
  caution <- if (length(nonpositive)) sprintf(
    "Use caution when interpreting results for %s because some participants had zero or negative follow-up times.",
    .cifreport_english_list(nonpositive))
  estimates <- if (any(is.finite(table$estimate))) sprintf(
    "Kaplan-Meier or Aalen-Johansen estimates and %s confidence intervals (CIs), where available, are shown in Table %d.",
    level, number)
  counts <- if (nrow(table))
    "The counts in the table include observed events through each time, including events at that time."
  curves <- if (length(figures)) sprintf("Kaplan-Meier or Aalen-Johansen curves are shown in %s.",
    .cifreport_english_list(paste("Figure", figures)))
  c(ranges, caution, estimates, counts, curves)
}

.cifreport_ratio_text <- function(table, number, level, reference = NULL) {
  if (!nrow(table)) {
    return(c(sprintf("No risk-ratio comparisons were available for Table %d.", number),
      if (length(reference)) sprintf("The specified reference group was %s.", reference)))
  }
  references <- unique(table$reference_label)
  leading <- sprintf("Table %d reports crude risk ratios using %s as the reference %s.",
    number, .cifreport_english_list(references), if (length(references) == 1L) "group" else "groups")
  interpretation <- "Risk ratios greater than 1 indicate a lower event risk in the reference group."
  available <- is.finite(table$conf.low) & is.finite(table$conf.high)
  outside <- table[available & (table$conf.high < 1 | table$conf.low > 1), , drop = FALSE]
  conclusion <- if (nrow(outside)) {
    comparisons <- sprintf("%s (%s) at time %s for %s versus %s", outside$outcome_id,
      tolower(outside$event_type), format(outside$timepoint, trim = TRUE),
      outside$comparison_label, outside$reference_label)
    sprintf("The %s confidence interval excluded 1 for %s.",
            level, .cifreport_english_list(comparisons))
  }
  c(leading, interpretation, conclusion)
}

.cifreport_test_text <- function(table, number) {
  available <- table[is.finite(table$p.value) & !is.na(table$method), , drop = FALSE]
  if (!nrow(available)) return(sprintf("No curve-comparison tests were available for Table %d.", number))
  methods <- unique(available$method)
  leading <- sprintf("Table %d compares survival or cumulative incidence curves using %s.",
                     number, .cifreport_english_list(methods))
  results <- vapply(methods, function(method) {
    outcomes <- unique(available$outcome_id[available$method == method & available$p.value < 0.05])
    if (length(outcomes)) {
      sprintf("For %s, p < 0.05 for %s.", method, .cifreport_english_list(outcomes))
    } else {
      sprintf("For %s, no outcomes had p < 0.05.", method)
    }
  }, character(1L))
  c(leading, results)
}

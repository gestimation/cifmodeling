.cifreport_plot <- function(curve, labels, scale, outcome.type, conf.int, risk_scale) {
  # ggsurvfit lays out the risk table on the active graphics device.
  grDevices::pdf(NULL, width = 8, height = 5.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  args <- list(formula_or_fit = curve, outcome.type = outcome.type, conf.int = conf.int,
               label.y = scale, type.y = if (risk_scale) "risk" else "surv")
  if (length(labels) > 1L) {
    strata <- paste0("group=", names(labels))
    args$level.strata <- strata
    args$order.strata <- strata
    args$label.strata <- stats::setNames(unname(labels), strata)
  }
  p <- as_ggplot(do.call(cifplot, args))
  if (length(labels) == 1L) p <- p + ggplot2::labs(subtitle = unname(labels))
  ggsurvfit::ggsurvfit_build(p)
}

.cifreport_render <- function(report) {
  work <- tempfile("cifreport-")
  if (!dir.create(work)) stop("Could not create the report working directory.", call. = FALSE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  template <- system.file("rmarkdown", "cifreport.Rmd", package = "cifmodeling", mustWork = TRUE)
  if (!file.copy(template, work)) stop("Could not stage the report template.", call. = FALSE)
  env <- list2env(list(report = report), parent = environment(.cifreport_render))
  staged <- rmarkdown::render(file.path(work, "cifreport.Rmd"), output_file = "report.html",
    output_dir = work, intermediates_dir = work, knit_root_dir = work, envir = env, quiet = TRUE)
  if (!file.copy(staged, report$output.file, overwrite = TRUE))
    stop("Could not write report to ", report$output.file, ".", call. = FALSE)
}

.cifreport_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  .cifreport_html_entities(gsub(">", "&gt;", x, fixed = TRUE))
}

.cifreport_baseline_html <- function(table) {
  # gtsummary's gt indentation inserts nonbreaking spaces that can become
  # literal <U+00A0> strings in non-UTF-8 locales. Indent with CSS instead.
  html_table <- gtsummary::as_gt(table, include = -"indent")
  rows <- which(table$table_body$row_type %in% c("level", "missing"))
  if (length(rows)) {
    html_table <- gt::tab_style(html_table,
      style = gt::cell_text(indent = gt::px(16)),
      locations = gt::cells_body(columns = "label", rows = rows))
  }
  # knitr captures cat() in the native locale. Character references preserve
  # Unicode labels when that locale cannot represent them.
  .cifreport_html_entities(gt::as_raw_html(html_table))
}

.cifreport_html_entities <- function(html) {
  html <- enc2utf8(html)
  non_ascii <- gregexpr("[^\\x00-\\x7f]+", html, perl = TRUE)
  regmatches(html, non_ascii) <- lapply(regmatches(html, non_ascii), function(parts) {
    vapply(parts, function(part) paste0(sprintf("&#%d;", utf8ToInt(part)), collapse = ""),
           character(1L))
  })
  html
}

.cifreport_print_table <- function(x, labels = NULL) {
  if (!nrow(x)) {
    cat("<p>No results available. See the analysis status log.</p>\n")
    return(invisible(NULL))
  }
  if (!is.null(labels)) x <- x[, names(labels), drop = FALSE]
  centered <- names(x)[vapply(x, is.numeric, logical(1L))]
  centered <- union(centered, intersect(names(x), c("code.event1", "code.event2", "code.censoring")))
  for (field in intersect(c("estimate", "conf.low", "conf.high"), names(x))) {
    value <- x[[field]]
    x[[field]] <- ifelse(is.na(value), "", if ("comparison_label" %in% names(x))
      sprintf("%.3f", value) else sprintf("%.2f%%", 100 * value))
  }
  for (field in intersect(c("statistic", "p.value"), names(x))) {
    value <- x[[field]]
    x[[field]] <- ifelse(is.na(value), "", sprintf("%.3f", value))
  }
  if ("timepoint" %in% names(x)) {
    times <- unique(x$timepoint[!is.na(x$timepoint)])
    labels_time <- as.character(times)
    ambiguous <- duplicated(labels_time) | duplicated(labels_time, fromLast = TRUE)
    labels_time[ambiguous] <- sprintf("%.17g", times[ambiguous])
    x$timepoint <- labels_time[match(x$timepoint, times)]
  }
  x[] <- lapply(x, function(value) replace(as.character(value), is.na(value), ""))
  table <- gt::gt(x)
  if (!is.null(labels)) table <- gt::cols_label(table, .list = as.list(labels))
  table <- gt::tab_options(table, column_labels.font.weight = "bold")
  table <- gt::cols_align(table, align = "left")
  if (length(centered)) table <- gt::cols_align(table, align = "center", columns = !!centered)
  cat(.cifreport_html_entities(gt::as_raw_html(table)), "\n\n", sep = "")
}

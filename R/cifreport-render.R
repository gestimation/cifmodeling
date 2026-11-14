.cifreport_plot <- function(curve, labels, scale, outcome.type, conf.int, risk_scale, cifplot.args) {
  # ggsurvfit lays out the risk table on the active graphics device.
  grDevices::pdf(NULL, width = 8, height = 5.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  args <- list(formula_or_fit = curve, outcome.type = outcome.type, conf.int = conf.int,
               label.y = scale, type.y = if (risk_scale) "risk" else "surv")
  if ("type.y" %in% names(cifplot.args)) args$label.y <- NULL
  args[names(cifplot.args)] <- cifplot.args
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

.cifreport_render <- function(report, formats, browser) {
  work <- tempfile("cifreport-")
  if (!dir.create(work)) stop("Could not create the report working directory.", call. = FALSE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  template <- system.file("rmarkdown", "cifreport.Rmd", package = "cifmodeling", mustWork = TRUE)
  if (!file.copy(template, work)) stop("Could not stage the report template.", call. = FALSE)
  env <- list2env(list(report = report), parent = environment(.cifreport_render))
  html <- rmarkdown::render(file.path(work, "cifreport.Rmd"), output_file = "report.html",
    output_dir = work, intermediates_dir = work, knit_root_dir = work, envir = env, quiet = TRUE)
  staged <- rep(html, length(formats))
  if ("docx" %in% formats) {
    word <- file.path(work, "report.docx")
    word_input <- file.path(work, "word.html")
    # Pandoc uses the HTML title as the Word title. Drop the repeated heading.
    source <- paste(readLines(html, encoding = "UTF-8", warn = FALSE), collapse = "\n")
    source <- sub("(?ms)^<h1>.*?</h1>", "", source, perl = TRUE)
    writeLines(source, word_input, useBytes = TRUE)
    rmarkdown::pandoc_convert(word_input, from = "html", to = "docx", output = word, verbose = FALSE)
    staged[formats == "docx"] <- .cifreport_format_word(word, work)
  }
  if ("pdf" %in% formats) {
    pdf <- file.path(work, "report.pdf")
    .cifreport_print_pdf(browser, html, pdf, work)
    staged[formats == "pdf"] <- pdf
  }
  for (i in seq_along(staged)) {
    if (!file.copy(staged[i], report$output.file[i], overwrite = TRUE))
      stop("Could not write report to ", report$output.file[i], ".", call. = FALSE)
  }
}

.cifreport_format_word <- function(word, work) {
  directory <- file.path(work, "word-package")
  entries <- utils::unzip(word, list = TRUE)$Name
  utils::unzip(word, exdir = directory)
  path <- file.path(directory, "word/document.xml")
  document <- xml2::read_xml(path)
  ns <- xml2::xml_ns(document)
  add <- function(parent, element, ...) xml2::xml_add_child(parent, paste0("w:", element), ...)

  # Letter with 0.7-inch side margins leaves 7.1 inches for the tables.
  section <- xml2::xml_find_first(document, "//w:sectPr", ns)
  add(section, "pgSz", `w:w` = "12240", `w:h` = "15840")
  add(section, "pgMar", `w:top` = "1080", `w:right` = "1008", `w:bottom` = "1080",
    `w:left` = "1008", `w:header` = "540", `w:footer` = "540", `w:gutter` = "0")

  for (table in xml2::xml_find_all(document, "//w:tbl", ns)) {
    columns <- xml2::xml_find_all(table, "w:tblGrid/w:gridCol", ns)
    properties <- xml2::xml_find_first(table, "w:tblPr", ns)
    margins <- add(properties, "tblCellMar", .where = 0)
    for (side in c("left", "right")) add(margins, side, `w:type` = "dxa", `w:w` = "55")
    rows <- xml2::xml_find_all(table, "w:tr", ns)
    sizes <- lapply(rows, function(row) {
      cells <- xml2::xml_find_all(row, "w:tc", ns)
      if (length(cells) == length(columns)) nchar(xml2::xml_text(cells), type = "width")
    })
    widths <- apply(do.call(rbind, sizes), 2L, function(x) max(5, min(30, max(x))))
    widths <- floor(10224 * widths / sum(widths))
    xml2::xml_set_attr(columns, "w:w", as.character(widths), ns = ns)
    layout <- xml2::xml_find_first(properties, "w:tblLayout", ns)
    if (inherits(layout, "xml_missing")) layout <- add(properties, "tblLayout", .where = 0)
    xml2::xml_set_attr(layout, "w:type", "fixed", ns = ns)
    for (run in xml2::xml_find_all(table, ".//w:r[w:t]", ns)) {
      formatting <- xml2::xml_find_first(run, "w:rPr", ns)
      if (inherits(formatting, "xml_missing")) formatting <- add(run, "rPr", .where = 0)
      add(formatting, "sz", `w:val` = "18")
      add(formatting, "szCs", `w:val` = "18")
    }
  }
  xml2::write_xml(document, path)
  output <- file.path(work, "formatted.docx")
  zip::zipr(output, entries[!endsWith(entries, "/")], root = directory, mode = "mirror")
  output
}

.cifreport_find_browser <- function() {
  configured <- Sys.getenv("CIFREPORT_BROWSER")
  candidates <- c(unname(Sys.which(c("google-chrome", "google-chrome-stable", "chromium",
    "chromium-browser", "microsoft-edge", "msedge", "chrome"))))
  if (.Platform$OS.type == "windows") {
    roots <- Sys.getenv(c("ProgramFiles", "ProgramFiles(x86)", "LOCALAPPDATA"))
    roots <- roots[nzchar(roots)]
    candidates <- c(candidates, file.path(roots, "Google/Chrome/Application/chrome.exe"),
      file.path(roots, "Microsoft/Edge/Application/msedge.exe"))
  } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
    candidates <- c(candidates, "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
      "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
      "/Applications/Chromium.app/Contents/MacOS/Chromium")
  }
  if (nzchar(configured)) candidates <- path.expand(configured)
  available <- file.exists(candidates) & !dir.exists(candidates)
  if (!any(available)) stop("PDF output requires Chrome, Chromium, or Microsoft Edge. ",
    "Set CIFREPORT_BROWSER to the browser executable.", call. = FALSE)
  candidates[which(available)[1L]]
}

.cifreport_print_pdf <- function(browser, html, pdf, work) {
  path <- normalizePath(html, winslash = "/", mustWork = TRUE)
  encoded <- gsub("%2F", "/", utils::URLencode(enc2utf8(path), reserved = TRUE, repeated = TRUE), fixed = TRUE)
  encoded <- sub("^([A-Za-z])%3A/", "\\1:/", encoded)
  url <- paste0(if (startsWith(path, "/")) "file://" else "file:///", encoded)
  args <- c("--headless", "--disable-gpu", "--no-first-run", "--no-pdf-header-footer",
    paste0("--user-data-dir=", file.path(work, "browser-profile")),
    paste0("--print-to-pdf=", pdf), url)
  log <- suppressWarnings(system2(browser, shQuote(args), stdout = TRUE, stderr = TRUE, timeout = 60))
  if (!is.null(attr(log, "status")) || !file.exists(pdf) || file.size(pdf) < 5L)
    stop("Could not generate the PDF report using ", basename(browser), ".\n",
      paste(utils::tail(log, 5L), collapse = "\n"), call. = FALSE)
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
      sprintf("%.3f", value) else sprintf("%.1f%%", 100 * value))
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

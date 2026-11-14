report_fixture <- function() {
  data.frame(time = rep(0:11, 3),
    status = c(1, 2, 1, 0, 2, 0, 1, 2, 0, 0, 0, 0,
               0, 1, 2, 0, 1, 0, 2, 1, 0, 0, 0, 0,
               2, 1, 2, 1, 0, 2, 1, 0, 0, 0, 0, 0),
    arm = factor(rep(c("A", "B", "C"), each = 12)),
    region = factor(rep(c("North", "South"), 18)),
    age = c(NA, 31:65), sex = factor(rep(c("F", "M"), 18)))
}

report_analyses <- function(d, outcomes = list(CIF = list(time = "time", event = "status")),
                            exposure = "arm", times = 4, subset.condition = NULL,
                            characteristics = c("age", "sex"), code.exposure.ref = NULL,
                            subgroup = NULL, outcome.type = NULL) {
  skip_if_not_installed("gtsummary", "2.6.1")
  skip_if_not_installed("survival")
  plan <- .cifreport_plan(subset.condition, exposure, outcomes, times, code.exposure.ref, subgroup,
    outcome.type = outcome.type)
  .cifreport_run(d, plan, characteristics, TRUE, 0.95, list(label = list(age ~ "Age in years")))
}

test_that("a combined report agrees with direct fixed-time and curve analyses", {
  d <- report_fixture()
  d$cohort <- ifelse(d$arm == "C", "C", "AB")
  d$survtime <- 2 * d$time
  d$death <- factor(ifelse(d$status == 0, "9", "2"))
  outcomes <- list(CIF = list(time = "time", event = "status", time.point = c(4, 7)),
    Survival = list(time = "survtime", event = "death",
      code.event1 = 2, code.censoring = 9, time.point = 8))
  result <- report_analyses(d, outcomes,
    exposure = c(Treatment = "arm", Region = "region"),
    code.exposure.ref = c(Treatment = "B"),
    subgroup = "cohort")
  expect_named(result$tables, c("All", "cohort = AB", "cohort = C"))
  expect_named(result$tables$All, c("Treatment", "Region"))
  tab <- result$tables$All$Treatment
  expect_s3_class(tab$table1, "tbl_summary")
  expect_true(any(tab$table1$table_body$var_label == "Age in years"))
  expect_true(any(tab$table1$table_body$row_type == "missing"))
  expect_equal(sum(tab$table2$total), 2 * nrow(d))
  expect_equal(tab$table2$time, rep(c("time", "survtime"), each = 3))
  expect_equal(tab$table2$event, rep(c("status", "death"), each = 3))
  expect_equal(tab$table2$code.event1, rep(c(1L, 2L), each = 3))
  expect_equal(tab$table2$code.event2, rep(c(2L, NA_integer_), each = 3))
  expect_equal(tab$table2$code.censoring, rep(c(0L, 9L), each = 3))
  expect_setequal(tab$table3$estimator, c("AJ", "KM"))
  expect_equal(tab$table2$primary + ifelse(is.na(tab$table2$competing), 0, tab$table2$competing) +
    tab$table2$censored, tab$table2$total)
  expect_named(result$figures$All$Treatment, names(outcomes))
  expect_true(all(vapply(result$figures$All$Treatment, inherits, logical(1), "ggplot")))
  expect_setequal(tab$table3$group_level, c("A", "B", "C"))
  expect_setequal(tab$table3$timepoint, c(4, 7, 8))
  expect_setequal(tab$table4$event_type, c("Primary event", "Competing event"))
  expect_setequal(tab$table4$comparison, c("A", "C"))
  expect_true(all(tab$table4$reference == "B"))
  expect_named(result$fits$All$Treatment$CIF$risk.ratios, c("4", "7"))
  for (id in names(outcomes)) {
    competing <- id == "CIF"
    formula <- if (competing) Event(time, status) ~ arm else Event(survtime, death) ~ arm
    type <- if (competing) "competing-risk" else "survival"
    codes <- if (competing) list() else list(code.event1 = 2, code.censoring = 9, code.event2 = 0)
    curve <- do.call(cifcurve, c(list(formula = formula, data = d, outcome.type = type), codes))
    x <- summary(curve, times = outcomes[[id]]$time.point, extend = TRUE)
    rows <- subset(tab$table3, outcome_id == id)
    expect_equal(rows$estimate, if (competing) 1 - x$surv else x$surv)
    expect_equal(rows$conf.low, if (competing) 1 - x$upper else x$lower)
    expect_equal(rows$conf.high, if (competing) 1 - x$lower else x$upper)
    test <- do.call(ciftest, c(list(formula = formula, data = d, outcome.type = type,
      test = if (competing) "gray" else "logrank"), codes))
    expect_equal(tab$table5$method[tab$table5$outcome_id == id], test$method)
    expect_equal(tab$table5$statistic[tab$table5$outcome_id == id], unname(test$statistic))
    expect_equal(tab$table5$parameter[tab$table5$outcome_id == id], unname(test$parameter))
    expect_equal(tab$table5$p.value[tab$table5$outcome_id == id], test$p.value)

    horizon <- outcomes[[id]]$time.point[1]
    model <- do.call(polyreg, c(list(nuisance.model = update(formula, . ~ 1),
      exposure = "arm", data = d, code.exposure.ref = "B", time.point = horizon,
      outcome.type = type), codes))
    x <- generics::tidy(model, event = "both", exponentiate = TRUE)
    ratios <- subset(tab$table4, outcome_id == id & timepoint == horizon)
    expect_equal(nrow(ratios), nrow(x))
    for (field in c("estimate", "conf.low", "conf.high", "p.value"))
      expect_equal(ratios[[field]], x[[field]])
  }
})

test_that("validation is scoped after population selection", {
  d <- report_fixture()
  d$arm[1L] <- NA
  d$cohort <- ifelse(is.na(d$arm), NA, as.character(d$arm))
  d$missing_time <- replace(d$time, 13L, NA)
  d$missing_event <- replace(d$status, 14L, NA)
  outcomes <- list(
    Valid = list(time = "time", event = "status"),
    Missing_time = list(time = "missing_time", event = "status"),
    Missing_event = list(time = "time", event = "missing_event"),
    Missing_column = list(time = "absent", event = "status"),
    Invalid_codes = list(time = "time", event = "status", code.event1 = 0),
    Invalid_times = list(time = "time", event = "status", time.point = Inf))
  result <- report_analyses(d, outcomes, exposure = c("arm", "region", "absent"), subgroup = "cohort")
  expect_named(result$tables$All, "region")
  expect_named(result$tables[["cohort = A"]], c("arm", "region"))
  expect_equal(unique(result$tables$All$region$table2$outcome_id), "Valid")
  a <- result$tables[["cohort = A"]]$arm
  expect_setequal(a$table2$outcome_id, c("Valid", "Missing_time", "Missing_event"))
  expect_equal(nrow(a$table4), 0)
  expect_equal(nrow(a$table5), 0)
  expect_s3_class(result$figures[["cohort = A"]]$arm$Valid, "ggplot")
  expect_true(all(result$status$status[result$status$analysis == "Outcome"] == "invalid_input"))
  expect_true(all(c("Grouping", "Subgroup", "Outcome", "Tables 4 and 5") %in% result$status$analysis))
})

test_that("subgroups stop after Table 2 for outcomes with zero events", {
  d <- report_fixture()
  d$cohort <- factor(ifelse(d$arm == "A", "A", "Other"), c("A", "Other", "Unused"))
  d$cohort[35] <- NA
  d$comp <- ifelse(d$arm == "A", 0L, d$status)
  d$surv <- ifelse(d$arm == "A", 0L, as.integer(d$status != 0L))
  d$alternate <- ifelse(d$status == 1L, 0L, d$status)
  outcomes <- list(CIF = list(time = "time", event = "comp"),
    Survival = list(time = "time", event = "surv", outcome.type = "survival"),
    Competing_only = list(time = "time", event = "alternate", outcome.type = "competing-risk"))
  result <- report_analyses(d, outcomes, subset.condition = "time <= 10",
    subgroup = "cohort")
  expect_named(result$tables, c("All", "cohort = A", "cohort = Other"))
  expect_equal(sum(result$tables$All$arm$table2$total), 3 * sum(d$time <= 10))
  a <- result$tables[["cohort = A"]]$arm
  expect_named(a, paste0("table", 1:5))
  expect_equal(a$table2$primary, c(0L, 0L, 0L))
  expect_equal(a$table2$competing, c(0L, NA_integer_, 3L))
  expect_equal(unique(result$tables[["cohort = A"]]$arm$table2$code.event2), c(2L, NA_integer_))
  expect_equal(unique(a$table3$outcome_id), "Competing_only")
  expect_null(result$fits[["cohort = A"]]$arm$CIF)
  expect_null(result$figures[["cohort = A"]]$arm$Survival)
  expect_s3_class(result$figures[["cohort = A"]]$arm$Competing_only, "ggplot")
  expect_setequal(result$status$outcome_id[result$status$analysis == "Tables 3-5 and Figure 1"],
                  c("CIF", "Survival"))
})

test_that("observed counts and estimates include exact ties without extrapolation", {
  d <- report_fixture()
  d$death <- d$status != 0
  d$time[d$arm == "A" & d$time == 11] <- 10
  times <- c(0, 1 - .Machine$double.eps, 1, 11, 12)
  r <- report_analyses(d, list(Survival = list(time = "time", event = "death",
    outcome.type = "survival")), times = times)
  tab <- r$tables$All$arm
  a <- subset(tab$table3, group_level == "A")
  b <- subset(tab$table3, group_level == "B")
  expect_equal(a$timepoint, times)
  expect_equal(a$primary[1], 1)
  expect_equal(b$primary[2:3], c(0, 1))
  expect_equal(b$estimate[2:3], c(1, 10 / 11))
  expect_equal(a$primary[4:5], rep(sum(d$death & d$arm == "A"), 2))
  expect_true(all(is.na(a$estimate[4:5])))
  expect_true(is.finite(b$estimate[4]))
  expect_true(all(is.na(b[b$timepoint == 12, c("estimate", "conf.low", "conf.high")])))
  expect_true(all(tab$table4$timepoint <= 10))
  expect_true(any(r$status$status == "unsupported" & r$status$timepoint == 12, na.rm = TRUE))
})

test_that("an outcome with no primary events retains its descriptive results", {
  d <- report_fixture()
  d$no_primary <- ifelse(d$status == 1, 0L, d$status)
  d$death <- d$status != 0
  r <- report_analyses(d, list(Zero = list(time = "time", event = "no_primary", outcome.type = "competing-risk"),
    Survival = list(time = "time", event = "death", outcome.type = "survival")), times = c(0, 4))
  tab <- r$tables$All$arm
  expect_true(all(subset(tab$table3, outcome_id == "Zero")$estimate == 0))
  expect_true(any(grepl("No primary events", r$status$message)))
  expect_true(nrow(subset(tab$table4, outcome_id == "Survival" & timepoint == 4)) > 0)
  expect_equal(tab$table5$outcome_id, "Survival")
})

test_that("the shared outcome type can be overridden or inferred per outcome", {
  d <- report_fixture()
  d$death <- d$status != 0L
  d <- d[d$status != 2L, ]
  common <- list(time = "time", event = "death")
  outcomes <- list(Shared = common, Survival = c(common, list(outcome.type = " SURV ")),
    Auto = c(common, list(outcome.type = "auto")))
  plan <- .cifreport_plan(NULL, "arm", outcomes, 4, outcome.type = " CR ")
  args <- lapply(plan$outcomes, function(x) .cifreport_outcome(x, d, d$arm)$arguments)
  expect_identical(vapply(args, `[[`, character(1L), "outcome.type"),
    c(Shared = "competing-risk", Survival = "survival", Auto = "survival"))
  expect_error(.cifreport_plan(NULL, "arm", outcomes, 4, outcome.type = "binomial"), "outcome.type")
})

test_that("a failed fit retains results at other horizons", {
  skip_if_not_installed("testthat", "3.2.0")
  original <- polyreg
  testthat::local_mocked_bindings(polyreg = function(..., time.point) {
    if (time.point == 4) stop("Example solver failure")
    original(..., time.point = time.point)
  }, .package = "cifmodeling")
  d <- report_fixture()
  d$status[d$status == 2 & d$time == 0] <- 0L
  r <- report_analyses(d, times = c(0, 4, 7))
  tab <- r$tables$All$arm
  expect_equal(unique(tab$table4$timepoint), 7)
  expect_equal(nrow(tab$table5), 1)
  expect_true(all(is.finite(tab$table3$estimate)))
  expect_true(any(r$status$status == "estimation_failed" & r$status$message == "Example solver failure"))
  expect_true(any(grepl("requires both event types", r$status$message)))
})

test_that("an absent reference affects only risk ratios", {
  r <- report_analyses(report_fixture(), exposure = c(Arm = "arm"), code.exposure.ref = "D")
  tab <- r$tables$All$Arm
  expect_true(nrow(tab$table3) > 0)
  expect_equal(nrow(tab$table4), 0)
  expect_equal(nrow(tab$table5), 1)
  expect_true(any(r$status$status == "invalid_input" & r$status$analysis == "Table 4"))
})

test_that("subset conditions use the package evaluation conventions", {
  d <- report_fixture()
  cutoff <- 4
  envir <- environment()
  for (condition in list("time >= cutoff", ~time >= cutoff, quote(time >= cutoff),
                         expression(time >= cutoff), d$time >= cutoff)) {
    expect_equal(.cifreport_population(condition, d, envir), d[d$time >= cutoff, ])
  }
  expect_equal(.cifreport_population(NULL, d, envir), d)
  expect_error(.cifreport_plan(list(All = NULL, Adult = ~age >= 18), "arm",
    list(CIF = list(time = "time", event = "status")), 4), "must define one population")
  expect_error(.cifreport_population(c("time > 0", "time < 9"), d, envir), "one expression")
  expect_equal(.cifreport_filename("trial"), paste0(format(Sys.Date(), "%Y%m%d"), "_trial_analysis_report.docx"))
})

test_that("complementary estimators reverse probability limits", {
  d <- report_fixture()
  d$death <- d$status != 0
  outcomes <- list(
    CIF = list(time = "time", event = "status", estimator = "1-AJ"),
    Survival = list(time = "time", event = "death", outcome.type = "survival", estimator = "1-KM"))
  r <- report_analyses(d, outcomes, exposure = "arm", code.exposure.ref = "absent")
  x <- r$tables$All$arm$table3
  aj <- summary(cifcurve(Event(time, status) ~ arm, d), times = 4)
  km <- summary(cifcurve(Event(time, death) ~ arm, d), times = 4)
  expect_equal(x$estimator, rep(c("1-AJ", "1-KM"), each = 3))
  expect_equal(x$estimate, c(aj$surv, 1 - km$surv))
  expect_equal(x$conf.low, c(aj$lower, 1 - km$upper))
  expect_equal(x$conf.high, c(aj$upper, 1 - km$lower))
  expect_true(all(vapply(r$figures$All$arm, inherits, logical(1), "ggplot")))
})

test_that("exposure display levels and order agree across tables and plots", {
  skip_if_not_installed("gtsummary", "2.6.1")
  d <- report_fixture()
  d$cohort <- ifelse(d$arm == "C", "C", "AB")
  plan <- .cifreport_plan(NULL, c(Treatment = "arm"),
    list(CIF = list(time = "time", event = "status")), 4,
    code.exposure.ref = "B", subgroup = "cohort",
    strata = list(level.strata = list(Treatment = c("A", "B", "C")),
                  order.strata = c("C", "A"), label.strata = c("Alpha", "Beta", "Gamma")))
  result <- .cifreport_run(d, plan, "age", FALSE, 0.95, list())
  tab <- result$tables$All$Treatment
  header <- tab$table1$table_styling$header
  expect_equal(as.character(header$modify_stat_level[grepl("^stat_", header$column)]),
               c("Gamma", "Alpha", "Beta"))
  expect_equal(tab$table2$group_level, c("C", "A", "B"))
  expect_equal(tab$table2$group_label, c("Gamma", "Alpha", "Beta"))
  expect_equal(tab$table3$group_label, tab$table2$group_label)
  expect_equal(unique(tab$table4$comparison_label), c("Gamma", "Alpha"))
  expect_equal(unique(tab$table4$reference), "B")
  expect_equal(unique(tab$table4$reference_label), "Beta")
  expect_equal(result$tables[["cohort = AB"]]$Treatment$table2$group_label, c("Alpha", "Beta"))
  expect_equal(result$tables[["cohort = C"]]$Treatment$table2$group_label, "Gamma")
  plot <- result$figures$All$Treatment$CIF
  expect_s3_class(plot, "patchwork")
  grob_labels <- function(x) {
    c(x$label, unlist(lapply(x$grobs, grob_labels)), unlist(lapply(x$children, grob_labels)))
  }
  display <- c("Gamma", "Alpha", "Beta")
  expect_identical(intersect(grob_labels(patchwork::patchworkGrob(plot)), display), display)
  labels <- .cifreport_strata(plan$exposure$Treatment, d, d)
  expect_equal(labels, c(C = "Gamma", A = "Alpha", B = "Beta"))
  expect_error(.cifreport_strata(list(variable = "arm", order.strata = "Z"), d, d), "unknown")
  expect_error(.cifreport_strata(list(variable = "arm", label.strata = c(A = "B")), d, d), "distinct")
})

test_that("report filenames can be omitted, supplied, or disabled", {
  skip_if_not_installed("testthat", "3.2.0")
  testthat::local_mocked_bindings(.cifreport_dependencies = function(...) NULL,
    .cifreport_run = function(...) list(tables = list()),
    .cifreport_render = function(report, formats, browser) {
      for (file in report$output.file) writeLines("report", file)
    })
  d <- report_fixture()
  outcomes <- list(CIF = list(time = "time", event = "status"))
  directory <- tempfile()
  dir.create(directory)
  old_directory <- setwd(directory)
  on.exit({
    setwd(old_directory)
    unlink(directory, recursive = TRUE)
  }, add = TRUE)
  default <- cifreport(d, outcomes, "arm", time.point = 4)
  expect_identical(basename(default$output.file), .cifreport_filename("d"))
  expect_true(file.exists(default$output.file))
  relative <- cifreport(d, outcomes, "arm", time.point = 4, filename.cifreport = "report.htm")
  expect_identical(relative$output.file, normalizePath("report.htm", winslash = "/"))
  legacy <- cifreport(d, outcomes, "arm", time.point = 4, output.file = "legacy.html")
  expect_true(file.exists(legacy$output.file))
  disabled <- cifreport(d, outcomes, "arm", time.point = 4, data.name = "disabled", filename.cifreport = NULL)
  expect_true("output.file" %in% names(disabled))
  expect_null(disabled$output.file)
  expect_false(file.exists(.cifreport_filename("disabled")))
  expect_error(cifreport(d, outcomes, "arm", filename.cifreport = "missing/report.html"), "does not exist")
  expect_error(cifreport(d, outcomes, "arm", filename.cifreport = "report.png"), "filenames ending")
  files <- c("second.docx", "first.html")
  multiple <- cifreport(d, outcomes, "arm", time.point = 4, filename.cifreport = files)
  expect_identical(multiple$output.file, normalizePath(files, winslash = "/"))
  expect_error(cifreport(d, outcomes, "arm", filename.cifreport = c("report.htm", "./report.htm")),
    "distinct paths")
  expect_error(cifreport(d, outcomes, "arm", filename.cifreport = NULL, output.file = "report.html"),
    "Supply only one")
})

test_that("cifreport writes self-contained HTML and returns an invisible result", {
  skip_if_not_installed("gtsummary", "2.6.1")
  skip_if_not_installed("gt")
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("xml2")
  skip_if_not_installed("zip")
  skip_if_not(rmarkdown::pandoc_available())
  file <- tempfile(fileext = ".html")
  word_file <- tempfile(fileext = ".docx")
  on.exit(unlink(c(file, word_file)))
  writeLines("previous report", file)
  d <- subset(report_fixture(), arm != "C")
  d$cohort <- ifelse(d$time == 11, "Zero", "Other")
  attr(d$cohort, "label") <- "Study subset"
  attr(d$cohort, "labels") <- c("With events" = "Other", "Without events" = "Zero")
  attr(d$arm, "label") <- "Treatment arm"
  attr(d$arm, "labels") <- c("Original Alpha" = "A", "Arm Beta" = "B")
  d$response <- structure(rep(c(0, 1, 2, NA), length.out = nrow(d)),
    labels = c(No = 0, Yes = 1), label = "Baseline response")
  d$marker <- as.integer(d$arm) * 100 + d$time
  d$negative_time <- d$time - 1
  plot_args <- list(label.x = "Follow-up in months", legend.position = "bottom",
    breaks.x = c(0, 4, 8), add.risktable = FALSE, add.conf = FALSE,
    add.competing.risk.mark = TRUE)
  result <- withVisible(cifreport(d, list(CIF = list(time = "time", event = "status"),
    Invalid = list(time = "negative_time", event = "status")),
    "arm", c("age", "response", "marker"), time.point = 4,
    filename.cifreport = c(file, word_file), title = "CIF <example>",
    table1.args = list(label = list(age ~ "\u5e74\u9f62", marker ~ "Example marker"),
      type = list(response ~ "continuous", marker ~ "continuous")),
    subgroup = "cohort",
    level.strata = c("A", "B"), order.strata = c("B", "A"),
    label.strata = c(A = "Arm Alpha"), cifplot.args = plot_args))
  expect_false(result$visible)
  expect_identical(result$value$specifications$cifplot.args, plot_args)
  plots <- Filter(Negate(is.null), lapply(result$value$figures, function(x) x$arm$CIF))
  expect_length(plots, 2L)
  for (plot in plots) {
    expect_false(inherits(plot, "patchwork"))
    expect_identical(plot$labels$x, plot_args$label.x)
    expect_equal(plot$scales$get_scales("x")$breaks, plot_args$breaks.x)
    expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomPoint") &&
      identical(layer$aes_params$shape, 16), logical(1L))))
  }
  expect_error(cifreport(d, list(CIF = list(time = "time", event = "status")), "arm",
    cifplot.args = list(label.strata = c("A", "B"))), "Unknown cifplot.args field: label.strata")
  expect_equal(result$value$output.file, normalizePath(c(file, word_file), winslash = "/"))
  expect_true("word/document.xml" %in% utils::unzip(word_file, list = TRUE)$Name)
  html <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  cells <- xml2::xml_text(xml2::xml_find_all(xml2::read_html(file), "//table//td"), trim = TRUE)
  percentages <- unlist(regmatches(cells, gregexpr("[0-9]+(?:\\.[0-9]+)?%(?=$|\\)$)", cells, perl = TRUE)))
  expect_true(length(percentages) > 0L && all(grepl("\\.[0-9]%$", percentages)))
  first_table <- strsplit(html, "Table 2.", fixed = TRUE)[[1L]][1L]
  expect_match(first_table, "p-value", fixed = TRUE)
  text <- gsub("\u00a0", " ", html, fixed = TRUE)
  text <- gsub("[[:space:]]+", " ", text)
  for (sentence in c(
    "Table 1 summarizes baseline characteristics of participants in d.",
    "The distribution of Treatment arm was Arm Beta (n=12) and Arm Alpha (n=12).",
    "Table 2 shows the definitions of the outcomes used in the analysis.",
    "Interpret results for CIF with caution because at least one Treatment arm group had fewer than five primary events.",
    "Observed follow-up for CIF ranged from 0 to 11.",
    "Observed follow-up for Invalid ranged from -1 to 10.",
    "Use caution when interpreting results for CIF and Invalid because some participants had zero or negative follow-up times.",
    "Kaplan-Meier or Aalen-Johansen estimates and 95.0% confidence intervals (CIs), where available, are shown in Table 3.",
    "Kaplan-Meier or Aalen-Johansen curves are shown in Figure 2.",
    "Table 4 reports crude risk ratios using Arm Alpha as the reference group.",
    "Risk ratios greater than 1 indicate a lower event risk in the reference group.",
    "Table 5 compares survival or cumulative incidence curves using",
    "no outcomes had p &lt; 0.05.")) {
    expect_true(grepl(sentence, text, fixed = TRUE), info = sentence)
  }
  for (removed in c("None of the available 95% confidence intervals excluded 1.",
    "Ratios compare the event probability by the specified time with the reference group.",
    "Tests compare curves across follow-up.")) {
    expect_false(grepl(removed, text, fixed = TRUE))
  }
  expect_true(grepl("</table> </div> <p> Table 1 summarizes", text, fixed = TRUE))
  expect_true(grepl("</table> </div> <p> Table 2 shows", text, fixed = TRUE))
  expect_true(grepl("</table> </div> <p> Observed follow-up", text, fixed = TRUE))
  baseline_text <- regmatches(text, regexpr("<p> Table 1 summarizes[^<]+</p>", text))
  expect_true(length(baseline_text) == 1L && grepl("Example marker", baseline_text, fixed = TRUE))
  subgroup_notes <- regmatches(text, gregexpr(
    "<p> [^<]*?\\(subgroup analysis of Study subset (?:With|Without) events\\)\\.",
    text, perl = TRUE))[[1L]]
  expect_length(subgroup_notes, 7L)
  for (label in c("Table 1.", "Table 2.", "Table 3.", "Table 4.", "Table 5.", "Figure 1.",
                 "CIF &lt;example&gt;", "Analysis status log", "data:image/png;base64,", "\u5e74\u9f62",
                 "Estimator", "Exposure", "Arm Alpha", "Arm Beta"))
    expect_true(grepl(label, html, fixed = TRUE), info = label)
  for (label in c("Outcome definitions", "Probability scale", "Risk scale"))
    expect_false(grepl(label, html, fixed = TRUE), info = label)
  headings <- regmatches(html, gregexpr("<h2>\\s*Table [0-9]+\\.[^<]+</h2>", html))[[1L]]
  headings <- trimws(gsub("[[:space:]]+", " ",
    sub("</h2>$", "", sub("^<h2>\\s*", "", headings))))
  exposure_title <- "Treatment arm (Arm Beta, Arm Alpha)"
  titles <- c(paste("Description of baseline characteristics according to", exposure_title),
    paste("Definition of outcomes and", exposure_title),
    paste("Description of outcomes at specified time-points according to", exposure_title),
    paste("Crude risk ratios for the primary and competing events and", exposure_title, "at specified time-points"),
    paste("Results of tests for the primary event and", exposure_title))
  expect_equal(headings, paste0("Table ", 1:12, ". ", c(titles,
    paste(titles, "in Study subset With events"),
    paste(titles[1:2], "in Study subset Without events"))))
  styled_tables <- gregexpr('<table class="gt_table"', html, fixed = TRUE)[[1L]]
  expect_length(styled_tables, length(headings) + 1L)
  gaps <- gregexpr('</h2>\\s*<div class="cifreport-table-gap-after" aria-hidden="true">', html, perl = TRUE)[[1L]]
  expect_length(gaps, length(headings))
  gaps_before <- gregexpr('<div class="cifreport-table-gap" aria-hidden="true">\\s*</div>\\s*<h2>\\s*Table [0-9]+\\.',
    html, perl = TRUE)[[1L]]
  expect_length(gaps_before, length(headings))
  baseline <- result$value$tables$All$arm$table1$table_body
  expect_equal(unname(baseline$p.value[baseline$variable == "age" & baseline$row_type == "label"]),
               stats::wilcox.test(age ~ arm, data = d)$p.value)
  response <- baseline[baseline$variable == "response", ]
  expect_true(all(response$var_type == "categorical"))
  expect_equal(response$label[response$row_type == "level"], c("No", "Yes", "2"))
  expect_true(any(response$row_type == "missing"))
  expect_equal(result$value$tables$All$arm$table2$group_label, c("Arm Beta", "Arm Alpha"))
  expect_false(grepl("Original Alpha", html, fixed = TRUE))
  expect_false(grepl("All / arm", html, fixed = TRUE))
  expect_true(grepl("No events were observed in this subgroup.", html, fixed = TRUE))
  expect_false(grepl("(<|&lt;)U\\+00A0", html))
  expect_false(dir.exists(sub("\\.html$", "_files", file)))
})

test_that("comparison narratives apply strict CI and p-value thresholds", {
  ratios <- data.frame(outcome_id = c("Primary result", "Competing result", "Touches one", "Missing"),
    event_type = c("Primary event", "Competing event", "Primary event", "Primary event"),
    timepoint = c(36, 60, 36, 60), comparison_label = "DES", reference_label = "Placebo",
    conf.low = c(0.5, 1.1, 1, NA_real_), conf.high = c(0.99, 1.8, 1.5, 0.9))
  text <- .cifreport_ratio_text(ratios, 9, "90%")
  expect_match(text[1L], "using Placebo as the reference group", fixed = TRUE)
  expect_match(text[2L], "lower event risk in the reference group", fixed = TRUE)
  expect_match(text[3L], "90% confidence interval excluded 1", fixed = TRUE)
  expect_match(text[3L], "Primary result (primary event) at time 36 for DES versus Placebo", fixed = TRUE)
  expect_match(text[3L], "Competing result (competing event) at time 60 for DES versus Placebo", fixed = TRUE)
  expect_false(grepl("Touches one|Missing", text[3L]))
  expect_length(.cifreport_ratio_text(ratios[3:4, ], 9, "90%"), 2L)
  expect_match(.cifreport_ratio_text(ratios[0L, ], 9, "95%"), "No risk-ratio comparisons", fixed = TRUE)

  tests <- data.frame(outcome_id = c("Detected", "At threshold", "Another outcome"),
    method = c("Log-rank test", "Log-rank test", "Gray's test"), p.value = c(0.0499, 0.05, 0.3))
  text <- .cifreport_test_text(tests, 10)
  expect_match(text[2L], "For Log-rank test, p < 0.05 for Detected.", fixed = TRUE)
  expect_match(text[3L], "For Gray's test, no outcomes had p < 0.05.", fixed = TRUE)
  expect_false(any(grepl("At threshold", text, fixed = TRUE)))
  expect_match(.cifreport_test_text(tests[0L, ], 10), "No curve-comparison tests", fixed = TRUE)
})

test_that("Table 1 preserves category labels and indentation in a C locale", {
  skip_if_not_installed("gtsummary", "2.6.1")
  skip_if_not_installed("gt")
  old <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", old))
  Sys.setlocale("LC_CTYPE", "C")

  d <- data.frame(arm = c("A", "A", "B", "B"),
    symptom = factor(c("\u75c7\u72b6\u3042\u308a", NA, "Other", "Other")))
  table <- gtsummary::tbl_summary(d, by = "arm", include = "symptom",
    label = list(symptom ~ "\u75c7\u72b6"), missing = "ifany")
  html <- .cifreport_baseline_html(table)
  label <- paste0(sprintf("&#%d;", utf8ToInt("\u75c7\u72b6\u3042\u308a")), collapse = "")
  expect_true(grepl(label, html, fixed = TRUE))
  expect_true(grepl("Unknown", html, fixed = TRUE))
  expect_true(grepl("text-indent: 16px", html, fixed = TRUE))
  expect_false(grepl("(<|&lt;)U\\+00A0", html))
})

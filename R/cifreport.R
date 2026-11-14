#' Create a survival and competing-risk report
#'
#' @param data A data frame with one row per subject.
#' @param outcomes A named list of outcome definitions. Each is a list with
#'   column names `time` and `event`, and optionally `outcome.type`
#'   to override the shared outcome type, `code.event1`
#'   (default 1), `code.event2` (default 2, for competing risks), `code.censoring`
#'   (default 0), and `time.point` to override the shared time points.
#'   Use `estimator = "KM"` or `"1-KM"` for survival and `"AJ"` or
#'   `"1-AJ"` for competing risks. Defaults are `"KM"` and `"AJ"`.
#'   The name of each definition is its outcome label.
#' @param exposure Character vector of grouping column names, or a named list
#'   whose entries are column names or lists with `variable` and optional
#'   `code.exposure.ref`, `level.strata`, `order.strata`, and `label.strata`.
#'   Vector or list names identify exposure definitions.
#' @param code.exposure.ref Reference level for the risk ratios. Supply one
#'   level for all exposures or a named vector/list keyed by exposure section
#'   labels. A reference in an exposure definition takes precedence. By
#'   default the first observed factor level is the reference.
#' @param level.strata Expected exposure codes. Unnamed labels are matched to
#'   these codes by position. If omitted, the levels in the data are used.
#' @param order.strata Exposure codes in display order. Codes omitted from this
#'   vector follow in their default order.
#' @param label.strata Display labels, either named by the original exposure
#'   codes or given in the order of `level.strata`. Supply any strata argument
#'   as a vector, or as a named list keyed by exposure definition for separate
#'   settings. Exposure-specific settings in `exposure` take precedence.
#'   The same labels and order are used in tables and figures. Use original
#'   exposure codes for `code.exposure.ref`.
#' @param characteristics Character vector of baseline columns for Table 1.
#' @param time.point Shared numeric evaluation times. May be omitted when each
#'   outcome supplies its own times.
#' @param outcome.type Shared outcome type, `"survival"` or `"competing-risk"`.
#'   The default `NULL` (also `"auto"`) uses the same rule as [cifcurve()]:
#'   more than two observed event codes indicate competing risks. Otherwise
#'   the outcome is treated as survival. Detection uses the population selected
#'   by `subset.condition` before splitting into subgroups. Specify
#'   `"competing-risk"` explicitly if an expected event code is absent.
#' @param cifplot.args Named list of [cifplot()] display arguments shared by all
#'   figures. Use the individual arguments for axis labels and scales, legends,
#'   confidence bands, risk and estimate tables, event marks, quantiles, and
#'   styling. For example, `list(label.x = "Months", add.risktable = FALSE,
#'   legend.position = "bottom")`. Set stratum labels and order through the
#'   report's strata arguments. Use `add.conf` to show or hide the bands and
#'   `conf.int` at report level to set their confidence level. By default the
#'   y-axis follows the outcome's estimator. `type.y` changes the figure's
#'   scale only. Competing-event marks use the observed event times for each
#'   competing-risk outcome. Settings apply to single-curve figures, with
#'   multiple exposure groups overlaid when present.
#' @param subset.condition One condition defining the analysis population.
#'   Use a logical vector, character condition, one-sided formula, or quoted
#'   expression. Missing condition values are excluded. `NULL` includes all rows.
#' @param subgroup Name of one categorical column. In addition to the analysis
#'   population as a whole, the report includes each observed level of this
#'   column within that population. Missing values contribute only to the
#'   analysis population as a whole. Numeric category codes are supported.
#' @param filename.cifreport One or more paths ending in `.docx`, `.pdf`,
#'   `.html`, or `.htm`. Relative paths are resolved from the working directory.
#'   Parent directories must exist. Existing files are overwritten after rendering.
#'   If omitted, saves `yyyymmdd_<data.name>_analysis_report.docx` in the working
#'   directory. Supply `NULL` to return the results without writing any files.
#'   Multiple filenames share one analysis. Word output uses editable tables
#'   with 9-point text and requires the `xml2` and `zip` packages.
#'   PDF requires Chrome, Chromium, or Microsoft Edge. Set the `CIFREPORT_BROWSER`
#'   environment variable to the browser executable if it is not found automatically.
#' @param output.file Older name for `filename.cifreport`. Supply only one.
#' @param data.name Dataset name used in the default filename and title. By default it
#'   uses the supplied data object's name, or `"data"` for an inline expression.
#' @param title Report title. The default is
#'   `Statistical Analysis Report of <data.name> on <yyyymmdd>`.
#' @param overall Include an overall column in Table 1 when there are multiple
#'   groups.
#' @param conf.int Pointwise confidence level for the curve and risk ratio.
#' @param table1.args Named list of additional arguments to [gtsummary::tbl_summary()],
#'   such as `label` and `type`. Column `label` attributes are also supported.
#'   Characteristics with a named `labels` attribute are summarized as
#'   categorical using those value labels.
#'   When comparison is possible, Table 1 includes p-values from the default
#'   tests of [gtsummary::add_p()] across the observed exposure groups.
#'
#' @details
#' Each population and grouping variable produces a section containing a
#' baseline table, full-follow-up event counts, cumulative counts and
#' probabilities at requested times, crude risk ratios, curve comparison
#' tests, and outcome-specific plots. Observed cumulative counts include
#' events at exactly the requested time. Table 3 and, by default, the figure
#' use the selected estimator scale. Percentages in the tables are displayed
#' to one decimal place by default. KM is survival, 1-KM is event probability,
#' AJ is primary-event cumulative incidence, and 1-AJ is its complement.
#' Probabilities and confidence limits beyond
#' the last observed follow-up in a group are omitted. Risk ratios compare
#' event probabilities (1 minus survival for survival outcomes), with one
#' unadjusted fit per time. Competing-risk fits jointly model both causes.
#' Log-rank and Gray tests compare curves over follow-up.
#' Table titles use the exposure's variable label (`label` attribute) when
#' available and include its displayed levels when it has value labels.
#' Subgroup titles use the variable and value labels when available.
#' For exposure levels, `label.strata` takes precedence over the data's
#' value labels. Unlabelled codes are displayed as given.
#' Within a subgroup level, an outcome with zero observed events over the
#' entire follow-up stops after Table 2. For competing risks, this requires
#' zero primary and zero competing events across the exposure groups.
#'
#' Missing baseline values are summarized by default. Missing group values
#' skip the population/group combination. Missing outcome values skip that
#' outcome within the population. Single-level groups retain descriptive
#' results. A missing reference level affects risk ratios only. The analysis
#' log records invalid, skipped, unsupported and failed analyses and warnings.
#'
#' The report supports unweighted, right-censored data with one row per
#' subject, one primary event code and, for competing risks, one competing
#' event code. File output requires optional reporting packages and Pandoc. See
#' `vignette("v5_reporting_results")` for examples and supported configurations.
#'
#' @return Invisibly, a list with `output.file` (absolute paths in the order
#'   requested, or `NULL` if no file was saved), `specifications`,
#'   `tables`, `figures`, `fits`, `status` and `versions`. Tables,
#'   figures and fits are nested by population and grouping definition.
#'   Table 1 is a gtsummary object. Tables 2--5 are numeric result data
#'   frames with outcome and group identifiers. Group codes and display
#'   labels are available in Tables 2--4. Figures and fits are
#'   indexed by outcome, with risk-ratio fits indexed by time.
#' @examples
#' if (requireNamespace("gtsummary", quietly = TRUE) &&
#'     utils::packageVersion("gtsummary") >= "2.6.1") {
#'   data(diabetes.complications)
#'   report <- cifreport(
#'     diabetes.complications,
#'     outcomes = list(Retinopathy = list(time = "t", event = "epsilon")),
#'     exposure = "fruitq1", characteristics = c("age", "sex"),
#'     time.point = 8, filename.cifreport = NULL
#'   )
#' }
#' @seealso [cifcurve()], [cifplot()], [polyreg()], [ciftest()]
#' @md
#' @export
cifreport <- function(data, outcomes, exposure, characteristics = character(),
                      time.point = NULL, subset.condition = NULL,
                      filename.cifreport = NULL, title = NULL,
                      overall = FALSE, conf.int = 0.95, table1.args = list(),
                      code.exposure.ref = NULL, data.name = NULL, subgroup = NULL,
                      label.strata = NULL, order.strata = NULL, level.strata = NULL,
                      outcome.type = NULL, cifplot.args = list(), output.file = NULL) {
  data_expr <- substitute(data)
  if (is.null(data.name)) data.name <- if (is.symbol(data_expr)) as.character(data_expr) else "data"
  if (!.cifreport_string(data.name)) stop("`data.name` must be a nonempty string.", call. = FALSE)
  date <- format(Sys.Date(), "%Y%m%d")
  if (is.null(title)) title <- paste("Statistical Analysis Report of", data.name, "on", date)
  if (!is.data.frame(data) || anyDuplicated(names(data)) ||
      anyNA(names(data)) || any(!nzchar(names(data))))
    stop("`data` must be a data frame with unique, nonempty column names.", call. = FALSE)
  if (!is.character(characteristics) || anyNA(characteristics) ||
      any(!nzchar(characteristics)) || anyDuplicated(characteristics))
    stop("`characteristics` must contain unique column names.", call. = FALSE)
  if (!is.logical(overall) || length(overall) != 1L || is.na(overall))
    stop("`overall` must be TRUE or FALSE.", call. = FALSE)
  output_check_level(conf.int)
  .cifreport_check_fields(table1.args, arg = "table1.args")
  if (any(names(table1.args) %in% c("data", "by", "include")))
    stop("`table1.args` cannot set data, by or include.", call. = FALSE)
  .cifreport_check_fields(cifplot.args, c(
    "type.y", "label.x", "label.y", "label.legend", "limits.x", "limits.y", "breaks.x", "breaks.y",
    "use.coord.cartesian", "add.conf", "add.risktable", "add.estimate.table",
    "symbol.risk.table", "font.size.risk.table", "add.censor.mark", "shape.censor.mark",
    "size.censor.mark", "add.competing.risk.mark", "shape.competing.risk.mark",
    "size.competing.risk.mark", "add.intercurrent.event.mark", "intercurrent.event.time",
    "shape.intercurrent.event.mark", "size.intercurrent.event.mark", "add.quantile", "level.quantile",
    "style", "palette", "linewidth", "linetype", "font.family", "font.size", "legend.position"
  ), "cifplot.args")
  if (!.cifreport_string(title)) stop("`title` must be a nonempty string.", call. = FALSE)
  if (missing(filename.cifreport)) {
    filename.cifreport <- if (missing(output.file)) .cifreport_filename(data.name, date) else
      output.file %||% .cifreport_filename(data.name, date, extension = "html")
  } else if (!missing(output.file)) {
    stop("Supply only one of `filename.cifreport` and `output.file`.", call. = FALSE)
  }
  filename.cifreport <- .cifreport_paths(filename.cifreport)
  formats <- tolower(tools::file_ext(filename.cifreport))
  save_files <- length(formats) > 0L
  .cifreport_dependencies(save_files, word = "docx" %in% formats)
  browser <- if ("pdf" %in% formats) .cifreport_find_browser()

  strata <- list(label.strata = label.strata, order.strata = order.strata, level.strata = level.strata)
  plan <- .cifreport_plan(subset.condition, exposure, outcomes, time.point,
                          code.exposure.ref, subgroup, strata, outcome.type)
  result <- .cifreport_run(as.data.frame(data), plan, characteristics, overall, conf.int, table1.args,
                           subset.env = parent.frame(), cifplot.args = cifplot.args)
  result["output.file"] <- list(filename.cifreport)
  result$specifications <- c(plan, list(data.name = data.name, characteristics = characteristics, overall = overall,
                                       conf.int = conf.int, table1.args = table1.args,
                                       cifplot.args = cifplot.args, title = title))
  packages <- c("cifmodeling", "survival", "ggplot2", "ggsurvfit", "gtsummary",
                "cardx", "broom", "broom.helpers", "gt", "knitr", "rmarkdown")
  result$versions <- c(R = as.character(getRversion()),
    vapply(stats::setNames(packages, packages), function(package) {
      if (requireNamespace(package, quietly = TRUE))
        as.character(getNamespaceVersion(package)) else NA_character_
    }, character(1L)), pandoc = if (save_files) as.character(rmarkdown::pandoc_version()) else NA_character_)
  if (save_files) .cifreport_render(result, formats, browser)
  invisible(result)
}

.cifreport_dependencies <- function(save_files, word = FALSE) {
  for (package in c("gtsummary", if (save_files) c("gt", "knitr", "rmarkdown"),
                   if (word) c("xml2", "zip"))) {
    if (!requireNamespace(package, quietly = TRUE))
      stop("Reporting requires the ", package, " package.", call. = FALSE)
  }
  if (package_version(getNamespaceVersion("gtsummary")) < "2.6.1")
    stop("Reporting requires gtsummary >= 2.6.1.", call. = FALSE)
  if (save_files && !rmarkdown::pandoc_available())
    stop("Reporting requires Pandoc. See rmarkdown::find_pandoc().", call. = FALSE)
}

.cifreport_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.cifreport_filename <- function(data.name, date = format(Sys.Date(), "%Y%m%d"), extension = "docx") {
  safe_name <- gsub('[<>:"/\\\\|?*[:cntrl:]]', "_", data.name)
  paste0(date, "_", safe_name, "_analysis_report.", extension)
}

.cifreport_paths <- function(paths) {
  if (is.null(paths)) return(NULL)
  if (!is.character(paths) || !is.null(dim(paths)) || !length(paths) || anyNA(paths) ||
      any(!nzchar(paths)) || any(!tolower(tools::file_ext(paths)) %in% c("docx", "pdf", "html", "htm")))
    stop("`filename.cifreport` must contain filenames ending in .docx, .pdf, .html, or .htm.", call. = FALSE)
  paths <- path.expand(paths)
  directories <- dirname(paths)
  missing <- !dir.exists(directories)
  if (any(missing)) stop("The directory for `filename.cifreport` does not exist: ",
    paste(unique(directories[missing]), collapse = ", "), call. = FALSE)
  paths <- file.path(normalizePath(directories, winslash = "/", mustWork = TRUE), basename(paths))
  if (any(dir.exists(paths))) stop("`filename.cifreport` points to a directory.", call. = FALSE)
  if (anyDuplicated(if (.Platform$OS.type == "windows") tolower(paths) else paths))
    stop("`filename.cifreport` must contain distinct paths.", call. = FALSE)
  paths
}

.cifreport_check_fields <- function(x, allowed = NULL, arg = "definition") {
  nm <- names(x)
  if (!is.list(x) || (length(x) && (is.null(nm) || anyNA(nm) ||
      any(!nzchar(nm)) || anyDuplicated(nm))))
    stop("`", arg, "` must be a named list with unique names.", call. = FALSE)
  if (!is.null(allowed) && any(!nm %in% allowed))
    stop("Unknown ", arg, " field: ", paste(setdiff(nm, allowed), collapse = ", "), ".", call. = FALSE)
}

.cifreport_plan <- function(subset.condition, exposure, outcomes, time.point, code.exposure.ref = NULL,
                            subgroup = NULL, strata = list(), outcome.type = NULL) {
  if (is.list(subset.condition))
    stop("`subset.condition` must define one population. Use `subgroup` for categorical subsets.",
         call. = FALSE)
  if (!is.null(subgroup) && !.cifreport_string(subgroup))
    stop("`subgroup` must name one categorical column.", call. = FALSE)
  if (is.character(exposure)) exposure <- as.list(exposure)
  if (!is.list(exposure) || !length(exposure)) stop("Supply at least one exposure definition.", call. = FALSE)
  if (is.null(names(exposure))) names(exposure) <- rep("", length(exposure))
  unnamed <- !is.na(names(exposure)) & names(exposure) == ""
  names(exposure)[unnamed] <- vapply(exposure[unnamed], function(group) {
    variable <- if (is.list(group)) group[["variable"]] else group
    if (.cifreport_string(variable)) variable else ""
  }, character(1L))
  .cifreport_check_fields(exposure, arg = "exposure")
  .cifreport_check_fields(outcomes, arg = "outcomes")
  if (!length(outcomes)) stop("Supply at least one outcome.", call. = FALSE)
  outcome.type <- .cifreport_normalize_outcome_type(outcome.type)
  for (field in names(strata)) {
    if (is.list(strata[[field]]))
      .cifreport_check_fields(strata[[field]], names(exposure), field)
  }
  if (!is.null(code.exposure.ref) && (length(code.exposure.ref) != 1L ||
      !is.null(names(code.exposure.ref)) || is.list(code.exposure.ref))) {
    code.exposure.ref <- as.list(code.exposure.ref)
    .cifreport_check_fields(code.exposure.ref, names(exposure), "code.exposure.ref")
  }
  for (id in names(exposure)) {
    group <- exposure[[id]]
    if (.cifreport_string(group)) group <- list(variable = group)
    .cifreport_check_fields(group,
      c("variable", "code.exposure.ref", "label.strata", "order.strata", "level.strata"), "group")
    if (!.cifreport_string(group[["variable"]]))
      stop("Each group needs a `variable` column name.", call. = FALSE)
    reference <- if (is.list(code.exposure.ref)) code.exposure.ref[[id]] else code.exposure.ref
    group$code.exposure.ref <- group$code.exposure.ref %||% reference
    for (field in names(strata)) {
      value <- strata[[field]]
      if (is.list(value)) value <- value[[id]]
      group[[field]] <- group[[field]] %||% value
    }
    exposure[[id]] <- group
  }
  outcomes <- lapply(outcomes, function(outcome) {
    .cifreport_check_fields(outcome,
      c("time", "event", "outcome.type", "code.event1", "code.event2", "code.censoring", "time.point",
        "estimator"),
      "outcome")
    defaults <- list(outcome.type = outcome.type, code.event1 = 1L, code.event2 = 2L,
                     code.censoring = 0L, time.point = time.point)
    defaults[names(outcome)] <- outcome
    defaults
  })
  list(subset.condition = subset.condition, subgroup = subgroup, exposure = exposure, outcomes = outcomes)
}

.cifreport_population <- function(condition, data, envir) {
  if (is.character(condition)) condition <- parse(text = condition)
  if (is.expression(condition) && length(condition) != 1L)
    stop("`subset.condition` must contain one expression.")
  if (inherits(condition, "formula")) envir <- environment(condition) %||% envir
  do.call(createAnalysisDataset, list(formula = ~1, data = data,
    other.variables.analyzed = names(data), subset.condition = condition), quote = TRUE, envir = envir)
}

.cifreport_group <- function(data, variable) {
  if (!variable %in% names(data)) stop("Grouping column is absent: ", variable, ".")
  value <- data[[variable]]
  if ((!is.factor(value) && !is.character(value) && !is.logical(value) &&
       (!is.numeric(value) || is.complex(value))) || !is.null(dim(value)))
    stop("The grouping column must be categorical or an atomic vector of group codes.")
  if (anyNA(value)) stop("The grouping column contains missing values.")
  if (is.numeric(value) && any(!is.finite(value))) stop("Group codes must be finite.")
  factor(value, ordered = FALSE)
}

.cifreport_normalize_outcome_type <- function(x) {
  if (is.null(x) || (.cifreport_string(x) && identical(tolower(trimws(x)), "auto")))
    return(NULL)
  if (!.cifreport_string(x) || !util_normalize_outcome_type(x) %in% c("survival", "competing-risk"))
    stop("`outcome.type` must be \"survival\", \"competing-risk\", or \"auto\".", call. = FALSE)
  util_check_outcome_type(x)
}

.cifreport_outcome <- function(outcome, data, group, population = data) {
  if (!.cifreport_string(outcome$time) || !.cifreport_string(outcome$event) ||
      !all(c(outcome$time, outcome$event) %in% names(data)))
    stop("Outcome time and event must name columns in the data.")
  outcome.type <- .cifreport_normalize_outcome_type(outcome$outcome.type)
  times <- outcome$time.point
  if (!is.numeric(times) || is.complex(times) || !is.null(dim(times)) || !length(times) ||
      any(!is.finite(times)) || any(times < 0))
    stop("`time.point` must supply finite, non-negative evaluation times.")
  time <- data[[outcome$time]]
  status <- data[[outcome$event]]
  if (anyNA(time) || anyNA(status)) stop("The outcome time or event column contains missing values.")
  if (!is.numeric(time) || is.complex(time) || !is.null(dim(time)) ||
      is.complex(status) || !is.null(dim(status)))
    stop("Follow-up must be numeric and event codes must form a vector.")
  if (is.null(outcome.type)) {
    population_status <- population[[outcome$event]]
    observed <- tryCatch(
      util_normalize_time_event(rep.int(0, length(population_status)),
        population_status, allowed = NULL)$event,
      # Keep valid subgroups analyzable if another subgroup has malformed codes.
      error = function(e) util_normalize_time_event(time, status, allowed = NULL)$event
    )
    outcome.type <- util_detect_outcome_type_from_status(observed, auto_message = FALSE)
  }
  choices <- if (outcome.type == "survival") c("KM", "1-KM") else c("AJ", "1-AJ")
  estimator <- outcome$estimator %||% choices[1L]
  if (!.cifreport_string(estimator) || !estimator %in% choices)
    stop("For this outcome, `estimator` must be ", paste(choices, collapse = " or "), ".")
  codes <- list(code.event1 = outcome$code.event1, code.censoring = outcome$code.censoring)
  if (outcome.type == "competing-risk") codes$code.event2 <- outcome$code.event2
  codes <- as.list(do.call(util_validate_event_codes, codes))
  response <- util_normalize_time_event(time, status, unlist(codes, use.names = FALSE))
  if (anyNA(response$event)) stop("Event codes are outside the supported integer range.")
  list(data = data.frame(followup = response$time, status = response$event, group = group),
       arguments = c(list(outcome.type = outcome.type), codes),
       time.variable = outcome$time, event.variable = outcome$event, estimator = estimator,
       times = sort(unique(as.numeric(times))))
}

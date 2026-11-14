#' Create a survival and competing-risk report
#'
#' @param data A data frame with one row per subject.
#' @param outcomes A named list of outcome definitions. Each is a list with
#'   column names `time` and `event`, and optionally `outcome.type`
#'   (`"competing-risk"`, the default, or `"survival"`), `code.event1`
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
#' @param subset.condition One condition defining the analysis population.
#'   Use a logical vector, character condition, one-sided formula, or quoted
#'   expression. Missing condition values are excluded. `NULL` includes all rows.
#' @param subgroup Name of one categorical column. In addition to the analysis
#'   population as a whole, the report includes each observed level of this
#'   column within that population. Missing values contribute only to the
#'   analysis population as a whole. Numeric category codes are supported.
#' @param output.file HTML filename, relative to the current directory or
#'   absolute. Its parent directory must exist. The default is
#'   `yyyymmdd_<data.name>_analysis_report.html` in the working directory.
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
#' events at exactly the requested time. Table 3 and the figure use the
#' selected estimator scale. KM is survival, 1-KM is event probability,
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
#' event code. It requires optional reporting packages and Pandoc. See
#' `vignette("v10_cifreport")` for examples and supported configurations.
#'
#' @return Invisibly, a list with `output.file`, `specifications`,
#'   `tables`, `figures`, `fits`, `status` and `versions`. Tables,
#'   figures and fits are nested by population and grouping definition.
#'   Table 1 is a gtsummary object. Tables 2--5 are numeric result data
#'   frames with outcome and group identifiers. Group codes and display
#'   labels are available in Tables 2--4. Figures and fits are
#'   indexed by outcome, with risk-ratio fits indexed by time.
#' @examples
#' if (requireNamespace("gtsummary", quietly = TRUE) &&
#'     utils::packageVersion("gtsummary") >= "2.6.1" &&
#'     requireNamespace("rmarkdown", quietly = TRUE) &&
#'     requireNamespace("gt", quietly = TRUE) && rmarkdown::pandoc_available()) {
#'   data(diabetes.complications)
#'   report <- cifreport(
#'     diabetes.complications,
#'     outcomes = list(Retinopathy = list(time = "t", event = "epsilon")),
#'     exposure = "fruitq1", characteristics = c("age", "sex"),
#'     time.point = 8, output.file = tempfile(fileext = ".html")
#'   )
#' }
#' @seealso [cifcurve()], [cifplot()], [polyreg()], [ciftest()]
#' @md
#' @export
cifreport <- function(data, outcomes, exposure, characteristics = character(),
                      time.point = NULL, subset.condition = NULL,
                      output.file = NULL, title = NULL,
                      overall = FALSE, conf.int = 0.95, table1.args = list(),
                      code.exposure.ref = NULL, data.name = NULL, subgroup = NULL,
                      label.strata = NULL, order.strata = NULL, level.strata = NULL) {
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
  if (!.cifreport_string(title)) stop("`title` must be a nonempty string.", call. = FALSE)
  if (is.null(output.file)) output.file <- .cifreport_filename(data.name, date)
  if (!.cifreport_string(output.file) || !grepl("\\.html?$", output.file, ignore.case = TRUE))
    stop("`output.file` must be an HTML filename.", call. = FALSE)
  output.file <- path.expand(output.file)
  output.file <- file.path(normalizePath(dirname(output.file), winslash = "/", mustWork = TRUE),
                           basename(output.file))
  if (dir.exists(output.file)) stop("`output.file` points to a directory.", call. = FALSE)
  .cifreport_dependencies()

  strata <- list(label.strata = label.strata, order.strata = order.strata, level.strata = level.strata)
  plan <- .cifreport_plan(subset.condition, exposure, outcomes, time.point,
                          code.exposure.ref, subgroup, strata)
  result <- .cifreport_run(as.data.frame(data), plan, characteristics, overall, conf.int, table1.args,
                           subset.env = parent.frame())
  result$output.file <- output.file
  result$specifications <- c(plan, list(data.name = data.name, characteristics = characteristics, overall = overall,
                                       conf.int = conf.int, table1.args = table1.args, title = title))
  packages <- c("cifmodeling", "survival", "ggplot2", "ggsurvfit", "gtsummary",
                "cardx", "broom", "broom.helpers", "gt", "knitr", "rmarkdown")
  result$versions <- c(R = as.character(getRversion()),
    vapply(stats::setNames(packages, packages), function(package) {
      if (requireNamespace(package, quietly = TRUE))
        as.character(getNamespaceVersion(package)) else NA_character_
    }, character(1L)), pandoc = as.character(rmarkdown::pandoc_version()))
  .cifreport_render(result)
  invisible(result)
}

.cifreport_dependencies <- function() {
  for (package in c("gtsummary", "gt", "knitr", "rmarkdown")) {
    if (!requireNamespace(package, quietly = TRUE))
      stop("Reporting requires the ", package, " package.", call. = FALSE)
  }
  if (package_version(getNamespaceVersion("gtsummary")) < "2.6.1")
    stop("Reporting requires gtsummary >= 2.6.1.", call. = FALSE)
  if (!rmarkdown::pandoc_available())
    stop("Reporting requires Pandoc. See rmarkdown::find_pandoc().", call. = FALSE)
}

.cifreport_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.cifreport_filename <- function(data.name, date = format(Sys.Date(), "%Y%m%d")) {
  safe_name <- gsub('[<>:"/\\\\|?*[:cntrl:]]', "_", data.name)
  paste0(date, "_", safe_name, "_analysis_report.html")
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
                            subgroup = NULL, strata = list()) {
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
    defaults <- list(outcome.type = "competing-risk", code.event1 = 1L, code.event2 = 2L,
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

.cifreport_outcome <- function(outcome, data, group) {
  if (!.cifreport_string(outcome$time) || !.cifreport_string(outcome$event) ||
      !all(c(outcome$time, outcome$event) %in% names(data)))
    stop("Outcome time and event must name columns in the data.")
  if (!.cifreport_string(outcome$outcome.type) ||
      !outcome$outcome.type %in% c("survival", "competing-risk"))
    stop("`outcome.type` must be \"survival\" or \"competing-risk\".")
  choices <- if (outcome$outcome.type == "survival") c("KM", "1-KM") else c("AJ", "1-AJ")
  estimator <- outcome$estimator %||% choices[1L]
  if (!.cifreport_string(estimator) || !estimator %in% choices)
    stop("For this outcome, `estimator` must be ", paste(choices, collapse = " or "), ".")
  times <- outcome$time.point
  if (!is.numeric(times) || is.complex(times) || !is.null(dim(times)) || !length(times) ||
      any(!is.finite(times)) || any(times < 0))
    stop("`time.point` must supply finite, non-negative evaluation times.")
  codes <- list(code.event1 = outcome$code.event1, code.censoring = outcome$code.censoring)
  if (outcome$outcome.type == "competing-risk") codes$code.event2 <- outcome$code.event2
  codes <- as.list(do.call(util_validate_event_codes, codes))
  time <- data[[outcome$time]]
  status <- data[[outcome$event]]
  if (anyNA(time) || anyNA(status)) stop("The outcome time or event column contains missing values.")
  if (!is.numeric(time) || is.complex(time) || !is.null(dim(time)) ||
      is.complex(status) || !is.null(dim(status)))
    stop("Follow-up must be numeric and event codes must form a vector.")
  response <- util_normalize_time_event(time, status, unlist(codes, use.names = FALSE))
  if (anyNA(response$event)) stop("Event codes are outside the supported integer range.")
  list(data = data.frame(followup = response$time, status = response$event, group = group),
       arguments = c(list(outcome.type = outcome$outcome.type), codes),
       time.variable = outcome$time, event.variable = outcome$event, estimator = estimator,
       times = sort(unique(as.numeric(times))))
}

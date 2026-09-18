#!/usr/bin/env Rscript

# Optional developer audit. cmprsk is deliberately not a package Suggests.
if (!requireNamespace("cmprsk", quietly = TRUE)) {
  stop("Install cmprsk to run this optional validation script.", call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools to load the package source tree.", call. = FALSE)
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(".github/scripts/validate-gray-cmprsk.R", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
devtools::load_all(root, quiet = TRUE)

fixture_environment <- new.env(parent = baseenv())
sys.source(
  file.path(root, "tests", "testthat", "fixtures", "gray_cmprsk_fixtures.R"),
  envir = fixture_environment
)

results <- lapply(fixture_environment$gray_cmprsk_fixtures, function(fixture) {
  df <- data.frame(
    time = fixture$time,
    status = fixture$status,
    group = factor(fixture$group, levels = unique(fixture$group))
  )
  package_result <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = FALSE,
    rho = fixture$rho
  )
  reference_result <- cmprsk::cuminc(
    ftime = df$time,
    fstatus = df$status,
    group = df$group,
    cencode = 0,
    rho = fixture$rho
  )
  current <- unname(reference_result$Tests[1L, c("stat", "pv", "df")])
  names(current) <- c("statistic", "p.value", "df")
  package_values <- c(
    statistic = unname(package_result$statistic),
    p.value = package_result$p.value,
    df = unname(package_result$parameter)
  )
  frozen <- c(
    statistic = fixture$statistic,
    p.value = fixture$p.value,
    df = fixture$df
  )
  data.frame(
    id = fixture$id,
    source = c("cifmodeling", "cmprsk-current", "cmprsk-frozen"),
    rbind(package_values, current, frozen),
    row.names = NULL,
    check.names = FALSE
  )
})

results <- do.call(rbind, results)
print(results, digits = 16, row.names = FALSE)

by_fixture <- split(results, results$id)
ok <- vapply(by_fixture, function(x) {
  max(abs(x$statistic - x$statistic[1L])) < 2e-12 &&
    max(abs(x$p.value - x$p.value[1L])) < 2e-12 &&
    length(unique(x$df)) == 1L
}, logical(1))
if (!all(ok)) {
  stop("At least one Gray validation fixture differs from cmprsk.", call. = FALSE)
}

message(
  "All fixtures agree with cmprsk ",
  as.character(utils::packageVersion("cmprsk")),
  "."
)

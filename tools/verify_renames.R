#!/usr/bin/env Rscript
# Verify no old symbol remains in code regions.

suppressPackageStartupMessages({
  library(jsonlite); library(fs)
})

map <- jsonlite::fromJSON("tools/rename_map.json")
old_names <- names(map)

scan_targets <- c("R", "tests", "vignettes", "inst", "NAMESPACE", "DESCRIPTION", "README.Rmd", "pkgdown")
paths <- dir_ls(".", recurse = TRUE, type = "file")
paths <- paths[file_exists(paths)]
paths <- paths[grepl(paste0("^\\./(", paste(scan_targets, collapse="|"), ")"), paths)]

hit <- list()
for (nm in old_names) {
  pat <- sprintf("(?<![A-Za-z0-9_.])%s(?![A-Za-z0-9_.])", nm)
  found <- lapply(paths, function(p) {
    lines <- readLines(p, warn = FALSE)
    which(grepl(pat, lines, perl = TRUE))
  })
  indices <- which(vapply(found, length, integer(1)) > 0)
  if (length(indices)) {
    hit[[nm]] <- Map(function(p, idx) list(file = p, lines = idx), paths[indices], found[indices])
  }
}

if (length(hit)) {
  cat("Old names still found:\n")
  print(hit)
  quit(status = 1)
} else {
  cat("All clear: no old names found in target paths.\n")
}

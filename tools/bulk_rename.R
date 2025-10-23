#!/usr/bin/env Rscript
# Safe bulk rename of internal R functions based on a JSON map.
# Usage: Rscript tools/bulk_rename.R

suppressPackageStartupMessages({
  library(jsonlite); library(fs); library(stringi)
})

root <- path_abs(".")
map  <- fromJSON("tools/rename_map.json")

# ---------- helpers ----------
word_boundary <- function(x) sprintf("(?<![A-Za-z0-9_.])%s(?![A-Za-z0-9_.])", stringi::stri_escape_regex(x))

read_txt  <- function(p) paste0(readLines(p, warn = FALSE), collapse = "\n")
write_txt <- function(p, s) writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], p, useBytes = TRUE)

should_touch <- function(p) {
  ext <- tolower(path_ext(p))
  bn  <- path_file(p)
  if (grepl("^\\.git|^\\.Rproj\\.user", p)) return(FALSE)
  if (bn %in% c("DESCRIPTION")) return(TRUE)
  ext %in% c("r","rmd","rd","qmd","md","yaml","yml","rproj","txt","rds","json")
}

# ---------- files ----------
files <- dir_ls(path = root, recurse = TRUE, type = "file")
files <- files[vapply(files, should_touch, logical(1))]

# Exclude some logs and pkg cache
files <- files[!grepl("(^|/)\\.(git|github|venv|renv|quarto|cache)(/|$)", files)]
files <- files[!grepl("(^|/)(README|NEWS|cran-comments).*\\.md$", files)] # 説明文は壊さない（必要なら後で）

# ---------- replacement rules ----------
replace_code_regions <- function(txt, replacements) {
  # R/ *.R, NAMESPACE, *.Rmd code chunks、Rdの\\code{} などを優先
  # 1) 全体を対象にするが、言語的にワード境界で限定
  for (old in names(replacements)) {
    new <- replacements[[old]]
    # bare tokens & backticked
    pat1 <- word_boundary(old)
    pat2 <- sprintf("`%s`", stringi::stri_escape_regex(old))
    txt  <- gsub(pat1, new, txt, perl = TRUE)
    txt  <- gsub(pat2, sprintf("`%s`", new), txt, perl = TRUE)
  }
  txt
}

# ---------- run ----------
cat(sprintf("Renaming %d symbols across %d files...\n", length(map), length(files)))
for (f in files) {
  old <- read_txt(f)
  new <- replace_code_regions(old, map)
  if (!identical(old, new)) {
    write_txt(f, new)
    cat("updated:", f, "\n")
  }
}

cat("Done. Now run: devtools::document(); devtools::check()\n")

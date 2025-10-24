suppressPackageStartupMessages({
  library(jsonlite); library(fs); library(stringi)
})

escape_regex <- function(x) {
  # NꕶGXP[v・
  gsub("([][{}()^$.*+?|\\-\\\\])", "\\\\\\1", x)
}

word_boundary <- function(x) sprintf(
  "(?<![A-Za-z0-9_.])%s(?![A-Za-z0-9_.])",
  escape_regex(x)
)

root <- path_abs(".")
map  <- fromJSON("tools/rename_map.json")

read_txt  <- function(p) paste0(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
write_txt <- function(p, s) writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], p, useBytes = TRUE)

should_touch <- function(p) {
  ext <- tolower(path_ext(p))
  bn  <- path_file(p)
  if (grepl("^\\.git|^\\.Rproj\\.user", p)) return(FALSE)
  if (bn %in% c("DESCRIPTION")) return(TRUE)
  ext %in% c("r","rmd","rd","qmd","md","yaml","yml","rproj","txt","rds","json")
}

files <- dir_ls(path = root, recurse = TRUE, type = "file")
files <- files[vapply(files, should_touch, logical(1))]

files <- files[!grepl("(^|/)\\.(git|github|venv|renv|quarto|cache)(/|$)", files)]
files <- files[!grepl("(^|/)(README|NEWS|cran-comments).*\\.md$", files)] # ־͉ȂiKvȂ迪・ﾁ・j

replace_code_regions <- function(txt, replacements) {
  for (old in names(replacements)) {
    new <- replacements[[old]]
    # bare tokens & backticked
    pat1 <- word_boundary(old)
    pat2 <- sprintf("`%s`", escape_regex(old))
    txt  <- gsub(pat1, new, txt, perl = TRUE)
    txt  <- gsub(pat2, sprintf("`%s`", new), txt, perl = TRUE)
  }
  txt
}

# ---------- run ----------
args <- commandArgs(trailingOnly = TRUE)

map_path <- if (length(args) >= 1) args[1] else "tools/rename_map.json"
root_dir <- if (length(args) >= 2) args[2] else "."

if (!file.exists(map_path)) {
  stop(sprintf("rename map JSON not found: %s", map_path))
}

map <- jsonlite::fromJSON(map_path, simplifyVector = TRUE)
if (!is.list(map) || !length(map)) stop("map is empty or malformed")

all_files <- fs::dir_ls(root_dir, recurse = TRUE, type = "file", fail = FALSE)
files <- all_files[
  grepl("\\.(R|r)$", all_files) |
  grepl("\\.(R|r)(md|MD)$", all_files) |
  grepl("\\.(q|Q)md$", all_files) |
  grepl("\\.(Rd|rd)$", all_files)
]

skip_dirs <- c("^\\.git/", "^\\.Rproj\\.user/", "^renv/", "^packrat/")
if (length(files)) {
  keep <- !Reduce(`|`, lapply(skip_dirs, function(p) grepl(p, files)))
  files <- files[keep]
}

if (!length(files)) stop("no target files found under: ", root_dir)

cat(sprintf("Renaming %d symbols across %d files...\n", length(map), length(files)))

n_updated <- 0L
for (f in files) {
  if (!file.exists(f)) {
    warning(sprintf("skip (not found): %s", f), call. = FALSE); next
  }
  info <- tryCatch(fs::file_info(f), error = function(e) NULL)
  if (!is.null(info) && is.finite(info$size) && info$size > 50e6) {
    warning(sprintf("skip (too large >50MB): %s", f), call. = FALSE); next
  }

  old <- tryCatch(read_txt(f), error = function(e) {
    warning(sprintf("skip (cannot open): %s ; %s", f, e$message), call. = FALSE)
    return(NULL)
  })
  if (is.null(old)) next

  new <- replace_code_regions(old, map)

  if (!identical(old, new)) {
    bak <- paste0(f, "~")
    tryCatch(fs::file_copy(f, bak, overwrite = TRUE), error = function(e) {})
    tryCatch({
      write_txt(f, new)
      cat("updated:", f, "\n")
      n_updated <- n_updated + 1L
    }, error = function(e) {
      warning(sprintf("write failed: %s ; %s", f, e$message), call. = FALSE)
    })
  }
}
cat(sprintf("Done. Updated %d/%d files.\n", n_updated, length(files)))
cat("Now run: devtools::document(); devtools::check()\n")
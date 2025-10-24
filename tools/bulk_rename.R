#!/usr/bin/env Rscript
# Safe bulk rename of internal R functions based on a JSON map.
# Usage: Rscript tools/bulk_rename.R

# 置換：stringi 依存をやめる
# suppressPackageStartupMessages({ library(jsonlite); library(fs) })
# library(stringi) は削除

suppressPackageStartupMessages({
  library(jsonlite); library(fs)
})

# 追加：正規表現エスケープ（最小限で十分）
escape_regex <- function(x) {
  # [, ], {, }, (, ), ^, $, ., *, +, ?, |, -, \ をエスケープ
  gsub("([][{}()^$.*+?|\\-\\\\])", "\\\\\\1", x)
}

word_boundary <- function(x) sprintf(
  "(?<![A-Za-z0-9_.])%s(?![A-Za-z0-9_.])",
  escape_regex(x)
)

read_txt  <- function(p) paste0(readLines(p, warn = FALSE), collapse = "\n")
write_txt <- function(p, s) writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], p, useBytes = TRUE)

# 追加: ファイル種別判定
is_r_file    <- function(f) grepl("\\.[rR]$", f)
is_rmd_file  <- function(f) grepl("\\.(R|r)(md|MD)$", f) || grepl("\\.(q|Q)md$", f)
is_rd_file   <- function(f) grepl("\\.(Rd|rd)$", f)

# 置換ロジックをファイルごとに適用
# --- 追加: Rmd/qmd のコードチャンクを安全に置換するユーティリティ ---
process_rmd_chunks <- function(s, fun) {
  # ```{r ...}\n ... \n``` をマッチ
  pat <- "(?s)```\\{r[^}]*\\}.*?\\n```"
  m <- gregexpr(pat, s, perl = TRUE)
  chunks <- regmatches(s, m)[[1]]
  if (length(chunks) == 0) return(s)
  # 各チャンクに処理を適用して差し戻す
  newchunks <- vapply(chunks, fun, character(1))
  regmatches(s, m) <- list(newchunks)
  s
}

# --- 既存: エスケープ & ワード境界 ---
escape_regex <- function(x) gsub("([][{}()^$.*+?|\\-\\\\])", "\\\\\\1", x)
word_boundary <- function(x) sprintf("(?<![A-Za-z0-9_.])%s(?![A-Za-z0-9_.])", escape_regex(x))

# --- 既存: ファイル種別判定 ---
is_r_file   <- function(f) grepl("\\.[rR]$", f)
is_rmd_file <- function(f) grepl("\\.(R|r)(md|MD)$", f) || grepl("\\.(q|Q)md$", f)
is_rd_file  <- function(f) grepl("\\.(Rd|rd)$", f)

# --- 置換ヘルパ ---
replace_token_and_backtick <- function(s, replacements) {
  for (old in names(replacements)) {
    new <- replacements[[old]]
    s <- gsub(word_boundary(old), new, s, perl = TRUE)                                   # 素のトークン
    s <- gsub(sprintf("`%s`", escape_regex(old)), sprintf("`%s`", new), s, perl = TRUE)  # バッククォート
  }
  s
}

replace_quoted_in_R <- function(s, replacements) {
  for (old in names(replacements)) {
    new <- replacements[[old]]
    # "OLD" / 'OLD' → "NEW" / 'NEW'
    s <- gsub(sprintf('"%s"', escape_regex(old)), sprintf('"%s"', new), s, perl = TRUE)
    s <- gsub(sprintf("'\\Q%s\\E'", old),         sprintf("'%s'", new), s, perl = TRUE)
  }
  s
}

# --- NEW: 関数定義 (foo <- function(...) の左辺) を確実に置換 ---
replace_function_defs <- function(s, replacements) {
  # ^ or \n の直後で、空白→旧名→空白→<-→空白→function(
    # 例:   reg_normalize_covariate <- function( ...
  for (old in names(replacements)) {
    new <- replacements[[old]]
    pat <- sprintf("(?m)(^|\\n)([ \\t]*)%s([ \\t]*)<-(\\s*)function\\s*\\(",
                   escape_regex(old))
    repl <- sprintf("\\1\\2%s\\3<-\\4function(", new)
    s <- gsub(pat, repl, s, perl = TRUE)
  }
  s
}

# --- 差し替え: ファイルに応じて置換する本体 ---
replace_code_regions <- function(txt, replacements, filename = NULL) {
  # 1) まず関数定義を確実に置換
  out <- replace_function_defs(txt, replacements)

  # 2) それから従来のトークン置換
  out <- replace_token_and_backtick(out, replacements)

  # 3) ファイル種別によって文字列も置換
  if (!is.null(filename) && is_r_file(filename)) {
    out <- replace_quoted_in_R(out, replacements)
  } else if (!is.null(filename) && is_rmd_file(filename)) {
    out <- process_rmd_chunks(
      out,
      function(chunk) replace_quoted_in_R(chunk, replacements)
    )
  }
  out
}


# ---------- run ----------
cat(sprintf("Renaming %d symbols across %d files...\n", length(map), length(files)))
for (f in files) {
  old <- read_txt(f)
  new <- replace_code_regions(old, map, filename = f)
  if (!identical(old, new)) {
    write_txt(f, new)
    cat("updated:", f, "\n")
  }
}
cat("Done. Now run: devtools::document(); devtools::check()\n")

## Test environments
- Local:
  - Windows 11, R 4.4.2: R CMD check --as-cran
  - macOS 14 (Apple Silicon), R 4.4.2: R CMD check --as-cran
  - Ubuntu 22.04, R 4.4.2: R CMD check --as-cran
- GitHub Actions:
  - Ubuntu-latest, R release and devel — no ERROR/WARNING
  - Windows-latest, R release — no ERROR/WARNING
  - macOS-latest, R release — no ERROR/WARNING
- win-builder:
  - R release — no ERROR/WARNING
  - R devel   — no ERROR/WARNING

## R CMD check results
0 errors | 0 warnings | 2 notes

* NOTE: New submission.
  This is the first CRAN submission of the package.

* NOTE: Help files contain \Sexpr{} expressions but no prebuilt PDF manual.
  The package does not ship a prebuilt manual; CRAN builds the manual.
  The \Sexpr{} calls are self-contained and do not depend on external resources.

* NOTE: checking for future file timestamps ... unable to verify current time.
  This is environment-specific and not indicative of a package issue.

## Downstream dependencies
None (new package).

## NOTE
Possibly misspelled words in DESCRIPTION found in win-builder are proper names or package names: Aalen, Johansen, Kaplan, ggsurvfit, modelsummary, subdistribution, survfit.

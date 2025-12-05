# cifmodeling 0.9.5

* First CRAN release.
* Documentation and tests were polished for CRAN submission.

# cifmodeling 0.9.4

* Sixth CRAN submission.

Invalid DESCRIPTION file, software names were single quoted. Rebuilt and resubmitted.

# cifmodeling 0.9.3

* Fifth CRAN submission.

Addressed CRAN feedback regarding file-system and environment management.
Removed development-time scripts under tools/ that performed bulk renaming
and wrote to the package directory or the current working directory.
These files are no longer included in the source tarball.
Eliminated the remaining use of .GlobalEnv in internal helpers.
panel_as_formula_global() now evaluates formulas in the calling environment
(parent.frame()) instead of .GlobalEnv, fully complying with the
“no modification of the .GlobalEnv” CRAN policy.
No user-visible changes; only internal cleanup for CRAN compliance.
>>>>>>> d7a80ad7267e1786d6a8dbb7001f1915d56255a7

# cifmodeling 0.9.2

* Fourth CRAN submission.
* Invalid file URLs  
  (articles/gallery.html, articles/polyreg.html, articles/formulas.html)
  were replaced by URLs in the web and verified using urlchecker::url_check(). 
  Rebuilt and resubmitted.

# cifmodeling 0.9.1

* Third CRAN submission.
* The DOI reported by CRAN (`10.21105/joss.00510`) does not appear to be present
  in the package source. Verified using recursive search and a freshly created
  tarball. Rebuilt and resubmitted.

# cifmodeling 0.9.0

* Second CRAN submission.

# cifmodeling 0.8.2

* Bug fixes
* README.Rmd, Vignettes and site were updated. 
* Cumulative hazard and log-log plots were implemented

# cifmodeling 0.8.1

* README.Rmd, Vignettes and site were updated. 
* polyreg S3 methods were implemented

# cifmodeling 0.8.0

* Breaking change: `cifplot()` and `cifpanel()` now return structured objects
  with explicit `plot`/`list.plot` and `patchwork` elements, metadata, and
  printing handled via new S3 methods. Automatic printing is limited to
  interactive sessions, and panel outputs store engine metadata in
  `print.info$engine`.

* Standardized outcome.type canonical values to lowercase; legacy uppercase and short aliases remain accepted. Internals now use tolower().

# cifmodeling 0.7.0

* Bug fixes
* README.Rmd, Vignettes and site were updated. 

# cifmodeling 0.6.0

* Initial CRAN submission.

# cifmodeling 0.5.0

* Bug fixes

# cifmodeling 0.4.0

* Bug fixes

# cifmodeling 0.3.0

* `cifpanel()` and panel mode of * `cifplot()` were developed.

# cifmodeling 0.2.0

* Bug fixes

# cifmodeling 0.1.0

* `cifcurve()`, `cifplot()` and `polyreg()` were developed.

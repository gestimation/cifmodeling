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

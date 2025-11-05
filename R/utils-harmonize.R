#' @keywords internal
harmonize_engine_output <- function(out) {
  n <- length(out$time %||% numeric())
  out$lower               <- out$lower               %||% out$low
  out$upper               <- out$upper               %||% out$high
  out$`std.err.cif`        <- out$`std.err.cif`        %||% rep(NA_real_, n)
  out$`influence.function` <- out$`influence.function` %||% list()
  out$strata               <- out$strata               %||% integer()
  out$`strata.levels`      <- out$`strata.levels`      %||% integer()
  out$type                 <- out$type                 %||% "kaplan-meier"
  out$method               <- out$method               %||% "Kaplan–Meier"
  out$`conf.type`          <- out$`conf.type`          %||% "arcsine-square root"
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

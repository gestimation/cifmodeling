#' Create a survival or competing-risks response
#'
#' A lightweight response constructor used in \code{cifcurve()} and \code{polyreg()}
#' to pass survival and competing-risks data via a model formula.
#'
#' @param time Numeric vector of follow-up times (non-negative).
#' @param event Integer (0=censor, 1,2,...) or a character/factor vector whose levels
#'   are numeric codes "0","1","2",... for competing events.
#'
#' @return An object of class \code{"Event"} (a 2-column matrix) with columns \code{time}, \code{event}.
#' @examples
#' ## event: 0=censor, 1=primary, 2=competing
#' # polyreg(nuisance.model = Event(t, epsilon) ~ 1, data = df, outcome.type="COMPETING-RISK")
#'
#' @export
Event <- function(time, event) {
  te <- normalize_time_event(time, event)
  ss <- cbind(time = te$time, event = te$event)
  dimnames(ss) <- list(NULL, c("time","event"))
  attr(ss, "type") <- "right"
  class(ss) <- c("Event", class(ss))
  ss
}

normalize_time_event <- function(time, event, allowed = NULL) {
  if (missing(time))  .err("req", arg = "time")
  if (missing(event)) .err("req", arg = "event")
  if (!is.numeric(time)) .err("numeric", arg = "time")
  if (any(time < 0, na.rm = TRUE)) .err("nonneg", arg = "time")

  if (is.numeric(event)) {
    if (any(event < 0, na.rm = TRUE) || any(event != floor(event), na.rm = TRUE)) {
      .err("ev_codes", allowed = "{0,1,2,...}",
           found = paste(unique(event[!is.na(event)]), collapse = ", "))
    }
    status <- suppressWarnings(as.integer(event))
  } else if (is.logical(event)) {
    status <- ifelse(is.na(event), NA_integer_, as.integer(event))
  } else if (is.factor(event) || is.character(event)) {
    ev_chr <- as.character(event)
    ok <- !is.na(ev_chr)
    if (!all(grepl("^[0-9]+$", ev_chr[ok]))) {
      .err("ev_codes", allowed = "'0','1','2',...",
           found = paste(unique(ev_chr[ok]), collapse = ", "))
    }
    status <- suppressWarnings(as.integer(ev_chr))
  } else {
    .err("ev_type")
  }

  if (length(status) != length(time)) {
    .err("len_mismatch", x = "time", y = "event",
         nx = length(time), ny = length(status))
  }

  if (!is.null(allowed)) {
    ok <- is.na(status) | status %in% allowed
    if (!all(ok)) {
      .err("ev_codes",
           allowed = paste0("{", paste(allowed, collapse = ","), "}"),
           found   = paste(sort(unique(status[!ok])), collapse = ", "))
    }
  }
  list(time = as.numeric(time), event = status)
}

untangle.specials <- function(tt, special, order = 1) {
  spc <- attr(tt, "specials")[[special]]
  if (length(spc) == 0)
    return(list(vars = character(0), terms = numeric(0)))
  facs <- attr(tt, "factors")
  fname <- dimnames(facs)
  ff <- apply(facs[spc, , drop = FALSE], 2, sum)
  list(vars = (fname[[1]])[spc], tvar = spc - attr(tt, "response"),
       terms = seq(ff)[ff & match(attr(tt, "order"), order, nomatch = 0)])
}

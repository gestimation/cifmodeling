.cifcurve_interval_fit <- function(input, call) {
  event <- .cifcurve_event_codes(input)
  fit <- calculateCIF_IC(input$t, input$right, event, input$w, input$strata)
  codes <- names(fit$support)
  labels <- if (is.null(input$strata_name)) "All" else paste0(input$strata_name, "=", codes)
  strata <- if (!is.null(input$strata_name)) stats::setNames(unname(fit$strata), labels)
  positive <- input$w > 0
  left <- split(input$t[positive], input$strata[positive], drop = TRUE)
  censored <- positive & event == 0L
  censor <- split(input$t[censored], input$strata[censored])
  curves <- Map(function(support, time, left, censor) {
    list(support = support, time = time, end = max(time), left = sort(left),
         censor = sort(unique(censor)))
  }, fit$support, split(fit$time, rep(seq_along(codes), fit$strata)), left[codes], censor[codes])
  names(curves) <- labels
  grid <- lapply(curves, function(curve) .cifcurve_interval_at(curve, curve$time))
  estimate <- unlist(lapply(grid, `[[`, "cif"), use.names = FALSE)
  display <- input$display
  keep <- match(codes, levels(input$strata))
  for (field in c("levels", "codes")) {
    if (length(display[[field]]) == nlevels(input$strata)) {
      display[[field]] <- stats::setNames(display[[field]][keep], labels)
    }
  }
  rownames(fit$diagnostics) <- labels
  unavailable <- rep(NA_real_, length(fit$time))
  out <- list(
    time = fit$time, surv = 1 - estimate, n = unname(fit$n), strata = strata,
    n.risk = unavailable, n.event = unavailable, n.censor = unavailable,
    n.risk.confirmed = unlist(lapply(grid, `[[`, "n.risk.confirmed"), use.names = FALSE),
    estimate.status = unlist(lapply(grid, `[[`, "estimate.status"), use.names = FALSE),
    aj = estimate, type = "interval", method = "Interval-censored NPMLE (uniform)",
    conf.type = "none", call = call, diagnostics = fit$diagnostics,
    reporting = list(curves = curves, outcome.type = input$outcome.type,
                     event.code = input$prepared$code.event1, display = display)
  )
  class(out) <- c("cifcurve_interval", "cifcurve", "survfit")
  out
}

.cifcurve_interval_at <- function(curve, times) {
  support <- curve$support
  target <- support[which(support$cause == 1L & support$mass > 0), , drop = FALSE]
  fraction <- outer(times, target$right, `>=`) * 1
  finite.interval <- target$left < target$right
  status <- rep("estimated", length(times))
  if (any(finite.interval)) {
    left <- target$left[finite.interval]
    width <- target$right[finite.interval] - left
    progress <- sweep(outer(times, left, `-`), 2L, width, `/`)
    inside <- outer(times, left, `>`) & outer(times, target$right[finite.interval], `<`)
    status[rowSums(inside) > 0L] <- "interpolated"
    fraction[, finite.interval] <- pmax(pmin(progress, 1), 0)
  }
  cif <- pmin(drop(fraction %*% target$mass), 1)
  tail.start <- support$left[is.na(support$cause) & support$mass > 0]
  unresolved <- times > min(c(Inf, tail.start))
  cif[unresolved] <- NA_real_
  status[unresolved] <- "unresolved_tail"
  confirmed <- if (is.null(curve$left)) rep(NA_integer_, length(times)) else
    length(curve$left) - findInterval(times, curve$left, left.open = TRUE)
  list(cif = cif, estimate.status = status, n.risk.confirmed = confirmed)
}

.cifcurve_interval_summary <- function(object, times, censored = FALSE, scale = 1,
                                       extend = FALSE, rmean = "none",
                                       data.frame = FALSE, dosum, ...) {
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("scale must be a positive finite number.", call. = FALSE)
  }
  if (!is.null(rmean) && !identical(rmean, "none") && !identical(rmean, FALSE)) {
    stop("Restricted mean summaries are unavailable for Interval().", call. = FALSE)
  }
  use.grid <- missing(times)
  if (!use.grid) {
    if (!is.numeric(times) || !length(times) || any(!is.finite(times) | times < 0)) {
      stop("times must be finite and nonnegative.", call. = FALSE)
    }
    times <- sort(times)
  }
  curves <- object$reporting$curves
  rows <- Map(function(curve, label) {
    t <- if (use.grid) curve$time else times[extend | times <= curve$end]
    at <- .cifcurve_interval_at(curve, t)
    unavailable <- rep(NA_real_, length(t))
    data.frame(time = t / scale, n.risk = unavailable, n.event = unavailable,
      n.risk.confirmed = at$n.risk.confirmed,
      n.censor = unavailable, surv = 1 - at$cif, cif = at$cif,
      std.err = unavailable, lower = unavailable, upper = unavailable,
      estimate.status = at$estimate.status,
      strata = rep(label, length(t)))
  }, curves, names(curves))
  out <- do.call(rbind, unname(rows))
  out$strata <- if (length(object$strata)) factor(out$strata, levels = names(curves))
  if (data.frame) return(out)
  structure(c(as.list(out), list(n = object$n, type = object$type,
    call = object$call, method = object$method)),
    class = c("summary.cifcurve_interval", "summary.survfit"))
}

#' @export
print.cifcurve_interval <- function(x, ...) {
  cat(x$method, "\n", sep = "")
  print(x$call)
  print(stats::setNames(x$n, names(x$reporting$curves)))
  invisible(x)
}

#' @export
print.summary.cifcurve_interval <- function(x, ...) {
  cols <- c("strata", "time", "cif", "surv", "n.risk.confirmed", "estimate.status")
  print(as.data.frame(unclass(x)[intersect(cols, names(x))]), row.names = FALSE, ...)
  invisible(x)
}

#' @export
plot.cifcurve_interval <- function(x, ...) {
  .cifcurve_interval_plot_error()
}

.cifcurve_interval_plot_error <- function() {
  stop("Plotting Interval() data requires passing the fitted cifcurve object to cifplot().",
       call. = FALSE)
}

.cifcurve_check_plot_interval <- function(x) {
  if (inherits(x, "cifcurve_interval")) .cifcurve_interval_plot_error()
  if (!inherits(x, "formula") || length(x) != 3L || !is.call(x[[2L]])) {
    return(invisible(NULL))
  }
  constructor <- x[[2L]][[1L]]
  if (is.call(constructor) &&
      (identical(constructor[[1L]], quote(`::`)) ||
       identical(constructor[[1L]], quote(`:::`)))) {
    constructor <- constructor[[3L]]
  }
  if (identical(constructor, quote(Interval))) .cifcurve_interval_plot_error()
  invisible(NULL)
}

.cifplot_interval_curve <- function(fit, type.y, linewidth, linetype, visual.info) {
  if (!type.y %in% c("risk", "surv")) {
    stop("Interval() plots support type.y = 'risk' or 'surv'.", call. = FALSE)
  }
  curves <- fit$reporting$curves
  rows <- Map(function(curve, label) {
    time <- sort(unique(c(0, curve$time)))
    cif <- .cifcurve_interval_at(curve, time)$cif
    path <- data.frame(time = time, cif = cif)
    support <- curve$support
    exact <- which(support$cause == 1L & support$left == support$right & support$mass > 0)
    jump <- match(support$right[exact], time)
    before <- data.frame(time = time[jump], cif = pmax(0, cif[jump] - support$mass[exact]))
    path <- rbind(before, path)
    path <- path[order(path$time), , drop = FALSE]
    path$estimate <- if (type.y == "risk") path$cif else 1 - path$cif
    path$strata <- label
    path[!is.na(path$estimate), , drop = FALSE]
  }, curves, names(curves))
  path <- do.call(rbind, unname(rows))
  path$strata <- factor(path$strata, levels = names(curves))
  .data <- NULL
  mapping <- if (!length(fit$strata)) NULL else if (isTRUE(linetype)) {
    ggplot2::aes(colour = .data$strata, linetype = .data$strata)
  } else {
    ggplot2::aes(colour = .data$strata)
  }
  p <- ggplot2::ggplot(path, ggplot2::aes(x = .data$time, y = .data$estimate)) +
    ggplot2::geom_path(mapping = mapping, linewidth = linewidth) +
    ggsurvfit::theme_ggsurvfit_default() +
    ggplot2::labs(colour = NULL, linetype = NULL, fill = NULL)
  if (isTRUE(visual.info$add.censor.mark)) {
    marks <- do.call(rbind, lapply(names(curves), function(label) {
      curve <- curves[[label]]
      time <- curve$censor
      if (!length(time)) return(NULL)
      cif <- .cifcurve_interval_at(curve, time)$cif
      data.frame(time = time, estimate = if (type.y == "risk") cif else 1 - cif,
                 strata = factor(label, levels = names(curves)))
    }))
    if (!is.null(marks)) {
      p <- p + ggplot2::geom_point(data = marks[!is.na(marks$estimate), ],
        mapping = if (length(fit$strata)) ggplot2::aes(colour = .data$strata),
        shape = visual.info$shape.censor.mark, size = visual.info$size.censor.mark,
        show.legend = FALSE)
    }
  }
  attr(p, "cifmodeling_interval_scale") <- type.y
  list(out_survfit_object = p, type.y = type.y,
       label.y = if (type.y == "risk") "Cumulative incidence" else "Survival")
}

.cifplot_interval_table <- function(plot, fit, visual.info, style.info) {
  if (!isTRUE(visual.info$add.risktable) && !isTRUE(visual.info$add.estimate.table)) {
    return(NULL)
  }
  built <- suppressWarnings(ggplot2::ggplot_build(plot))
  panel <- built$layout$panel_params[[1L]]
  times <- panel$x$get_breaks()
  times <- times[is.finite(times) & times >= 0]
  if (!length(times)) return(NULL)
  info <- summary(fit, times = times, extend = TRUE, data.frame = TRUE)
  keys <- names(fit$reporting$curves)
  labels <- .cif_display_levels(keys, attr(plot, "cifmodeling_display_labels"),
                               fallback = sub("^[^=]*=", "", keys))
  if (!length(fit$strata)) info$strata <- keys[1L]
  colour <- if (length(fit$strata)) built$plot$scales$get_scales("colour")
  order <- if (is.null(colour)) keys else intersect(colour$get_limits(), keys)
  info <- info[info$strata %in% order, , drop = FALSE]
  info$strata <- factor(info$strata, levels = rev(order))
  label.map <- stats::setNames(unname(labels), keys)
  symbol <- switch(visual.info$symbol.risk.table %||% "",
                   square = "\u25A0", circle = "\u25CF", triangle = "\u25B2", "")
  if (length(order) > 1L && nzchar(symbol)) label.map[] <- paste(symbol, label.map)
  size <- visual.info$font.size.risk.table %||% 3
  stats <- list()
  interpolated <- info$estimate.status == "interpolated"
  if (isTRUE(visual.info$add.risktable)) {
    stats[["Confirmed at risk"]] <- as.character(info$n.risk.confirmed)
  }
  if (isTRUE(visual.info$add.estimate.table)) {
    estimate <- if (identical(attr(plot, "cifmodeling_interval_scale"), "risk")) info$cif else info$surv
    stats[["Estimate"]] <- ifelse(is.na(estimate), "\u2014",
                                  paste0(sprintf("%.2f", estimate), ifelse(interpolated, "*", "")))
  }
  .data <- NULL
  tables <- lapply(names(stats), function(name) {
    info$value <- stats[[name]]
    caption <- if (name == "Estimate" && any(interpolated)) {
      "* Interpolated assuming uniform allocation within the support interval."
    }
    ggplot2::ggplot(info, ggplot2::aes(x = .data$time, y = .data$strata,
                                         label = .data$value)) +
      ggplot2::geom_text(size = size, family = style.info$font.family) +
      ggplot2::scale_y_discrete(labels = if (length(fit$strata)) label.map else "") +
      ggplot2::labs(title = name, caption = caption) +
      plot_theme_risktable_font(font.family = style.info$font.family,
                                plot.title.size = style.info$font.size)
  })
  aligned <- ggsurvfit::ggsurvfit_align_plots(c(list(plot), tables))
  patchwork::wrap_plots(aligned, ncol = 1L,
                        heights = c(1, rep(.12 + .07 * length(order), length(tables))))
}

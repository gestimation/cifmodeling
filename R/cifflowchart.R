#' Create a flowchart of exclusions, treatment groups, and outcome status
#'
#' `cifflowchart()` creates a lightweight flowchart of exclusions, treatment
#' groups, and outcome status with `DiagrammeR`. It can be used as a simple
#' CONSORT-like flowchart helper, but it is not intended to fully automate a
#' CONSORT-compliant diagram.
#'
#' @param formula A formula. Supported forms are `response ~ arm`,
#'   `response ~ 1`, and `Event(time, status) ~ arm`.
#' @param data A data frame.
#' @param time.point Optional time point for `Event(time, status)` formulas.
#' @param pre.exclude,post.exclude Optional exclusion variables. Logical vectors
#'   use `TRUE` as excluded. Character and factor vectors use non-missing,
#'   non-empty values as exclusion reasons.
#' @param subset.condition Optional expression evaluated in `data` before
#'   counting.
#' @param na.action Included for consistency with other formula interfaces.
#' @param outcome.type Reserved for future use.
#' @param code.event1,code.event2,code.censoring Status codes for event-history
#'   formulas.
#' @param label.strata Optional labels for treatment groups.
#' @param order.strata Optional order for treatment groups.
#' @param label.events Optional labels for outcome states. A named vector may be
#'   used to map existing labels to display labels.
#' @param title Optional graph title.
#' @param percent Logical; show percentages in nodes.
#' @param digits Number of digits for percentages.
#' @param ... Reserved for future use.
#'
#' @return A `DiagrammeR::grViz()` htmlwidget.
#' @export
cifflowchart <- function(formula,
                         data,
                         time.point = NULL,
                         pre.exclude = NULL,
                         post.exclude = NULL,
                         subset.condition = NULL,
                         na.action = na.omit,
                         outcome.type = NULL,
                         code.event1 = 1,
                         code.event2 = 2,
                         code.censoring = 0,
                         label.strata = NULL,
                         order.strata = NULL,
                         label.events = NULL,
                         title = NULL,
                         percent = TRUE,
                         digits = 1,
                         ...) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' is required to use cifflowchart().", call. = FALSE)
  }
  prepared <- .flowchart_prepare_data(
    formula = formula, data = data, time.point = time.point,
    pre.exclude = substitute(pre.exclude), post.exclude = substitute(post.exclude),
    subset.condition = substitute(subset.condition), na.action = na.action,
    outcome.type = outcome.type, code.event1 = code.event1, code.event2 = code.event2,
    code.censoring = code.censoring, label.strata = label.strata,
    order.strata = order.strata, label.events = label.events,
    percent = percent, digits = digits, ...)
  DiagrammeR::grViz(.flowchart_make_dot(prepared, title = title))
}

.flowchart_prepare_data <- function(formula, data, time.point = NULL,
                                    pre.exclude = NULL, post.exclude = NULL,
                                    subset.condition = NULL, na.action = na.omit,
                                    outcome.type = NULL, code.event1 = 1,
                                    code.event2 = 2, code.censoring = 0,
                                    label.strata = NULL, order.strata = NULL,
                                    label.events = NULL, percent = TRUE,
                                    digits = 1, ...) {
  if (!inherits(formula, "formula")) stop("'formula' must be a formula.", call. = FALSE)
  data <- as.data.frame(data)
  if (!is.null(subset.condition) && !identical(subset.condition, quote(NULL))) {
    keep <- eval(subset.condition, data, parent.frame())
    if (!is.logical(keep) || length(keep) != nrow(data)) stop("'subset.condition' must evaluate to a logical vector with one value per row.", call. = FALSE)
    data <- data[keep & !is.na(keep), , drop = FALSE]
  }

  lhs <- formula[[2]]
  rhs <- formula[[3]]
  rhs_vars <- all.vars(rhs)
  if (length(rhs_vars) > 1) stop("The right-hand side of 'formula' must be ~ 1 or a single group variable.", call. = FALSE)
  has_group <- length(rhs_vars) == 1
  mode <- if (is.call(lhs) && identical(lhs[[1]], as.name("Event"))) "event" else "categorical"
  if (mode == "event" && length(lhs) != 3) stop("Event() formulas must be Event(time, status); Event(status) is not supported.", call. = FALSE)
  if (mode == "categorical" && !is.null(time.point)) stop("'time.point' can only be used with Event(time, status) formulas.", call. = FALSE)

  n <- nrow(data)
  pre <- .flowchart_exclude(eval(pre.exclude, data, parent.frame()), n, "pre.exclude")
  post_all <- .flowchart_exclude(eval(post.exclude, data, parent.frame()), n, "post.exclude")
  both <- pre$excluded & post_all$excluded
  if (any(both)) warning("Some rows are both pre.exclude and post.exclude; pre.exclude takes precedence.", call. = FALSE)

  included <- !pre$excluded
  group_raw <- if (has_group) data[[rhs_vars]] else factor(rep("All", n), levels = "All")
  group_chr <- .flowchart_levels(group_raw, order.strata)
  group_chr[is.na(group_chr)] <- "Missing group"
  if (has_group && any(is.na(group_raw))) warning("Missing group values are displayed as 'Missing group'.", call. = FALSE)
  group_label <- .flowchart_relabel(group_chr, label.strata)

  outcome <- if (mode == "categorical") {
    out <- data[[all.vars(lhs)[1]]]
    vals <- .flowchart_levels(out, NULL)
    vals[is.na(vals)] <- "Missing outcome"
    .flowchart_relabel(vals, label.events)
  } else {
    time <- data[[as.character(lhs[[2]])]]; status <- data[[as.character(lhs[[3]])]]
    .flowchart_event_status(time, status, time.point, code.event1, code.event2, code.censoring, label.events)
  }

  .flowchart_count_nodes(pre, post_all, included, group_chr, group_label, outcome,
                         has_group, percent, digits)
}

.flowchart_count_nodes <- function(pre, post_all, included, group_key, group_label,
                                   outcome, has_group, percent, digits) {
  total_n <- length(included)
  pre_counts <- .flowchart_tab(pre$reason[pre$excluded])
  groups <- if (has_group) {
    lev <- attr(group_key, "levels", exact = TRUE)
    c(lev[lev %in% group_key[included]], setdiff(unique(group_key[included]), lev))
  } else "All"
  group_labels <- vapply(groups, function(g) group_label[match(g, group_key)], character(1))
  group_n <- vapply(groups, function(g) sum(included & group_key == g), integer(1))
  post_counts <- lapply(groups, function(g) .flowchart_tab(post_all$reason[included & group_key == g & post_all$excluded]))
  analysis <- included & !post_all$excluded
  outcome_counts <- lapply(groups, function(g) .flowchart_tab(outcome[analysis & group_key == g]))
  analysis_n <- vapply(groups, function(g) sum(analysis & group_key == g), integer(1))
  list(total_n = total_n, pre_counts = pre_counts, groups = groups,
       group_labels = unname(group_labels), group_n = group_n,
       post_counts = post_counts, outcome_counts = outcome_counts,
       analysis_n = analysis_n, has_group = has_group, percent = percent,
       digits = digits)
}

.flowchart_make_dot <- function(x, title = NULL) {
  esc <- function(z) gsub('"', '\\"', z, fixed = TRUE)
  label <- function(name, n, denom) {
    paste0(esc(name), "\\nn = ", n, if (x$percent && denom > 0) sprintf(paste0(" (%.", x$digits, "f%%)"), 100 * n / denom) else "")
  }
  lines <- c("digraph cifflowchart {", "graph [rankdir = TB];", "node [shape = box, style = rounded, fontname = Helvetica];", "edge [fontname = Helvetica];")
  if (!is.null(title)) lines <- c(lines, paste0("labelloc = t; label = \"", esc(title), "\";"))
  lines <- c(lines, paste0("all [label = \"", label("All patients", x$total_n, x$total_n), "\"];"))
  parent <- "all"
  if (length(x$pre_counts)) {
    lines <- c(lines, paste0("pre [label = \"", label("Excluded before treatment", sum(x$pre_counts), x$total_n), "\n", esc(.flowchart_reason_lines(x$pre_counts)), "\"];"), "all -> pre;", paste0("included [label = \"", label("Included before treatment", sum(x$group_n), x$total_n), "\"];"), "all -> included;")
    parent <- "included"
  }
  for (i in seq_along(x$groups)) {
    gp <- paste0("group", i); lines <- c(lines, paste0(gp, " [label = \"", label(x$group_labels[i], x$group_n[i], sum(x$group_n)), "\"];"), paste0(parent, " -> ", gp, ";"))
    if (length(x$post_counts[[i]])) lines <- c(lines, paste0("post", i, " [label = \"", label("Excluded after treatment", sum(x$post_counts[[i]]), x$group_n[i]), "\n", esc(.flowchart_reason_lines(x$post_counts[[i]])), "\"];"), paste0(gp, " -> post", i, ";"))
    for (j in seq_along(x$outcome_counts[[i]])) lines <- c(lines, paste0("out", i, "_", j, " [label = \"", label(names(x$outcome_counts[[i]])[j], unname(x$outcome_counts[[i]])[j], x$analysis_n[i]), "\"];"), paste0(gp, " -> out", i, "_", j, ";"))
  }
  paste(c(lines, "}"), collapse = "\n")
}

.flowchart_exclude <- function(x, n, arg) {
  if (is.null(x)) return(list(excluded = rep(FALSE, n), reason = rep(NA_character_, n)))
  if (length(x) != n) stop(sprintf("'%s' must have one value per row.", arg), call. = FALSE)
  if (is.logical(x)) {
    if (any(is.na(x))) stop(sprintf("'%s' must not contain NA when it is logical.", arg), call. = FALSE)
    return(list(excluded = x, reason = ifelse(x, "Excluded", NA_character_)))
  }
  if (is.factor(x) || is.character(x)) {
    y <- as.character(x); excluded <- !is.na(y) & y != ""
    return(list(excluded = excluded, reason = ifelse(excluded, y, NA_character_)))
  }
  if (is.numeric(x)) stop(sprintf("Numeric '%s' is not supported.", arg), call. = FALSE)
  stop(sprintf("'%s' must be NULL, logical, factor, or character.", arg), call. = FALSE)
}

.flowchart_event_status <- function(time, status, tau, e1, e2, cens, labels) {
  out <- rep(NA_character_, length(status))
  if (is.null(tau)) {
    out[status == e1] <- "Event 1"; out[status == e2] <- "Event 2"; out[status == cens] <- "Censored"
  } else {
    out[time <= tau & status == e1] <- "Event 1 by tau"; out[time <= tau & status == e2] <- "Event 2 by tau"; out[time < tau & status == cens] <- "Censored before tau"; out[time >= tau] <- "Event-free at tau"
  }
  out[is.na(out)] <- "Missing outcome"
  .flowchart_relabel(out, labels)
}

.flowchart_levels <- function(x, order = NULL) {
  y <- as.character(x)
  lev <- if (!is.null(order)) as.character(order) else if (is.factor(x)) levels(x) else unique(y[!is.na(y)])
  out <- as.character(factor(y, levels = lev))
  attr(out, "levels") <- lev
  out
}

.flowchart_relabel <- function(x, labels = NULL) {
  if (is.null(labels)) return(x)
  lab <- as.character(labels)
  if (!is.null(names(labels))) {
    lev <- attr(x, "levels", exact = TRUE)
    hit <- match(x, names(labels)); x[!is.na(hit)] <- lab[hit[!is.na(hit)]]
    if (!is.null(lev)) { lev_hit <- match(lev, names(labels)); lev[!is.na(lev_hit)] <- lab[lev_hit[!is.na(lev_hit)]]; attr(x, "levels") <- lev }
    x
  } else {
    ux <- unique(x[!is.na(x)]); hit <- match(x, ux); x[!is.na(hit)] <- lab[hit[!is.na(hit)]]
    lev <- attr(x, "levels", exact = TRUE); if (!is.null(lev)) { lev_hit <- match(lev, ux); lev[!is.na(lev_hit)] <- lab[lev_hit[!is.na(lev_hit)]]; attr(x, "levels") <- lev }
    x
  }
}
.flowchart_tab <- function(x) {
  lev <- attr(x, "levels", exact = TRUE)
  x <- x[!is.na(x)]
  if (!length(x)) return(integer(0))
  if (is.null(lev)) lev <- unique(x)
  table(factor(x, levels = lev[lev %in% x]))
}
.flowchart_reason_lines <- function(tab) paste(paste0(names(tab), ": ", as.integer(tab)), collapse = "\\n")

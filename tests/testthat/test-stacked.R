test_that("stacked invariants: upper equals cumulative CIF and top equals sum", {
  skip_if_not_installed("survival")
  df <- data.frame(
    time  = c(0, 1, 2, 2, 3, 3),
    cause = c("A", "A", "A", "B", "A", "B"),
    cif   = c(0, 0.10, 0.20, 0.05, 0.25, 0.10)
  )
  sf <- cifmodeling:::cifplot_build_stacked_survfit(df, cause_order = c("A", "B"), add_t0 = TRUE)
  bounds <- attr(sf, "stack_bounds")$data
  top_vals <- bounds[bounds$cause == "B", ]
  top_by_time <- tapply(top_vals$ymax, top_vals$time, mean)
  grid <- sort(unique(df$time))
  causes <- c("A", "B")
  step_matrix <- vapply(causes, function(cause) {
    idx <- df$cause == cause
    vals <- rep(NA_real_, length(grid))
    if (any(idx)) {
      pos <- match(df$time[idx], grid)
      vals[pos] <- df$cif[idx]
    }
    last <- 0
    for (i in seq_along(vals)) {
      if (!is.na(vals[i])) last <- vals[i]
      vals[i] <- last
    }
    vals[!is.finite(vals)] <- 0
    cummax(pmin(pmax(vals, 0), 1))
  }, numeric(length(grid)))
  sum_cif <- colSums(step_matrix)
  names(sum_cif) <- grid
  expect_equal(as.numeric(top_by_time), as.numeric(sum_cif[names(top_by_time)]), tolerance = 1e-12)
  expect_true(all(top_by_time >= 0 & top_by_time <= 1, na.rm = TRUE))
})

test_that("stacked fallback engine produces a ggplot", {
  df <- data.frame(time = c(0, 1, 2), cause = c("A", "A", "B"), cif = c(0, 0.2, 0.1))
  sf <- cifmodeling:::cifplot_build_stacked_survfit(df, add_t0 = TRUE)
  p <- cifmodeling:::cifplot_draw_stacked_from_survfit(sf, engine = "ggplot")
  expect_s3_class(p, "ggplot")
})

test_that("cifplot stacked rendering returns a plot", {
  skip_if_not_installed("survival")
  df <- data.frame(time = c(1, 2, 3, 4), epsilon = c(1, 2, 0, 1))
  fit <- cifmodeling::cifcurve(cifmodeling::Event(time, epsilon) ~ 1,
                               data = df,
                               outcome.type = "competing-risk")
  p <- cifmodeling::cifplot(fit,
                            type.y = "stacked",
                            engine = "ggplot",
                            add.conf = FALSE,
                            add.risktable = FALSE,
                            add.estimate.table = FALSE,
                            add.censor.mark = FALSE,
                            add.competing.risk.mark = FALSE,
                            add.intercurrent.event.mark = FALSE,
                            add.quantile = FALSE)
  expect_s3_class(p, "ggplot")
})

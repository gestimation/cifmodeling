test_that("default n.risk.type preserves weighted risk sets", {
  df <- data.frame(
    t = c(1, 2, 2, 3),
    epsilon = c(1, 0, 2, 0),   # <- 2 を入れて competing-risk を成立させる
    w = c(1, 2, 3, 4)
  )

  fit_default <- cifcurve(Event(t, epsilon) ~ 1,
                          data = df,
                          weights = "w",
                          outcome.type = "competing-risk",
                          conf.type = "none")
  fit_weighted <- cifcurve(Event(t, epsilon) ~ 1,
                           data = df,
                           weights = "w",
                           outcome.type = "competing-risk",
                           conf.type = "none",
                           n.risk.type = "weighted")

  expect_equal(fit_default$n.risk, fit_weighted$n.risk)
  expect_equal(fit_default$n.risk, fit_default$n.risk.weighted)
  expect_equal(fit_default$n.risk.type, "weighted")
})


test_that("n.risk.type switches among weighted, unweighted, and ess outputs", {
  df <- data.frame(
    t = c(1, 2, 2, 3, 4),
    epsilon = c(1, 0, 2, 0, 1),
    w = c(1, 2, 3, 4, 5)
  )

  fit_weighted <- cifcurve(Event(t, epsilon) ~ 1,
                           data = df,
                           weights = "w",
                           outcome.type = "competing-risk",
                           conf.type = "none",
                           n.risk.type = "weighted")
  fit_unweighted <- cifcurve(Event(t, epsilon) ~ 1,
                             data = df,
                             weights = "w",
                             outcome.type = "competing-risk",
                             conf.type = "none",
                             n.risk.type = "unweighted")
  fit_ess <- cifcurve(Event(t, epsilon) ~ 1,
                      data = df,
                      weights = "w",
                      outcome.type = "competing-risk",
                      conf.type = "none",
                      n.risk.type = "ess")

  expect_equal(fit_weighted$n.risk,   ceiling(fit_weighted$n.risk.weighted))
  expect_equal(fit_unweighted$n.risk, ceiling(fit_unweighted$n.risk.unweighted))
  expect_equal(fit_ess$n.risk,        ceiling(fit_ess$n.risk.ess))

  expect_equal(fit_weighted$n.risk.weighted, fit_unweighted$n.risk.weighted)
  expect_equal(fit_weighted$n.risk.unweighted, fit_unweighted$n.risk.unweighted)
  expect_equal(fit_weighted$n.risk.ess, fit_unweighted$n.risk.ess)

  expect_true(all(c("n.risk.weighted","n.risk.unweighted","n.risk.ess","n.risk.type") %in% names(fit_weighted)))
  expect_equal(fit_unweighted$n.risk.type, "unweighted")
  expect_equal(fit_ess$n.risk.type, "ess")
})

test_that("Kish ESS matches hand calculations with and without strata", {
  compute_ess <- function(time, weight) {
    uniq_t <- sort(unique(time))
    vapply(
      uniq_t,
      function(tt) {
        idx <- time >= tt
        sum_w <- sum(weight[idx])
        sum_w2 <- sum(weight[idx]^2)
        if (sum_w <= 0 || sum_w2 <= 0) 0 else (sum_w^2) / sum_w2
      },
      numeric(1)
    )
  }

  df <- data.frame(
    t = c(1, 2, 2, 3),
    epsilon = c(1, 0, 2, 0),   # <- 2 を入れる
    w = c(1, 2, 3, 4)
  )
  fit_no_strata <- cifcurve(Event(t, epsilon) ~ 1,
                            data = df,
                            weights = "w",
                            outcome.type = "competing-risk",
                            conf.type = "none",
                            n.risk.type = "ess")
  expect_equal(
    fit_no_strata$n.risk.ess,
    compute_ess(df$t, df$w),
    tolerance = 1e-12
  )

  df_strata <- transform(df, strata = factor(c("a", "a", "b", "b")))
  fit_strata <- cifcurve(Event(t, epsilon) ~ strata,
                         data = df_strata,
                         weights = "w",
                         outcome.type = "competing-risk",
                         conf.type = "none",
                         n.risk.type = "ess")
  expect_equal(
    fit_strata$n.risk.ess,
    c(
      compute_ess(df_strata$t[df_strata$strata == "a"], df_strata$w[df_strata$strata == "a"]),
      compute_ess(df_strata$t[df_strata$strata == "b"], df_strata$w[df_strata$strata == "b"])
    ),
    tolerance = 1e-12
  )
})

test_that("weights accepts quoted and unquoted column names (survfit-like)", {
  df <- data.frame(
    t = c(1, 2, 2, 3, 4),
    epsilon = c(1, 0, 2, 0, 1),
    w = c(1, 2, 3, 4, 5)
  )

  fit_quoted <- cifcurve(Event(t, epsilon) ~ 1,
                         data = df,
                         weights = "w",
                         outcome.type = "competing-risk",
                         conf.type = "none")

  fit_unquoted <- cifcurve(Event(t, epsilon) ~ 1,
                           data = df,
                           weights = w,
                           outcome.type = "competing-risk",
                           conf.type = "none")

  expect_equal(fit_unquoted$n.risk.weighted, fit_quoted$n.risk.weighted)
  expect_equal(fit_unquoted$n.risk, fit_quoted$n.risk)
})

test_that("cifplot forwards n.risk.type when fitting from a formula", {
  df <- data.frame(
    t = c(1, 2, 2, 3, 4),
    epsilon = c(1, 0, 2, 0, 1),
    w = c(1, 2, 3, 4, 5)
  )

  fit_ref <- cifcurve(
    Event(t, epsilon) ~ 1,
    data = df,
    weights = "w",
    outcome.type = "competing-risk",
    conf.type = "none",
    n.risk.type = "ess"
  )

  plot_out <- cifplot(
    Event(t, epsilon) ~ 1,
    data = df,
    weights = "w",
    outcome.type = "competing-risk",
    conf.type = "none",
    n.risk.type = "ess",
    add.risktable = FALSE,
    print.panel = FALSE
  )

  expect_equal(plot_out$survfit.info$survfit$n.risk.type, "ess")
  expect_equal(plot_out$survfit.info$survfit$n.risk, fit_ref$n.risk)
})

test_that("cifpanel allows panel-wise n.risk.type selection", {
  df <- data.frame(
    t = c(1, 2, 2, 3, 4),
    epsilon = c(1, 0, 2, 0, 1),
    w = c(1, 2, 3, 4, 5)
  )

  panel_out <- cifpanel(
    formulas = list(Event(t, epsilon) ~ 1, Event(t, epsilon) ~ 1),
    data = df,
    weights = "w",
    outcome.type = "competing-risk",
    conf.type = "none",
    n.risk.type = list("weighted", "e"),
    legend.position = "none",
    print.panel = FALSE
  )

  expect_equal(panel_out$survfit.info$n.risk.type[[1]], "weighted")
  expect_equal(panel_out$survfit.info$n.risk.type[[2]], "ess")
  expect_equal(panel_out$survfit.info$curves[[1]]$n.risk.type, "weighted")
  expect_equal(panel_out$survfit.info$curves[[2]]$n.risk.type, "ess")

  fit_panel2 <- cifcurve(
    Event(t, epsilon) ~ 1,
    data = df,
    weights = "w",
    outcome.type = "competing-risk",
    conf.type = "none",
    n.risk.type = "ess"
  )
  expect_equal(panel_out$survfit.info$curves[[2]]$n.risk, fit_panel2$n.risk)
})

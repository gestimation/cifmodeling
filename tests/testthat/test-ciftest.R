# tests/testthat/test-ciftest-logrank-components.R
testthat::local_edition(3)

testthat::test_that("calculate_log_rank() matches survdiff for unweighted log-rank (rho=0,gamma=0)", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,4,5,6,6,7),
    status = c(1,0,1,2,1,0,1,1,0,0), # 1=event1, 2=other (treated as censor here), 0=censor
    trt    = factor(c("A","A","B","B","A","C","B","C","C","A"), levels = c("A","B","C")),
    str    = factor(c("s1","s1","s1","s1","s2","s2","s2","s2","s1","s2"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)   # treat status==2 as censor for cause1
  w   <- rep(1, nrow(df))

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = 0, gamma = 0
  )

  # survdiff comparison (gamma not supported; gamma=0 here)
  ev1 <- as.integer(df$status == 1)
  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt, data = df, rho = 0)

  # chisq from our score/var (df=K-1)
  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)

  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-8)

  # also check O-E matches our score components for non-ref groups
  oe <- as.numeric(sdiff$obs - sdiff$exp)
  names(oe) <- names(sdiff$obs)
  # our score corresponds to groups except reference (first level)
  testthat::expect_equal(unname(comp$score), unname(oe[-1]), tolerance = 1e-8)
})

testthat::test_that("matches survdiff for FH rho>0 (gamma=0)", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,4,5,6,6,7),
    status = c(1,0,1,2,1,0,1,1,0,0),
    trt    = factor(c("A","A","B","B","A","C","B","C","C","A"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))
  ev1 <- as.integer(df$status == 1)

  rho <- 1
  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = rho, gamma = 0
  )
  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt, data = df, rho = rho)

  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)
  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-6)
})

testthat::test_that("stratified version matches survdiff with strata()", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,4,5,6,6,7),
    status = c(1,0,1,2,1,0,1,1,0,0),
    trt    = factor(c("A","A","B","B","A","C","B","C","C","A"), levels = c("A","B","C")),
    str    = factor(c("s1","s1","s1","s1","s2","s2","s2","s2","s1","s2"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))
  ev1 <- as.integer(df$status == 1)

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = df$str,
    data = df, exposure = "trt",
    rho = 0, gamma = 0
  )

  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt + survival::strata(str), data = df, rho = 0)
  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)

  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-8)
})

testthat::test_that("reference level via code.exposure.ref is respected", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,4,5,6,6,7),
    status = c(1,0,1,1,1,1,0,0),
    trt    = factor(c("A","A","B","A","C","B","C","A"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))
  ev1 <- as.integer(df$status == 1)

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    code.exposure.ref = "B",
    rho = 0, gamma = 0,
    prefix = "a"
  )

  testthat::expect_equal(comp$ref, "B")
  testthat::expect_true(all(grepl("^a_", names(comp$score))))

  # Compare with survdiff on relevelled trt
  df2 <- df
  df2$trt <- stats::relevel(df2$trt, ref = "B")
  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt, data = df2, rho = 0)
  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)

  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-8)
})

testthat::test_that("integer weights equal replication (frequency-weight sanity check)", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,5,6,6),
    status = c(1,0,1,2,1,1,1,0),
    trt    = factor(c("A","A","B","B","A","C","B","C"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)

  w_int <- c(1,2,1,1,2,1,1,2)   # integer weights
  comp_w <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w_int,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = 0, gamma = 0
  )

  df_rep <- df[rep(seq_len(nrow(df)), w_int), , drop = FALSE]
  eps_rep <- ifelse(df_rep$status == 1, 1L, 0L)
  comp_rep <- cifmodeling:::calculate_log_rank(
    t = df_rep$time, epsilon = eps_rep, weights = rep(1, nrow(df_rep)),
    strata = rep.int(1L, nrow(df_rep)),
    data = df_rep, exposure = "trt",
    rho = 0, gamma = 0
  )

  testthat::expect_equal(comp_w$score, comp_rep$score, tolerance = 1e-10)
  testthat::expect_equal(comp_w$var,   comp_rep$var,   tolerance = 1e-10)
})

testthat::test_that("gamma>0 runs and returns finite outputs (smoke test)", {
  df <- data.frame(
    time   = c(1,2,2,4,5,6,7,7),
    status = c(1,0,1,1,0,1,1,0),
    trt    = factor(c("A","A","B","A","C","B","C","A"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = 0, gamma = 1
  )

  testthat::expect_true(all(is.finite(comp$score)))
  testthat::expect_true(all(is.finite(comp$var)))
})

testthat::test_that("errors on one-level exposure", {
  df <- data.frame(
    time = c(1,2,3),
    status = c(1,0,1),
    trt = factor(c("A","A","A"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w <- rep(1, nrow(df))

  testthat::expect_error(
    cifmodeling:::calculate_log_rank(
      t = df$time, epsilon = eps, weights = w,
      strata = rep.int(1L, nrow(df)),
      data = df, exposure = "trt",
      rho = 0, gamma = 0
    ),
    "only one level|no valid levels|Effect estimation is not possible",
    ignore.case = TRUE
  )
})

testthat::test_that("A12 from calculate_A12_logrank_weightit matches direct finite-difference under linear weights", {
  # Skip gracefully if not yet available
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found in cifmodeling namespace.")
  }
  if (!exists("calculate_log_rank", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_log_rank not found in cifmodeling namespace.")
  }

  set.seed(1)

  n <- 60
  df <- data.frame(
    time = sample(1:8, n, replace = TRUE),
    status_raw = sample(c(0L, 1L, 2L), n, replace = TRUE, prob = c(0.55, 0.35, 0.10)),
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trt = factor(sample(c("A","B","C"), n, replace = TRUE), levels = c("A","B","C"))
  )

  # epsilon normalized for log-rank: 1=event1, others treated as censor (0)
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)

  # Ensure at least one event in each stratum (otherwise score in that stratum is zero and derivative can be degenerate)
  # If not, tweak a couple of records deterministically
  if (sum(epsilon[df$str == "s1"]) == 0L) epsilon[which(df$str == "s1")[1L]] <- 1L
  if (sum(epsilon[df$str == "s2"]) == 0L) epsilon[which(df$str == "s2")[1L]] <- 1L

  t <- df$time
  strata <- df$str

  # ----- Known (linear) weight mechanism -----
  # w(beta) = w0 + X %*% beta, with beta0=0 => w0
  # then dw/dB = X exactly
  q <- 2
  Xw <- matrix(runif(n*q, min = -0.2, max = 0.2), nrow = n, ncol = q)
  colnames(Xw) <- paste0("B", 1:q)

  w0 <- rep(1, n)  # positive baseline weights

  # Mock weightit object with required fields
  weightit_mock <- list(
    weights = w0,
    Mparts = list(dw_dBtreat = Xw)
  )

  rho <- 0
  gamma <- 0
  prob.bound <- 1e-7

  # ----- Run A12 function (under test) -----
  out <- cifmodeling:::calculate_A12_logrank_weightit(
    t = t,
    epsilon = as.integer(epsilon),
    strata = strata,
    data = df,
    exposure = "trt",
    weightit = weightit_mock,
    rho = rho,
    gamma = gamma,
    prob.bound = prob.bound,
    fd_rel_step = 1e-4   # make the step not too tiny for numerical stability
  )

  A12_hat <- out$A12
  steps <- out$fd_steps

  # ----- "Truth" via direct finite difference w.r.t beta using the known mapping w(beta)=w0+Xw beta -----
  # Use the SAME step sizes as calculate_A12... to avoid step-size confounding.
  base <- cifmodeling:::calculate_log_rank(
    t = t,
    epsilon = as.integer(epsilon),
    weights = w0,
    strata = strata,
    data = df,
    exposure = "trt",
    rho = rho,
    gamma = gamma,
    prob.bound = prob.bound
  )

  p <- length(base$score)   # p = K-1
  A12_true <- matrix(0, nrow = p, ncol = q)
  rownames(A12_true) <- names(base$score)
  colnames(A12_true) <- colnames(Xw)

  for (j in seq_len(q)) {
    h <- steps[j]
    if (h == 0) next

    # beta0 +/- h e_j => w +/- h * Xw[,j]
    w_plus  <- w0 + h * Xw[, j]
    w_minus <- w0 - h * Xw[, j]

    s_plus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_plus,
      strata = strata, data = df, exposure = "trt",
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    s_minus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_minus,
      strata = strata, data = df, exposure = "trt",
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    dU_dBj <- (as.numeric(s_plus) - as.numeric(s_minus)) / (2 * h)

    # A12 is on mean-score scale (1/n) dU_total/dB^T
    A12_true[, j] <- dU_dBj / n
  }

  # ----- Assertions -----
  testthat::expect_equal(dim(A12_hat), dim(A12_true))
  testthat::expect_equal(rownames(A12_hat), rownames(A12_true))
  testthat::expect_equal(colnames(A12_hat), colnames(A12_true))

  # Linear weights => should match essentially to machine precision
  testthat::expect_equal(A12_hat, A12_true, tolerance = 1e-10)
})

testthat::test_that("A12 is zero when dw/dB is zero", {
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found in cifmodeling namespace.")
  }

  set.seed(2)
  n <- 40
  df <- data.frame(
    time = sample(1:6, n, replace = TRUE),
    status_raw = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.6, 0.4)),
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trt = factor(sample(c("A","B","C"), n, replace = TRUE), levels = c("A","B","C"))
  )
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)

  w0 <- rep(1, n)
  Xw0 <- matrix(0, nrow = n, ncol = 3)
  colnames(Xw0) <- paste0("B", 1:3)

  weightit_mock <- list(weights = w0, Mparts = list(dw_dBtreat = Xw0))

  out <- cifmodeling:::calculate_A12_logrank_weightit(
    t = df$time,
    epsilon = as.integer(epsilon),
    strata = df$str,
    data = df,
    exposure = "trt",
    weightit = weightit_mock,
    rho = 0,
    gamma = 0,
    fd_rel_step = 1e-4
  )

  testthat::expect_true(all(abs(out$A12) < 1e-14))
})


# --- helper: "true" A12 via direct finite-difference in beta (exact w(beta±h e_j)) ---
A12_true_beta_fd <- function(t, epsilon, strata, data, exposure,
                             beta0, w_map, h_vec,
                             rho = 0, gamma = 0, prob.bound = 1e-7) {
  base <- cifmodeling:::calculate_log_rank(
    t = t, epsilon = as.integer(epsilon), weights = w_map(beta0),
    strata = strata, data = data, exposure = exposure,
    rho = rho, gamma = gamma, prob.bound = prob.bound
  )
  score0 <- as.numeric(base$score)
  p <- length(score0)
  q <- length(beta0)
  n <- length(t)

  A12 <- matrix(0, nrow = p, ncol = q)
  rownames(A12) <- names(base$score) %||% paste0("score", seq_len(p))
  colnames(A12) <- paste0("B", seq_len(q))

  for (j in seq_len(q)) {
    h <- h_vec[j]
    if (!is.finite(h) || h == 0) next

    bp <- beta0; bp[j] <- bp[j] + h
    bm <- beta0; bm[j] <- bm[j] - h

    s_plus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_map(bp),
      strata = strata, data = data, exposure = exposure,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    s_minus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_map(bm),
      strata = strata, data = data, exposure = exposure,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    dU <- (as.numeric(s_plus) - as.numeric(s_minus)) / (2 * h)
    A12[, j] <- dU / n  # mean-score scale
  }
  A12
}

# --- helper: run A12 function under test and compute error vs truth ---
run_and_error <- function(t, epsilon, strata, data, exposure,
                          weightit_mock, beta0, w_map,
                          fd_rel_step, rho = 0, gamma = 0, prob.bound = 1e-7) {

  out <- cifmodeling:::calculate_A12_logrank_weightit(
    t = t,
    epsilon = as.integer(epsilon),
    strata = strata,
    data = data,
    exposure = exposure,
    weightit = weightit_mock,
    rho = rho, gamma = gamma, prob.bound = prob.bound,
    fd_rel_step = fd_rel_step
  )
  A12_hat <- out$A12
  h_vec <- out$fd_steps

  A12_true <- A12_true_beta_fd(
    t = t, epsilon = epsilon, strata = strata, data = data, exposure = exposure,
    beta0 = beta0, w_map = w_map, h_vec = h_vec,
    rho = rho, gamma = gamma, prob.bound = prob.bound
  )

  err <- max(abs(A12_hat - A12_true))
  step_scale <- max(h_vec^2)

  list(err = err, step2 = step_scale, out = out, A12_true = A12_true)
}

testthat::test_that("Nonlinear exp weights: A12 approximation error scales like O(h^2)", {
  # skip if functions are not available yet
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found.")
  }
  if (!exists("calculate_log_rank", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_log_rank not found.")
  }

  set.seed(10)
  n <- 120
  df <- data.frame(
    time = sample(1:10, n, replace = TRUE),
    status_raw = rbinom(n, 1, 0.35),          # 1=event1, 0=censor
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trt = factor(sample(c("A","B","C"), n, replace = TRUE), levels = c("A","B","C"))
  )
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)
  # ensure at least one event per stratum
  if (sum(epsilon[df$str == "s1"]) == 0L) epsilon[which(df$str == "s1")[1L]] <- 1L
  if (sum(epsilon[df$str == "s2"]) == 0L) epsilon[which(df$str == "s2")[1L]] <- 1L

  t <- df$time
  strata <- df$str
  exposure <- "trt"

  # nonlinear exp weights: w(beta)=exp(X beta)
  q <- 2
  X <- matrix(runif(n*q, min = -0.08, max = 0.08), nrow = n, ncol = q)
  colnames(X) <- paste0("B", 1:q)
  beta0 <- c(0.15, -0.10)

  w_map <- function(beta) as.numeric(exp(X %*% beta))

  w0 <- w_map(beta0)
  dw <- w0 * X  # exact dw/dB at beta0

  weightit_mock <- list(
    weights = w0,
    Mparts = list(dw_dBtreat = dw)
  )

  # two step sizes (halve rel step)
  r1 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 1e-3, rho = 0, gamma = 0)
  r2 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 5e-4, rho = 0, gamma = 0)

  testthat::expect_true(is.finite(r1$err) && is.finite(r2$err))
  testthat::expect_true(r2$err > 0)  # should not be exactly 0 for nonlinear mapping

  # If error is O(h^2), err_ratio ~ step2_ratio
  err_ratio  <- r1$err / r2$err
  step_ratio <- r1$step2 / r2$step2

  # allow slack (numerical noise / piecewise behavior)
  testthat::expect_true(err_ratio / step_ratio > 0.4 && err_ratio / step_ratio < 2.5)

  # and error should decrease when step decreases
  testthat::expect_true(r2$err < r1$err)
})

testthat::test_that("Nonlinear logit-IPW weights (binary): A12 approximation error scales like O(h^2)", {
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found.")
  }
  if (!exists("calculate_log_rank", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_log_rank not found.")
  }

  set.seed(11)
  n <- 150
  df <- data.frame(
    time = sample(1:12, n, replace = TRUE),
    status_raw = rbinom(n, 1, 0.30),
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trtbin = factor(rbinom(n, 1, 0.5), levels = c(0,1))  # binary exposure
  )
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)
  if (sum(epsilon[df$str == "s1"]) == 0L) epsilon[which(df$str == "s1")[1L]] <- 1L
  if (sum(epsilon[df$str == "s2"]) == 0L) epsilon[which(df$str == "s2")[1L]] <- 1L

  t <- df$time
  strata <- df$str
  exposure <- "trtbin"
  A <- as.integer(df$trtbin == 1L)

  q <- 2
  X <- matrix(runif(n*q, min = -0.12, max = 0.12), nrow = n, ncol = q)
  colnames(X) <- paste0("B", 1:q)
  beta0 <- c(0.10, -0.05)

  # define propensity and IPW weights
  p_map <- function(beta) as.numeric(stats::plogis(X %*% beta))

  # keep away from 0/1 but choose beta0, X so it should not bind
  clip <- function(p, lo = 1e-3, hi = 1 - 1e-3) pmin(pmax(p, lo), hi)

  w_map <- function(beta) {
    p <- clip(p_map(beta))
    A / p + (1 - A) / (1 - p)
  }

  p0 <- clip(p_map(beta0))
  w0 <- w_map(beta0)

  # analytic dw/dB at beta0 (no clipping active assumed)
  # dp/dB = p(1-p) X
  # dw/dB = X * [ -(A)(1-p)/p + (1-A)p/(1-p) ]
  fac <- (-(A) * (1 - p0) / p0) + ((1 - A) * p0 / (1 - p0))
  dw <- X * fac  # rowwise multiply

  weightit_mock <- list(
    weights = w0,
    Mparts = list(dw_dBtreat = dw)
  )

  r1 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 1e-3, rho = 0, gamma = 0)
  r2 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 5e-4, rho = 0, gamma = 0)

  testthat::expect_true(is.finite(r1$err) && is.finite(r2$err))
  testthat::expect_true(r2$err > 0)

  err_ratio  <- r1$err / r2$err
  step_ratio <- r1$step2 / r2$step2

  testthat::expect_true(err_ratio / step_ratio > 0.4 && err_ratio / step_ratio < 2.5)
  testthat::expect_true(r2$err < r1$err)
})

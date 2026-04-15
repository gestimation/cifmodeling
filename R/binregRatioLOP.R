#' Fit direct binomial regression for restricted mean time lost
#' using a log-odds product parameterization
#'
#' Fits an inverse probability of censoring weighted (IPCW) ratio regression model
#' for the conditional proportion of outcome attributable to a specific cause in
#' competing risks data. The target quantity is
#' \deqn{E\{Z_1(\tau)\mid X\} / E\{Z(\tau)\mid X\},}
#' where \eqn{Z_1(\tau)} is the cause-specific restricted time lost up to the time
#' horizon \eqn{\tau} for cause 1 and \eqn{Z(\tau)} is the corresponding total
#' restricted time lost up to \eqn{\tau}.
#' The model is parameterized through a nuisance log-odds product model together
#' with a log percentage ratio parameter for a binary exposure. The final column
#' of the design matrix is assumed to be a binary exposure coded as 0/1, and the
#' remaining columns are treated as adjustment covariates.
#'
#' Right censoring is handled through IPCW. When `cens.weights` is not supplied,
#' censoring weights are estimated internally from `cens.model` using
#' [mets::phreg()] and predicted censoring survival probabilities.
#'
#' @param formula A model formula with an `Event()` response on the left-hand side.
#'   The right-hand side defines the regression design matrix. The final column of
#'   the resulting design matrix must correspond to a binary exposure coded as 0/1.
#' @param data A data frame containing the variables in `formula`.
#' @param cause Integer or vector of integers specifying the event code(s) treated
#'   as the primary cause of interest.
#' @param time Numeric scalar giving the time horizon \eqn{\tau}.
#' @param beta Optional numeric vector of starting values. Its length must equal
#'   the number of columns in the design matrix. By default, nuisance coefficients
#'   are initialized at `0.1` and the exposure log percentage ratio parameter at
#'   `log(1.2)`.
#' @param type Character string specifying the estimating equation. `"I"` gives
#'   the basic IPCW estimator, and `"II"` gives the augmented estimator that
#'   updates the estimating equation using censoring-related augmentation terms.
#' @param offset Optional numeric vector of offsets. Defaults to a vector of zeros.
#' @param weights Optional observation weights for the estimating equation.
#'   Defaults to 1 for all observations.
#' @param cens.weights Optional censoring survival probabilities used for IPCW.
#'   When supplied, internal fitting of the censoring model is skipped.
#' @param cens.model A right-hand side formula for the censoring model used when
#'   `cens.weights` is `NULL`. The default `~ + 1` gives marginal censoring weights.
#' @param se Logical; if `TRUE`, compute robust influence-function-based standard
#'   errors.
#' @param kaplan.meier Logical; if `TRUE`, allow Kaplan-Meier-type baseline
#'   estimation in the censoring model fit when applicable.
#' @param cens.code Integer code used for censoring in the event variable.
#'   Defaults to `0`.
#' @param no.opt Logical; if `TRUE`, skip numerical optimization and evaluate the
#'   estimating equations at the supplied `beta`.
#' @param augmentation Optional numeric vector used to augment the estimating
#'   equation. Defaults to a vector of zeros.
#' @param outcome Character string specifying the outcome scale. Currently only
#'   `"rmtl"`,  restricted mean time lost decomposition up to `time`, is implemented.
#' @param model Character string identifying the model family. Currently only
#'   `"log-odds"` is implemented.
#' @param Ydirect Optional user-supplied outcome matrix replacing the internally
#'   constructed IPCW outcome.
#' @param ... Additional arguments passed through model-frame construction.
#'
#' @details
#' Let \eqn{A} denote the binary exposure in the final column of the design matrix
#' and let \eqn{L} denote the remaining covariates. The function models the
#' conditional percentage for the primary cause under a log-odds product
#' parameterization, returning fitted percentages under `A = 0` and `A = 1`.
#'
#' The returned object includes coefficient estimates, naive and robust variance
#' estimates, estimated influence functions, fitted percentages, and quantities
#' related to the censoring model fit.
#'
#' @return
#' An object of class `"binreg"` containing at least the following components:
#' \describe{
#'   \item{coef}{Estimated regression coefficients. The last coefficient is the
#'   log percentage ratio parameter for the binary exposure.}
#'   \item{se.robust}{Robust standard errors based on the estimated influence
#'   functions.}
#'   \item{robvar}{Robust variance-covariance matrix.}
#'   \item{iid}{Estimated influence functions used for robust variance estimation.}
#'   \item{p0, p1}{Estimated conditional percentages under exposure levels 0 and 1.}
#'   \item{p}{Estimated conditional percentage at the observed exposure level.}
#'   \item{converged}{Logical indicating whether the numerical solver reported
#'   convergence.}
#' }
#'
#' @seealso [Event()], [calculatePercentageLOP()], [polyreg()], [cifcurve()]
#'
#' @examples
#' ## event: 0 = censoring, 1 = primary cause, 2 = competing cause
#' data(diabetes.complications)
#'
#' fit <- binregRatioLOP(
#'   Event(t, epsilon) ~ fruitq1,
#'   data = diabetes.complications,
#'   time = 8,
#'   cause = 1,
#'   type = "I"
#' )
#'
#' fit$coef
#' fit$se.robust
#'
#' ## Estimated percentages under A = 0 and A = 1 when there are no adjustment covariates
#' calculatePercentageLOP(
#'   fit$coef,
#'   X_L = matrix(numeric(0), nrow = 1),
#'   offset = 0
#' )
#'
#' @name binregRatioLOP
#' @export
#' @importFrom stats model.frame model.extract model.matrix terms update.formula predict
#' @importFrom survival Surv
#' @importFrom nleqslv nleqslv
#' @importFrom numDeriv jacobian
#' @importFrom mets phreg construct_id sumstrata revcumsumstrata cumsumstrata
binregRatioLOP <- function(formula,
                           data,
                           cause = 1,
                           time = NULL,
                           beta = NULL,
                           type = c("II", "I"),
                           offset = NULL,
                           weights = NULL,
                           cens.weights = NULL,
                           cens.model = ~ + 1,
                           se = TRUE,
                           kaplan.meier = TRUE,
                           cens.code = 0,
                           no.opt = FALSE,
                           augmentation = NULL,
                           outcome = "rmtl",
                           model = "log-odds",
                           Ydirect = NULL,
                           ...) {
  cl <- match.call()
  type <- match.arg(type)
  outcome <- match.arg(outcome)

  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())

  Yresp <- model.extract(m, "response")
  if (!inherits(Yresp, "Event")) stop("Expected an 'Event'-object")

  if (ncol(Yresp) == 2) {
    exit <- Yresp[, 1]
    status <- Yresp[, 2]
  } else {
    stop("Only right-censored data are currently supported.")
  }

  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- untangle.specials(Terms, "cluster")
    Terms <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }

  if (!is.null(attributes(Terms)$specials$strata)) {
    ts <- untangle.specials(Terms, "strata")
    Terms <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
  }

  if (!is.null(attributes(Terms)$specials$offset)) {
    ts <- untangle.specials(Terms, "offset")
    Terms <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }

  X <- model.matrix(Terms, m)
  if (ncol(X) == 0) {
    stop("X must contain at least one column; last column must be binary exposure A.")
  }
  X <- as.matrix(X)

  call.id <- id
  conid <- construct_id(id, nrow(X), as.data = TRUE)
  name.id <- conid$name.id
  id <- conid$id
  nid <- conid$nid

  if (is.null(offset)) offset <- rep(0, length(exit))
  if (is.null(weights)) weights <- rep(1, length(exit))
  if (is.null(time)) stop("Must give time for modelling.")
  if (is.null(beta)) {
    p <- ncol(X)
    beta <- c(rep(0.1, p - 1L), log(1.2))
  }
  if (is.null(augmentation)) augmentation <- rep(0, ncol(X))

  do_censor_fit <- is.null(cens.weights)

  statusC <- 1L * (status %in% cens.code)
  statusE <- (status %in% cause) & (exit <= time)
  if (sum(statusE) == 0) warning("No events of type 1 before time.")

  ucauses <- sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc) == 0) Causes <- ucauses else Causes <- ucauses[-ccc]

  data$id__ <- id
  data$exit <- exit
  data$statusC <- statusC
  cens.strata <- cens.nstrata <- NULL

  nevent <- sum((status %in% cause) * (exit <= time))
  obs <- (exit <= time & (!statusC)) | (exit >= time)

  if (do_censor_fit) {
    kmt <- kaplan.meier
    formC <- update.formula(cens.model, Surv(exit, statusC) ~ .)
    resC <- phreg(formC, data)
    if (resC$p > 0) kmt <- FALSE
    exittime <- pmin(exit, time)
    cens.weights <- suppressWarnings(
      predict(
        resC, data,
        times = exittime,
        individual.time = TRUE,
        se = FALSE,
        km = kmt,
        tminus = TRUE
      )$surv
    )
    cens.strata <- resC$strata[order(resC$ord)]
    cens.nstrata <- resC$nstrata
  } else {
    formC <- NULL
    resC <- NULL
  }

  if (!is.null(Ydirect)) {
    Y <- Ydirect * obs / cens.weights
  } else if (outcome == "cif") {
    Y <- cbind(
      c((status %in% Causes) * (exit <= time) / cens.weights),
      c((status %in% cause)  * (exit <= time) / cens.weights)
    )
  } else {
    Y <- cbind(
      c((time - pmin(exit, time)) * obs / cens.weights),
      c((status %in% cause) * (time - pmin(exit, time)) * obs / cens.weights)
    )
  }
  Yipcw <- Y

  make_lop_prob <- function(pp, Xmat, offset_vec) {
    X_A <- as.vector(Xmat[, ncol(Xmat), drop = FALSE])

    if (!all(is.na(X_A) | (abs(X_A) < 1e-12 | abs(X_A - 1) < 1e-12))) {
      warning("The last column of X is assumed to be binary exposure A coded as 0/1.")
    }

    X_L <- Xmat[, -ncol(Xmat), drop = FALSE]
    p0_p1 <- calculatePercentageLOP(pp, X_L, offset_vec)
    p <- p0_p1[, 1] * (1 - X_A) + p0_p1[, 2] * X_A

    list(
      X_A = X_A,
      X_L = X_L,
      p0 = p0_p1[, 1],
      p1 = p0_p1[, 2],
      p  = p
    )
  }

  score_vec <- function(pp, Xmat, Ymat, weights_vec, offset_vec, augmentation_vec) {
    fit <- make_lop_prob(pp, Xmat, offset_vec)
    resid <- Ymat[, 1] * fit$p - Ymat[, 2]
    colSums(weights_vec * Xmat * c(resid)) + augmentation_vec
  }

  fit_stage <- function(start, augmentation_vec) {
    gradient_fn <- function(pp) {
      score_vec(pp, X, Y, weights, offset, augmentation_vec)
    }

    if (no.opt) {
      pp <- start
      sol <- NULL
      converged_flag <- TRUE
      termcd_val <- NA_integer_
      message_val <- NA_character_
    } else {
      sol <- nleqslv::nleqslv(start, gradient_fn, control = list(trace = 0))
      pp <- sol$x
      converged_flag <- sol$termcd %in% c(1, 2)
      termcd_val <- sol$termcd
      message_val <- sol$message

      if (!converged_flag) {
        warning(sprintf(
          "nleqslv may not have fully converged: termcd=%s, message=%s",
          termcd_val, message_val
        ))
      }
    }

    grad <- gradient_fn(pp)
    hess <- numDeriv::jacobian(gradient_fn, pp)

    solve_error <- NA_character_
    ihess <- tryCatch(
      solve(hess),
      error = function(e) {
        solve_error <<- conditionMessage(e)
        NULL
      }
    )

    fit <- make_lop_prob(pp, X, offset)
    resid <- Y[, 1] * fit$p - Y[, 2]
    Dlogl <- weights * X * c(resid)

    if (is.null(ihess)) {
      iid <- matrix(NA_real_, nrow = nid, ncol = ncol(X))
      robvar <- matrix(NA_real_, ncol = ncol(X), nrow = ncol(X))
      modelvar <- NULL
      se_robust <- rep(NA_real_, ncol(X))
      se_model <- rep(NA_real_, ncol(X))
      converged_flag <- FALSE
    } else {
      iid <- apply(Dlogl %*% t(ihess), 2, sumstrata, id, max(id) + 1)
      robvar <- crossprod(iid)
      modelvar <- tryCatch(solve(hess), error = function(e) NULL)
      se_robust <- sqrt(diag(robvar))
      se_model <- if (!is.null(modelvar)) sqrt(diag(modelvar)) else rep(NA_real_, length(pp))
    }

    list(
      par = pp,
      coef = pp,
      ploglik = sum(weights * resid^2),
      gradient = grad,
      hessian = hess,
      ihessian = ihess,
      id = id,
      Dlogl = Dlogl,
      iid = iid,
      iid.naive = iid,
      robvar = robvar,
      var = robvar,
      se.robust = se_robust,
      modelvar = modelvar,
      se.model = se_model,
      p = fit$p,
      p0 = fit$p0,
      p1 = fit$p1,
      resid = resid,
      solution = sol,
      augmentation = augmentation_vec,
      converged = converged_flag,
      termcd = termcd_val,
      solver_message = message_val,
      solve_error = solve_error
    )
  }

  fit1 <- fit_stage(beta, augmentation)

  coefI <- fit1$coef
  gradientI <- fit1$gradient
  hessianI <- fit1$hessian
  ihessianI <- fit1$ihessian
  iidI_naive <- fit1$iid

  need_censor_work <- do_censor_fit && (isTRUE(se) || identical(type, "II"))

  MGCiidI <- matrix(0, nrow = nid, ncol = ncol(X))
  MGCiid  <- matrix(0, nrow = nid, ncol = ncol(X))
  augmentationII <- augmentation
  id_out <- id

  if (need_censor_work) {
    ord <- resC$ord

    X_ord <- X[ord, , drop = FALSE]
    Y_ord <- Y[ord, , drop = FALSE]
    exit_ord <- exit[ord]
    offset_ord <- offset[ord]

    fit_ord <- make_lop_prob(fit1$coef, X_ord, offset_ord)
    Yo_ord <- Y_ord[, 1] * fit_ord$p - Y_ord[, 2]

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0, length(xx$strata))
    S0i[xx$jumps + 1]  <- 1 / resC$S0
    S0i2[xx$jumps + 1] <- 1 / (resC$S0^2)

    h <- apply(X_ord * Yo_ord, 2, revcumsumstrata, xx$strata, xx$nstrata)
    btime <- 1 * (exit_ord < time)

    IhdLam0 <- apply(h * S0i2 * btime, 2, cumsumstrata, xx$strata, xx$nstrata)
    U <- matrix(0, nrow(xx$X), ncol(X_ord))
    U[xx$jumps + 1, ] <- (resC$jumptimes < time) * h[xx$jumps + 1, ] / c(resC$S0)
    MGt <- (U - IhdLam0) * c(xx$weights)

    mid <- max(xx$id)
    MGCiidI <- apply(MGt, 2, sumstrata, xx$id, mid + 1)
    MGCiid  <- MGCiidI

    if (type == "II") {
      hYt <- revcumsumstrata(Yo_ord, xx$strata, xx$nstrata)
      IhdLam0_y <- cumsumstrata(hYt * S0i2 * btime, xx$strata, xx$nstrata)

      U_y <- rep(0, length(xx$strata))
      U_y[xx$jumps + 1] <- (resC$jumptimes < time) * hYt[xx$jumps + 1] / c(resC$S0)

      MGt_y <- X_ord * c(U_y - IhdLam0_y) * c(xx$weights)
      MGtiid <- apply(MGt_y, 2, sumstrata, xx$id, mid + 1)

      augmentationII <- augmentation + colSums(MGtiid)

      EXt <- apply(X_ord, 2, revcumsumstrata, xx$strata, xx$nstrata)
      IEXhYtdLam0 <- apply(EXt * c(hYt) * S0i * S0i2 * btime,
                           2, cumsumstrata, xx$strata, xx$nstrata)

      U2 <- matrix(0, nrow(xx$X), ncol(X_ord))
      U2[xx$jumps + 1, ] <- (resC$jumptimes < time) *
        hYt[xx$jumps + 1] * EXt[xx$jumps + 1, ] / c(resC$S0)^2

      MGt2 <- (U2 - IEXhYtdLam0) * c(xx$weights)
      MGCiid2 <- apply(MGt2, 2, sumstrata, xx$id, mid + 1)

      MGCiid <- MGCiid + (MGtiid - MGCiid2)
    }

    id_out <- xx$id
  }

  if (type == "II") {
    fit_final <- fit_stage(fit1$coef, augmentationII)
  } else {
    fit_final <- fit1
  }

  val <- c(list(coef = fit_final$coef), fit_final)
  if (length(val$coef) == length(colnames(X))) {
    names(val$coef) <- colnames(X)
  }

  val <- c(
    val,
    list(
      time = time,
      formula = formula,
      formC = formC,
      exit = exit,
      cens.weights = cens.weights,
      cens.strata = cens.strata,
      cens.nstrata = cens.nstrata,
      model.frame = m,
      n = length(exit),
      nevent = nevent,
      ncluster = nid
    )
  )

  val$call <- cl
  val$MGciid <- MGCiid
  val$id <- id_out
  val$call.id <- call.id
  val$name.id <- name.id
  val$nid <- nid

  val$iidI <- iidI_naive
  val$iid.naive <- fit_final$iid
  val$naive.robvar <- crossprod(val$iid.naive)
  val$se.naive.robust <- sqrt(diag(val$naive.robvar))

  if (isTRUE(se) && do_censor_fit) {
    val$iid  <- fit_final$iid + MGCiid  %*% t(fit_final$ihessian)
    val$iidI <- iidI_naive  + MGCiidI %*% t(ihessianI)
  } else {
    val$iid <- fit_final$iid
  }

  if (is.matrix(val$iid) && length(name.id) == nrow(val$iid)) {
    rownames(val$iid) <- name.id
  }
  if (is.matrix(val$iidI) && length(name.id) == nrow(val$iidI)) {
    rownames(val$iidI) <- name.id
  }
  if (is.matrix(val$iid.naive) && length(name.id) == nrow(val$iid.naive)) {
    rownames(val$iid.naive) <- name.id
  }

  val$robvar <- crossprod(val$iid)
  val$se.robust <- sqrt(diag(val$robvar))
  val$var <- val$robvar
  val$se.coef <- val$se.robust
  val$naive.var <- val$naive.robvar

  val$cause <- cause
  val$cens.code <- cens.code
  val$augmentation <- augmentationII
  val$model <- model[1]
  val$outcome <- outcome[1]
  val$Yipcw <- Yipcw
  val$Causes <- Causes
  val$nevent <- nevent

  val$coefI <- coefI
  val$varI <- crossprod(val$iidI)
  val$se.coefI <- sqrt(diag(val$varI))
  val$gradientI <- gradientI
  val$hessianI <- hessianI
  val$ihessianI <- ihessianI

  val$converged <- fit_final$converged
  val$termcd <- fit_final$termcd
  val$solver_message <- fit_final$solver_message
  val$solve_error <- fit_final$solve_error
  val$n_event1_tau <- sum((status %in% cause) & (exit <= time), na.rm = TRUE)
  val$n_censored_tau <- sum(statusC == 1 & exit <= time, na.rm = TRUE)
  val$n_obs_tau <- sum(obs, na.rm = TRUE)
  val$n_A1 <- sum(X[, ncol(X)] == 1, na.rm = TRUE)
  val$n_A0 <- sum(X[, ncol(X)] == 0, na.rm = TRUE)

  class(val) <- "binreg"
  val
}

#' Calculate fitted percentages from log-odds product regression coefficients
#'
#' Converts regression coefficients from `binregRatioLOP()` into fitted
#' percentages under exposure levels 0 and 1.
#'
#' The coefficient vector is assumed to consist of nuisance coefficients for the
#' covariates `X_L`, followed by a final coefficient representing the log
#' percentage ratio parameter for the binary exposure.
#'
#' @param beta Numeric vector of regression coefficients. Its length must equal
#'   `ncol(X_L) + 1`, where the final element is the exposure log percentage ratio
#'   parameter.
#' @param X_L Numeric matrix of adjustment covariates excluding the binary
#'   exposure. A vector is coerced to a one-column matrix.
#' @param offset Numeric vector of offsets with one value per row of `X_L`.
#' @param tol Numeric tolerance used to switch to a numerically stable expression
#'   when the nuisance linear predictor is close to zero.
#' @param eps Small positive constant used to truncate fitted percentages away
#'   from 0 and 1.
#'
#' @details
#' For each observation, the function returns two fitted percentages:
#' one for exposure level 0 and one for exposure level 1. These are obtained from
#' the nuisance linear predictor and the log percentage ratio parameter under the
#' log-odds product parameterization.
#'
#' @return
#' A numeric matrix with two columns:
#' \describe{
#'   \item{p_0}{Estimated percentage under exposure level 0.}
#'   \item{p_1}{Estimated percentage under exposure level 1.}
#' }
#'
#' @examples
#' ## One binary exposure and no adjustment covariates
#' beta <- c(log(0.75))
#' X_L <- matrix(numeric(0), nrow = 1)
#' calculatePercentageLOP(beta, X_L = X_L, offset = 0)
#'
#' ## One adjustment covariate plus the exposure parameter
#' beta <- c(0.2, log(0.75))
#' X_L <- matrix(c(-1, 1), ncol = 1)
#' calculatePercentageLOP(beta, X_L = X_L, offset = c(0, 0))
#'
#' @name calculatePercentageLOP
#' @export
calculatePercentageLOP <- function(beta, X_L, offset, tol = 1e-8, eps = 1e-10) {
  n <- length(offset)

  if (is.null(dim(X_L))) {
    X_L <- matrix(X_L, nrow = n)
  }

  nL <- ncol(X_L)

  if (length(beta) != nL + 1L) {
    stop("length(beta) must be ncol(X_L) + 1 (last element is theta).")
  }

  if (nL == 0L) {
    phi <- offset
  } else {
    phi <- as.vector(X_L %*% beta[seq_len(nL)] + offset)
  }

  theta <- rep(beta[nL + 1L], n)
  expphi <- exp(phi)
  exptheta <- exp(theta)

  p_0 <- p_1 <- rep(NA_real_, n)
  small <- abs(phi) < tol

  if (any(small)) {
    p_0[small] <- 1 / (1 + exptheta[small])
    p_1[small] <- exptheta[small] * p_0[small]
  }

  if (any(!small)) {
    idx <- which(!small)

    denomi_1 <- -(exptheta[idx] + 1) * expphi[idx]
    inside <- exp(2 * phi[idx]) * (exptheta[idx] + 1)^2 +
      4 * exp(theta[idx] + phi[idx]) * (1 - expphi[idx])
    inside <- pmax(inside, 0)

    denomi_2 <- sqrt(inside)
    denomi <- denomi_1 + denomi_2
    numera <- 2 * exptheta[idx] * (1 - expphi[idx])

    p_0[idx] <- denomi / numera
    p_1[idx] <- exptheta[idx] * p_0[idx]
  }

  p_0 <- pmin(pmax(p_0, eps), 1 - eps)
  p_1 <- pmin(pmax(p_1, eps), 1 - eps)

  cbind(p_0, p_1)
}

#' @title Fit coherent regression models of CIFs using polytomous log odds products
#'
#' @description
#' `polyreg()` fits regression models of CIFs, targeting familiar effect measures
#' (risk ratios, odds ratios and subdistribution hazard ratios).
#' Modeling the nuisance structure using polytomous log odds products ensures that
#' the sum of cause-specific CIFs does not exceed one, and enables coherent modelling
#' of the multiplicative effects.
#'
#' This function follows a familiar formula–data workflow: the outcome and
#' covariates other than the exposure are specified through a formula in `nuisance.model`
#' (with `Event()` or `survival::Surv()` on the LHS), and the exposure of interest
#' is given by a separate variable name in `exposure`. The fitted object contains
#' tidy summaries of exposure effects (point estimates, SEs, CIs, and p-values)
#' and can be summarised with `summary.polyreg()` or formatted with external tools
#' such as `modelsummary::modelsummary()`.
#'
#' @param nuisance.model A `formula` describing the outcome and
#'   nuisance covariates, excluding the exposure of interest.
#'   The LHS must be `Event(time, status)` or `survival::Surv(time, status)`.
#' @param exposure A character string giving the name of the categorical exposure
#'   variable in `data`.
#' @param strata Optional character string with the name of the stratification
#'   variable used to adjust for dependent censoring (default `NULL`).
#' @param data A data frame containing the outcome, exposure and nuisance
#'   covariates referenced by `nuisance.model`.
#' @param subset.condition Optional character string giving a logical condition to subset
#' `data` (default `NULL`).
#' @param na.action A function specifying the action to take on missing values (default `na.omit`).
#' @param code.event1 Integer code of the event of interest (default `1`).
#' @param code.event2 Integer code of the competing event (default `2`).
#' @param code.censoring  Integer code of censoring (default `0`).
#' @param code.exposure.ref Integer code identifying the reference exposure
#'   category (default `0`).
#' @param effect.measure1 Character string specifying the effect measure for the
#'   primary event. Supported values are `"RR"`, `"OR"` and
#'   `"SHR"`.
#' @param effect.measure2 Character string specifying the effect measure for the
#'   competing event. Supported values are `"RR"`, `"OR"` and
#'   `"SHR"`.
#' @param time.point
#'   Numeric time point at which the exposure effect is evaluated for
#'   time-point models. Required for `"competing-risk"` and `"survival"`
#'   outcomes.
#' @param outcome.type Character string selecting the outcome type. Valid values are
#'   `"competing-risk"`, `"survival"`, `"binomial"`, `"proportional-survival"`,
#'   and `"proportional-competing-risk"`. The default is `"competing-risk"`.
#'   If explicitly set to `NULL`, `polyreg()` attempts to infer the outcome type from the data: if the
#'   event variable has more than two distinct levels, `"competing-risk"`
#'   is assumed; otherwise, `"survival"` is assumed. Abbreviations such as
#'   `"S"` or `"C"` are accepted; mixed or ambiguous inputs trigger
#'   automatic detection from the event coding in `data`.
#' @param conf.level Confidence level for Wald-type intervals (default `0.95`).
#' @param report.nuisance.parameter Logical; if `TRUE`, the returned object
#'   includes estimates of the nuisance model parameters (default `FALSE`).
#' @param report.optim.convergence Logical; if `TRUE`, optimization
#'   convergence summaries are returned (default `FALSE`).
#' @param report.sandwich.conf Logical or `NULL`. When `TRUE`,
#' confidence intervals based on sandwich variance are computed.
#' When `FALSE`, they are omitted (default `TRUE`).
#' This confidence interval is default for time-point models
#' (`"outcome.type=competing-risk"`, `"survival"` or `"binomial"`) and
#' is not available otherwise.
#' @param report.boot.conf Logical or `NULL`. When `TRUE`, bootstrap
#' confidence intervals are computed. When `FALSE`, they are omitted.
#' If `NULL`, the function chooses based on `outcome.type` (default `NULL`).
#' This confidence interval is default for proportional models
#' (`outcome.type="proportional-competing-risk"` or `"proportional-survival"`).
#' @param boot.bca Logical indicating the bootstrap confidence interval method.
#'   Use `TRUE` for bias-corrected and accelerated intervals or `FALSE`
#'   for the normal approximation (default `FALSE`).
#' @param boot.multiplier Character string specifying the wild bootstrap weight distribution.
#' One of `"rademacher"`, `"mammen"` or `"gaussian"` (default `"rademacher"`).
#' @param boot.replications Integer giving the number of bootstrap replications
#'   (default `200`).
#' @param boot.seed Numeric seed used for resampling of bootstrap.
#' @param nleqslv.method Character string specifying the solver used in
#'   \pkg{nleqslv()}. Available choices are `"Broyden"` and `"Newton"`.
#' @param optim.parameter1 Numeric tolerance for convergence of the outer loop
#'    (default `1e-6`).
#' @param optim.parameter2 Numeric tolerance for convergence of the inner loop
#'    (default `1e-6`).
#' @param optim.parameter3 Numeric constraint on the absolute value of
#'   parameters (default `100`).
#' @param optim.parameter4 Integer maximum number of outer loop iterations
#'   (default `50`).
#' @param optim.parameter5 Integer maximum number of `nleqslv`
#'   iterations per outer iteration (default `50`).
#' @param optim.parameter6 Integer maximum number of iterations for the
#'   Levenberg-Marquardt routine (default `50`).
#' @param optim.parameter7 Numeric convergence tolerance for the
#'   Levenberg-Marquardt routine (default `1e-10`).
#' @param optim.parameter8 Numeric tolerance for updating the Hessian in the
#'   Levenberg-Marquardt routine (default `1e-6`).
#' @param optim.parameter9 Numeric starting value for the Levenberg-Marquardt
#'   damping parameter lambda (default `1e-6`).
#' @param optim.parameter10 Numeric upper bound for lambda in the
#'   Levenberg-Marquardt routine (default `40`).
#' @param optim.parameter11 Numeric lower bound for lambda in the
#'   Levenberg-Marquardt routine (default `0.025`).
#' @param optim.parameter12 Numeric multiplicative increment applied to lambda
#'   when the Levenberg-Marquardt step is successful (default `2`).
#' @param optim.parameter13 Numeric multiplicative decrement applied to lambda
#'   when the Levenberg-Marquardt step is unsuccessful (default `0.5`).
#' @param data.initial.values Optional data frame providing starting values for
#'   the optimization (default `NULL`).
#' @param normalize.covariate Logical indicating whether covariates should
#'   be centered and scaled prior to optimization (default `TRUE`).
#' @param terminate.time.point Logical indicating whether time points
#'   that contribute estimation are terminated by min of max follow-up times
#'   of each exposure level (default `TRUE`).
#' @param prob.bound Numeric lower bound used to internally truncate probabilities away
#'   from 0 and 1 (default `1e-5`).
#'
#' @details
#'
#' ### Overview
#' `polyreg()` implements **log odds product modeling** for CIFs at user-specified
#' time points, focusing on multiplicative effects of a categorical exposure, or
#' constant effects over time like Cox regression and Fine-Gray models. It estimates
#' multiplicative effects such as **risk ratios**, **odds ratios**, or
#' **subdistribution hazard ratios**, while ensuring that the probabilities across
#' competing events sum to one. This is achieved through
#' **reparameterization using polytomous log odds products**, which fits so-called
#' effect-measure models and nuisance models on multiple competing events
#' simultaneously. Additionally, `polyreg()` supports direct binomial regression
#' for survival outcomes and the Richardson model for binomial outcomes,
#' both of which use log odds products.
#'
#' ### Key arguments
#' -   `nuisance.model`: a formula with `Event()` or `survivai::Surv()`
#' describing the outcome and nuisance covariates, excluding the exposure of interest.
#' -   `exposure`: name of the categorical exposure variable
#' -   `effect.measure1` and `effect.measure2`: the effect measures
#' for event1 and event2 (`"RR"`, `"OR"` or `"SHR"`).
#' -   `outcome.type`: type of the outcome variable (`"competing-risk"`, `"survival"`,
#' `"binomial"`, `"proportional-survival"` or `"proportional-competing-risk"`).
#' -   `time.point`: time point(s) at which the exposure effect is evaluated.
#' Required for `"competing-risk"` and `"survival"` outcomes.
#' -   `strata`: name of the stratification variable used for IPCW adjustment for dependent censoring.
#'
#' ### Outcome type and event status coding
#'
#' The `outcome.type` argument must be set to:
#'
#' - Effects on cumulative incidence probabilities at a specific time:
#'   `"competing-risk"`.
#' - Effects on a risk at a specific time: `"survival"`.
#' - Common effects on cumulative incidence probabilities over time:
#'   `"proportional-competing-risk"`.
#' - Common effects on a risk over time: `"proportional-survival"`.
#' - Effects on a risk of a binomial outcome: `"binomial"`.
#'
#' | Setting                         | Codes                                        | Meaning                                       |
#' |---------------------------------|----------------------------------------------|-----------------------------------------------|
#' | competing-risk                  | `code.event1`, `code.event2`, `code.censoring` | event of interest / competing event / censoring |
#' | competing-risk (default)        | `code.event1 = 1`, `code.event2 = 2`, `code.censoring = 0` | event of interest / competing event / censoring |
#' | survival                        | `code.event1`, `code.censoring`             | event / censoring                             |
#' | survival (default)             | `code.event1 = 1`, `code.censoring = 0`     | event / censoring                             |
#' | survival (ADaM-ADTTE)           | `code.event1 = 0`, `code.censoring = 1`     | set to match ADaM convention                  |
#' | proportional-survival           | `code.event1`, `code.censoring`             | event / censoring                             |
#' | proportional-survival (default) | `code.event1 = 1`, `code.censoring = 0`     | event / censoring                             |
#' | proportional-survival (ADaM)    | `code.event1 = 0`, `code.censoring = 1`     | set to match ADaM convention                  |
#' | proportional-competing-risk     | `code.event1`, `code.event2`, `code.censoring` | event of interest / competing event / censoring |
#' | proportional-competing-risk (default) | `code.event1 = 1`, `code.event2 = 2`, `code.censoring = 0` | event of interest / competing event / censoring |
#'
#' ### Effect measures for categorical exposure
#'
#' Choose the effect scale for event 1 and (optionally) event 2:
#'
#' | Argument          | Applies to       | Choices              | Default |
#' |-------------------|------------------|----------------------|---------|
#' | `effect.measure1` | event of interest | `"RR"`, `"OR"`, `"SHR"` | `"RR"`  |
#' | `effect.measure2` | competing event   | `"RR"`, `"OR"`, `"SHR"` | `"RR"`  |
#'
#' - `RR`: risk ratio at `time.point` or common over time.
#' - `OR`: odds ratio at `time.point` or common over time.
#' - `SHR`: subdistribution hazard ratio or common over time.
#'
#' ### Inference and intervals (advanced)
#'
#' | Argument             | Meaning                                         | Default      |
#' |----------------------|-------------------------------------------------|--------------|
#' | `conf.level`         | Wald-type CI level                              | `0.95`       |
#' | `report.sandwich.conf` | Sandwich variance CIs                        | `TRUE`       |
#' | `report.boot.conf`   | Bootstrap CIs (used by `"proportional-*"` types) | `NULL`       |
#' | `boot.bca`           | Use BCa intervals (else normal approximation)   | `FALSE`      |
#' | `boot.multiplier`    | Method for wild bootstrap                       | `"rademacher"` |
#' | `boot.replications`  | Bootstrap replications                          | `200`        |
#' | `boot.seed`          | Seed for resampling                             | `46`         |
#'
#' ### Optimization & solver controls (advanced)
#'
#' `polyreg()` solves estimating equations with optional inner routines.
#'
#' | Argument          | Role                            | Default   |
#' |-------------------|---------------------------------|-----------|
#' | `nleqslv.method`  | Root solver                     | `"Newton"` |
#' | `optim.parameter1`, `optim.parameter2` | Outer / inner convergence tolerances | `1e-6`, `1e-6` |
#' | `optim.parameter3`| Parameter absolute bound        | `100`     |
#' | `optim.parameter4`| Max outer iterations            | `50`      |
#' | `optim.parameter5`| Max `nleqslv` iterations per outer | `50` |
#' | `optim.parameter6:13` | Levenberg–Marquardt controls (iterations, tolerances, lambda) | see defaults |
#'
#' ### Data handling and stability
#'
#' | Argument             | Meaning                                            | Default          |
#' |----------------------|----------------------------------------------------|------------------|
#' | `subset.condition`   | Expression (as character) to subset `data`         | `NULL`           |
#' | `na.action`          | NA handling function                               | `stats::na.omit` |
#' | `normalize.covariate`| Center/scale nuisance covariates                   | `TRUE`           |
#' | `terminate.time.point` | Truncate support by exposure-wise follow-up maxima | `TRUE`        |
#' | `prob.bound`         | Truncate probabilities away from 0/1 (numerical guard) | `1e-5`      |
#' | `data.initial.values`| Optional starting values data frame                | `NULL`           |
#'
#' ### Downstream use
#'
#' `polyreg()` returns an object of class `"polyreg"` that contains
#' regression coefficients (`coef`), variance-covariance matrix (`vcov`)
#' and a list of event-wise \emph{tidy} and \emph{glance} tables (`summary`).
#' Users should typically access results via the S3 methods:
#'
#' - `coef()` — extract regression coefficients.
#' - `vcov()` — extract the variance–covariance matrix
#'   (sandwich or bootstrap, depending on `outcome.type` and the
#'   `report.*` arguments).
#' - `nobs()` — number of observations used in the fit.
#' - `summary()` — print an event-wise, modelsummary-like table of estimates,
#'   CIs and p-values, and return the underlying list of tidy/glance tables invisibly.
#'
#' For backward compatibility, components named `coefficient` and `cov`
#' may also be present and mirror `coef` and `vcov`, respectively.
#' The `summary` component can be passed to external functions such as
#' `modelsummary()` for further formatting, if desired.
#'
#' ### Reproducibility and conventions
#'
#' - If convergence warnings appear, relax/tighten tolerances or cap the parameter
#'   bound (`optim.parameter1`–`3`) and inspect the output with
#'   `report.optim.convergence = TRUE`.
#' - If necessary, modify other `optim.parameter`, provide user-specified
#'   initial values, or reduce the number of nuisance parameters (e.g., provide
#'   a small set of time points contributing to estimation when using
#'   `"proportional-survival"` or `"proportional-competing-risk"`).
#' - Set `boot.seed` for reproducible bootstrap results.
#' - Match CDISC ADaM conventions via `code.event1 = 0`, `code.censoring = 1`
#'   (and, if applicable, `code.event2` for competing events).
#'
#' @importFrom nleqslv nleqslv
#' @importFrom boot boot boot.ci
#' @importFrom Rcpp sourceCpp
#' @importFrom stats IQR as.formula binomial coef glm mad median
#' @importFrom stats model.extract model.frame model.matrix update na.omit na.pass
#' @importFrom stats pnorm qnorm rbinom reformulate rexp sd setNames terms time var
#' @importFrom stats runif rnorm
#'
#' @useDynLib cifmodeling, .registration = TRUE
#'
#' @return
#' A list of class `"polyreg"` containing the fitted exposure effects and
#' supporting results. Key components and methods include:
#'
#' - `coef`: regression coefficients on the chosen effect-measure scale
#' - `vcov`: variance–covariance matrix of the regression coefficients
#' - `diagnostic.statistics`: a data frame with inverse probability weights,
#'   influence function contributions, and predicted potential outcomes
#' - `summary`: event-wise tidy/glance summaries used by
#'   `summary.polyreg()` or `modelsummary::modelsummary()`
#' - additional elements storing convergence information and internal
#'   tuning parameters.
#'
#' Standard S3 methods are available: `coef.polyreg()`, `vcov.polyreg()`,
#' `nobs.polyreg()`, and `summary.polyreg()`.
#'
#' @examples
#' data(diabetes.complications)
#' output <- polyreg(
#'   nuisance.model = Event(t, epsilon) ~ +1,
#'   exposure = "fruitq1",
#'   data = diabetes.complications,
#'   effect.measure1 = "RR",
#'   effect.measure2 = "RR",
#'   time.point = 8,
#'   outcome.type = "competing-risk"
#' )
#'
#' coef(fit)
#' vcov(fit)
#' nobs(fit)
#' summary(fit)
#'
#' @name polyreg
#' @section Lifecycle:
#' \lifecycle{experimental}
#'
#' @seealso [cifcurve()] for KM/AJ estimators; [cifplot()] for display of a CIF; [cifpanel()] for display of multiple CIFs; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
polyreg <- function(
    nuisance.model,
    exposure,
    strata = NULL,
    data,
    subset.condition = NULL,
    na.action = na.omit,
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.exposure.ref = 0,
    effect.measure1 = "RR",
    effect.measure2 = "RR",
    time.point = NULL,
    outcome.type = "competing-risk",
    conf.level = 0.95,
    report.nuisance.parameter = FALSE,
    report.optim.convergence = FALSE,
    report.sandwich.conf = TRUE,
    report.boot.conf = NULL,
    boot.bca = FALSE,
    boot.multiplier = "rademacher",
    boot.replications = 200,
    boot.seed = 46,
    nleqslv.method = "Newton",
    optim.parameter1 = 1e-6,
    optim.parameter2 = 1e-6,
    optim.parameter3 = 100,
    optim.parameter4 = 50,
    optim.parameter5 = 50,
    optim.parameter6 = 50,
    optim.parameter7 = 1e-10,
    optim.parameter8 = 1e-6,
    optim.parameter9 = 1e-6,
    optim.parameter10 = 40,
    optim.parameter11 = 0.025,
    optim.parameter12 = 2,
    optim.parameter13 = 0.5,
    data.initial.values = NULL,
    normalize.covariate = TRUE,
    terminate.time.point = TRUE,
    prob.bound = 1e-7
) {

  #######################################################################################################
  # 1. Pre-processing (function: checkSpell, checkInput, reg_normalize_covariate, sortByCovariate)
  #######################################################################################################
  computation.time0 <- proc.time()
  outcome.type  <- util_check_outcome_type(outcome.type, formula=formula, data=data)
  ce <- reg_check_effect.measure(effect.measure1, effect.measure2)
  ci <- reg_check_input(data, nuisance.model, exposure, code.event1, code.event2, code.censoring, code.exposure.ref, outcome.type, conf.level, report.sandwich.conf, report.boot.conf, nleqslv.method, normalize.covariate)
  normalize.covariate <- ci$normalize.covariate
  report.sandwich.conf <- ci$report.sandwich.conf
  report.boot.conf <- ci$report.boot.conf

  data <- createAnalysisDataset(formula=nuisance.model, data=data, other.variables.analyzed=c(exposure, strata), subset.condition=subset.condition, na.action=na.action)
  out_normalizeCovariate <- reg_normalize_covariate(nuisance.model, data, normalize.covariate, outcome.type, ci$out_readExposureDesign$exposure.levels)
  normalized_data <- out_normalizeCovariate$normalized_data
  tp <- reg_read_time.point(nuisance.model, normalized_data, ci$out_readExposureDesign$x_a, outcome.type, code.censoring, terminate.time.point, time.point)
  index.vector <- reg_index_for_parameter(NA, ci$x_l, ci$x_a, length(tp))

  estimand <- list(
    effect.measure1      = ce$effect.measure1,
    effect.measure2      = ce$effect.measure2,
    time.point           = tp,
    code.event1          = code.event1,
    code.event2          = code.event2,
    code.censoring       = code.censoring,
    code.exposure.ref    = code.exposure.ref,
    exposure.levels      = ci$out_readExposureDesign$exposure.levels,
    index.vector         = index.vector
  )

  boot.method <- list(
    report.sandwich.conf = report.sandwich.conf,
    report.boot.conf     = report.boot.conf,
    boot.bca             = boot.bca,
    boot.multiplier      = boot.multiplier,
    boot.replications    = boot.replications,
    boot.seed            = boot.seed
  )

  optim.method <- list(
    nleqslv.method    = nleqslv.method,
    optim.parameter1  = optim.parameter1,
    optim.parameter2  = optim.parameter2,
    optim.parameter3  = optim.parameter3,
    optim.parameter4  = optim.parameter4,
    optim.parameter5  = optim.parameter5,
    optim.parameter6  = optim.parameter6,
    optim.parameter7  = optim.parameter7,
    optim.parameter8  = optim.parameter8,
    optim.parameter9  = optim.parameter9,
    optim.parameter10 = optim.parameter10,
    optim.parameter11 = optim.parameter11,
    optim.parameter12 = optim.parameter12,
    optim.parameter13 = optim.parameter13
  )

  #######################################################################################################
  # 2. Pre-processing and Calculating initial values alpha_beta_0 (function: calculateInitialValues)
  #######################################################################################################
  if (outcome.type == "competing-risk" || outcome.type == "survival" || outcome.type == "binomial") {
    alpha_beta_0 <- getInitialValues(
      formula = nuisance.model,
      data = normalized_data,
      exposure = exposure,
      data.initial.values = data.initial.values,
      estimand = estimand,
      specific.time = estimand$time.point,
      outcome.type = outcome.type,
      prob.bound = prob.bound
    )
  } else if (outcome.type == "proportional-survival" || outcome.type == "proportional-competing-risk") {
    alpha_beta_0 <- getInitialValuesProportional(
      formula = nuisance.model,
      data = normalized_data,
      outcome.type = outcome.type,
      exposure = exposure,
      estimand = estimand,
      data.initial.values = data.initial.values,
      prob.bound = prob.bound,
      out_normalizeCovariate = out_normalizeCovariate
    )
  }

  #######################################################################################################
  # 3. Calculating IPCW (function: calculateIPCW, calculateIPCWMatrix)
  #######################################################################################################
  if (outcome.type == "competing-risk" || outcome.type == "survival") {
    ip.weight.matrix <- calculateIPCW(nuisance.model, normalized_data, code.censoring, strata, estimand$time.point)
  } else if (outcome.type == "binomial") {
    ip.weight.matrix <- matrix(1,nrow(normalized_data),1)
  } else if (outcome.type == "proportional-survival" || outcome.type == "proportional-competing-risk") {
    ip.weight.matrix <- calculateIPCWMatrix(nuisance.model, normalized_data, code.censoring, strata, estimand, out_normalizeCovariate)
  }

  #######################################################################################################
  # 4. Parameter estimation (functions: estimating_equation_ipcw, _survival, _proportional)
  #######################################################################################################
  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
      call_and_capture <- function(fun, ...) {
      out_ipcw <<- do.call(fun, list(...))
      out_ipcw$ret
      }
      estimating_equation_i <- function(p) call_and_capture(
        estimating_equation_ipcw,
        formula = nuisance.model, data = normalized_data, exposure = exposure, outcome.type=outcome.type,
        ip.weight.matrix = ip.weight.matrix, alpha_beta = p, estimand = estimand,
        optim.method = optim.method, prob.bound = prob.bound,
        initial.CIFs = initial.CIFs
      )
    setInitialCIFs <- function(new.CIFs) initial.CIFs <<- new.CIFs
    getResults     <- function() out_ipcw
      list(
      estimating_equation_i = estimating_equation_i,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
    )
  }

  assessConvergence <- function(new_params, current_params, current_obj_value, optim.parameter1, optim.parameter2, optim.parameter3) {
    assessRelativeDifference <- function(new, old) {
      max(abs(new - old) / pmax(1, abs(old)))
    }
    is_stalled <- function(x, stall_patience = 3, stall_eps = 1e-3) {
      n <- length(x)
      if (n < stall_patience) return(FALSE)
      recent <- x[(n - stall_patience + 1):n]
      rng <- range(recent)
      rel_diff <- (diff(rng) / max(1e-12, mean(recent)))
      rel_diff <= stall_eps
    }
    if (any(abs(new_params) > optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max.absolute.difference <- max(param_diff)
    relative.difference <- assessRelativeDifference(new_params, current_params)
    obj_value <- drop(crossprod(obj$estimating_equation_i(new_params)))

    converged <- (relative.difference <= optim.parameter1) || (obj_value <= optim.parameter2) || is_stalled(c(current_obj_value, obj_value))
    criteria1 <- (relative.difference <= optim.parameter1)
    criteria2 <- (obj_value <= optim.parameter2)
    criteria3 <- is_stalled(c(current_obj_value, obj_value))
    converged  <- (criteria1 || criteria2 || criteria3)
    converged.by <- if (!converged) NA_character_
    else if (criteria1) "Converged in relative difference"
    else if (criteria2) "Converged in objective function"
    else "Stalled"

    list(converged = converged, converged.by=converged.by, relative.difference = relative.difference, max.absolute.difference = max.absolute.difference, obj_value = obj_value)
  }

  obj <- makeObjectiveFunction()
  iteration <- 0L
  max.absolute.difference <- Inf
  out_nleqslv <- NULL
  current_params <- alpha_beta_0
  current_obj_value <- numeric(0)
  trace_df  <- NULL
  store_params <- TRUE

  while ((iteration < optim.parameter4) & (max.absolute.difference > optim.parameter1)) {
    iteration <- iteration + 1
    prev_params <- current_params

    out_nleqslv <- nleqslv(
      prev_params,
      obj$estimating_equation_i,
      method  = reg_choose_nleqslv_method(nleqslv.method),
      control = list(maxit = optim.parameter5, allowSingular = FALSE)
    )
    new_params <- out_nleqslv$x
    current_obj_value <- drop(crossprod(obj$estimating_equation_i(new_params)))

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    ac <- assessConvergence(new_params, prev_params, current_obj_value, optim.parameter1, optim.parameter2, optim.parameter3)

    nleqslv.info <- extractOptimizationInfo(out_nleqslv, nleqslv.method)
    computation.time.second <- as.numeric((proc.time() - computation.time0)[3])

    trace_df <- append_trace(
      trace_df,
      iteration = iteration,
      computation.time.second = computation.time.second,
      nleqslv.method = nleqslv.method,
      nleqslv.info = nleqslv.info,
      objective.function = ac$obj_value,
      relative.difference = ac$relative.difference,
      max.absolute.difference = ac$max.absolute.difference,
      converged.by = if (ac$converged) ac$converged.by else FALSE,
      coefficient = if (store_params) new_params else NULL
    )

    current_params <- new_params
    converged.by <- ac$converged.by
    objective.function = ac$obj_value
    max.absolute.difference <- ac$max.absolute.difference
    relative.difference <- ac$relative.difference
    if (ac$converged) break
  }
  out_getResults <- obj$getResults()

  #######################################################################################################
  # 5. Calculating variance (functions: calculateCov, calculateCovSurvival)
  #######################################################################################################

  report_var <- FALSE
  if (!isTRUE(report.sandwich.conf) && !isTRUE(report.boot.conf)) {
    report_var <- FALSE
  } else if (isTRUE(report.sandwich.conf) && !isTRUE(report.boot.conf)) {
    report_var <- TRUE
  } else if (isTRUE(report.boot.conf)) {
    if (outcome.type %in% c("proportional-survival","proportional-competing-risk")) {
      report_var <- FALSE
    } else if (outcome.type %in% c("survival","binomial","competing-risk")) {
      report_var <- TRUE
    } else {
      stop(sprintf("Unsupported outcome.type for boot rule: %s", outcome.type))
    }
  }

  out_calculateCov <- switch(
    outcome.type,
    "competing-risk"   = calculateCov(out_getResults, estimand, boot.method, prob.bound),
    "survival"         = calculateCovSurvival(out_getResults, estimand, boot.method, prob.bound),
    "binomial"         = calculateCovSurvival(out_getResults, estimand, boot.method, prob.bound),
    "proportional-survival"     = NULL,
    "proportional-competing-risk"= NULL,
    stop(sprintf("Unsupported outcome.type for covariance: %s", outcome.type))
  )

  if (!report_var) {
    out_calculateCov <- NULL
  }

  out_normalizeEstimate <- reg_normalize_estimate(
    outcome.type               = outcome.type,
    normalize.covariate        = normalize.covariate,
    current_params             = current_params,
    out_getResults             = out_getResults,
    estimand                   = estimand,
    prob.bound                 = prob.bound,
    out_normalizeCovariate     = out_normalizeCovariate,
    out_calculateCov           = out_calculateCov
  )

  alpha_beta_estimated <- out_normalizeEstimate$alpha_beta_estimated
  cov_estimated        <- out_normalizeEstimate$cov_estimated
  cov_bootstrap        <- out_normalizeEstimate$cov_bootstrap

  #######################################################################################################
  # 6. Calculating bootstrap confidence interval (functions: boot, solveEstimatingEquation)
  #######################################################################################################
  if (isTRUE(report.boot.conf) && (outcome.type=="proportional-survival" || outcome.type=="proportional-competing-risk")) {
    set.seed(boot.seed)
    boot.coef     <- rep(NA,2)
    boot.coef_se  <- rep(NA,2)
    boot.p_value  <- rep(NA,2)
    boot.conf_low <- rep(NA,2)
    boot.conf_high<- rep(NA,2)

    if (outcome.type=="proportional-competing-risk") {
      index_coef    <- c(length(estimand$time.point) + 1, 2*length(estimand$time.point) + 2)
    } else if (outcome.type=="proportional-survival") {
      index_coef    <- c(length(estimand$time.point) + 1)
    }

    boot_function <- function(data, indices) {
      res <- try({
        coef <- solveEstimatingEquation(nuisance.model=nuisance.model, exposure=exposure, strata=strata, normalized_data = data[indices, , drop = FALSE], outcome.type=outcome.type, estimand=estimand, optim.method=optim.method, out_normalizeCovariate=out_normalizeCovariate, prob.bound=prob.bound, alpha_beta_0=alpha_beta_0)
      }, silent = TRUE)
      if (inherits(res, "try-error")) NA_real_ else as.numeric(res)
    }

    out_boot <- boot(normalized_data, boot_function, R = boot.replications)
    for (j in index_coef) {
      if (isTRUE(boot.bca)) {
        out_boot.ci       <- boot.ci(out_boot, conf = conf.level, index = index_coef[j], type = c("norm", "bca"))
        boot.coef[j]      <- (out_boot.ci$normal[2] + out_boot.ci$normal[3])/2
        ci_range          <- out_boot.ci$normal[3] - out_boot.ci$normal[2]
        boot.coef_se[j]   <- ci_range/2/qnorm(1 - (1-conf.level)/2)
        boot.p_value[j]   <- 2 * (1 - pnorm(abs(boot.coef[j]) / boot.coef_se[j]))
        boot.conf_low[j]  <- out_boot.ci$bca[4]
        boot.conf_high[j] <- out_boot.ci$bca[5]
      } else {
        ok      <- apply(out_boot$t, 1, function(x) all(is.finite(x)))
        t_ok    <- out_boot$t[ok, , drop = FALSE]
        mean_ok <- colMeans(t_ok)
        sd_ok   <- apply(t_ok, 2, sd)
        ci_ok   <- cbind(
          lower = mean_ok - qnorm(0.975) * sd_ok,
          upper = mean_ok + qnorm(0.975) * sd_ok
        )
        boot.coef[j]      <- mean_ok[j]
        boot.coef_se[j]   <- sd_ok[j]
        boot.p_value[j]   <- 2 * (1 - pnorm(abs(boot.coef[j]) / boot.coef_se[j]))
        boot.conf_low[j]  <- ci_ok[j,1]
        boot.conf_high[j] <- ci_ok[j,2]
      }
    }
    out_bootstrap <- list(boot.coef=boot.coef, boot.coef_se=boot.coef_se, boot.p_value=boot.p_value, boot.conf_low=boot.conf_low, boot.conf_high=boot.conf_high)
  } else {out_bootstrap <- NULL}

  #######################################################################################################
  # 7. Output (functions: reportSurvival, reportCOMPETING-RISK, reportPrediction)
  #######################################################################################################
  out_summary <- reportEffects(
    outcome.type, report.nuisance.parameter, report.optim.convergence, report.sandwich.conf, report.boot.conf,
    nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated, cov_bootstrap,
    out_bootstrap, out_getResults, iteration, converged.by, objective.function, max.absolute.difference, relative.difference,
    out_nleqslv, conf.level, optim.method$nleqslv.method
  )
  if (outcome.type == "competing-risk" || outcome.type == "survival" || outcome.type == "binomial") {
    data$influence.function <- out_calculateCov$influence.function
    data$ip.weight <- out_getResults$ip.weight
    data$potential.CIFs <- out_getResults$potential.CIFs
  }
  out_data <- data

  out <- list(
    coef               = alpha_beta_estimated,
    vcov               = cov_estimated,
    vcov_bootstrap     = cov_bootstrap,
    bootstrap          = out_bootstrap,
    summary            = out_summary,
    diagnostics        = out_data,
    outcome.type       = outcome.type,
    exposure           = exposure,
    conf.level         = conf.level,
    estimand           = estimand,
    boot.method        = boot.method,
    optim.method       = optim.method,
    optim.inofo        = trace_df,
    call               = match.call()
  )
  class(out) <- "polyreg"
  return(out)
}

#' @title Fits regression models of cumulative incidence functions based on polytomous
#' log-odds products and stratified IPCW estimator
#' @description The direct polytomous regression enables coherent modeling and
#' estimation of a variety of multiplicative effects of a categorical exposure under
#' several outcome types, including competing risks, survival and binomial outcomes.
#' @param nuisance.model A \code{\link[stats]{formula}} describing the outcome and
#'   nuisance covariates, excluding the exposure of interest.
#' @param exposure A character string giving the name of the categorical exposure
#'   variable in \code{data}.
#' @param strata Optional character string with the name of the stratification
#'   variable used to adjust for dependent censoring. Defaults to \code{NULL}.
#' @param data A data frame containing the outcome, exposure and nuisance
#'   covariates referenced by \code{nuisance.model}.
#' @param subset.condition Optional expression (as a character string) defining a
#'   subset of \code{data} to analyse. Defaults to \code{NULL}.
#' @param na.action A function specifying the action to take on missing values.
#'   The default is \code{\link[stats]{na.omit}}.
#' @param code.event1 Integer code corresponding to the first event of interest.
#'   Defaults to \code{1}.
#' @param code.event2 Integer code corresponding to the competing event. Defaults
#'   to \code{2}.
#' @param code.censoring Integer code representing censoring. Defaults to
#'   \code{0}.
#' @param code.exposure.ref Integer code identifying the reference exposure
#'   category. Defaults to \code{0}.
#' @param effect.measure1 Character string specifying the effect measure for the
#'   primary event. Supported values are \code{"RR"}, \code{"OR"} and
#'   \code{"SHR"}.
#' @param effect.measure2 Character string specifying the effect measure for the
#'   competing event. Supported values are \code{"RR"}, \code{"OR"} and
#'   \code{"SHR"}.
#' @param time.point Numeric time point at which the exposure effect is
#'   evaluated. Required for survival and competing risk analyses.
#' @param outcome.type Character string selecting the outcome type. Valid values
#'   are \code{"COMPETING-RISK"}, \code{"SURVIVAL"}, \code{"BINOMIAL"},
#'   \code{"PROPORTIONAL"} and \code{"POLY-PROPORTIONAL"}. Defaults to
#'   \code{"COMPETING-RISK"}.
#' If \code{NULL} (default), the function automatically infers the outcome type
#' from the data: if the event variable has more than two unique levels,
#' \code{"COMPETING-RISK"} is assumed; otherwise, \code{"SURVIVAL"} is used.
#' You can also use abbreviations such as \code{"S"} or \code{"C"}.
#' Mixed or ambiguous inputs (e.g., \code{c("S", "C")}) trigger automatic
#' detection based on the event coding in \code{data}.
#' @param conf.level Confidence level for Wald-type intervals. Defaults to
#'   \code{0.95}.
#' @param report.nuisance.parameter Logical; if \code{TRUE}, the returned object
#'   includes estimates of the nuisance model parameters. Defaults to
#'   \code{FALSE}.
#' @param report.optim.convergence Logical; if \code{TRUE}, optimisation
#'   convergence summaries are returned. Defaults to \code{FALSE}.
#' @param report.sandwich.conf Logical or \code{NULL}. When \code{TRUE}, confidence
#'   intervals based on sandwich variance are computed. When \code{FALSE}, they are
#'   omitted. Defaults to \code{TRUE}.
#' @param report.boot.conf Logical or \code{NULL}. When \code{TRUE}, bootstrap
#'   confidence intervals are computed. When \code{FALSE}, they are omitted. If
#'   \code{NULL}, the function chooses based on \code{outcome.type}.
#' @param boot.bca Logical indicating the bootstrap confidence interval method.
#'   Use \code{TRUE} for bias-corrected and accelerated intervals or \code{FALSE}
#'   for the normal approximation. Defaults to \code{FALSE}.
#' @param boot.parameter1 Integer giving the number of bootstrap replications.
#'   Defaults to \code{200}.
#' @param boot.parameter2 Numeric seed used for resampling of bootstrap.
#' @param nleqslv.method Character string defining the solver used by
#'   \code{\link[nleqslv]{nleqslv}}. Available choices include \code{"nleqslv"},
#'   \code{"Broyden"}, \code{"Newton"}, \code{"optim"}, \code{"BFGS"} and
#'   \code{"SANN"}.
#' @param optim.parameter1 Numeric tolerance for convergence of the outer loop.
#'   Defaults to \code{1e-6}.
#' @param optim.parameter2 Numeric tolerance for convergence of the inner loop.
#'   Defaults to \code{1e-6}.
#' @param optim.parameter3 Numeric constraint on the absolute value of
#'   parameters. Defaults to \code{100}.
#' @param optim.parameter4 Integer maximum number of outer loop iterations.
#'   Defaults to \code{50}.
#' @param optim.parameter5 Integer maximum number of \code{nleqslv}
#'   iterations per outer iteration. Defaults to \code{50}.
#' @param optim.parameter6 Integer maximum number of iterations for the
#'   Levenberg-Marquardt routine. Defaults to \code{50}.
#' @param optim.parameter7 Numeric convergence tolerance for the
#'   Levenberg-Marquardt routine. Defaults to \code{1e-10}.
#' @param optim.parameter8 Numeric tolerance for updating the Hessian in the
#'   Levenberg-Marquardt routine. Defaults to \code{1e-6}.
#' @param optim.parameter9 Numeric starting value for the Levenberg-Marquardt
#'   damping parameter lambda. Defaults to \code{1e-6}.
#' @param optim.parameter10 Numeric upper bound for lambda in the
#'   Levenberg-Marquardt routine. Defaults to \code{40}.
#' @param optim.parameter11 Numeric lower bound for lambda in the
#'   Levenberg-Marquardt routine. Defaults to \code{0.025}.
#' @param optim.parameter12 Numeric multiplicative increment applied to lambda
#'   when the Levenberg-Marquardt step is successful. Defaults to \code{2}.
#' @param optim.parameter13 Numeric multiplicative decrement applied to lambda
#'   when the Levenberg-Marquardt step is unsuccessful. Defaults to \code{0.5}.
#' @param data.initial.values Optional data frame providing starting values for
#'   the optimisation. Defaults to \code{NULL}.
#' @param should.normalize.covariate Logical indicating whether covariates should
#'   be centred and scaled prior to optimization. Defaults to \code{TRUE}.
#' @param should.terminate.time.point Logical indicating whether time points
#'   that contribute estimation are terminated by min of max follow-up times
#'   of each exposure level. Defaults to \code{TRUE}.
#' @param prob.bound Numeric lower bound used to truncate probabilities away
#'   from 0 and 1. Defaults to \code{1e-5}.
#'
#' @details
#'
#' ### Overview
#' `polyreg()` implements **log-odds product modeling** for CIFs at user-specified
#' time points, focusing on multiplicative effects of a categorical exposure, or
#' constant effects over time like Cox regression and Fine-Gray models. It estimates
#' multiplicative effects such as **risk ratios**, **odds ratios**, or
#' **subdistribution hazard ratios**, while ensuring that the probabilities across
#' competing events sum to one. This is achieved through
#' **reparameterization using polytomous log-odds products**, which fits so-called
#' effect-measure models and nuisance models on multiple competing events
#' simultaneously. Additionally, `polyreg()` supports direct binomial regression
#' for survival outcomes and the Richardson model for binomial outcomes,
#' both of which use log-odds products.
#'
#' The function follows the familiar **formula + data** syntax with `Event()` or
#' `Surv()` and outputs tidy results, including point estimates, standard errors,
#' confidence intervals, and p-values. Its results can be easily summarized with
#' `summary()` or combined with tools such as **modelsummary** or **broom** for reporting.
#'
#' ### Key arguments
#' -   `nuisance.model` — a formula describing the outcome and nuisance covariates,
#' excluding the exposure of interest.
#' -   `exposure` — specifies the categorical exposure variable
#' -   `effect.measure1` and `effect.measure2` — specifies the effect measures
#' for event1 and event2 (`"RR"`, `"OR"` or `"SHR"`).
#' -   `outcome.type` selects the outcome type (`"COMPETING-RISK"`, `"SURVIVAL"`,
#' `"BINOMIAL"`, `"PROPORTIONAL"` or `"POLY-PROPORTIONAL"`).
#' -   `time.point` — specifies time point at which the exposure effect is evaluated.
#' Required for `"COMPETING-RISK"` and `"SURVIVAL"` outcomes.
#' -   `strata` — specifies a stratification variable used to adjust for dependent censoring.
#'
#' ### Outcome type and event status coding
#'
#' The `outcome.type` argument must be set to:
#' -   Effects on cumulative incidence probabilities at a specific time: `"COMPETING-RISK"`
#' -   Effects on a risk at a specific time: `"SURVIVAL"`
#' -   Common effects on cumulative incidence probabilities over time: `"POLY-PROPORTIONAL"`
#' -   Common effects on a risk over time: `"PROPORTIONAL"`
#' -   Effects on a risk of a binomial outcome: `"BINOMIAL"`
#'
#' | Setting | Codes | Meaning |
#' |---|---|---|
#' | COMPETING-RISK | `code.event1`, `code.event2`, `code.censoring` | event of interest / competing event / censoring |
#' | COMPETING-RISK (default)| `code.event1=1`, `code.event2=2`, `code.censoring1` | event of interest / competing event / censoring |
#' | SURVIVAL | `code.event1`, `code.censoring`                   | event / censoring |
#' | SURVIVAL  (default)| `code.event1=1`, `code.censoring=0`     | event / censoring |
#' | SURVIVAL (ADaM-ADTTE) | `code.event1=0`, `code.censoring=1` | set to match ADaM convention |
#' | PROPORTIONAL | `code.event1`, `code.censoring`                   | event / censoring |
#' | PROPORTIONAL  (default)| `code.event1=1`, `code.censoring=0`     | event / censoring |
#' | PROPORTIONAL (ADaM-ADTTE) | `code.event1=0`, `code.censoring=1` | set to match ADaM convention |
#' | POLY-PROPORTIONAL | `code.event1`, `code.event2`, `code.censoring` | event of interest / competing event / censoring |
#' | POLY-PROPORTIONAL (default)| `code.event1=1`, `code.event2=2`, `code.censoring1` | event of interest / competing event / censoring |
#'
#' ### Effect measures for categorical exposure
#'
#' Choose the effect scale for event 1 and (optionally) event 2:
#'
#' | Argument | Applies to | Choices | Default |
#' |---|---|---|---|
#' | `effect.measure1` |event of interest | `"RR"`, `"OR"`, `"SHR"` | `"RR"` |
#' | `effect.measure2` | competing event | `"RR"`, `"OR"`, `"SHR"` | `"RR"` |
#'
#' - `RR`: risk ratio at `time.point` or common over time.
#' - `OR`: odds ratio at `time.point` or common over time.
#' - `SHR`: subdistribution hazard ratio  or common over time.
#'
#' ### Inference and intervals (advanced)
#'
#' | Argument | Meaning | Default |
#' |---|---|---|
#' | `conf.level` | Wald-type CI level | `0.95` |
#' | `report.sandwich.conf` | Sandwich variance CIs | `TRUE` |
#' | `report.boot.conf` | Bootstrap CIs (use if `"PROPORTIONAL"` or `"POLY-PROPORTIONAL"`) | `NULL` |
#' | `boot.bca` | Use BCa intervals (else normal approximation) | `FALSE` |
#' | `boot.parameter1` | Bootstrap reps | `200` |
#' | `boot.parameter2` | Seed for resampling | `46` |
#'
#' ### Optimization & solver controls (advanced)
#'
#' `polyreg()` solves estimating equations with optional inner routines.
#'
#' | Argument | Role | Default |
#' |---|---|---|
#' | `nleqslv.method` | Root solver | `"nleqslv"` |
#' | `optim.parameter1` / `optim.parameter2` | Outer/inner convergence tolerances | `1e-6`, `1e-6` |
#' | `optim.parameter3` | Parameter absolute bound | `100` |
#' | `optim.parameter4` | Max outer iterations | `50` |
#' | `optim.parameter5` | Max `nleqslv` iters per outer | `50` |
#' | `optim.parameter6:13` | Levenberg–Marquardt controls (iters/tolerances/lambda) | see defaults |
#'
#' Tips:
#' - If convergence warnings appear, relax/tighten tolerances or cap the parameter
#'   bound (`optim.parameter3`). Inspect `report.optim.convergence = TRUE`.
#'
#' ### Data handling and stability
#'
#' | Argument | Meaning | Default |
#' |---|---|---|
#' | `subset.condition` | An expression (as character) to subset `data` | `NULL` |
#' | `na.action` | NA handling function | `stats::na.omit` |
#' | `should.normalize.covariate` | Center/scale nuisance covariates | `TRUE` |
#' | `should.terminate.time.point` | Truncate support by exposure-wise follow-up maxima | `TRUE` |
#' | `prob.bound` | Truncate probabilities away from 0/1 (numerical guard) | `1e-5` |
#' | `data.initial.values` | Optional starting values data frame | `NULL` |
#'
#' ### Returned object and downstream use
#'
#' This function returns a list object that includes:
#' -   `coefficient` — regression coefficients
#' -   `cov` — variance-covariance matrix for regression coefficients
#' -   `diagnostic.statistics` — a data frame containing inverse probability weights,
#' influence functions, and predicted potential outcomes
#' -   `summary` — a summary of estimated exposure effects
#'
#' Use `summary` output with `msummary()` to display formatted results. The regression
#' coefficients and their variance-covariance matrix are provided as `coefficient`
#' and `cov`, respectively, with the first element corresponding to the intercept term,
#' subsequent elements to the covariates in `nuisance.model`, and the last element
#' to the variable specified by `exposure=`. Finally, `diagnostic.statistics` is
#' a data frame containing inverse probability weights, influence functions, and
#' predicted values of the potential outcomes of individual observations.
#'
#' ### Reproducibility and conventions
#'
#' - Set `boot.parameter2` for reproducible bootstrap results.
#' - Match CDISC ADaM conventions via `code.event1 = 0`, `code.censoring = 1`
#'   (and, if applicable, `code.event2` for competing events).
#' - Use `strata` when censoring may depend on baseline covariates (IPCW stratification).
#'
#' @importFrom nleqslv nleqslv
#' @importFrom boot boot boot.ci
#' @importFrom Rcpp sourceCpp
#' @importFrom stats IQR as.formula binomial coef glm mad median
#' @importFrom stats model.extract model.frame model.matrix na.omit na.pass
#' @importFrom stats pnorm qnorm rbinom reformulate rexp sd setNames terms time var
#' @useDynLib cifmodeling, .registration = TRUE
#'
#' @return A list containing fitted exposure effects and supporting results. The
#'   main components include \code{coefficient} (estimated exposure and
#'   covariate effects), \code{cov} (their variance-covariance matrix),
#'   \code{summary} (a tidy summary table compatible with
#'   \code{\link[modelsummary]{msummary}}) and \code{diagnostic.statistics}
#'   (inverse probability weights, influence functions and predicted potential
#'   outcomes).
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
#'   outcome.type = "COMPETING-RISK"
#' )
#' if (requireNamespace("modelsummary", quietly = TRUE)) {
#' modelsummary::msummary(output$summary, statistic = c("conf.int", "p.value"), exponentiate = TRUE)
#' }
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
    outcome.type = "COMPETING-RISK",
    conf.level = 0.95,
    report.nuisance.parameter = FALSE,
    report.optim.convergence = FALSE,
    report.sandwich.conf = TRUE,
    report.boot.conf = NULL,
    boot.bca = FALSE,
    boot.parameter1 = 200,
    boot.parameter2 = 46,
    nleqslv.method = "nleqslv",
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
    should.normalize.covariate = TRUE,
    should.terminate.time.point = TRUE,
    prob.bound = 1e-5
) {

  #######################################################################################################
  # 1. Pre-processing (function: checkSpell, checkInput, normalizeCovariate, sortByCovariate)
  #######################################################################################################
#  computation.time0 <- proc.time()
  outcome.type  <- check_outcome.type(outcome.type, formula=formula, data=data)
  ce <- check_effect.measure(effect.measure1, effect.measure2)
  ci <- check_input_polyreg(data, nuisance.model, exposure, code.event1, code.event2, code.censoring, code.exposure.ref, outcome.type, conf.level, report.sandwich.conf, report.boot.conf, nleqslv.method, should.normalize.covariate)
  should.normalize.covariate <- ci$should.normalize.covariate
  report.sandwich.conf <- ci$report.sandwich.conf
  report.boot.conf <- ci$report.boot.conf

  data <- createAnalysisDataset(formula=nuisance.model, data=data, other.variables.analyzed=c(exposure, strata), subset.condition=subset.condition, na.action=na.action)
  out_normalizeCovariate <- normalizeCovariate(nuisance.model, data, should.normalize.covariate, outcome.type, ci$out_readExposureDesign$exposure.levels)
  normalized_data <- out_normalizeCovariate$normalized_data
  tp <- read_time.point(nuisance.model, normalized_data, ci$out_readExposureDesign$x_a, outcome.type, code.censoring, should.terminate.time.point, time.point)
  index.vector <- calculateIndexForParameter(NA, ci$x_l, ci$x_a, length(tp))

  estimand <- list(
    effect.measure1=ce$effect.measure1,
    effect.measure2=ce$effect.measure2,
    time.point=tp,
    code.event1=code.event1,
    code.event2=code.event2,
    code.censoring=code.censoring,
    code.exposure.ref=code.exposure.ref,
    exposure.levels=ci$out_readExposureDesign$exposure.levels,
    index.vector=index.vector
  )

  optim.method <- list(
    nleqslv.method = nleqslv.method,
    optim.parameter1 = optim.parameter1,
    optim.parameter2 = optim.parameter2,
    optim.parameter3 = optim.parameter3,
    optim.parameter4 = optim.parameter4,
    optim.parameter5 = optim.parameter5,
    optim.parameter6 = optim.parameter6,
    optim.parameter7 = optim.parameter7,
    optim.parameter8 = optim.parameter8,
    optim.parameter9 = optim.parameter9,
    optim.parameter10 = optim.parameter10,
    optim.parameter11 = optim.parameter11,
    optim.parameter12 = optim.parameter12,
    optim.parameter13 = optim.parameter13
  )

  #######################################################################################################
  # 2. Pre-processing and Calculating initial values alpha_beta_0 (function: calculateInitialValues)
  #######################################################################################################
  if (outcome.type == "COMPETING-RISK" || outcome.type == "SURVIVAL" || outcome.type == "BINOMIAL") {
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
  } else if (outcome.type == "PROPORTIONAL" || outcome.type == "POLY-PROPORTIONAL") {
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
  if (outcome.type == "COMPETING-RISK" || outcome.type == "SURVIVAL") {
    ip.weight.matrix <- calculateIPCW(nuisance.model, normalized_data, code.censoring, strata, estimand$time.point)
  } else if (outcome.type == "BINOMIAL") {
    ip.weight.matrix <- matrix(1,nrow(normalized_data),1)
  } else if (outcome.type == "PROPORTIONAL" || outcome.type == "POLY-PROPORTIONAL") {
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
      recent <- x[(n - stall_patience + 1L):n]
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
  nleqslv_method  <- choose_nleqslv_method(nleqslv.method)
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
      method  = nleqslv_method,
      control = list(maxit = optim.parameter5, allowSingular = FALSE)
    )
    new_params <- out_nleqslv$x
    current_obj_value <- drop(crossprod(obj$estimating_equation_i(new_params)))

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    ac <- assessConvergence(new_params, prev_params, current_obj_value, optim.parameter1, optim.parameter2, optim.parameter3)

    nleqslv.info <- extractOptimizationInfo(out_nleqslv, nleqslv.method)
#    computation.time.second <- as.numeric((proc.time() - computation.time0)[3])

    trace_df <- append_trace(
      trace_df,
      iteration = iteration,
      computation.time.second = 46,
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
  out_calculateCov <- switch(
    outcome.type,
    "COMPETING-RISK"   = calculateCov(out_getResults, estimand, prob.bound),
    "SURVIVAL"         = calculateCovSurvival(out_getResults, estimand, prob.bound),
    "BINOMIAL"         = calculateCovSurvival(out_getResults, estimand, prob.bound),
    "PROPORTIONAL"     = NULL,
    "POLY-PROPORTIONAL"= NULL,
    stop(sprintf("Unsupported outcome.type for covariance: %s", outcome.type))
  )

  out_normalizeEstimate <- normalizeEstimate(
    outcome.type               = outcome.type,
    report.sandwich.conf       = report.sandwich.conf,
    should.normalize.covariate = should.normalize.covariate,
    current_params             = current_params,
    out_getResults             = out_getResults,
    estimand                   = estimand,
    prob.bound                 = prob.bound,
    out_normalizeCovariate     = out_normalizeCovariate,
    out_calculateCov           = out_calculateCov
  )

  alpha_beta_estimated <- out_normalizeEstimate$alpha_beta_estimated
  cov_estimated        <- out_normalizeEstimate$cov_estimated

  #######################################################################################################
  # 6. Calculating bootstrap confidence interval (functions: boot, solveEstimatingEquation)
  #######################################################################################################
  if (isTRUE(report.boot.conf)) {
    set.seed(boot.parameter2)
    boot.coef     <- rep(NA,2)
    boot.coef_se  <- rep(NA,2)
    boot.p_value  <- rep(NA,2)
    boot.conf_low <- rep(NA,2)
    boot.conf_high<- rep(NA,2)

    if (outcome.type=="POLY-PROPORTIONAL") {
      index_coef    <- c(length(estimand$time.point) + 1, 2*length(estimand$time.point) + 2)
    } else if (outcome.type=="PROPORTIONAL") {
      index_coef    <- c(length(estimand$time.point) + 1)
    } else if (outcome.type=="COMPETING-RISK") {
      index_coef <- seq_len(2*out_normalizeCovariate$n_covariate+4)
    } else if (outcome.type=="BINOMIAL" | outcome.type=="SURVIVAL") {
      index_coef <- seq_len(out_normalizeCovariate$n_covariate+2)
    }

    boot_function <- function(data, indices) {
      res <- try({
        coef <- solveEstimatingEquation(nuisance.model=nuisance.model, exposure=exposure, strata=strata, normalized_data = data[indices, , drop = FALSE], outcome.type=outcome.type, estimand=estimand, optim.method=optim.method, out_normalizeCovariate=out_normalizeCovariate, prob.bound=prob.bound, alpha_beta_0=alpha_beta_0)
      }, silent = TRUE)
      if (inherits(res, "try-error")) NA_real_ else as.numeric(res)
    }

    out_boot <- boot(normalized_data, boot_function, R = boot.parameter1)
    for (j in index_coef) {
      if (isTRUE(boot.bca)) {
        out_boot.ci <- boot.ci(out_boot, conf = conf.level, index = index_coef[j], type = c("norm", "bca"))
        boot.coef[j] <- (out_boot.ci$normal[2] + out_boot.ci$normal[3])/2
        ci_range <- out_boot.ci$normal[3] - out_boot.ci$normal[2]
        boot.coef_se[j] <- ci_range/2/qnorm(1 - (1-conf.level)/2)
        boot.p_value[j] <- 2 * (1 - pnorm(abs(boot.coef[j]) / boot.coef_se[j]))
        boot.conf_low[j] <- out_boot.ci$bca[4]
        boot.conf_high[j] <- out_boot.ci$bca[5]
      } else {
        ok <- apply(out_boot$t, 1L, function(x) all(is.finite(x)))
        t_ok <- out_boot$t[ok, , drop = FALSE]
#        cat("Skipped:", sum(!ok), " / ", length(ok), "replicates\n")

        mean_ok <- colMeans(t_ok)
        sd_ok   <- apply(t_ok, 2, sd)
        ci_ok <- cbind(
          lower = mean_ok - qnorm(0.975) * sd_ok,
          upper = mean_ok + qnorm(0.975) * sd_ok
        )
        boot.coef[j] <- mean_ok[j]
        boot.coef_se[j] <- sd_ok[j]
        boot.p_value[j] <- 2 * (1 - pnorm(abs(boot.coef[j]) / boot.coef_se[j]))
        boot.conf_low[j] <- ci_ok[j,1]
        boot.conf_high[j] <- ci_ok[j,2]
      }
    }
    out_bootstrap <- list(
      boot.coef=boot.coef, boot.coef_se=boot.coef_se, boot.p_value=boot.p_value, boot.conf_low=boot.conf_low, boot.conf_high=boot.conf_high
    )
  } else {out_bootstrap <- NULL}

  #######################################################################################################
  # 7. Output (functions: reportSurvival, reportCOMPETING-RISK, reportPrediction)
  #######################################################################################################
  out_summary <- reportEffects (
    outcome.type, report.nuisance.parameter, report.optim.convergence, report.sandwich.conf, report.boot.conf,
    nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated,
    out_bootstrap, out_getResults, iteration, converged.by, objective.function, max.absolute.difference, relative.difference,
    out_nleqslv, conf.level, optim.method$nleqslv.method
  )
  if (outcome.type == "COMPETING-RISK" || outcome.type == "SURVIVAL" || outcome.type == "BINOMIAL") {
    data$influence.function <- out_calculateCov$influence.function
    data$ip.weight <- out_getResults$ip.weight
    data$potential.CIFs <- out_getResults$potential.CIFs
  }
  out_data <- data
  out <- list(summary=out_summary, coefficient=alpha_beta_estimated, cov=cov_estimated, bootstrap=out_bootstrap, diagnostic.statistics=out_data, optimization.info=trace_df)
  return(out)
}

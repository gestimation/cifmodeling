#' Simulate 2-arm competing risks data
#'
#' @param n              Total sample size.
#' @param cif1_control   Numeric matrix with two columns: time and baseline CIF for cause 1 in the control group.
#' @param cif2_control   Same as cif1_control, for cause 2.
#' @param rr_cif1        Numeric scalar, multiplicative effect of treatment on cause 1.
#' @param rr_cif2        Numeric scalar, multiplicative effect of treatment on cause 2.
#' @param rr_L           Numeric scalar, multiplicative effect shared by covariates (L1, L2).
#' @param censoring_rate Optional numeric, rate parameter for exponential censoring if necessary.
#' @param tau            Numeric maximum follow-up time.
#' @param link_cif1      Character, link function for cause 1: c("cloglog","logistic","log").
#' @param link_cif2      Character, link function for cause 2: c("cloglog","logistic","log").
#' @param estimated_ps   Logical. TRUE to use estimated PS in IPTW, FALSE to use true PS.
#' @param overlap_ps     Character, "strong", "partial", or "weak" overlap in the propensity scores.
#' @param binary_L1      Logical, TRUE if L1 is binary, FALSE if L1 is normal.
#' @param binary_L2      Logical, TRUE if L2 is binary, FALSE if L2 is normal.
#' @param seed           Optional integer for reproducibility.
#'
#' @return A data.frame with columns:
#'   - id        : subject ID
#'   - A         : treatment indicator (0/1)
#'   - L1, L2    : baseline covariates
#'   - ps_true   : true propensity score P(A=1 | L)
#'   - ps_ipw    : PS actually used in IPTW (true or estimated)
#'   - s_ipw     : stabilised IPTW
#'   - time      : observed time (event or censoring)
#'   - status    : 0 = censored, 1 = event of interest, 2 = competing event
generate_competing_risk_data <- function(
    n,
    cif1_control,
    cif2_control,
    rr_cif1,
    rr_cif2,
    rr_L,
    censoring_rate = NULL,
    tau = NULL,
    link_cif1 = c("cloglog", "logistic", "log"),
    link_cif2 = c("cloglog", "logistic", "log"),
    estimated_ps = TRUE,
    overlap_ps   = c("strong", "partial", "weak"),
    binary_L1    = FALSE,
    binary_L2    = FALSE,
    seed         = NULL
) {
  link_cif1  <- match.arg(link_cif1)
  link_cif2  <- match.arg(link_cif2)
  overlap_ps <- match.arg(overlap_ps)

  if (!is.null(seed)) set.seed(seed)

  L1 <- if (binary_L1) rbinom(n, 1, 0.5) - 0.5 else rnorm(n)
  L2 <- if (binary_L2) rbinom(n, 1, 0.5) - 0.5 else rnorm(n)

  if (overlap_ps == "strong") {
    beta1 <- 0.5; beta2 <- 0
  } else if (overlap_ps == "partial") {
    beta1 <- 1.0; beta2 <- 0
  } else { # "weak"
    beta1 <- 1.0; beta2 <- 1.0
  }
  beta0   <- 0
  ps_true <- plogis(beta0 + beta1 * L1 + beta2 * L2)
  A       <- rbinom(n, 1, ps_true)

  if (isTRUE(estimated_ps)) {
    df_ps  <- data.frame(A = A, L1 = L1, L2 = L2)
    fit_ps <- glm(A ~ L1 + L2, data = df_ps, family = binomial())
    ps_hat <- drop(plogis(model.matrix(fit_ps) %*% coef(fit_ps)))
    ps_ipw <- ps_hat
  } else {
    ps_ipw <- ps_true
  }

  pA <- mean(A)
  numer <- ifelse(A == 1, pA, 1 - pA)
  denom <- ifelse(A == 1, ps_ipw, 1 - ps_ipw)
  s_ipw   <- numer / denom

  X    <- cbind(A, L1, L2)
  lrr1 <- c(log(rr_cif1), log(rr_L), log(rr_L))
  lrr2 <- c(log(rr_cif2), log(rr_L), log(rr_L))

  rr1  <- exp(drop(X %*% lrr1))
  rr2  <- exp(drop(X %*% lrr2))

  max_f1 <- generate_cif_treated(cif1_control, rr = max(rr1), type = link_cif1)
  max_f2 <- generate_cif_treated(cif2_control, rr = max(rr2), type = link_cif2)
  if (max_f1 + max_f2 > 1) {
    warning("Models do not satisfy CIF1 + CIF2 <= 1 at the maximum time.")
  }

  T1 <- generate_subdistribution(cif1_control,rr1,type=link_cif1)
  T2 <- generate_subdistribution(cif2_control,rr2,type=link_cif2)

  cif1_treated <- generate_cif_treated(cif1_control,rr1,type=link_cif1)
  cif2_treated <- generate_cif_treated(cif2_control,rr2,type=link_cif2)
  dies <- rbinom(n,1,cif1_treated+cif2_treated)
  sel1 <- rbinom(n,1,cif2_treated/(cif1_treated+cif2_treated))+1
  epsilon  <- dies*(sel1)
  T1$epsilon <- epsilon

  T1$time <- T1$timecause
  T1$time2 <- T2$timecause
  T1$status <- epsilon

  if (is.null(tau))  {
    tau <- tail(cif1_control, 1)
    tau <- tau[1,1]
  }
  T1 <- dtransform(T1, time=tau,   epsilon==0)
  T1 <- dtransform(T1, status=0,   epsilon==0)
  T1 <- dtransform(T1, time=time2, epsilon==2)
  T1 <- dtransform(T1, status=2,   epsilon==2)
output1 <- survfit(Surv(time, status)~1, data=data)

  if (!is.null(censoring_rate))  {
    cc <- rexp(n)/censoring_rate
    T1$status <- ifelse(T1$time<cc,T1$status,0)
    T1$time <- pmin(T1$time, cc)
  }

  T1$A <- A
  T1$L1 <- L1
  T1$L2 <- L2
  T1$ps_true <- ps_true
  T1$ps_ipw <- ps_ipw
  T1$s_ipw <- s_ipw
  return(T1)
}

generate_cif_treated <- function(cif,rr,type=c("log","cloglog","logistic")) {
  mcif <- max(cif[,2])
  if (type[1]=="log") mcif <- mcif*rr
  if (type[1]=="cloglog") mcif <- 1- exp(-mcif*rr)
  if (type[1]=="logistic") mcif <- mcif* rr/(1 + mcif * rr)
  return(mcif)
}

lin.approx <- function (x2, xfx, x = 1) {
  breaks <- xfx[, x]
  fx <- xfx[, -x]
  ri <- fast.approx(breaks, x2, type = "left")
  maxindex <- which(ri == length(breaks))
  rip1 <- ri + 1
  rip1[maxindex] <- length(breaks)
  rrr <- (x2 - breaks[ri])/(breaks[rip1] - breaks[ri])
  rrr[maxindex] <- 0
  res <- rrr * (fx[rip1] - fx[ri]) + fx[ri]
  res[is.na(res)] <- tail(fx, 1)
  res[is.na(rrr)] <- fx[ri][is.na(rrr)]
  return(res)
}

dtransform <- function (data, ...) {
  if (is.vector(data))
    data <- data.frame(data)
  e <- eval(substitute(list(...)), data, parent.frame())
  tags <- names(e)
  condn <- match("", tags)
  if (!is.na(condn)) {
    condition <- TRUE
    cond <- e[[condn[1]]]
    whereT <- which(cond)
    e[[condn]] <- NULL
    tags <- tags[-condn]
  }
  else condition <- FALSE
  inx <- match(tags, names(data))
  matched <- !is.na(inx)
  matchedtags <- seq(length(e))[matched]
  if (any(matched)) {
    k <- 1
    if (condition == TRUE) {
      for (i in inx[matched]) {
        mk <- matchedtags[k]
        if (length(e[[mk]]) == 1)
          data[whereT, i] <- e[[mk]]
        else data[whereT, i] <- e[[mk]][whereT]
        k <- k + 1
      }
    }
    else data[inx[matched]] <- e[matched]
    data <- data.frame(data)
  }
  if (!all(matched)) {
    if (condition)
      for (i in seq(length(e))[!matched]) {
        if (length(e[[i]]) == 1)
          e[[i]] <- rep(e[[i]], nrow(data))
        e[[i]][!cond] <- NA
      }
    data <- cbind(data, data.frame(e[!matched]))
  }
  return(data)
}

generate_subdistribution <-function (cumhazard, rr, n = NULL, entry = NULL, type = "cloglog", startcum = c(0, 0), ...)
{
  if (!is.null(n))
    rr <- rep(1, n)
  entry = NULL
  logit <- function(p) log(p/(1 - p))
  if (cumhazard[1, 2] > 0)
    cumhazard <- rbind(startcum, cumhazard)
  breaks <- cumhazard[, 1]
  cumh <- cumhazard[, 2]
  mm <- tail(breaks, 1)
  n <- length(rr)
  if (type == "cloglog") {
    F1tau <- 1 - exp(-tail(cumh, 1) * rr)
    ttt <- -log(1 - runif(n) * F1tau)/rr
  }
  else if (type == "logistic") {
    F1tau <- tail(cumh, 1) * rr/(1 + tail(cumh, 1) * rr)
    v <- runif(n) * F1tau
    ttt <- exp(logit(v))/rr
  }
  else if (type == "rr" | type == "cif") {
    F1tau <- tail(cumh, 1)
    ttt <- runif(n) * F1tau
    F1tau <- F1tau * rr
  }
  else stop(" cloglog or logistic or give function (fun=) \n")
  entry <- cumentry <- rep(0, n)
  ttte <- ttt + cumentry
  rrx <- lin.approx(ttt, cumhazard, x = -1)
  timecause <- rrx
  rrx <- ifelse(rrx > mm, mm, rrx)
  status <- rbinom(n, 1, F1tau)
  rrx[status == 0] <- mm
  dt <- data.frame(entry = entry, time = rrx, status = status,
                   rr = rr, F1tau = F1tau, timecause = timecause)
  attr(dt, "cumhaz") <- cumhazard
  return(dt)
}

cif1 <- cbind(c(0,10,20,100),c(0,0.15,0.2,0.22))
cif2 <- cbind(c(0,10,20,100),c(0,0.7,0.75,0.78))


data <- generate_competing_risk_data(20, cif1, cif2, 1, 1, 1, link_cif1="logistic", link_cif2="logistic", overlap_ps="strong")
data$w <- 2

library(survival)
library(mstate)
data$event <- factor(data$status, 0:2, c("censoring", "event1", "event2"))
output1 <- survfit(Surv(time, event)~1, data=data)
print(output1$pstate)
print(output1$std.err)
output2 <- survfit(Surv(time, event)~1, data=data, w=w)
print(output2$pstate)
print(output2$std.err)
#print(output1$lower)
#print(output1$upper)

output3 <- cifcurve(Surv(time, event)~1, data=data, outcome.type="competing-risk", error="delta")
print(output3$std.err)
output4 <- cifcurve(Surv(time, event)~1, data=data, outcome.type="competing-risk", error="delta", weights="w")
print(output4$std.err)

print(output2$surv)
print(output1$std.err)
print(output2$std.err)











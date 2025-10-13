calculateRMST <-function(out_survfit, tau, alpha = 0.05)
{
  idx = out_survfit$time <= tau
  wk.time = sort(c(out_survfit$time[idx], tau))
  wk.surv = out_survfit$surv[idx]
  wk.n.risk = out_survfit$n.risk[idx]
  wk.n.event = out_survfit$n.event[idx]
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst
  wk.var <- ifelse((wk.n.risk - wk.n.event) == 0, 0, wk.n.event/(wk.n.risk*(wk.n.risk - wk.n.event)))
  wk.var = c(wk.var, 0)
  rmst.var = sum(cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se = sqrt(rmst.var)
  out = matrix(0, 2, 4)
  out[1, ] = c(rmst, rmst.se, rmst - qnorm(1 - alpha/2) * rmst.se,
               rmst + qnorm(1 - alpha/2) * rmst.se)
  out[2, ] = c(tau - out[1, 1], rmst.se, tau - out[1, 4], tau -
                 out[1, 3])
  rownames(out) = c("RMST", "RMTL")
  colnames(out) = c("Est.", "se", paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""), paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""))
  Z = list()
  Z$result = out
  Z$rmst = out[1, ]
  Z$rmtl = out[2, ]
  Z$tau = tau
  Z$rmst.var = rmst.var
  Z$fit = out_survfit
  return(Z)
}


createAtRiskMatrix <- function(t) {
  atrisk <- outer(t, t, "<=")
  return(atrisk)
}

calculateKaplanMeier <- function(t, d){
  n <- length(t)
  data <- data.frame(t = t, d = d, id = 1:n)
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_d <- sorted_data$d
  sorted_id <- sorted_data$id
  atrisk <- createAtRiskMatrix(sorted_t)
  n_atrisk <- rowSums(atrisk)
  s <- 1 - sorted_d / n_atrisk
  km <- cumprod(s)
  data <- data.frame(id=sorted_id, km=km)
  sorted_data <- data[order(data$id), ]
  km <- sorted_data$km
  return(km)
}

calculateNelsonAalen <- function(t, d) {
  atrisk <- createAtRiskMatrix(t)
  n_atrisk <- rowSums(atrisk)
  na <- d / n_atrisk
  return(na)
}

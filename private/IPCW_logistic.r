#' ---
#' title: "IPCWロジスティック回帰"
#' author: "Shigetaka Kobari"
#' date: "2025-01-02"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE---------------------------------------------------------------------------
library(MASS)
#library(matrixStats)

library(survival)
#library(ggsurvfit)
#library(survminer)
set.seed(42)

#' 
#' ## セッション1. 回帰分析
#' 
#' ### 乱数の発生
#' 
## ---------------------------------------------------------------------------------------------------
n = 10
mean_x = 0
std_dev_x = 1

covariate = matrix(rnorm(n, mean = mean_x, sd = std_dev_x), ncol = 1)
print(covariate)

#' 
#' ### 列ベクトルの結合
#' 
## ---------------------------------------------------------------------------------------------------
intercept =  matrix(1, n, 1)
x = cbind(intercept, covariate)
print(x)

#' 
#' ### ベクトルの積
#' 
## ---------------------------------------------------------------------------------------------------
alpha_true = matrix(c(7, 7), 2, 1)
conditional_ey = x %*% alpha_true
ey = mean(conditional_ey)

#' 
#' ### 単回帰に伴う乱数の発生
## ---------------------------------------------------------------------------------------------------
error = matrix(rnorm(n, mean = 0, sd = 1), ncol = 1)
y = conditional_ey+error
print(y)

#' 
#' ### 逆行列と一般化逆行列
#' 
## ---------------------------------------------------------------------------------------------------
# 逆行列
inv_xx = solve(t(x) %*% x)
# alpha_hat = (x^Tx)^(-1)x^Ty
alpha_hat = inv_xx %*% t(x) %*% y
print(alpha_hat)

# 一般化逆行列
sx = cbind(intercept, intercept, covariate)
# inv_sxx = solve(t(sx) %*% sx)
ginv_sxx =ginv(t(sx) %*% sx) 
alpha_hat = ginv_sxx %*% t(sx) %*% y
print(alpha_hat)

#' 
#' ### 正規方程式の解
#' 
## ---------------------------------------------------------------------------------------------------
normal_equation = function(alpha){
  residuals = y-x %*% alpha
  gradient = t(x) %*% residuals
  return(sum(gradient^2))
}

alpha_0 = matrix(c(-1, -1), 2, 1)

# BFGS: 準ニュートン法
# fn: 最小化する関数
result = optim(par=alpha_0, fn=normal_equation, method="BFGS")
print(result)

#' 
#' ### 線形回帰を実行する関数
#' 
## ---------------------------------------------------------------------------------------------------
fit_linear_regression <- function(formula, data){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  
  alpha_0 <- rep(0, ncol(x))
  
  normal_equation <- function(alpha){
    residuals <- y-x %*% alpha
    gradient <- t(x) %*% residuals
    return(sum(gradient^2))
  }
  
  result <- optim(par=alpha_0, fn=normal_equation, method="BFGS")
  beta <- result$par
  residuals <- y - x %*% beta
  fitted.values <- x %*% beta
  
  # モデルの自由度を計算
  n <- length(y)        # 観測値の数
  p <- ncol(x)          # 説明変数の数（係数の数）
  df_residual <- n - p  # 残差の自由度
  
  # QR分解を追加
  qr_decomp <- qr(x)
  
  # lm型オブジェクトを模倣
  z <- list(
    coefficients = setNames(beta, colnames(x)),
    residuals = residuals,
    fitted.values = fitted.values,
    terms = mt,
    call = cl,
    rank = qr_decomp$rank,
    qr = qr_decomp,
    xlevels = .getXlevels(mt, mf),
    model = mf,
    df.residual = df_residual # 残差自由度を追加
  )
  
  # lm型クラスを設定
  class(z) <- "lm"
  return(z)
}

df = data.frame(y=y, x=covariate)
fit = fit_linear_regression(y~x, df)
print(fit$coefficients)
summary(fit)

## ---------------------------------------------------------------------------------------------------
fit = lm(y~x, df)
summary(fit)

#' 
#' ## セッション2. ロジスティック回帰
#' 
#' ### ロジスティック回帰に従う乱数の発生
#' 
## ---------------------------------------------------------------------------------------------------
n_0 = 100
n_1 = 100
n = n_0+n_1
one = matrix(1, n, 1)

y_0 = matrix(rbinom(n_0, size = 1, prob = 0.5), ncol=1)
y_1 = matrix(rbinom(n_1, size = 1, prob = 0.2689414), ncol=1)

y = rbind(y_0, y_1)

#' 
#' ### アルゴリズムのアウトライン: モデルの特定→方程式→数値計算
#' 
## ---------------------------------------------------------------------------------------------------
intercept = matrix(1, n, 1)
a = matrix(c(rep(0, n_0), rep(1, n_1)), ncol = 1)
x = cbind(intercept, a)

beta_true = matrix(c(0, -1), 2, 1)
xbeta = x %*% beta_true
conditional_p = exp(xbeta) / (1 + exp(xbeta))

p0_hat = sum(y[a == 0]) / sum(a == 0) 
p1_hat = sum(y[a == 1]) / sum(a == 1) 
print(p0_hat)
print(p1_hat)

#' 
#' ### スコア方程式の解
#' 
## ---------------------------------------------------------------------------------------------------
# スコア方程式
score_equation = function(beta){
  xbeta = x %*% beta
  conditional_p = exp(xbeta) / (1 + exp(xbeta))
  w = 1 / (conditional_p * (1 - conditional_p))
  #w = 1
  zero = t(x) %*% (w * (y - conditional_p)) 
  return(sum(zero^2))
}


beta_0 = matrix(c(0, 0), 2, 1)
result = optim(par=beta_0, fn=score_equation, method="BFGS")

beta = result$par
print(beta)

x_tmp = matrix(c(1, 1, 0, 1), nrow = 2, ncol = 2)
# 列優先で埋まる
#      [,1] [,2]
#[1,]    1    0
#[2,]    1    1
xbeta = x_tmp %*% beta 

p0_hat_p1_hat = exp(xbeta) / (1 + exp(xbeta))
print(p0_hat_p1_hat)


#' 
#' ### ロジスティック回帰を実行する関数
#' 
## ---------------------------------------------------------------------------------------------------
fit_logistic <- function(formula, data) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  
  # スコア方程式（最尤推定）
  score_equation <- function(beta) {
    xbeta <- x %*% beta
    conditional_p <- exp(xbeta) / (1 + exp(xbeta))
    zero <- t(x) %*% (y - conditional_p)
    return(sum(zero^2))
  }
  
  # 初期値
  beta_0 <- rep(0, ncol(x))
  
  # 最適化
  result <- optim(par = beta_0, fn = score_equation, method = "BFGS")
  
  # 結果の格納
  beta <- result$par
  fitted_values <- 1 / (1 + exp(-x %*% beta))  # フィッティングされた値
  residuals <- y - fitted_values              # 残差
  
  # Null deviance（全ての予測値が平均の場合の偏差）
  y_mean <- mean(y)
  null_deviance <- sum(-2 * (y * log(y_mean) + (1 - y) * log(1 - y_mean)))
  
  # Residual deviance
  residual_deviance <- sum(-2 * (y * log(fitted_values) + (1 - y) * log(1 - fitted_values)))
  
  # 自由度
  df_null <- length(y) - 1
  df_residual <- length(y) - length(beta)
  
  # AICの計算
  log_likelihood <- sum(y * log(fitted_values) + (1 - y) * log(1 - fitted_values))
  aic <- -2 * log_likelihood + 2 * length(beta)
  
  # 返り値の構造
  result_list <- list(
    coefficients = setNames(beta, colnames(x)),
    fitted.values = fitted_values,
    residuals = residuals,
    formula = formula,
    terms = mt,
    data = data,
    call = cl,
    null.deviance = null_deviance,
    residual.deviance = residual_deviance,
    deviance = residual_deviance,  # devianceを追加
    df.null = df_null,
    df.residual = df_residual,
    aic = aic,  # AICを追加
    x = x,
    y = y
  )
  
  # クラス付与
  class(result_list) <- c("glm", "lm")
  return(result_list)
}

df = data.frame(y=y, a=a)
fit = fit_logistic(y~a, df)
print(fit)

#' 
## ---------------------------------------------------------------------------------------------------
# Rのglm関数の結果と比較
fit = glm(y~a, df, family = binomial(link = "logit"))
print(fit)

#' 
#' ## セッション3. Kaplan-Meier推定量
#' 
#' ### 3年生存確率0.5となるハザードの計算
#' 
## ---------------------------------------------------------------------------------------------------
tau = 3
survival_at_tau = 0.5
hazard = -log(survival_at_tau) / tau
print(hazard)

#' 
#' ### 指数乱数の発生
#' 
## ---------------------------------------------------------------------------------------------------
n = 10
one = matrix(1, n, 1)
true_time = matrix(rexp(n), ncol=1) / hazard
print(true_time)

# 確認のため、SASで作成した乱数を読み込む
#df3 = read.csv('sandwich_test_for_km.csv')
#true_time = df3$true_time

mean_t = mean(true_time)
hazard_hat = 1 / mean_t
print(hazard_hat)


#' 
#' ### 打ち切り変数の発生
## ---------------------------------------------------------------------------------------------------
# tau = 3より大きい場合は打ち切り
censoring_time = matrix(rep(tau, n), ncol = 1)
d = as.integer(true_time <= censoring_time) 
t = ifelse(d == 1, true_time, censoring_time)
print(d)
print(t)

#' 
#' ### アットリスク行列
#' 
## ---------------------------------------------------------------------------------------------------
t_matrix = matrix(rep(t, each = n), nrow = n, byrow = TRUE) 
atrisk = t(t_matrix) >= t_matrix # atrisk [i, j]: 観測t_iの時点で個体jがリスク集団に含まれる
n_atrisk = rowSums(atrisk) 
print(n_atrisk)

#' 
#' ### 行列のソート
#' 
## ---------------------------------------------------------------------------------------------------
# dfを経由すればソートしやすい
data = data.frame(t = t, d = d)
# 第1列（t）でソート
sorted_data = data[order(data$t), ]  # t を基準にソート
sorted_t = sorted_data$t            # ソートされた観測時間
sorted_d = sorted_data$d

#' 
#' ### 累積和関数を用いたKaplan-Meier推定量
#' 
## ---------------------------------------------------------------------------------------------------
t_matrix = matrix(rep(sorted_t, each = n), nrow = n, byrow = TRUE) 
atrisk = t(t_matrix) >= t_matrix
n_atrisk = rowSums(atrisk) 

s = 1 - sorted_d / n_atrisk # 生存確率の更新
log_s = log(s) # 対数変換
km = exp(cumsum(log_s)) # 累積和をとったあと指数変換
print(sorted_t)
print(km)
print(n_atrisk)





#' 
## ---------------------------------------------------------------------------------------------------
# 生存時間のパッケージと一致するか確認
data = data.frame(sorted_t = sorted_t, sorted_d = sorted_d)
km_fit = survfit(Surv(sorted_t, sorted_d) ~ 1, data=data)
ggsurvplot(km_fit, data = data, 
           conf.int = FALSE,               
           xlab = "Time",                 
           ylab = "Survival Probability", 
           title = "Kaplan-Meier Curve", 
           risk.table = TRUE) 

#' 
#' ### Kaplan-Meier推定量の関数
## ---------------------------------------------------------------------------------------------------
# Kaplan-Meier推定量の関数
kaplan_meier = function(t, d){
  # dfを作成
  n = length(t)
  data = data.frame(t = t, d = d, id = 1:n)
  
  # 観察時間でsort
  sorted_data = data[order(data$t), ]
  sorted_t = sorted_data$t
  sorted_d = sorted_data$d
  sorted_id = sorted_data$id
  
  # リスク集合のサイズを計算
  t_matrix = matrix(rep(sorted_t, each = n), nrow = n, byrow = TRUE) 
  atrisk = t(t_matrix) >= t_matrix
  n_atrisk = rowSums(atrisk) 
  
  s = 1 - sorted_d / n_atrisk 
  log_s = log(s)
  km = exp(cumsum(log_s)) 
  
  # sortされたものを戻す
  data = data.frame(id=sorted_id, km=km)
  sorted_data = data[order(data$id), ]
  km = sorted_data$km
  
  return(km)
}

# Nelson-Aalen推定量の関数
nelson_aalen <- function(t, d){
  # dfを作成
  n = length(t)
  data = data.frame(t = t, d = d, id = 1:n)
  
  # 観察時間でsort
  sorted_data = data[order(data$t), ]
  sorted_t = sorted_data$t
  sorted_d = sorted_data$d
  sorted_id = sorted_data$id
  
  # リスク集合のサイズを計算
  t_matrix = matrix(rep(sorted_t, each = n), nrow = n, byrow = TRUE) 
  atrisk = t(t_matrix) >= t_matrix
  n_atrisk = rowSums(atrisk) 
  
  na = sorted_d / n_atrisk
  
  # sortされたものを戻す
  data = data.frame(id=sorted_id, na=na)
  sorted_data = data[order(data$id), ]
  na = sorted_data$na
  return(na)
}

s_hat = kaplan_meier(t, d)
print(t)
print(d)
print(s_hat)

#' 
#' 
#' ## セッション4. IPCWロジスティック回帰 
#' 
#' ### 3年生存確率0.5/0.2689414となるハザードの計算
#' 
## ---------------------------------------------------------------------------------------------------
tau = 3
hazard0 = -log(0.5) / tau
hazard1 = -log(0.2689414) / tau

#' 
## ---------------------------------------------------------------------------------------------------
# 準備
n_0 = 100
n_1 = 100
n = n_0 + n_1

#' 
#' ### 指数乱数の発生
#' 
## ---------------------------------------------------------------------------------------------------
# 指数乱数の発生
true_time0 = matrix(rexp(n_0, rate=hazard0), ncol=1) 
true_time1 = matrix(rexp(n_1, rate=hazard1), ncol=1)
true_time = rbind(true_time0, true_time1)
# 打ち切り分布の生成
censoring_time = matrix(rexp(n, rate=1/5), ncol=1) 

# 確認のため、SASで作成した乱数を読み込む
#df4 = read.csv('sandwich_test.csv')
#true_time = df4$true_time
#censoring_time = df4$censoring_time

# 打ち切りフラグ
d = as.integer(true_time <= censoring_time)
# 観測時間
t = ifelse(d == 1, true_time, censoring_time)
print(d)
print(t)

#' 
#' ### ロジスティク回帰
## ---------------------------------------------------------------------------------------------------
intercept = matrix(1, n, 1)
a = matrix(c(rep(0, n_0), rep(1, n_1)), ncol = 1)
x = cbind(intercept, a)
y = as.integer(true_time <= tau)  

beta_0 = matrix(c(0, 0), 2, 1)
result = optim(par=beta_0, fn=score_equation, method="BFGS")

beta = result$par

x_tmp = matrix(c(1, 1, 0, 1), nrow = 2, ncol = 2)
xbeta = x_tmp %*% beta 

p0_hat_p1_hat = exp(xbeta) / (1 + exp(xbeta)) 
print(p0_hat_p1_hat)

# 確認
p0_hat = sum(y[a == 0]) / sum(a == 0) 
p1_hat = sum(y[a == 1]) / sum(a == 1) 
print(p0_hat)
print(p1_hat)

#' 
#' ### IPCW推定量
## ---------------------------------------------------------------------------------------------------
IP_WEIGHT = function(t, d, tau){
  n = length(t)
  one = rep(1, n)
  
  # 閾値条件と打ち切りフラグ
  use = as.integer(t <= tau)
  censoring = one-d
  
  # Kaplan-Meier推定量
  km_0 = kaplan_meier(t, censoring)
  #print(km_0)
  
  # 最小値の計算
  usekm_0 = use*km_0 + 2*(one-use)
  #print(usekm_0)
  minkm_0 = min(usekm_0)
  #print(minkm_0)
  
  # Kaplan-Meier推定量の逆数の計算
  tmp1 = ifelse(km_0 > 0, 1 / km_0, 0)
  tmp2 = if (minkm_0 == 0) 0 else 1 / minkm_0
  
  # 重みの計算
  tmp3 = (t <= tau) * (1 - censoring) * tmp1
  tmp4 = (t > tau) * tmp2
  ip_weight = tmp3 + tmp4
  
  return(ip_weight)
}

ip_weight = IP_WEIGHT(t, d, tau)
print(ip_weight)

#' 
## ---------------------------------------------------------------------------------------------------
# 簡単な例で確認する
t0 = c(2, 3, 6, 8, 10) # 観測時間
d0 = c(1, 0, 1, 0, 1) # イベント発生フラグ
tau0 =  5 # 閾値
weight = IP_WEIGHT(t0, d0, tau0)
# km_0 = c(1.000 0.750 0.750 0.375 0.375)
# minkm_0 = 0.750, 1/minkm_0  = 1.333
print(weight)

# これはSASと一致

#' 
#' 
#' ### 重み付きスコア方程式
## ---------------------------------------------------------------------------------------------------
# x, y, ip_weights
ipcw_equation = function(beta){
  xbeta = x %*% beta
  conditional_p = exp(xbeta) / (1 + exp(xbeta))
  w = 1 / (conditional_p * (1 - conditional_p)) 
  wy = ip_weight * y_censored  
  zero = t(x) %*% (w * (wy - conditional_p)) 
  return(sum(zero^2))
}

beta_0 = matrix(c(0, -1), 2, 1)
y_censored = ifelse(t <= tau & d == 1, 1, 0)
ip_weight = IP_WEIGHT(t, d, tau)
result = optim(par=beta_0, fn=ipcw_equation, method="BFGS")

beta = result$par

x_tmp = matrix(c(1, 1, 0, 1), nrow = 2, ncol = 2)
xbeta = x_tmp %*% beta 

p0_hat_p1_hat = exp(xbeta) / (1 + exp(xbeta)) 
print(p0_hat_p1_hat)

#' 
#' ### IPCW推定を実行する関数
#' 
## ---------------------------------------------------------------------------------------------------
fit_ipcw = function(t, d, tau, x){
  beta_0 = matrix(c(0, -1), 2, 1)
  y_censored = ifelse(t <= tau & d == 1, 1, 0)
  ip_weight = IP_WEIGHT(t, d, tau)
  # x, y, ip_weights
  ipcw_equation = function(beta){
    xbeta = x %*% beta
    conditional_p = exp(xbeta) / (1 + exp(xbeta))
    w = 1 / (conditional_p * (1 - conditional_p)) 
    wy = ip_weight * y_censored  
    zero = t(x) %*% (w * (wy - conditional_p)) 
    return(sum(zero^2))
  }
  result = optim(par=beta_0, fn=ipcw_equation, method="BFGS")
  beta = result$par

  x_tmp = matrix(c(1, 1, 0, 1), nrow = 2, ncol = 2)
  xbeta = x_tmp %*% beta 
  
  p0_hat_p1_hat = exp(xbeta) / (1 + exp(xbeta)) 
  return(p0_hat_p1_hat)
}

fit = fit_ipcw(t, d, tau, x)
print(fit)

#' 
#' ### サンドイッチ分散
#' 
#' #### スコア分散
#' 
## ---------------------------------------------------------------------------------------------------
computing_score = function(beta, x, y_censored, ip_weight){
  xbeta = x %*% beta
  conditional_p = exp(xbeta) / (1 + exp(xbeta))
  w = 1 / (conditional_p * (1 - conditional_p)) 
  wy = ip_weight * y_censored  
  # スコア行列の計算
  score = matrix(0, nrow = nrow(x), ncol = ncol(x)) 
  for (i in seq_len(ncol(x))) {
    score[, i] <- x[, i] * (w * (wy - conditional_p))
  }
  return(score) 
}

score = computing_score(beta, x, y_censored, ip_weight)
print(score)

#' 
#' #### 打ち切り過程のマルチンゲール第1項
#' #### 打ち切り過程のマルチンゲール第2項
#' 
## ---------------------------------------------------------------------------------------------------
computing_cov_score <- function(score, t, d, x, tau){
  n = length(t)
  one = rep(1, n)
  # 打ち切りフラグ
  censoring = one-d
  
  # リスクセット行列
  t = as.vector(t)
  t1 = matrix(t, nrow = n, ncol = n, byrow = TRUE)    # 行列形式のt, tが列・行いずれであっても行方向に並べる
  t2 = t(t1)                            # 転置行列
  atrisk = (t2 >= t1)                   # リスクセット
  atrisk_ = (t2 <= t1)                  # リスクセットの反対
  
  # マルチンゲール第1項
  cm_1 = diag(censoring)
  
  # マルチンゲール第2項
  censoring_dna = nelson_aalen(t, censoring) # ネルソンアーレン推定量
  cm_2 = atrisk * matrix(censoring_dna, nrow = n, ncol = n, byrow = TRUE)
  censoring_martingale = cm_1 - cm_2
  
  # Kaplan-Meier推定量 (censoringに対する)
  censoring_km = kaplan_meier(t, censoring)
  censoring_km[censoring_km == 0] = 0.0001 # # 0回避
  
  censoring_mkm = censoring_martingale / matrix(censoring_km, nrow = n, ncol = n, byrow = TRUE)
  
  # 生存関数(イベントに対する)
  survival_km = kaplan_meier(t, d)
  survival_km[survival_km == 0] <- 0.0001 # 0回避
  
  si = matrix(1, nrow = n, ncol = ncol(x))
  
  for (i in seq_len(ncol(x))) {
    tmp1 = t(x[, i]) %*% (atrisk * d)
    tmp2 = tmp1 / survival_km / n
    tmp2 = matrix(tmp2, nrow = n, ncol = n, byrow = TRUE) 
    integrand = tmp2 * censoring_mkm
    
    for (j in seq_len(n)) {
      use = (t <= tau)
      cusum = cumsum(integrand[j, ])
      usecusum = use * cusum
      integral = max(usecusum)
      si[j, i] = score[j, i] + integral
    }
  }
  
  cov_score = t(si) %*% si / n
  return(cov_score)
}

#' 
#' #### ヘシアン行列
#' 
## ---------------------------------------------------------------------------------------------------
computing_hesse = function(beta, x, ip_weight){
  n = nrow(x)
  xbeta = x %*% beta
  conditional_p = exp(xbeta) / (1 + exp(xbeta))
  w = ip_weight / (conditional_p * (1 - conditional_p))
  w = matrix(w, nrow = n, ncol = ncol(x)) 
  hesse = t(x) %*% (w * x) / n
  return(hesse)
}

#' 
## ---------------------------------------------------------------------------------------------------
sandwich_variance = function(beta, x, y_censored, ip_weight, t, d, tau){
  # スコア行列の計算
  score = computing_score(beta, x, y_censored, ip_weight)
  # スコアの共分散行列の計算
  cov_score = computing_cov_score(score, t, d, x, tau)
  # ヘシアン行列の計算
  hesse = computing_hesse(beta, x, ip_weight)
  # サンドイッチ分散推定量の計算
  hesse_inv = solve(hesse)                   # ヘシアン行列の逆行列
  n = nrow(x)                                # サンプル数
  print(dim(hesse_inv))
  print(dim(cov_score))
  print(dim(hesse_inv))
  cov_beta = hesse_inv %*% cov_score %*% hesse_inv / n
  
  # 結果をリストとして返す
  list(
    cov_score = cov_score,
    hesse = hesse,
    cov_beta = cov_beta
  )
}

print(sandwich_variance(beta, x, y_censored, ip_weight, t, d, tau))

## ---------------------------------------------------------------------------------------------------
# 関数にまとめる
fit_ipcw = function(t, d, tau, x){
  n = nrow(x)
  beta_0 = matrix(c(0, -1), 2, 1)
  y_censored = ifelse(t <= tau & d == 1, 1, 0)
  
  IP_WEIGHT = function(t, d, tau){
    n = length(t)
    one = rep(1, n)
    
    # 閾値条件と打ち切りフラグ
    use = as.integer(t <= tau)
    censoring = one-d
    
    # Kaplan-Meier推定量
    km_0 = kaplan_meier(t, censoring)
    #print(km_0)
    
    # 最小値の計算
    usekm_0 = use*km_0 + 2*(one-use)
    #print(usekm_0)
    minkm_0 = min(usekm_0)
    #print(minkm_0)
    
    # Kaplan-Meier推定量の逆数の計算
    tmp1 = ifelse(km_0 > 0, 1 / km_0, 0)
    tmp2 = if (minkm_0 == 0) 0 else 1 / minkm_0
    
    # 重みの計算
    tmp3 = (t <= tau) * (1 - censoring) * tmp1
    tmp4 = (t > tau) * tmp2
    ip_weight = tmp3 + tmp4
    
    return(ip_weight)
  }
  ip_weight = IP_WEIGHT(t, d, tau)
  
  # x, y, ip_weights
  ipcw_equation = function(beta){
    xbeta = x %*% beta
    conditional_p = exp(xbeta) / (1 + exp(xbeta))
    w = 1 / (conditional_p * (1 - conditional_p)) 
    wy = ip_weight * y_censored  
    zero = t(x) %*% (w * (wy - conditional_p)) 
    return(sum(zero^2))
  }
  result = optim(par=beta_0, fn=ipcw_equation, method="BFGS")
  beta = result$par
  
  x_tmp = matrix(c(1, 1, 0, 1), nrow = 2, ncol = 2)
  xbeta = x_tmp %*% beta 
  
  p0_hat_p1_hat = exp(xbeta) / (1 + exp(xbeta)) 
  
  computing_score = function(beta, x, y_censored, ip_weight){
    xbeta = x %*% beta
    conditional_p = exp(xbeta) / (1 + exp(xbeta))
    w = 1 / (conditional_p * (1 - conditional_p)) 
    wy = ip_weight * y_censored  
    # スコア行列の計算
    score = matrix(0, nrow = nrow(x), ncol = ncol(x)) 
    for (i in seq_len(ncol(x))) {
      score[, i] <- x[, i] * (w * (wy - conditional_p))
    }
    return(score) 
  }
  
  computing_cov_score <- function(score, t, d, x, tau){
    n = length(t)
    one = rep(1, n)
    # 打ち切りフラグ
    censoring = one-d
    
    # リスクセット行列
    t = as.vector(t)
    t1 = matrix(t, nrow = n, ncol = n, byrow = TRUE)    # 行列形式のt, tが列・行いずれであっても行方向に並べる
    t2 = t(t1)                            # 転置行列
    atrisk = (t2 >= t1)                   # リスクセット
    atrisk_ = (t2 <= t1)                  # リスクセットの反対
    
    # マルチンゲール第1項
    cm_1 = diag(censoring)
    
    # マルチンゲール第2項
    censoring_dna = nelson_aalen(t, censoring) # ネルソンアーレン推定量
    cm_2 = atrisk * matrix(censoring_dna, nrow = n, ncol = n, byrow = TRUE)
    censoring_martingale = cm_1 - cm_2
    
    # Kaplan-Meier推定量 (censoringに対する)
    censoring_km = kaplan_meier(t, censoring)
    censoring_km[censoring_km == 0] = 0.0001 # # 0回避
    
    censoring_mkm = censoring_martingale / matrix(censoring_km, nrow = n, ncol = n, byrow = TRUE)
    
    # 生存関数(イベントに対する)
    survival_km = kaplan_meier(t, d)
    survival_km[survival_km == 0] <- 0.0001 # 0回避
    
    si = matrix(1, nrow = n, ncol = ncol(x))
    
    for (i in seq_len(ncol(x))) {
      tmp1 = t(x[, i]) %*% (atrisk * d)
      tmp2 = tmp1 / survival_km / n
      tmp2 = matrix(tmp2, nrow = n, ncol = n, byrow = TRUE) 
      integrand = tmp2 * censoring_mkm
      
      for (j in seq_len(n)) {
        use = (t <= tau)
        cusum = cumsum(integrand[j, ])
        usecusum = use * cusum
        integral = max(usecusum)
        si[j, i] = score[j, i] + integral
      }
    }
    
    cov_score = t(si) %*% si / n
    return(cov_score)
  }
  
  computing_hesse = function(beta, x, ip_weight){
    n = nrow(x)
    xbeta = x %*% beta
    conditional_p = exp(xbeta) / (1 + exp(xbeta))
    w = ip_weight / (conditional_p * (1 - conditional_p))
    w = matrix(w, nrow = n, ncol = ncol(x)) 
    hesse = t(x) %*% (w * x) / n
    return(hesse)
  }
  # スコア行列の計算
  score = computing_score(beta, x, y_censored, ip_weight)
  # スコアの共分散行列の計算
  cov_score = computing_cov_score(score, t, d, x, tau)
  # ヘシアン行列の計算
  hesse = computing_hesse(beta, x, ip_weight)
  # サンドイッチ分散推定量の計算
  hesse_inv = solve(hesse) 
  cov_beta = hesse_inv %*% cov_score %*% hesse_inv / n
  std_err = sqrt(diag(cov_beta))
  
  result_list <- list(
    coefficients = setNames(beta, colnames(x)),
    std_error = setNames(std_err, colnames(x)),
    n = n
  )
  class(result_list) <- "ipcw_fit"
  return(result_list)
}

summary.ipcw_fit <- function(object){
  cat("Inverse Probability of Censoring Weighted (IPCW) Regression\n")
  cat("--------------------------------------------------------\n")
  cat("Number of observations:", object$n, "\n\n")
  
  coef_table <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = object$std_error,
    row.names = names(object$coefficients)
  )
  print(coef_table)
}

fit = fit_ipcw(t, d, tau, x)
summary(fit)

## ---------------------------------------------------------------------------------------------------
B = 1000
p0 = rep(0, B)
p1 = rep(0, B)

tau = 3

for (i in 1:B){
  # 指数乱数の発生
  true_time0 = matrix(rexp(n_0, rate=hazard0), ncol=1) 
  true_time1 = matrix(rexp(n_1, rate=hazard1), ncol=1)
  true_time = rbind(true_time0, true_time1)
  # 打ち切り分布の生成
  censoring_time = matrix(rexp(n, rate=1/5), ncol=1) 
  # 打ち切りフラグ
  d = as.integer(true_time <= censoring_time)
  # 観測時間
  t = ifelse(d == 1, true_time, censoring_time)
  
  fit = fit_ipcw(t, d, tau, x)
  beta = fit$coefficients
  x_tmp = matrix(c(1, 1, 0, 1), nrow = 2, ncol = 2)
  xbeta = x_tmp %*% beta 
  
  p0_hat_p1_hat = exp(xbeta) / (1 + exp(xbeta)) 
  p0[i] = p0_hat_p1_hat[1, 1]
  p1[i] = p0_hat_p1_hat[2, 1]
}

#' 
## ---------------------------------------------------------------------------------------------------
print(mean(p0))
print(mean(p1))
hist(p0)
hist(p1)


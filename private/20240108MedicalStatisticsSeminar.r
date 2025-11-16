library(mets)
data(melanoma)
c1 <- phreg(Surv(days,status!=2)~+ulc+thick+sex,melanoma)
summary(c1)
gof(c1)

library(prodlim)
library(mets)
data(bmt)
a <- prodlim(Hist(time,cause)~+tcell+platelet,data=bmt)
par(mfrow=c(1,2))
plot(a,ylim=c(0,0.5),confint=FALSE,atrisk=FALSE,main="TRM")



library(mets)
data(bmt)
cif1 <- cifreg(Event(time,cause)~tcell+platelet+age,bmt,propodds=NULL,cause=1)
summary(cif1)


set.seed(2)
ss <- rexp(100)
gg <- factor(sample(1:3,100,replace=TRUE),1:3,c('a','b','c'))
cc <- sample(0:2,100,replace=TRUE)
strt <- sample(1:2,100,replace=TRUE)
print(xx <- cuminc(ss,cc,gg,strt))
plot(xx,lty=1,color=1:6)


# 必要なパッケージの読み込み
library(survival)
library(survminer)

# データの読み込み
data(lung)

# サバイバルオブジェクトの作成
fit <- survfit(Surv(time, status) ~ sex, data = lung)

# Kaplan-Meier曲線のプロット
ggsurvplot(
  fit, 
  data = lung,
  pval = TRUE,             # p値を表示
  conf.int = TRUE,         # 信頼区間を表示
  risk.table = TRUE,       # リスクテーブルを表示
  risk.table.col = "strata", # リスクテーブルの色を設定
  linetype = "strata",     # 線の種類を設定
  surv.median.line = "hv", # 中央生存時間の線を表示
  ggtheme = theme_minimal(), # テーマを設定
  palette = c("#E7B800", "#2E9FDF") # カラーパレットを設定
)

fit2 <- cif(Surv(time, status) ~ sex, data = lung)


# 必要なパッケージの読み込み
library(cmprsk)
library(survminer)

# データの読み込み
data(lung)

# イベントの種類を設定（ここでは1をイベント、0を検閲とします）
lung$status <- ifelse(lung$status == 2, 1, 0)

# CIFの計算
cif_fit <- cuminc(ftime = lung$time, fstatus = lung$status, group = lung$sex)

# 累積発生率曲線のプロット
ggcompetingrisks(
  cif_fit,
  conf.int = TRUE,          # 信頼区間を表示
  ggtheme = theme_minimal(), # テーマを設定
  palette = c("#E7B800", "#2E9FDF"), # カラーパレットを設定
  xlab = "Time",            # x軸のラベルを設定
  ylab = "Cumulative Incidence", # y軸のラベルを設定
  title = "Cumulative Incidence Function by Sex" # タイトルを設定
)



# イベントの種類を設定（ここでは1をイベント、0を検閲とします）
lung$status <- ifelse(lung$status == 2, 1, 0)

# CIFの計算
cif_fit <- cuminc(ftime = lung$time, fstatus = lung$status, group = lung$sex)

# 累積発生率曲線のプロット
ggcompetingrisks(
  cif_fit,
  conf.int = TRUE,          # 信頼区間を表示
  ggtheme = theme_minimal(), # テーマを設定
  palette = c("#E7B800", "#2E9FDF"), # カラーパレットを設定
  xlab = "Time",            # x軸のラベルを設定
  ylab = "Cumulative Incidence", # y軸のラベルを設定
  title = "Cumulative Incidence Function by Sex" # タイトルを設定
)



# 必要なパッケージの読み込み
library(cmprsk)
library(survminer)

# データの読み込み
data(lung)

# イベントの種類を設定（ここでは1をイベント、0を検閲とします）
lung$status <- ifelse(lung$status == 2, 1, 0)

# CIFの計算
cif_fit <- cuminc(ftime = lung$time, fstatus = lung$status, group = lung$sex)

# 累積発生率曲線のプロット
ggcompetingrisks(
  cif_fit,
  conf.int = TRUE,          # 信頼区間を表示
  ggtheme = theme_minimal(), # テーマを設定
  palette = c("#E7B800", "#2E9FDF"), # カラーパレットを設定
  xlab = "Time",            # x軸のラベルを設定
  ylab = "Cumulative Incidence", # y軸のラベルを設定
  title = "Cumulative Incidence Function by Sex" # タイトルを設定
)

library(survival)
data(lung, package = "survival")

library(cmprsk)
library(ggplot2)
library(dplyr)

# CIFの計算
cif_fit <- cuminc(ftime = lung$time, fstatus = lung$status, group = lung$sex)

# データフレームに変換
cif_data <- as.data.frame(cif_fit)

# データの整形
cif_data <- cif_data %>% mutate(group = factor(ifelse(grepl("1", rownames(cif_data)), "Male", "Female"))) %>% rename(time = time, cif = est)


# 累積発生率曲線のプロット
ggplot(cif_data, aes(x = time, y = cif, color = group)) +
  geom_step() +
  labs(
    title = "Cumulative Incidence Function by Sex",
    x = "Time",
    y = "Cumulative Incidence",
    color = "Sex"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#E7B800", "#2E9FDF"))
###########################################################################################
# このRスクリプトは生物統計家がRパッケージ開発を学ぶための資料である
# RStudioの環境設定については, RiffomonasプロジェクトのYouTube動画を参考にしている
# https://www.youtube.com/watch?v=XjolVT16YNw


###########################################################################################
# 【1】RStudioの環境設定

#------------------------------------------------------------------------------------------
#-ステップ1. 事前準備----------------------------------------------------------------------

# GitHubの個人アカウント作成（アカウント名とメールアドレス必須）
# 詳細は以下のURL参照
# https://docs.github.com/ja/get-started/start-your-journey/creating-an-account-on-github

#------------------------------------------------------------------------------------------
#-ステップ2. 環境設定----------------------------------------------------------------------

# RStudioへのdevtoolsのインストール
#-- 以下のコードを実行
install.packages("devtools")
library(devtools)

#------------------------------------------------------------------------------------------
# Rパッケージ開発のためのフォルダ作成<path of R project>・ドキュメントの自動生成
#-- 以下のコードを実行
create_package("<path of R project>")

#------------------------------------------------------------------------------------------
# RStudioにおけるGit環境設定
#-- 以下のコードを実行
#-- エラーが生じた場合はGitインストールやGitHubメードアドレス・アカウント名登録を行う
library(devtools)
use_git()

#-- GitHubメードアドレス・アカウント名の登録
#system("git config --global user.name 'Your Name'")
#system("git config --global user.email 'youremail@example.com'")
#system("git config user.name 'Your Name'")
#system("git config user.email 'youremail@example.com'")

#------------------------------------------------------------------------------------------
# Rパッケージのドキュメント編集
# DESCRIPTIONファイルの編集（RStudio右下Filesから開ける）
# MIT + file LICENSEの設定
#-- 以下のコードを実行
use_mit_license()

#------------------------------------------------------------------------------------------
#-ステップ3. プログラミング----------------------------------------------------------------

# これ以降は, Rスクリプトファイルに関数をプログラミングする段階
# 関数のプログラミング→test関数による動作確認→関数の更新の繰り返し

#------------------------------------------------------------------------------------------
# Rスクリプト全体の実行
#-- 以下のコードを実行
load_all()

#------------------------------------------------------------------------------------------
# Roxygenの追加（RStudio右上CodeメニューInsert Roxygen Skeletonから実行）

#------------------------------------------------------------------------------------------
# NAMESPACEやRdファイルなどドキュメントの作成・更新（roxygenなどから情報を取得）
#-- 以下のコードを実行
document()

#------------------------------------------------------------------------------------------
# テスト環境のフォルダ作成・ドキュメントの自動生成
use_testthat()

# test_that関数の作成（RStudio右下Filesのフォルダtestsにスクリプトが生成される）
use_test("<name of R function>")

# test_that関数を用いた自作関数のテスト
test()

#------------------------------------------------------------------------------------------
# GitHubとの連携
#-- 以下のコードを実行
use_github()

#------------------------------------------------------------------------------------------
# バージョン管理・バックアップ・パッケージ公開はGitHubが利用可能
# 必要に応じてRStudio右上Gitタブでcommit, push, pullを行う

#------------------------------------------------------------------------------------------
#-ステップ4. ドキュメント作成・パッケージビルド・公開--------------------------------------

#------------------------------------------------------------------------------------------
# README.Rmd作成のための環境設定
#-- 以下のコードを実行
use_readme_rmd()

# RマークダウンファイルREADME.Rmdの作成・更新
#-- 以下のコードを実行
build_readme()

#------------------------------------------------------------------------------------------
# パッケージのバージョン管理はDESCRIPTIONファイル編集により行う

#------------------------------------------------------------------------------------------
# チェック
check()


###########################################################################################
# 【2】Rの行列・データフレーム・関数の使い方

#------------------------------------------------------------------------------------------
# ステップ1. 行列の次元と積

a <- 1:10
b <- 1:10
print(a)          # 1から10までの値を持つベクトルができた
print(a*b)        # aとbを掛け算すると1から100までの値を持つベクトルに
print(length(a))  # aの長さは10

A <- as.matrix(a) # ベクトルaを行列Aに変換した
B <- as.matrix(b)
print(dim(A))     # 関数dimで次元を確認すると10×1行列

print(A*b)        # 行列Aとベクトルbの「要素ごとの積」
print(A*B)        # 行列Aと行列Bの「要素ごとの積」

AA <- cbind(a, a)
print(AA*b)        # 行列AAとベクトルbの「要素ごとの積」
c <- 1:2
C <- t(as.matrix(c))
print(AA*c)        # 行列AAとベクトルCの「要素ごとの積」

print(A %*% B)    # 行列の積は%*%で計算する（右の式は行列の次元があわないから計算できない）
print(t(A) %*% B) # 行列Aの転置行列と行列Bの積

#------------------------------------------------------------------------------------------
# ステップ2. ベクトル・行列の特殊計算

print(cumsum(a))  # 関数cumsumでベクトルaの累積の和を計算
print(cumprod(a)) # 関数cumprodでベクトルaの累積の積を計算
?cumsum           # 知らない関数がでてきたら

C <- matrix(1:9, nrow = 3, byrow = TRUE)
print(C)
print(C[2,2])     # 行列の要素
print(sum(C))     # 行列の和
print(rowSums(C)) # 行列の行の和
print(colSums(C)) # 行列の列の和

#------------------------------------------------------------------------------------------
# ステップ3. データフレーム
# データフレームは行列と異なり, 数値とテキストを混在することができ, 各行・各列がラベルを持っている

example <- data.frame(age = c(45, 32, 20, 25), sex = c("M", "F", "F", "M"), id = 1:4)
print(example)
print(as.matrix(example))
print(as.numeric(as.matrix(example)))

# データフレームXの列iはX[i]で参照できるが, ラベル（変数名）xということがわかればX[["x"]]やX$xと書くこともできる
print(example$age)
print(example[1])
print(example[["age"]])

# データフレームの行数や列数はnrow関数やncol関数で取得できる
print(nrow(example))
print(ncol(example))

# データフレームの行名や列名はrownames関数やcolnames関数で取得できる
print(rownames(example))
print(colnames(example))

# データフレームの行名や列名を変更するときはrownames関数やcolnames関数を使って変更する
colnames(example) <- c("age_in_year", "sex_in_MF", "id")
print(example)

# order(): 数値またはテキストから順位情報を読み取る関数
print(order(example$age))
print(order(example$sex))

# ある要素の順位情報に従ってデータフレームの順番を変えたいときは, order関数またはsort関数が利用できる
sorted_data <- example[order(example$age), ]
print(sorted_data)

#------------------------------------------------------------------------------------------
# ステップ4. 関数の自作
# 以下のコードを実行すると与えられた生存時間データからアットリスク数を計算する関数が作成される

calculateNumberAtRisk <- function(time) {
  atrisk <- outer(time, time, "<=")          # なぜこの関数でアットリスク数が得られるかは時間があるときに考えよ
  n_atrisk <- rowSums(atrisk)
  return(n_atrisk)
}

lambda <- 1                                  # ハザードを設定
n <- 10
survival_time <- rexp(n, rate = lambda)      # 指数分布から生存時間を発生
print(survival_time)
n_atrisk <- calculateNumberAtRisk(survival_time)
print(n_atrisk)                              # アットリスク数が正しいか確認すること

# 【参考】関数について情報を得るにはgetAnywhere関数と?を用いる
getAnywhere(outer)
?outer

#------------------------------------------------------------------------------------------
# order関数, outer関数, cumprod関数を組み合わせると生存時間データからのKaplan-Meier推定量を計算できる

calculateKaplanMeier1 <- function(t, d){
  n <- length(t)
  data <- data.frame(t = t, d = d, id = 1:n)
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_d <- sorted_data$d
  sorted_id <- sorted_data$id
  atrisk <- outer(sorted_t, sorted_t, "<=") # なぜこの関数でアットリスク数が得られるかは時間があるときに考えよ
  n_atrisk <- rowSums(atrisk)
  s <- 1 - sorted_d / n_atrisk
  km <- cumprod(s)
  data <- data.frame(id=sorted_id, km=km)
  sorted_data <- data[order(data$id), ]
  km <- sorted_data$km
  return(km)
}

survival_time <- rexp(N, rate = 0.2)     # 指数分布から生存時間を発生
censoring_time <- rexp(N, rate = 0.5)       # 指数分布から打ち切り時間を発生
event <- as.numeric(survival_time<=censoring_time)
observed_time <- event*survival_time + (1-event)*censoring_time

km <- calculateKaplanMeier1(observed_time, event)
print(observed_time)
print(event)
print(km)


###########################################################################################
# 【3】関数の高速化

#------------------------------------------------------------------------------------------
# Rで行われる計算のほとんどは行列演算であり, プログラミング上も,
# ベクトルや行列の線型演算として書いた方が, 可読性・速度が高い
# 一方で, べき乗や微分積分など線型ではない計算をプログラミングするときは
# forループなどを用いる必要がある
# ベクトルXの要素iはX[i]により参照できる
# 以下のように, iについてループを回すと, ベクトルの要素ごとの計算ができる

x <- 1:10
n <- length(x)
cum1 <- 0
for(i in 1:n){
  cum1 <- cum1 + x[i]
}
print(cum1)

# forループを書くときには, 最低限でもループの最初・最後について, 挙動を確認しておくべきである
# 上のような簡単なアルゴリズムですら, たとえばcumが定義されていないと, 最初の計算がうまくいかない

cum2 <- NULL
for(i in 1:n){
  cum2 <- cum2 + x[i]
}
print(cum2)

# 最後のループも要確認である. 以下の例では, インデックスの終点がベクトルの範囲を超えてしまっている

n <- length(x)+1
cum3 <- 0
for(i in 1:n){
  cum3 <- cum3 + x[i]
}
print(cum3)

# 【参考】同じように計算がうまいいかない場合でも, 得られる結果はケースバイケースである
# 上の例では, cum1とcum2ではクラスも保持する値も異なる
class(cum2)
class(cum3)
print(is.na(cum2))
print(is.na(cum3))


#------------------------------------------------------------------------------------------
# すでに述べたように, Rのforループは計算速度が遅く, 既存の関数を利用した方が速い
# 関数の速度を比較するためにはmicrobenchmark関数を利用できる

calculateCum <- function(x) {
  cum <- 0
  for(i in 1:length(x)){
    cum <- cum + x[i]
  }
  return(cum)
}

library(microbenchmark)
x <- 1:100
microbenchmark(calculateCum(x), cumsum(x), times=200)


#------------------------------------------------------------------------------------------
# 例題として, calculateKaplanMeier1関数を高速化してみよう
# sapply関数を用いると可読性が下がるがコンパクトにプログラムを書ける

calculateKaplanMeier2 <- function(t, d){
  unique_t <- sort(unique(t))
  n_risk <- sapply(unique_t, function(s) sum(t >= s))
  n_event <- sapply(unique_t, function(s) sum(d[t == s]))
  km <- cumprod(1 - n_event / n_risk)
  return(km)
}

#------------------------------------------------------------------------------------------
# 似たようなコードでも使用する関数によって速度は変わる
# 以下のコードではパッケージRfastのsort_uniqueを利用

library(Rfast)
calculateKaplanMeier3 <- function(t, d){
  unique_t <- sort_unique(t)
  n_risk <- sapply(unique_t, function(s) sum(t >= s))
  n_event <- sapply(unique_t, function(s) sum(d[t == s]))
  km <- cumprod(1 - n_event / n_risk)
  return(km)
}

#------------------------------------------------------------------------------------------
# Rで計算するよりC++で計算する方が速い
# Rでは複数の方法でC++を利用可能であり, パッケージRccpもそのひとつ
# パッケージRccpでは, use_rcpp()を用いて生成されるcppファイルに関数を書くか方法と
# Rファイルに直接cppFunctionを定義する方法の2通りがある
# Rcpp::compileAttributes() により必要に応じてRcppExports.* を再生成


library(Rcpp)
cppFunction('
NumericVector calculateKaplanMeier4(NumericVector t, NumericVector d) {
  int n = t.size();
  NumericVector km(n);
  IntegerVector indices = seq(0, n-1);
  std::sort(indices.begin(), indices.end(), [&](int i, int j) {
    return t[i] < t[j];
  });

  NumericVector sorted_t(n);
  NumericVector sorted_d(n);
  IntegerVector sorted_id(n);
  for (int i = 0; i < n; i++) {
    sorted_t[i] = t[indices[i]];
    sorted_d[i] = d[indices[i]];
    sorted_id[i] = indices[i] + 1;
  }

  int n_atrisk = n;
  double km_i = 1.0;
  for (int i = 0; i < n; i++) {
    if (i > 0 && sorted_t[i] != sorted_t[i-1]) {
      n_atrisk = n - i;
    }
    km_i *= (1 - sorted_d[i] / (double)n_atrisk);
    km[i] = km_i;
  }
  return km;
}
')

# 【参考】C++の使用
# C++の文法にはRとは異なる点がいくつかある.
# たとえば変数を宣言するタイミングを厳密に考える必要があるし
# 配列のインデックスもゼロスタートである（X[0]）.

# パッケージRcpp特有のC++使用ルールについては以下のサイトが参考になる
# https://teuder.github.io/rcpp4everyone_ja/

# Rスクリプトではなくcppファイルに関数を書く方法については以下のYouTube動画を参照のこと
# https://www.youtube.com/watch?v=ZzNwRNbD1gU
# 以下のコードを実行するとcppファイルが生成される
# use_rcpp("<name of function>")

#------------------------------------------------------------------------------------------
# microbenchmark関数を用いて計算速度を比較すると以下の結果になる
# このケースではsurvfit関数が遅いように見えるが, 気のせいであり, 調べるとわかるがsurvfit関数の計算は高度に最適化されている

createTestData <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
  one <- rep(1, n)
  t <- c(1:(n/2), 1:(n/2))
  epsilon <- rep(1, n)
  epsilon[2] <- 2
  epsilon[3] <- 2
  if (first_zero==TRUE) {
    epsilon[1] <- 0
    epsilon[n/2+1] <- 0
  }
  if (last_zero==TRUE) {
    epsilon[n/2] <- 0
    epsilon[n] <- 0
  }
  w <- rep(w, n)
  if (logical_strata==TRUE) {
    strata <- (t %% 2 == 0)
  } else {
    strata <- as.factor((t %% 2 == 0))
  }
  if (na_strata==TRUE) {
    strata[1] <- NA
  }
  subset <- rep(1, n)
  if (subset_present==TRUE) {
    subset[1] <- 0
  }
  d <- as.numeric(epsilon>0)
  return(data.frame(id = 1:n, t = t, epsilon = epsilon, d = d, w = w, strata = strata, subset=subset))
}

library(microbenchmark)
library(survival)
library(ggsurvfit)
df_test <- createTestData(10, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=TRUE)
microbenchmark(calculateKaplanMeier1(df_test$t, df_test$d),
               calculateKaplanMeier2(df_test$t, df_test$d),
               calculateKaplanMeier3(df_test$t, df_test$d),
               calculateKaplanMeier4(df_test$t, df_test$d),
               survfit(Surv(t, d)~1, data=df_test),
               km.curve(Surv(t, d)~1, data=df_test, use.ggsurvfit = FALSE, use.polyreg = FALSE),
               times = 20)

#a <- km.curve(Surv(t, d)~1, data=df_test, use.ggsurvfit = FALSE, use.polyreg = FALSE, error="siatis")
#summary(a)
#df_test <- createTestData(10, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
#df_test$strata_ <- as.numeric(df_test$strata)
#df_test$strata_
#calculateKM_rcpp(t=df_test$t, d=df_test$d, error = "greenwood", strata=df_test$strata_)

#------------------------------------------------------------------------------------------
# 【参考】フォーミュラ型
# survfit関数で用いられているフォーミュラ型は, 「~」の左側と右側に変数名を記述するR特有の表現
# 多くの場合, 左辺はアウトカム, 右辺は共変量
# データフレームのどの変数を用いるかを指定するときに有用

example1 <- y ~ x1+x2
class(example1)
print(example1[1])
print(example1[2])
print(example1[3])

# survfit関数では, フォーミュラ（Surv(observed_time, event)~1）から変数名を読み取っている
# Surv(, )またはEvent(, ): 生存時間アウトカムを表すための関数（パッケージsurvivalまたはmetsに依存するので注意）
# terms(): フォーミュラから変数名を読み取るための関数
# all.vars(): terms()の出力から変数名を読み取るための関数
# ~+1: 定数項

nuisance.model <- Event(t_test, epsilon_test) ~ +1
data <- df_test

lhs <- terms(nuisance.model)
time_var <- all.vars(lhs)[1] # 時間の変数名
event_var <- all.vars(lhs)[2] # イベントの変数名
time <- data[[time_var]] # 時間の情報
event_all_ <- data[[event_var]] # イベントコードの情報
print(time)
print(event_all_)

# 別の読み取り方もある
# 以下のコードではフォーミュラ（nuisance.model）とデータフレーム（df_test）を
# モデルフレーム型（example2）に変換し, データを取得している
# モデルフレームはデータそのもので, モデルフレームを扱うために用意された関数が利用できることも特長
# 一方でフォーミュラは変数名の情報にすぎず, 実在するデータと結びついていない
# model.frame(): フォーミュラとデータフレームを結びつけるための関数
# model.extract(): モデルフレームからデータを読み取るための関数

nuisance.model <- Event(t_test, epsilon_test) ~ +1
special <- NULL
out_terms <- terms(nuisance.model, special, data = df_test)
example2 <- match.call()
example2$formula <- out_terms
example2[[1]] <- as.name("model.frame")
example2 <- eval(example2, parent.frame())
y <- model.extract(example2, "response")
time <- y[, 1] # 時間の情報
event_all_ <- y[, 2] # イベントコードの情報

class(example2)
print(y)
print(time)
print(event_all_)

# model.matrix(): モデルフレームから回帰モデルのデザイン行列を作るための関数
# デザイン行列が特異行列にならないように自動的にダミー変数を作成してくれるが
# デザイン行列のコーディングによって, 回帰係数の解釈が異なるため, 仕様を理解する必要がある
# コーディングは, 対比（contrasts.arg）によって指定できる
# デフォルト: contr.treatment（切片項は1, ダミー変数は0と1でコーディング）

stratified.model <- Event(t_test, epsilon_test) ~ strata_test
special <- NULL
out_terms <- terms(stratified.model, special, data = df_test)
example3 <- match.call()
example3$formula <- out_terms
example3[[1]] <- as.name("model.frame")
example3 <- eval(example3, parent.frame())
design_matrix_t <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.treatment"))
design_matrix_s <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.sum"))
design_matrix_p <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.poly"))
design_matrix_S <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.SAS"))
design_matrix_h <- model.matrix(out_terms, example3, contrasts.arg = list(strata_test = "contr.helmert"))
print(design_matrix_S)

# フォーミュラから変数xを減らすときは, マイナスを用いて-xと指定する
# 以下のフォーミュラは切片項を含まないデザイン行列に対応する
# このデザイン行列を用いて線型モデルを当てはめると, 共変量の水準ごとの平均が推定される

no.intercept.model <- Event(t_test, epsilon_test) ~ -1+strata_test
special <- NULL
out_terms <- terms(no.intercept.model, special, data = df_test)
example4 <- match.call()
example4$formula <- out_terms
example4[[1]] <- as.name("model.frame")
example4 <- eval(example4, parent.frame())
design_matrix_n <- model.matrix(out_terms, example4)
print(design_matrix_n)

#------------------------------------------------------------------------------------------
# 【参考】ggsurvfit関数とkm.curve関数
# ggplot系の生存曲線を描くための関数のひとつにggsurvfit関数がある
# ggsurvfit関数は, km.curve関数やsurvfit関数のアウトプット（survfit型オブジェクト）を自動で図示してくれる

data(diabetes.complications)
diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
out_km <- km.curve(Surv(t, d)~fruitq1, data=diabetes.complications, use.ggsurvfit = FALSE, label.strata=c("High fruit intake", "Low fruit intake"))
#out_km <- survfit(Surv(t, d)~fruitq1, data=diabetes.complications)

# ggplot系の関数は"+"を用いてグラフに要素を追加できる
library(ggsurvfit);
out_ggsurvfit <- ggsurvplot(out_km)
print(out_ggsurvfit)

# テーマを変更
out_ggsurvfit <- ggsurvfit(out_km) +
  theme_classic()
print(out_ggsurvfit)

# フォントを変更
out_ggsurvfit <- ggsurvfit(out_km) +
  theme_classic() +
  labs(x = "Time (months)",
     y = "Risk of diabetic retinopathy") +
  lims(x = c(0, 10),
       y = c(0, 1)) +
  theme(axis.title = element_text(size = 16, family = "serif"),
      axis.text = element_text(size = 14, family = "serif"),
      legend.position = "top",
      legend.text = element_text(size = 14, family = "serif"))
print(out_ggsurvfit)

# 信頼区間と打ち切りを表示
out_ggsurvfit <- ggsurvfit(out_km) +
  add_confidence_interval() +
  add_censor_mark() +
  labs(x = "Time (months)",
       y = "Risk of diabetic retinopathy") +
  lims(x = c(0, 10),
       y = c(0, 1)) +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 10, by = 2)),
                  y_scales = list(breaks = seq(0, 1, by = 0.2)))+
  theme_classic() +
  scale_color_manual(values = c('#54738E', '#82AC7C')) +
  scale_fill_manual(values = c('#54738E', '#82AC7C')) +
  theme(axis.title = element_text(size = 16, family = "serif"),
        axis.text = element_text(size = 14, family = "serif"),
        legend.position = "top",
        legend.text = element_text(size = 14, family = "serif"))
print(out_ggsurvfit)


# km.curve関数は自動で生存曲線を描くことができるが, 内部ではggsurvfit関数にsurvfit型オブジェクトを渡している
out_km <- km.curve(Surv(t, d)~fruitq1, data=diabetes.complications, label.strata=c("High fruit intake", "Low fruit intake"), time.point = c(4,8))
out_km <- km.curve(Surv(t, d)~fruitq1, data=diabetes.complications, label.strata=c("High fruit intake", "Low fruit intake"), conf.type = "none")

out_km <- km.curve(Surv(t, d)~fruitq1, data=diabetes.complications, label.strata=c("High fruit intake", "Low fruit intake"), lims.x = c(0, 9))

#------------------------------------------------------------------------------------------
# 【参考】gtsummaryパッケージ

if (!require("gtsummary")) install.packages("gtsummary")
if (!require("flextable")) install.packages("flextable")
library(gtsummary)
library(flextable)
library(broom)
data(diabetes.complications)

# ベースラインデータの表
table1 <- diabetes.complications |>
  select(fruitq1, age, sex) |>
  tbl_summary(by = fruitq1,
              type = list(c(age) ~ "continuous",
                          c(sex) ~ "categorical"),
              statistic = list(all_continuous() ~ "{mean} ({sd})")) |>
  add_n()
table1

# ベースラインデータの表+p値
get_kruskal.fisher.test <- function(data, variable, by, ...) {
  variable <- data[[variable]]
  by <- data[[by]]
  if (is.numeric(variable)|is.integer(variable)) {
    p <- kruskal.test(variable ~ as.factor(by))$p.value
  } else if (is.factor(variable)|is.logical(variable)) {
    table <- table(variable, by)
    p <- fisher.test(table)$p.value
  }
  if (p<0.05) {
    p_text <- c("<0.05")
  } else
  {
    p_text <- paste0(round(p, digit=2))
  }
  return(p_text)
}

table1_p_value <- diabetes.complications |>
  select(fruitq1, age, sex, sbp) |>
  tbl_summary(by = fruitq1,
              type = list(c(age) ~ "continuous",
                          c(sbp) ~ "continuous",
                          c(sex) ~ "categorical"),
              statistic = list(all_continuous() ~ "{mean} ({sd})")) |>
  add_stat(fns = everything() ~ get_kruskal.fisher.test) |>
  modify_header(add_stat_1 = "**p**", all_stat_cols() ~ "**{level}**") |>
  as_gt() |>
  gt::tab_source_note(gt::md("*Kruskal-Wallis test or Fisher exact test*"))
table1_p_value


# ベースラインデータの表+平均の差
table1_difference <- diabetes.complications |>
  select(fruitq1, age, sex, sbp) |>
  tbl_summary(by = fruitq1,
              type = list(c(age) ~ "continuous",
                          c(sex) ~ "dichotomous",
                          c(sbp) ~ "continuous"),
              statistic = list(all_continuous() ~ "{mean} ({sd})")) |>
  add_difference()
table1_difference

# ベースラインデータの表+平均の差
get_glm <- function(data, variable, by, ...) {
  glm(data[[variable]] ~ as.factor(data[[by]])) |>
    broom::tidy() |>
    tail(1) |>
    select(estimate, p.value)
}
table1_difference <- diabetes.complications |>
  select(fruitq1, age, sbp) |>
  tbl_summary(by = fruitq1,
              type = list(c(age) ~ "continuous",
                          c(sbp) ~ "continuous"),
              statistic = list(all_continuous() ~ "{mean} ({sd})")) |>
  add_stat(fns = everything() ~ get_glm) |>
  modify_header(estimate = "**Difference**", p.value = "**p**") |>
  as_gt() |>
  gt::tab_source_note(gt::md("*Differences are estimated by linear models*"))
table1_difference

# イベント数の集計
library(labelled)
attr(diabetes.complications$t, "label") <- "Follow-up years"
diabetes.complications$d1 <- as.integer(diabetes.complications$epsilon==1)
attr(diabetes.complications$d1, "label") <- "Diabetic retinopathy"
diabetes.complications$d2 <- as.integer(diabetes.complications$epsilon==2)
attr(diabetes.complications$d2, "label") <- "Macrovascular diseases"
table2 <- diabetes.complications |>
  select(fruitq1, d1, d2, t) |>
  tbl_summary(by = fruitq1,
              type = list(c(d1) ~ "dichotomous",
                          c(d2) ~ "dichotomous",
                          c(t) ~ "continuous"),
              statistic = list(all_continuous() ~ "{min}-{max}")) |>
  add_stat(fns = everything() ~ get_kruskal.fisher.test) |>
  modify_header(add_stat_1 = "**p**", all_stat_cols() ~ "**{level}**") |>
  as_gt() |>
  gt::tab_source_note(gt::md("*Kruskal-Wallis test or Fisher exact test*"))
table2

# Kaplan-Meier推定量
#install.packages("cardx")
#library(survival)
out_survfit1 <- survfit(Surv(t, d1)~1, diabetes.complications)
out_survfit2 <- survfit(Surv(t, d1)~sex, diabetes.complications)
out_survfit3 <- survfit(Surv(t, d1)~fruitq1, diabetes.complications)
out_survfit <- list(out_survfit1, out_survfit2, out_survfit3)
table3 <- tbl_survfit(out_survfit, times = c(4, 8), label_header = "**{time} years**") |>
  modify_spanning_header(all_stat_cols() ~ "**Free from diabetic retinopathy**")
table3

# Cox回帰
attr(diabetes.complications$sex, "label") <- "Sex (woman v.s. man)"
attr(diabetes.complications$fruitq1, "label") <- "Fruit intake (low v.s. high)"
out_coxph <- coxph(Surv(t, d1) ~ sex+fruitq1, data=diabetes.complications)
tbl_regression(out_coxph, exponentiate = TRUE)



# Jack knife計算速度の評価

library(jackknifeR)
library(boot)
jknife <- function(x, fxn, ci=0.95) {
  theta <- fxn(x)
  n <- nrow(x)
  partials <- rep(0,n)
  for (i in 1:n){
    partials[i] <- fxn(x[-i,])
  }
  pseudos <- (n*theta) - (n-1)*partials
  jack.est <- mean(pseudos)
  jack.se <- sqrt(var(pseudos)/n)
  alpha = 1-ci
  CI <- qt(alpha/2,n-1,lower.tail=FALSE)*jack.se
  jack.ci <- c(jack.est - CI, jack.est + CI)
  list(est=jack.est, se=jack.se, ci=jack.ci, pseudos = pseudos, partials=partials)
}

fn12 <- function(data) {
  out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
  return(out_km$surv[30])
}

fn3 <- function(data, i) {
  data <- data[i,]
  out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
  return(out_km$surv[30])
}

fn3w <- function(data, w) {
  out_km <- calculateKM_rcpp(data$t, data$d, w, as.integer(data$strata), "none")
  return(out_km$surv[30])
}

#----------------------------------------
data(diabetes.complications)
diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
diabetes.complications$w <- rep(1,nrow(diabetes.complications))
df_test_ <- diabetes.complications[1:200,]
jk1 <- jackknife(statistic = fn12, d = 1, data = df_test, numCores= 2)
print(jk1)

jk2 <- jknife(df_test_, fn12, ci=0.95)
print(jk2$se)

jk3 <- empinf(data=df_test_, statistic=fn3, type="jack", stype="i")
se3 <- sd(jk3)/sqrt(length(jk3))
print(se3)

jk3w <- empinf(data=df_test_, statistic=fn3w, type="jack", stype="w")
se3w <- sd(jk3w)/sqrt(length(jk3w))
print(se3w)

library(ggsurvfit)
out_km <- km.curve(Surv(t, d)~1, data=df_test_)
print(out_km$std.err)
print(out_km$time)
print(df_test_$t[30])



#----------------------------------------
createTestData <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
  one <- rep(1, n)
  t <- c(1:(n/2), 1:(n/2))
  epsilon <- rep(1, n)
  epsilon[2] <- 2
  epsilon[3] <- 2
  if (first_zero==TRUE) {
    epsilon[1] <- 0
    epsilon[n/2+1] <- 0
  }
  if (last_zero==TRUE) {
    epsilon[n/2] <- 0
    epsilon[n] <- 0
  }
  w <- rep(w, n)
  if (logical_strata==TRUE) {
    strata <- (t %% 2 == 0)
  } else {
    strata <- as.factor((t %% 2 == 0))
  }
  if (na_strata==TRUE) {
    strata[1] <- NA
  }
  subset <- rep(1, n)
  if (subset_present==TRUE) {
    subset[1] <- 0
  }
  d <- as.numeric(epsilon>0)
  return(data.frame(id = 1:n, t = t, epsilon = epsilon, d = d, w = w, strata = strata, subset=subset))
}

library(microbenchmark)
df_test <- createTestData(100, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
microbenchmark(jk1 <- jackknife(statistic = fn, d = 1, data = df_test, numCores= 2),
               jk2 <- jknife(df_test, fn, ci=0.95),
               jk3 <- empinf(data=df_test, statistic=fn_, type="jack", stype="i"),
               times = 20)



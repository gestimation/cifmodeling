// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {
  constexpr double EPS = 1e-12;

  inline bool feq(double a, double b, double eps=EPS){
    return std::fabs(a-b) <= eps;
  }

  struct Rec {
    double t;   // time
    int cause;  // 0=censor, 1,2=events (cause 1 used for CIF)
    int g;      // group 1..K
    int id;     // original row index (0..N-1)
  };

  // lower_bound + tolerance match; return index or -1
  inline int find_time_index(const std::vector<double>& grid, double target){
    auto it = std::lower_bound(grid.begin(), grid.end(), target);
    if (it == grid.end()) return -1;
    if (feq(*it, target)) return (int)std::distance(grid.begin(), it);
    if (it != grid.begin()){
      auto itp = it - 1;
      if (feq(*itp, target)) return (int)std::distance(grid.begin(), itp);
    }
    return -1;
  }
}

//' K-group Aalen–Johansen for cause 1 with individual influence functions
//'
//' @param t    numeric vector of event/censoring times (>=0)
//' @param epsilon integer vector: 0=censor, 1=cause1, 2=cause2 (>=2 treated as "any event" but CIF is for cause1)
//' @param strata integer vector of group labels in {1,2,...,K}
//'
//' @return list with
//'   - time:     numeric vector of cause-1 jump times (union across stratas), length m+1 with leading 0
//'   - n.risk:   K x (m+1) matrix of at-risk counts on the time grid
//'   - aj1:      K x (m+1) matrix of Aalen–Johansen CIF for cause 1
//'   - var_aj1:  K x (m+1) matrix of variances (sum of IF^2)
//'   - if_aj1:   list of length K; each is a matrix (n_g x (m+1)) of individual IF trajectories
//'
//' @details
//' Internally builds two time grids (unions across stratas):
//'   TTJP: all event times (cause>=1) with leading 0
//'   TJP:  cause-1 times (cause==1) with leading 0  <-- output grid
//' CIF/IF recursion is done on TTJP per strata, then mapped to TJP via index matching.
//' @examples
//' # res <- calculateIFofAJ(t, epsilon, strata)
// [[Rcpp::export]]
List calculateIFofAJ(NumericVector t, IntegerVector epsilon, IntegerVector strata){
  const int N = t.size();
  if (epsilon.size()!=N || strata.size()!=N) stop("Length mismatch.");
  for (int i=0;i<N;i++){
    if (!R_finite(t[i]) || t[i] < 0.0) stop("t must be finite and nonnegative.");
    if (epsilon[i] < 0) stop("epsilon must be >=0.");
    if (strata[i] < 1) stop("strata must be >=1.");
  }

  // ---- 0) Collect stratas and K
  IntegerVector guniq = sort_unique(strata);
  const int K = guniq.size();
  // Build a compact mapping (label -> 0..K-1)
  std::vector<int> lab2k( *std::max_element(guniq.begin(), guniq.end()) + 1, -1 );
  for (int k=0;k<K;k++) lab2k[guniq[k]] = k;

  // ---- 1) Build records and sort (time asc, event before censor at ties)
  std::vector<Rec> x; x.reserve(N);
  for (int i=0;i<N;i++){
    int gl = strata[i];
    int kg = lab2k[gl];
    if (kg < 0) stop("Internal: strata label mapping failed.");
    x.push_back({(double)t[i], (int)epsilon[i], gl, i});
  }
  std::sort(x.begin(), x.end(), [](const Rec& a, const Rec& b){
    if (!feq(a.t,b.t)) return a.t < b.t;
    // events first (epsilon>=1) then censor (0)
    return a.epsilon > b.epsilon;
  });

  // ---- 2) Build global tJP (any event) and TJP (epsilon1), both with leading 0
  std::vector<double> tJP; tJP.reserve(N+1); tJP.push_back(0.0);
  std::vector<double> TJP;  TJP.reserve(N+1);  TJP.push_back(0.0);
  for (int i=0;i<N;i++){
    if (x[i].epsilon >= 1){
      if (!feq(tJP.back(), x[i].t)) tJP.push_back(x[i].t);
      if (x[i].epsilon == 1){
        if (!feq(TJP.back(), x[i].t)) TJP.push_back(x[i].t);
      }
    }
  }
  const int TNJP = (int)tJP.size() - 1;  // index 1..TNJP
  const int NJP  = (int)TJP.size()  - 1;  // index 1..NJP

  // ---- 3) Per-strata counts on tJP: ny, nd_any, nd_c1
  NumericMatrix ny(K, TNJP+1);         // at-risk per strata at tJP
  NumericMatrix nd_any(K, TNJP+1);
  NumericMatrix nd_c1(K, TNJP+1);

  // at-risk: Y_g(t_j) = #{ i in g : T_i >= t_j }
  for (int j=1;j<=TNJP;j++){
    double tj = tJP[j];
    for (int i=0;i<N;i++){
      if (x[i].t + EPS >= tj){
        int kg = lab2k[x[i].g];
        ny(kg, j) += 1.0;
      }
    }
  }
  // jumps: counts at exact tJP[j]
  for (int j=1;j<=TNJP;j++){
    double tj = tJP[j];
    for (int i=0;i<N;i++){
      if (feq(x[i].t, tj)){
        int kg = lab2k[x[i].g];
        if (x[i].epsilon >= 1) nd_any(kg, j) += 1.0;
        if (x[i].epsilon == 1) nd_c1(kg, j)  += 1.0;
      }
    }
  }

  // ---- 4) Aalen–Johansen on tJP per strata: dlambda, surv, aj1, var via IF
  NumericMatrix dl_any(K, TNJP+1), dl_c1(K, TNJP+1);
  NumericMatrix surv(K, TNJP+1), aj1_tjp(K, TNJP+1), var_tjp(K, TNJP+1);

  // For IF storage per strata: first count strata sizes and collect indices in original order
  std::vector< std::vector<int> > idx_of_strata(K);
  for (int i=0;i<N;i++){
    int kg = lab2k[strata[i]];
    idx_of_strata[kg].push_back(i);
  }
  std::vector<int> n_g(K,0);
  for (int kg=0; kg<K; ++kg) n_g[kg] = (int)idx_of_strata[kg].size();

  // Recursions
  for (int kg=0; kg<K; ++kg){
    surv(kg,0) = 1.0;
    aj1_tjp(kg,0) = 0.0;
    var_tjp(kg,0) = 0.0;

    for (int j=1;j<=TNJP;j++){
      if (ny(kg,j) > 0.0){
        dl_any(kg,j) = nd_any(kg,j) / ny(kg,j);
        dl_c1(kg,j)  = nd_c1(kg,j)  / ny(kg,j);
      } else {
        dl_any(kg,j) = 0.0;
        dl_c1(kg,j)  = 0.0;
      }
      surv(kg,j)     = surv(kg,j-1) * (1.0 - dl_any(kg,j));
      aj1_tjp(kg,j) = aj1_tjp(kg,j-1) + surv(kg,j-1) * dl_c1(kg,j);
    }
  }

  // ---- 5) Influence functions: per individual per strata on tJP, and accumulate var
  // We'll store IF matrices per strata on tJP first (n_g x (TNJP+1)), then map to TJP.
  std::vector< NumericMatrix > IF_tjp(K);
  for (int kg=0; kg<K; ++kg){
    IF_tjp[kg] = NumericMatrix(n_g[kg], TNJP+1);
  }

  for (int kg=0; kg<K; ++kg){
    // For each individual in this strata, compute xs/xf over tJP using strata's processes
    for (int row=0; row<n_g[kg]; ++row){
      int i = idx_of_strata[kg][row];              // original individual index in data order
      // find that individual's record from sorted x: we need t_i and epsilon_i
      double Ti = t[i];
      int Ci = epsilon[i];

      std::vector<double> xs(TNJP+1, 0.0), xf(TNJP+1, 0.0);
      IF_tjp[kg](row,0) = 0.0;

      for (int j=1;j<=TNJP;j++){
        double dm0=0.0, dm1=0.0;

        if (Ti + EPS >= tJP[j]){
          if (feq(Ti, tJP[j]) && (Ci >= 1)) dm0 = 1.0 - dl_any(kg,j);
          else                                dm0 = 0.0 - dl_any(kg,j);

          if (feq(Ti, tJP[j]) && (Ci == 1)) dm1 = 1.0 - dl_c1(kg,j);
          else                                dm1 = 0.0 - dl_c1(kg,j);
        }
        if (ny(kg,j) > 0.0){
          xs[j] = xs[j-1] + dm0 / ny(kg,j);
          xf[j] = xf[j-1] + surv(kg,j-1) * ((dm1/ny(kg,j)) - (xs[j-1]*dl_c1(kg,j)));
        } else {
          xs[j] = xs[j-1];
          xf[j] = xf[j-1];
        }
        IF_tjp[kg](row, j) = xf[j];
        var_tjp(kg, j)    += xf[j] * xf[j];
      }
    }
  }

  // ---- 6) Map tJP -> TJP (output grid)
  // Build index map: for each TJP[j], find tJP index
  std::vector<int> IDX(NJP+1, 0); // 0 maps to 0
  for (int j=1;j<=NJP;j++){
    int jj = find_time_index(tJP, TJP[j]);
    if (jj < 0) stop("Internal mapping error: TJP not found in tJP.");
    IDX[j] = jj;
  }

  // Allocate outputs on TJP
  NumericVector time(NJP+1);
  for (int j=0;j<=NJP;j++) time[j] = TJP[j];

  NumericMatrix n_risk(K, NJP+1), aj1(K, NJP+1), var_aj1(K, NJP+1);
  for (int kg=0; kg<K; ++kg){
    n_risk(kg,0)  = 0.0;
    aj1(kg,0)     = 0.0;
    var_aj1(kg,0) = 0.0;
    for (int j=1;j<=NJP;j++){
      int jj = IDX[j];
      n_risk(kg,j)  = ny(kg, jj);
      aj1(kg,j)     = aj1_tjp(kg, jj);
      var_aj1(kg,j) = var_tjp(kg, jj);
    }
  }

  // Map IF matrices to TJP: for each strata, build (n_g x (NJP+1))
  List if_aj1(K);
  for (int kg=0; kg<K; ++kg){
    NumericMatrix M(n_g[kg], NJP+1);
    for (int r=0; r<n_g[kg]; ++r){
      M(r,0) = 0.0;
      for (int j=1;j<=NJP;j++){
        int jj = IDX[j];
        M(r,j) = IF_tjp[kg](r, jj);
      }
    }
    if_aj1[kg] = M;
  }

  // Optional: return strata sizes as atribute
  IntegerVector strata_sizes(K);
  for (int kg=0; kg<K; ++kg) strata_sizes[kg] = n_g[kg];

  return List::create(
    _["time"]     = time,
    _["n.risk"]   = n_risk,
    _["aj1"]      = aj1,
    _["var_aj1"]  = var_aj1,
    _["if_aj1"]   = if_aj1,
    _["strata_sizes"] = strata_sizes
  );
}



// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {
  inline double qnorm_from_conf(double conf) {
    return R::qnorm(0.5 + conf/2.0, 0.0, 1.0, 1, 0);
  }
  inline double clamp01(double x, double eps=1e-12){
    if (x < eps) return eps;
    if (x > 1.0 - eps) return 1.0 - eps;
    return x;
  }
  // 正規化: sum_j w[j] * dt[j] = 1 (j=1..m)
  std::vector<double> normalize_weight(const NumericVector& time, const NumericVector& w_in){
    const int m = time.size()-1;
    std::vector<double> w(m+1, 0.0);
    if (w_in.size() != (m+1)) stop("weight length mismatch.");
    double denom = 0.0;
    for (int j=1;j<=m;++j) denom += w_in[j] * (time[j]-time[j-1]);
    if (denom <= 0) stop("weight normalization failed (sum w*dt <= 0).");
    for (int j=1;j<=m;++j) w[j] = w_in[j]/denom;
    return w;
  }
  // Beta(P,Q)重み生成（基準は baseF、j=1..m で使用）
  std::vector<double> beta_weight(const NumericVector& time, const NumericVector& baseF, double P, double Q){
    const int m = time.size()-1;
    if (baseF.size() != (m+1)) stop("baseF length mismatch.");
    std::vector<double> w(m+1, 0.0);
    const double Fmax = std::max(baseF[m], 1e-12);
    for (int j=1;j<=m;++j){
      double r = baseF[j]/Fmax; if (r<0) r=0; if (r>1) r=1;
      w[j] = std::pow(1.0-r, P) * std::pow(r, Q);
    }
    // normalize
    double denom = 0.0;
    for (int j=1;j<=m;++j) denom += w[j]*(time[j]-time[j-1]);
    if (denom <= 0) stop("beta weight normalization failed.");
    for (int j=1;j<=m;++j) w[j] /= denom;
    return w;
  }
}

//' Weighted averages (K groups) with IF-based (JK) variance
//'
//' @param time   NumericVector length m+1 (0 included). From calculateIFofAJ()$time
//' @param aj1    K x (m+1) matrix of CIF (cause 1). From calculateIFofAJ()$aj1
//' @param if_aj1 List length K; each is matrix n_g x (m+1) of individual IF. From calculateIFofAJ()$if_aj1
//' @param contrasts NumericMatrix L x K. Each row sums to 0. Defines L contrasts for CIF差（diff_aj1）
//' @param weights List of weight vectors (each length m+1). If empty and P,Q given, Beta weights are created.
//' @param P      NumericVector of length R (optional). Beta-weight P for each weight grid.
//' @param Q      NumericVector of length R (optional). Beta-weight Q for each weight grid.
//' @param weight_base integer (0..K). 0=overall mean CIF as base; 1..K=use that group's CIF as base for Beta weights.
//' @param conf_int double (0,1), e.g., 0.95
//' @param ref    integer in 1..K. Reference group for rr_aj1 / or_aj1
//'
//' @return List with
//'   - weight: list of normalized weights (length R)
//'   - diff_aj1: data.frame (contrast_id, weight_id, est, se, lwr, upr, z, p)
//'   - rr_aj1:   data.frame (group,      weight_id, est, se, lwr, upr, z, p)  # group vs ref
//'   - or_aj1:   data.frame (group,      weight_id, est, se, lwr, upr, z, p)  # group vs ref
//'   - cov_diff: covariance matrix of IF projections for diff (size (L*R) x (L*R))
//'   - cov_rr:   covariance matrix for log(RR) projections ((K-1)*R x (K-1)*R)
//'   - cov_or:   covariance matrix for log(OR) projections ((K-1)*R x (K-1)*R)
//'
// [[Rcpp::export]]
List calculateWeightedAverage(
  NumericVector time,
  NumericMatrix aj1,
  List          if_aj1,
  NumericMatrix contrasts,
  List          weights = R_NilValue,
  Nullable<NumericVector> P = R_NilValue,
  Nullable<NumericVector> Q = R_NilValue,
  int           weight_base = 0,
  double        conf_int = 0.95,
  int           ref = 1
){
  const int m = time.size()-1;
  if (m <= 0) stop("time must have length >= 2.");
  const int K = aj1.nrow();
  if (aj1.ncol() != (m+1)) stop("aj1 must be K x (m+1).");
  if (if_aj1.size() != K)  stop("if_aj1 must be list of length K.");
  // quick check IF dims
  for (int g=0; g<K; ++g){
    NumericMatrix M = if_aj1[g];
    if (M.ncol() != (m+1)) stop("if_aj1[[g]] must have (m+1) columns.");
  }
  if (ref < 1 || ref > K) stop("ref out of range.");
  const int refg = ref-1;

  // contrasts L x K
  const int L = contrasts.nrow();
  if (contrasts.ncol() != K) stop("contrasts must be L x K.");
  // optional: check each row sums to ~0
  for (int l=0;l<L;++l){
    double s = 0.0; for (int g=0; g<K; ++g) s += contrasts(l,g);
    if (std::fabs(s) > 1e-8) stop("each contrast row must sum to 0.");
  }

  // ---- Build weight grids ----
  std::vector< std::vector<double> > W; // R x (m+1)
  int Rw = weights.size();
  NumericVector Pv, Qv;
  if (Rw == 0){
    // generate from P,Q
    if (P.isNull() || Q.isNull()) stop("weights empty: provide P and Q to generate Beta weights.");
    Pv = as<NumericVector>(P); Qv = as<NumericVector>(Q);
    if (Pv.size() != Qv.size() || Pv.size()==0) stop("P and Q must have same positive length.");
    Rw = Pv.size();
    // baseF
    NumericVector baseF(m+1);
    if (weight_base == 0){
      // overall mean CIF at each time
      for (int j=0;j<=m;++j){
        double s=0.0; for (int g=0; g<K; ++g) s += aj1(g,j);
        baseF[j] = s / K;
      }
    } else {
      int bg = weight_base-1;
      if (bg<0 || bg>=K) stop("weight_base out of range.");
      for (int j=0;j<=m;++j) baseF[j] = aj1(bg, j);
    }
    for (int r=0; r<Rw; ++r){
      W.push_back( beta_weight(time, baseF, Pv[r], Qv[r]) );
    }
  } else {
    // normalize given weights
    for (int r=0; r<Rw; ++r){
      NumericVector wr = weights[r];
      W.push_back( normalize_weight(time, wr) );
    }
  }

  const double zq = qnorm_from_conf(conf_int);
  const double EPS = 1e-12;

  // 便利なΔt
  std::vector<double> dt(m+1, 0.0);
  for (int j=1;j<=m;++j) dt[j] = time[j]-time[j-1];

  // =======================
  // 1) diff_aj1（任意コントラスト）: est, IF 射影, 分散
  // =======================
  const int NR_diff = L * Rw;
  NumericVector diff_est(NR_diff), diff_se(NR_diff), diff_lwr(NR_diff), diff_upr(NR_diff), diff_z(NR_diff), diff_p(NR_diff);
  IntegerVector diff_contrast_id(NR_diff), diff_weight_id(NR_diff);

  // 各コンボの IF 射影（個体別）を保持して共分散を作る
  // 先に全個体総数を把握
  int Ntot = 0;
  std::vector<int> offset_g(K,0);
  for (int g=0; g<K; ++g){ NumericMatrix M = if_aj1[g]; offset_g[g] = Ntot; Ntot += M.nrow(); }
  NumericMatrix IFproj_diff(Ntot, NR_diff); // (個体) x (L*R)

  for (int l=0; l<L; ++l){
    // 事前に contrast l の時点ごとの係数 c_g をコピー
    std::vector<double> c(K);
    for (int g=0; g<K; ++g) c[g] = contrasts(l,g);

    for (int r=0; r<Rw; ++r){
      const int idx = l*Rw + r;
      diff_contrast_id[idx] = l+1;
      diff_weight_id[idx]   = r+1;

      // 推定量（時間平均）
      double est = 0.0;
      for (int j=1;j<=m;++j){
        double f = 0.0;
        for (int g=0; g<K; ++g) f += c[g] * aj1(g,j);
        est += W[r][j] * dt[j] * f;
      }
      diff_est[idx] = est;

      // IF 射影と分散
      double var = 0.0;
      int rowOfs = 0;
      for (int g=0; g<K; ++g){
        NumericMatrix IFg = if_aj1[g];
        const int ng = IFg.nrow();
        for (int i=0; i<ng; ++i){
          double prj = 0.0;
          for (int j=1;j<=m;++j){
            prj += W[r][j] * dt[j] * (c[g] * IFg(i,j));
          }
          IFproj_diff(rowOfs+i, idx) = prj;
          var += prj*prj;
        }
        rowOfs += ng;
      }
      double se = std::sqrt(var);
      diff_se[idx]  = se;
      diff_lwr[idx] = est - zq*se;
      diff_upr[idx] = est + zq*se;
      if (se>0){
        double z = est/se;
        diff_z[idx] = z;
        diff_p[idx] = R::pchisq(z*z, 1.0, false, false);
      } else {
        diff_z[idx] = NA_REAL; diff_p[idx] = NA_REAL;
      }
    }
  }
  // 共分散（(L*R)x(L*R)）
  NumericMatrix cov_diff(NR_diff, NR_diff);
  for (int a=0; a<NR_diff; ++a){
    for (int b=a; b<NR_diff; ++b){
      double s=0.0; for (int i=0;i<Ntot;++i) s += IFproj_diff(i,a)*IFproj_diff(i,b);
      cov_diff(a,b)=cov_diff(b,a)=s;
    }
  }

  // =======================
  // 2) rr_aj1（参照群 vs 他群）：自然スケール平均 + logスケールIF
  // =======================
  const int Ncmp = K-1;
  const int NR_rr = Ncmp * Rw;
  NumericVector rr_est(NR_rr), rr_se(NR_rr), rr_lwr(NR_rr), rr_upr(NR_rr), rr_z(NR_rr), rr_p(NR_rr);
  IntegerVector rr_group(NR_rr), rr_weight_id(NR_rr);

  NumericMatrix IFproj_rr(Ntot, NR_rr); // log(RR)上の IF 射影

  // 参照群の IF 行列
  NumericMatrix IFref = if_aj1[refg];

  for (int g=0, pos=0; g<K; ++g){
    if (g==refg) continue;
    for (int r=0; r<Rw; ++r){
      const int idx = pos*Rw + r;
      rr_group[idx] = g+1;
      rr_weight_id[idx] = r+1;

      // 推定（自然スケールで平均）
      double est = 0.0;
      for (int j=1;j<=m;++j){
        double u = std::max(aj1(g,j),   EPS);
        double v = std::max(aj1(refg,j),EPS);
        est += W[r][j] * dt[j] * (u/v);
      }
      rr_est[idx] = est;

      // 分散：logスケール IF の時間積分
      double var = 0.0;
      int rowOfs = 0;
      for (int gg=0; gg<K; ++gg){
        NumericMatrix IFg = if_aj1[gg];
        const int ng = IFg.nrow();
        for (int i=0;i<ng;++i){
          double prj = 0.0;
          for (int j=1;j<=m;++j){
            double u = clamp01(aj1(g,j), EPS);
            double v = clamp01(aj1(refg,j), EPS);
            if (gg == g){
              prj += W[r][j] * dt[j] * (IFg(i,j)/u);
            } else if (gg == refg){
              prj += W[r][j] * dt[j] * (- IFg(i,j)/v);
            } // 他群は寄与0
          }
          // 個体のグローバル row index
          int gi_ofs = offset_g[gg];
          IFproj_rr(gi_ofs+i, idx) = prj;
          var += prj*prj;
        }
      }
      double se = std::sqrt(var);
      rr_se[idx] = se;
      if (se>0 && est>EPS){
        double var_log = var; // logスケール
        double se_log  = std::sqrt(var_log);
        double log_est = std::log(est);
        rr_lwr[idx] = std::exp(log_est - zq*se_log);
        rr_upr[idx] = std::exp(log_est + zq*se_log);
        double z    = log_est / se_log;
        rr_z[idx]   = z;
        rr_p[idx]   = R::pchisq(z*z, 1.0, false, false);
      } else {
        rr_lwr[idx]=rr_upr[idx]=rr_z[idx]=rr_p[idx]=NA_REAL;
      }
    }
    ++pos;
  }
  // 共分散（(K-1)*R x (K-1)*R）: log(RR) IF 射影
  NumericMatrix cov_rr(NR_rr, NR_rr);
  for (int a=0; a<NR_rr; ++a){
    for (int b=a; b<NR_rr; ++b){
      double s=0.0; for (int i=0;i<Ntot;++i) s += IFproj_rr(i,a)*IFproj_rr(i,b);
      cov_rr(a,b)=cov_rr(b,a)=s;
    }
  }

  // =======================
  // 3) or_aj1（参照群 vs 他群）：自然スケール平均 + logスケールIF
  // =======================
  NumericVector or_est(NR_rr), or_se(NR_rr), or_lwr(NR_rr), or_upr(NR_rr), or_z(NR_rr), or_p(NR_rr);
  IntegerVector or_group(NR_rr), or_weight_id(NR_rr);

  NumericMatrix IFproj_or(Ntot, NR_rr);

  for (int g=0, pos=0; g<K; ++g){
    if (g==refg) continue;
    for (int r=0; r<Rw; ++r){
      const int idx = pos*Rw + r;
      or_group[idx] = g+1;
      or_weight_id[idx] = r+1;

      // 推定（自然スケールで平均）
      double est = 0.0;
      for (int j=1;j<=m;++j){
        double u = clamp01(aj1(g,j),   EPS);
        double v = clamp01(aj1(refg,j),EPS);
        double ou = u/(1.0-u), ov=v/(1.0-v);
        est += W[r][j] * dt[j] * (ou/ov);
      }
      or_est[idx] = est;

      // 分散：log(OR) = logit(u) - logit(v)
      double var = 0.0;
      for (int gg=0; gg<K; ++gg){
        NumericMatrix IFg = if_aj1[gg];
        const int ng = IFg.nrow();
        for (int i=0;i<ng;++i){
          double prj = 0.0;
          for (int j=1;j<=m;++j){
            double u = clamp01(aj1(g,j),   EPS);
            double v = clamp01(aj1(refg,j),EPS);
            if (gg == g){
              prj += W[r][j] * dt[j] * (IFg(i,j) / (u*(1.0-u)));
            } else if (gg == refg){
              prj += W[r][j] * dt[j] * (- IFg(i,j) / (v*(1.0-v)));
            }
          }
          int gi_ofs = offset_g[gg];
          IFproj_or(gi_ofs+i, idx) = prj;
          var += prj*prj;
        }
      }
      double se = std::sqrt(var);
      or_se[idx] = se;
      if (se>0 && est>EPS){
        double var_log = var; // すでにlogスケール
        double se_log  = std::sqrt(var_log);
        double log_est = std::log(est);
        or_lwr[idx] = std::exp(log_est - zq*se_log);
        or_upr[idx] = std::exp(log_est + zq*se_log);
        double z    = log_est / se_log;
        or_z[idx]   = z;
        or_p[idx]   = R::pchisq(z*z, 1.0, false, false);
      } else {
        or_lwr[idx]=or_upr[idx]=or_z[idx]=or_p[idx]=NA_REAL;
      }
    }
    ++pos;
  }
  // 共分散（(K-1)*R x (K-1)*R）: log(OR) IF 射影
  NumericMatrix cov_or(NR_rr, NR_rr);
  for (int a=0; a<NR_rr; ++a){
    for (int b=a; b<NR_rr; ++b){
      double s=0.0; for (int i=0;i<Ntot;++i) s += IFproj_or(i,a)*IFproj_or(i,b);
      cov_or(a,b)=cov_or(b,a)=s;
    }
  }

  // 出力 DataFrame
  DataFrame diff_df = DataFrame::create(
    _["contrast_id"] = diff_contrast_id,
    _["weight_id"]   = diff_weight_id,
    _["est"] = diff_est,
    _["se"]  = diff_se,
    _["lwr"] = diff_lwr,
    _["upr"] = diff_upr,
    _["z"]   = diff_z,
    _["p"]   = diff_p
  );
  DataFrame rr_df = DataFrame::create(
    _["group"]     = rr_group,      // vs ref
    _["weight_id"] = rr_weight_id,
    _["est"] = rr_est,
    _["se"]  = rr_se,
    _["lwr"] = rr_lwr,
    _["upr"] = rr_upr,
    _["z"]   = rr_z,
    _["p"]   = rr_p
  );
  DataFrame or_df = DataFrame::create(
    _["group"]     = or_group,      // vs ref
    _["weight_id"] = or_weight_id,
    _["est"] = or_est,
    _["se"]  = or_se,
    _["lwr"] = or_lwr,
    _["upr"] = or_upr,
    _["z"]   = or_z,
    _["p"]   = or_p
  );

  // 重みリスト（R側で再利用しやすいよう返却）
  List Wout(Rw);
  for (int r=0;r<Rw;++r){
    NumericVector w(m+1); w[0]=0.0;
    for (int j=1;j<=m;++j) w[j]=W[r][j];
    Wout[r]=w;
  }

  return List::create(
    _["weight"]   = Wout,
    _["diff_aj1"] = diff_df,
    _["rr_aj1"]   = rr_df,
    _["or_aj1"]   = or_df,
    _["cov_diff"] = cov_diff,
    _["cov_rr"]   = cov_rr,
    _["cov_or"]   = cov_or
  );
}


// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

namespace {
  // 単純Cholesky（非特異を想定）。失敗したら false。
  bool chol_lower(const NumericMatrix& A, NumericMatrix& L) {
    const int p = A.nrow();
    if (A.ncol()!=p) return false;
    for (int i=0;i<p;++i){
      for (int j=0;j<=i;++j){
        double s = A(i,j);
        for (int k=0;k<j;++k) s -= L(i,k)*L(j,k);
        if (i==j){
          if (s <= 0.0) return false;
          L(i,j) = std::sqrt(s);
        } else {
          L(i,j) = s / L(j,j);
        }
      }
      for (int j=i+1;j<p;++j) L(i,j)=0.0;
    }
    return true;
  }

  // 共分散 -> 相関
  NumericMatrix cov_to_corr(const NumericMatrix& S) {
    const int p = S.nrow();
    NumericMatrix R(p,p);
    for (int i=0;i<p;++i){
      double sii = S(i,i);
      double di  = (sii>0.0) ? std::sqrt(sii) : NA_REAL;
      for (int j=0;j<p;++j){
        double sjj = S(j,j);
        double dj  = (sjj>0.0) ? std::sqrt(sjj) : NA_REAL;
        double denom = (R_finite(di) && R_finite(dj) && di>0.0 && dj>0.0) ? (di*dj) : NA_REAL;
        double v = S(i,j);
        R(i,j) = (R_finite(denom) ? (v/denom) : (i==j?1.0:NA_REAL));
      }
    }
    return R;
  }
}

//' Core maxCombo calculator from z & covariance
//'
//' @param z   NumericVector of standardized statistics (per weight × per contrast/comparison).
//'            Use NA to drop entries.
//' @param cov NumericMatrix covariance matrix aligned with z (same length).
//' @param n_simuation integer, Monte Carlo draws for MVN-based p-value (default 20000).
//'
//' @return list(stat, p, R, corr)
//'   - stat : observed max |z|
//'   - p    : Monte Carlo p-value under MVN(0, corr)
//'   - R    : number of z used (after NA removal)
//'   - corr : correlation matrix used
//'
//' @details
//' This implements the maxCombo test (max |Z|) using IF-based covariance (via input cov).
//' It draws X ~ MVN(0, corr) with Cholesky and computes Pr(max |X_r| >= stat).
// [[Rcpp::export]]
List calculateMaxCombo(NumericVector z, NumericMatrix cov, int n_simuation = 20000) {
  const int P = z.size();
  if (cov.nrow()!=P || cov.ncol()!=P) stop("cov shape must match length(z).");

  // 有効な（finite）要素だけ抽出
  std::vector<int> keep;
  keep.reserve(P);
  for (int i=0;i<P;++i) {
    if (R_finite(z[i])) keep.push_back(i);
  }
  const int R = (int)keep.size();
  if (R==0) return List::create(
    _["stat"] = NA_REAL, _["p"] = NA_REAL, _["R"] = 0, _["corr"] = NumericMatrix(0,0)
  );

  NumericVector z_use(R);
  NumericMatrix cov_use(R,R);
  for (int a=0;a<R;++a){
    z_use[a] = z[ keep[a] ];
    for (int b=0;b<R;++b) cov_use(a,b) = cov( keep[a], keep[b] );
  }

  // 観測統計量
  double zmax_obs = 0.0;
  for (int r=0;r<R;++r){
    double a = std::fabs(z_use[r]);
    if (a > zmax_obs) zmax_obs = a;
  }

  // 相関行列
  NumericMatrix corr = cov_to_corr(cov_use);

  // コレスキー分解（失敗ならBonferroni）
  NumericMatrix L(R,R);
  bool ok = chol_lower(corr, L);

  double p_mc = NA_REAL;
  if (ok) {
    int exceed = 0;
    for (int b=0;b<n_simuation;++b){
      // y ~ N(0,I)
      NumericVector y(R);
      for (int r=0;r<R;++r) y[r] = R::rnorm(0.0, 1.0);
      // x = L y
      NumericVector x(R);
      for (int i=0;i<R;++i){
        double s=0.0;
        for (int k=0;k<=i;++k) s += L(i,k)*y[k];
        x[i] = s;
      }
      double zmax_sim = 0.0;
      for (int r=0;r<R;++r){
        double a = std::fabs(x[r]);
        if (a > zmax_sim) zmax_sim = a;
      }
      if (zmax_sim >= zmax_obs - 1e-12) exceed++;
    }
    p_mc = (exceed + 1.0) / (n_simuation + 1.0);
  } else {
    // Choleskyが失敗したら Bonferroni で保守的に
    double pmin = 1.0;
    for (int r=0;r<R;++r){
      // 単独p（両側）: 2*(1-Phi(|z|))
      double pr = 2.0 * R::pnorm(-std::fabs(z_use[r]), 0.0, 1.0, /*lower_tail*/true, /*log_p*/false);
      if (R_finite(pr)) pmin = std::min(pmin, pr);
    }
    p_mc = std::min(1.0, R * pmin);
  }

  return List::create(
    _["stat"] = zmax_obs,
    _["p"]    = p_mc,
    _["R"]    = R,
    _["corr"] = corr
  );
}


// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Convenience wrapper: build maxCombo directly from calculateWeightedAverage() result
//'
//' @param out_cwa      List output from calculateWeightedAverage()
//' @param measure "diff", "rr", or "or"  （diff_aj1 / rr_aj1 / or_aj1 に対応）
//' @param n_simuation   Monte Carlo draws (default 20000)
//'
//' @return list(per_weight = data.frame, maxcombo = list(stat,p,R,corr))
//'
//' @details
//' - それぞれの measure で、per-weight（× per-contrast or per-group）の Z を連結して maxCombo を計算します。
//' - 共分散行列は out_cwa$cov_diff / cov_rr / cov_or を使用します。
//' - NAのZは自動で落とします（共分散も対応成分で抜粋）。
// [[Rcpp::export]]
List callMaxCombo(List out_cwa, std::string measure = "diff", int n_simuation = 20000) {
  DataFrame df;
  NumericMatrix cov;
  if (measure == "diff") {
    df  = as<DataFrame>(out_cwa["diff_aj1"]);
    cov = as<NumericMatrix>(out_cwa["cov_diff"]);
  } else if (measure == "rr") {
    df  = as<DataFrame>(out_cwa["rr_aj1"]);
    cov = as<NumericMatrix>(out_cwa["cov_rr"]);
  } else if (measure == "or") {
    df  = as<DataFrame>(out_cwa["or_aj1"]);
    cov = as<NumericMatrix>(out_cwa["cov_or"]);
  } else {
    stop("measure must be one of 'diff', 'rr', or 'or'.");
  }

  if (!df.containsElementNamed("z"))
    stop("The input list must contain a data.frame with a 'z' column for the chosen measure.");

  NumericVector z = df["z"];
  // コア関数へ
  Function core("calculateMaxCombo");
  List mc = core(z, cov, n_simuation);

  return List::create(
    _["per_weight"] = df,
    _["maxcombo"]   = mc
  );
}

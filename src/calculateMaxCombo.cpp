// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {
constexpr double EPS = 1e-12;

inline bool feq(double a, double b, double eps=EPS){
  double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
  return std::fabs(a-b) <= eps * scale;
}

struct Rec {
  double t;
  int    epsilon;
  int    g;
  int    id;
};

inline double clamp01(double x, double eps=1e-12){
  if (!R_finite(x)) return eps;
  if (x < eps) return eps;
  if (x > 1.0 - eps) return 1.0 - eps;
  return x;
}

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
 //' @param t numeric vector of event/censoring times (>=0)
 //' @param epsilon integer vector: 0=censor, 1=event1, 2=event2 (>=2 treated as "any event" but CIF is for event1)
 //' @param strata integer vector of strata indices in {1,2,...,K}
 //'
 //' @return list with
 //'   - time:     numeric vector of cause-1 jump times (union across strata), length m+1 with leading 0
 //'   - n.risk:   K x (m+1) matrix of at-risk counts on the time grid
 //'   - aj1:      K x (m+1) matrix of Aalen–Johansen CIF for event1
 //'   - var_aj1:  K x (m+1) matrix of variances (sum of IF^2)
 //'   - if_aj1:   list of length K; each is a matrix (n_g x (m+1)) of individual IF trajectories
 //'
 //' @details
 //' Internally builds two time grids (unions across strata):
 //'   TTJP: all event times (epsilon>=1) with leading 0
 //'   TJP:  cause-1 times (epsilon==1) with leading 0  <-- output grid
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

   IntegerVector guniq = sort_unique(strata);
   const int K = guniq.size();
   std::vector<int> lab2k( *std::max_element(guniq.begin(), guniq.end()) + 1, -1 );
   for (int k=0;k<K;k++) lab2k[guniq[k]] = k;

   std::vector<Rec> x; x.reserve(N);
   for (int i=0;i<N;i++){
     int gl = strata[i];
     int kg = lab2k[gl];
     if (kg < 0) stop("Internal: strata label mapping failed.");
     x.push_back({(double)t[i], (int)epsilon[i], gl, i});
   }
   std::sort(x.begin(), x.end(), [](const Rec& a, const Rec& b){
     if (!feq(a.t,b.t)) return a.t < b.t;
     return a.epsilon > b.epsilon;
   });

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
   const int TNJP = (int)tJP.size() - 1;
   const int NJP  = (int)TJP.size() - 1;

   NumericMatrix ny(K, TNJP+1);
   NumericMatrix nd_any(K, TNJP+1);
   NumericMatrix nd_c1(K, TNJP+1);

   std::vector< std::vector<const Rec*> > perG(K);
   for (const auto& r : x) perG[ lab2k[r.g] ].push_back(&r);
   for (int kg=0; kg<K; ++kg){
     auto& v = perG[kg];
     size_t pointer_lt = 0;
     size_t pointer_eq = 0;
     const size_t n_g = v.size();
     for (int j=1;j<=TNJP;++j){
       const double tj = tJP[j];
       while (pointer_lt < n_g && v[pointer_lt]->t + EPS < tj) ++pointer_lt;
       ny(kg, j) = (double)(n_g - pointer_lt);

       size_t tmp = pointer_eq;
       while (tmp < n_g && feq(v[tmp]->t, tj)){
         if (v[tmp]->epsilon >= 1) nd_any(kg, j) += 1.0;
         if (v[tmp]->epsilon == 1) nd_c1(kg, j)  += 1.0;
         ++tmp;
       }
       pointer_eq = tmp;
     }
   }

   NumericMatrix dl_any(K, TNJP+1), dl_c1(K, TNJP+1);
   NumericMatrix surv(K, TNJP+1), aj1_tjp(K, TNJP+1), var_tjp(K, TNJP+1);

   std::vector< std::vector<int> > idx_of_strata(K);
   for (int i=0;i<N;i++){
     int kg = lab2k[strata[i]];
     idx_of_strata[kg].push_back(i);
   }
   std::vector<int> n_g(K,0);
   for (int kg=0; kg<K; ++kg) n_g[kg] = (int)idx_of_strata[kg].size();

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

   std::vector< NumericMatrix > IF_tjp(K, NumericMatrix(0,0));
   for (int kg=0; kg<K; ++kg){
     IF_tjp[kg] = NumericMatrix(n_g[kg], TNJP+1);
   }

   for (int kg=0; kg<K; ++kg){
     for (int row=0; row<n_g[kg]; ++row){
       int i = idx_of_strata[kg][row];              // original individual index in data order
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

   std::vector<int> IDX(NJP+1, 0);
   for (int j=1;j<=NJP;j++){
     int jj = find_time_index(tJP, TJP[j]);
     if (jj < 0) stop("Internal mapping error: TJP not found in tJP.");
     IDX[j] = jj;
   }

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


//' Cause-specific Nelson–Aalen for cause 1 with individual influence functions (K groups)
 //'
 //' @param t numeric vector of observed times (>=0)
 //' @param epsilon integer vector: 0=censoring, 1=event1 (event of interest), 2+=event2 (competing risks, treated as censoring)
 //' @param strata integer vector of strata labels in {1,2,...,K}
 //'
 //' @return list with
 //'   - time:         numeric vector of cause-1 jump times (union across groups), length m+1 with leading 0
 //'   - n.risk:       K x (m+1) matrix of at-risk counts Y_g(t_j) on this grid
 //'   - cumhaz1:      K x (m+1) matrix of Nelson–Aalen cumulative hazard for cause 1
 //'   - var_cumhaz1:  K x (m+1) matrix of variances (sum of IF^2)
 //'   - if_cumhaz1:   list length K; each is matrix (n_g x (m+1)) of individual IF trajectories
 //'
 //' @details
 //' Competing events (cause >= 2) are treated as censoring for the cause-specific hazard of event1,
 //' i.e., they exit the risk set at their time. The time grid uses only cause-1 event times (with leading 0).
 //'
 //' This function is intended for log-rank–type testing pipelines (Fleming–Harrington weights, maxCombo(log-rank), etc.).
 //' @examples
 //' # res <- calculateIFofNA(t, epsilon, strata)
 // [[Rcpp::export]]

 List calculateIFofNA(NumericVector t, IntegerVector epsilon, IntegerVector strata){
   const int N = t.size();
   if (epsilon.size()!=N || strata.size()!=N) stop("Length mismatch.");
   for (int i=0;i<N;i++){
     if (!R_finite(t[i]) || t[i] < 0.0) stop("t must be finite and nonnegative.");
     if (strata[i] < 1) stop("strata must be >=1.");
   }

   // --- unique stratas & map to 0..K-1
   IntegerVector guniq = sort_unique(strata);
   const int K = guniq.size();
   std::vector<int> lab2k( *std::max_element(guniq.begin(), guniq.end()) + 1, -1 );
   for (int k=0;k<K;k++) lab2k[guniq[k]] = k;

   // --- build records & sort by time (events before censors at ties; for risk set counting, >= is used)
   std::vector<Rec> x; x.reserve(N);
   for (int i=0;i<N;i++){
     x.push_back({(double)t[i], (int)epsilon[i], (int)strata[i], i});
   }
   std::sort(x.begin(), x.end(), [](const Rec& a, const Rec& b){
     if (!feq(a.t,b.t)) return a.t < b.t;
     // event (epsilon==1) first, then others (including censor/competing)
     // (tie rule does not affect Y(t_j) since we use >= for risk counting)
     return (a.epsilon==1) > (b.epsilon==1);
   });

   // --- time grid: union of epsilon-1 event times (leading 0)
   std::vector<double> TJP; TJP.reserve(N+1); TJP.push_back(0.0);
   for (int i=0;i<N;i++){
     if (x[i].epsilon == 1){
       if (!feq(TJP.back(), x[i].t)) TJP.push_back(x[i].t);
     }
   }
   const int NJP = (int)TJP.size() - 1; // indices 1..NJP
   NumericVector time(NJP+1);
   for (int j=0;j<=NJP;j++) time[j] = TJP[j];

   // --- per-strata counts on TJP: Y_g(j), d1_g(j)
   NumericMatrix Y(K, NJP+1);      // risk set
   NumericMatrix d1(K, NJP+1);     // epsilon-1 events at the jump time

   // risk set: Y_g(t_j) = #{ i in g : T_i >= t_j } (competing acts as censoring: they still satisfy >= before jump)
   for (int j=1;j<=NJP;j++){
     double tj = TJP[j];
     for (int i=0;i<N;i++){
       if (x[i].t + EPS >= tj){
         int kg = lab2k[x[i].g];
         Y(kg, j) += 1.0;
       }
     }
   }
   // epsilon-1 event counts at exact jump
   for (int j=1;j<=NJP;j++){
     double tj = TJP[j];
     for (int i=0;i<N;i++){
       if (feq(x[i].t, tj) && x[i].epsilon==1){
         int kg = lab2k[x[i].g];
         d1(kg, j) += 1.0;
       }
     }
   }

   NumericMatrix cumhaz1(K, NJP+1);
   NumericMatrix var_cumhaz1(K, NJP+1);
   std::vector< std::vector<int> > idx_of_strata(K);
   for (int i=0;i<N;i++) idx_of_strata[ lab2k[strata[i]] ].push_back(i);
   std::vector<int> n_g(K);
   for (int g=0; g<K; ++g) n_g[g] = (int)idx_of_strata[g].size();

   std::vector< NumericMatrix > IF_store(K);
   for (int gk=0; gk<K; ++gk){
     cumhaz1(gk,0) = 0.0; var_cumhaz1(gk,0)=0.0;
     NumericMatrix IFg(n_g[gk], NJP+1);
     for (int r=0;r<n_g[gk];++r) IFg(r,0)=0.0;
     for (int j=1;j<=NJP;j++){
       if (Y(gk,j) > 0.0){
         cumhaz1(gk,j) = cumhaz1(gk,j-1) + d1(gk,j) / Y(gk,j);
       } else {
         cumhaz1(gk,j) = cumhaz1(gk,j-1);
       }
       const double Yj = Y(gk,j);
       const double d1j = d1(gk,j);
       for (int r=0; r<n_g[gk]; ++r){
         int i = idx_of_strata[gk][r];
         double Ti = t[i];
         int Ci = epsilon[i];
         double incr = 0.0;
         if (Yj > 0.0){
           double a = (feq(Ti, TJP[j]) && (Ci==1)) ? 1.0 : 0.0;
           double y = (Ti + EPS >= TJP[j]) ? 1.0 : 0.0;
           incr = (a - y * (d1j / Yj)) / Yj;
         }
         IFg(r, j) = IFg(r, j-1) + incr;
         var_cumhaz1(gk, j) += IFg(r, j) * IFg(r, j);
       }
     }
     IF_store[gk] = IFg;
   }

   List if_list(K);
   for (int gk=0; gk<K; ++gk) if_list[gk] = IF_store[gk];

   NumericMatrix n_risk = Y;
   return List::create(
     _["time"]        = time,
     _["n.risk"]      = n_risk,
     _["cumhaz1"]     = cumhaz1,
     _["var_cumhaz1"] = var_cumhaz1,
     _["if_cumhaz1"]  = if_list
   );
 }

namespace {
inline double qnorm_from_conf(double conf) {
  return R::qnorm(0.5 + conf/2.0, 0.0, 1.0, 1, 0);
}
std::vector<double> normalize_weight(const NumericVector& time, const NumericVector& w_in){
  const int m = time.size()-1;
  if (w_in.size() != (m+1)) stop("weight length mismatch.");
  std::vector<double> w(m+1, 0.0);
  double denom = 0.0;
  for (int j=1;j<=m;++j) denom += w_in[j] * (time[j]-time[j-1]);
  if (denom <= 0) stop("weight normalization failed (sum w*dt <= 0).");
  for (int j=1;j<=m;++j) w[j] = w_in[j]/denom;
  return w;
}

std::vector<double> beta_weight_from_base(const NumericVector& time, const NumericVector& base, double rho, double gamma){
  const int m = time.size()-1;
  if (base.size() != (m+1)) stop("base length mismatch.");
  std::vector<double> w(m+1, 0.0);

  double xmax = 1e-12;
  for (int j=0;j<=m;++j) if (R_finite(base[j])) xmax = std::max(xmax, (double)base[j]);
  for (int j=1;j<=m;++j){
    double r = base[j]/xmax; if (r<0) r=0; if (r>1) r=1;
    w[j] = std::pow(1.0-r, rho) * std::pow(r, gamma);
  }
  double denom = 0.0;
  for (int j=1;j<=m;++j) denom += w[j]*(time[j]-time[j-1]);
  if (denom <= 0) stop("beta weight normalization failed.");
  for (int j=1;j<=m;++j) w[j] /= denom;
  return w;
}
}

//' Linear weighted average with IF-based (JK) variance (generic: difference only)
 //'
 //' @param time   NumericVector length m+1 (0 included). analysis grid
 //' @param curve  K x (m+1) matrix of target trajectories (e.g., cumhaz1 or CIF)
 //' @param if_list List length K; each is matrix (n_g x (m+1)) of individual IF for the same target
 //' @param contrasts NumericMatrix L x K (row sums = 0)
 //' @param weights List of numeric vectors (each length m+1). If empty, rho/gamma/base supplied to build Beta-like weights.
 //' @param rho,gamma NumericVector (same length R) used when weights is empty.
 //' @param base NumericVector length m+1 used for Beta-like weights when weights is empty (e.g., overall CIF or pooled S)
 //' @param conf_int double (0,1)
 //'
 //' @return list(weight=list, table=data.frame(contrast_id, weight_id, est, se, lwr, upr, z, p),
 //'              cov = covariance matrix of IF-projections (size (L*R) x (L*R)))
 // [[Rcpp::export]]
 Rcpp::List calculateWeightedAverageLinear(
     Rcpp::NumericVector time,
     Rcpp::NumericMatrix curve,
     Rcpp::List          if_list,
     Rcpp::NumericMatrix contrasts,
     Rcpp::List          weights    = R_NilValue,
     Rcpp::Nullable<Rcpp::NumericVector> rho = R_NilValue,
     Rcpp::Nullable<Rcpp::NumericVector> gamma = R_NilValue,
     Rcpp::Nullable<Rcpp::NumericVector> base = R_NilValue,
     double              conf_int   = 0.95
 ){
   const int m = time.size()-1;
   if (m <= 0) Rcpp::stop("time must have length >= 2.");
   const int K = curve.nrow();
   if (curve.ncol() != (m+1)) Rcpp::stop("curve must be K x (m+1).");
   if ((int)if_list.size() != K) Rcpp::stop("if_list must be list of length K.");

   const int L = contrasts.nrow();
   if (contrasts.ncol()!=K) Rcpp::stop("contrasts must be L x K.");
   for (int l=0;l<L;++l){
     double s=0.0; for (int g=0; g<K; ++g) s+=contrasts(l,g);
     if (std::fabs(s) > 1e-8) Rcpp::stop("each contrast row must sum to 0.");
   }

   std::vector< std::vector<double> > W;
   int Rw = weights.size();
   if (Rw == 0){
     if (rho.isNull() || gamma.isNull() || base.isNull()) Rcpp::stop("When weights empty, supply rho, gamma, and base.");
     Rcpp::NumericVector rv = Rcpp::as<Rcpp::NumericVector>(rho);
     Rcpp::NumericVector gv = Rcpp::as<Rcpp::NumericVector>(gamma);
     if (rv.size()!=gv.size() || rv.size()==0) Rcpp::stop("rho and gamma must have same positive length.");
     Rcpp::NumericVector basev = Rcpp::as<Rcpp::NumericVector>(base);
     if (basev.size() != (m+1)) Rcpp::stop("base must be length m+1.");
     Rw = rv.size();
     for (int r=0;r<Rw;++r){
       W.push_back( beta_weight_from_base(time, basev, rv[r], gv[r]) );
     }
   } else {
     for (int r=0;r<Rw;++r){
       Rcpp::NumericVector wr = weights[r];
       W.push_back( normalize_weight(time, wr) );
     }
   }

   std::vector<double> dt(m+1, 0.0);
   for (int j=1;j<=m;++j) dt[j] = time[j]-time[j-1];

   int Ntot=0; std::vector<int> ofs(K,0);
   for (int g=0; g<K; ++g){ Rcpp::NumericMatrix M = if_list[g]; ofs[g]=Ntot; Ntot += M.nrow(); }

   const double zq = qnorm_from_conf(conf_int);
   const int NR = L * Rw;

   Rcpp::NumericVector est(NR), se(NR), lwr(NR), upr(NR), z(NR), p(NR);
   Rcpp::IntegerVector cid(NR), wid(NR);
   Rcpp::NumericMatrix cov(NR, NR);

   for (int l=0;l<L;++l){
     std::vector<double> c(K); for (int g=0; g<K; ++g) c[g]=contrasts(l,g);
     for (int r=0;r<Rw;++r){
       int idx = l*Rw + r;
       cid[idx]=l+1; wid[idx]=r+1;

       double estv=0.0;
       for (int j=1;j<=m;++j){
         double f=0.0; for (int g=0; g<K; ++g) f += c[g] * curve(g,j);
         estv += W[r][j] * dt[j] * f;
       }
       est[idx]=estv;

       double var=0.0;
       int ro=0;
       for (int g=0; g<K; ++g){
         Rcpp::NumericMatrix IFg = if_list[g];
         int ng = IFg.nrow();
         for (int i=0;i<ng;++i){
           double prj=0.0;
           for (int j=1;j<=m;++j) prj += W[r][j] * dt[j] * (c[g] * IFg(i,j));
           var += prj*prj;
         }
         ro += ng;
       }
       double sev = std::sqrt(var);
       se[idx]=sev;
       lwr[idx]=estv - zq*sev;
       upr[idx]=estv + zq*sev;
       if (sev>0){ double zv=estv/sev; z[idx]=zv; p[idx]=R::pchisq(zv*zv,1.0,false,false); }
       else { z[idx]=NA_REAL; p[idx]=NA_REAL; }
     }
   }

   std::fill(cov.begin(), cov.end(), 0.0);
   {
     int rowOfs = 0;
     for (int g=0; g<K; ++g){
       Rcpp::NumericMatrix IFg = if_list[g];
       const int ng = IFg.nrow();
       for (int i=0;i<ng;++i){
         std::vector<double> v(NR, 0.0);
         for (int l=0;l<L;++l){
           std::vector<double> c(K); for (int gg=0; gg<K; ++gg) c[gg]=contrasts(l,gg);
           for (int r=0;r<Rw;++r){
             int idx = l*Rw + r;
             double prj=0.0;
             for (int j=1;j<=m;++j){
               prj += W[r][j] * dt[j] * (c[g] * IFg(i,j));
             }
             v[idx] = prj;
           }
         }
         for (int a=0;a<NR;++a){
           if (v[a]==0.0) continue;
           for (int b=a;b<NR;++b){
             double s = v[a]*v[b];
             if (s!=0.0){
               cov(a,b) += s;
               if (b!=a) cov(b,a) += s;
             }
           }
         }
       }
       rowOfs += ng;
     }
   }

   Rcpp::List Wout(Rw);
   for (int r=0;r<Rw;++r){
     Rcpp::NumericVector w(m+1); w[0]=0.0;
     for (int j=1;j<=m;++j) w[j]=W[r][j];
     Wout[r]=w;
   }

   Rcpp::DataFrame tab = Rcpp::DataFrame::create(
     Rcpp::Named("contrast_id") = cid,
     Rcpp::Named("weight_id")   = wid,
     Rcpp::Named("est") = est,
     Rcpp::Named("se")  = se,
     Rcpp::Named("lwr") = lwr,
     Rcpp::Named("upr") = upr,
     Rcpp::Named("z")   = z,
     Rcpp::Named("p")   = p
   );

   return Rcpp::List::create(
     Rcpp::Named("weight") = Wout,
     Rcpp::Named("table")  = tab,
     Rcpp::Named("cov")    = cov
   );
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
     Nullable<NumericVector> rho = R_NilValue,
     Nullable<NumericVector> gamma = R_NilValue,
     int           weight_base = 0,
     double        conf_int = 0.95,
     int           ref = 1
 ){
   const int m = time.size()-1;
   if (m <= 0) stop("time must have length >= 2.");
   const int K = aj1.nrow();
   if (aj1.ncol() != (m+1)) stop("aj1 must be K x (m+1).");
   if (if_aj1.size() != K)  stop("if_aj1 must be list of length K.");

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

   std::vector< std::vector<double> > W;
   int Rw = weights.size();
   NumericVector Pv, Qv;
   if (Rw == 0){
     // generate from P,Q
     if (rho.isNull() || gamma.isNull()) stop("weights empty: provide rho and gamma to generate Beta weights.");
     Pv = as<NumericVector>(rho); Qv = as<NumericVector>(gamma);
     if (Pv.size() != Qv.size() || Pv.size()==0) stop("rho and gamma must have same positive length.");
     Rw = Pv.size();
     NumericVector baseF(m+1);
     if (weight_base == 0){
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
       W.push_back( beta_weight_from_base(time, baseF, Pv[r], Qv[r]) );
     }
   } else {
     for (int r=0; r<Rw; ++r){
       NumericVector wr = weights[r];
       W.push_back( normalize_weight(time, wr) );
     }
   }

   const double zq = qnorm_from_conf(conf_int);
   const double EPS = 1e-12;

   std::vector<double> dt(m+1, 0.0);
   for (int j=1;j<=m;++j) dt[j] = time[j]-time[j-1];

   // =======================
   // 1) diff_aj1（任意コントラスト）: est, IF 射影, 分散
   // =======================
   const int NR_diff = L * Rw;
   NumericVector diff_est(NR_diff), diff_se(NR_diff), diff_lwr(NR_diff), diff_upr(NR_diff), diff_z(NR_diff), diff_p(NR_diff);
   IntegerVector diff_contrast_id(NR_diff), diff_weight_id(NR_diff);

   int Ntot = 0;
   std::vector<int> offset_g(K,0);
   for (int g=0; g<K; ++g){ NumericMatrix M = if_aj1[g]; offset_g[g] = Ntot; Ntot += M.nrow(); }
   NumericMatrix IFproj_diff(Ntot, NR_diff); // (個体) x (L*R)

   for (int l=0; l<L; ++l){
     std::vector<double> c(K);
     for (int g=0; g<K; ++g) c[g] = contrasts(l,g);

     for (int r=0; r<Rw; ++r){
       const int idx = l*Rw + r;
       diff_contrast_id[idx] = l+1;
       diff_weight_id[idx]   = r+1;

       double est = 0.0;
       for (int j=1;j<=m;++j){
         double f = 0.0;
         for (int g=0; g<K; ++g) f += c[g] * aj1(g,j);
         est += W[r][j] * dt[j] * f;
       }
       diff_est[idx] = est;

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

   NumericMatrix IFproj_rr(Ntot, NR_rr);
   NumericMatrix IFref = if_aj1[refg];

   for (int g=0, pos=0; g<K; ++g){
     if (g==refg) continue;
     for (int r=0; r<Rw; ++r){
       const int idx = pos*Rw + r;
       rr_group[idx] = g+1;
       rr_weight_id[idx] = r+1;

       double est = 0.0;
       for (int j=1;j<=m;++j){
         double u = std::max(aj1(g,j),   EPS);
         double v = std::max(aj1(refg,j),EPS);
         est += W[r][j] * dt[j] * (u/v);
       }
       rr_est[idx] = est;

       double var = 0.0;
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
             }
           }
           int gi_ofs = offset_g[gg];
           IFproj_rr(gi_ofs+i, idx) = prj;
           var += prj*prj;
         }
       }
       double se = std::sqrt(var);
       rr_se[idx] = se;
       if (se>0 && est>EPS){
         double var_log = var;
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

   NumericMatrix cov_rr(NR_rr, NR_rr);
   {
     for (int gg=0; gg<K; ++gg){
       NumericMatrix IFg = if_aj1[gg];
       const int ng = IFg.nrow();
       for (int i=0;i<ng;++i){
         std::vector<double> v(NR_rr, 0.0);
         int pos = 0;
         for (int g=0; g<K; ++g){
           if (g==refg) continue;
           for (int r=0;r<Rw;++r){
             int idx = pos*Rw + r;
             double prj=0.0;
             for (int j=1;j<=m;++j){
               double u = clamp01(aj1(g,j),   EPS);
               double vref = clamp01(aj1(refg,j), EPS);
               if (gg == g)     prj += W[r][j]*dt[j]*( IFg(i,j)/u );
               else if (gg==refg) prj += W[r][j]*dt[j]*( -IFg(i,j)/vref );
             }
             v[idx] = prj;
           }
           ++pos;
         }
         for (int a=0;a<NR_rr;++a){
           if (v[a]==0.0) continue;
           for (int b=a;b<NR_rr;++b){
             double s = v[a]*v[b];
             if (s!=0.0){
               cov_rr(a,b)+=s;
               if (b!=a) cov_rr(b,a)+=s;
             }
           }
         }
       }
     }
   }

   // =======================
   // 3) or_aj1（参照群 vs 他群）：**幾何平均**（logit差の平均）に統一
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

       double theta = 0.0;
       for (int j=1;j<=m;++j){
         double u = clamp01(aj1(g,j),   EPS);
         double v = clamp01(aj1(refg,j),EPS);
         double logit_u = std::log(u) - std::log(1.0-u);
         double logit_v = std::log(v) - std::log(1.0-v);
         theta += W[r][j] * dt[j] * (logit_u - logit_v);
       }
       double est = std::exp(theta);
       or_est[idx] = est;

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
       double se_log = std::sqrt(var);
       or_se[idx] = se_log;
       if (se_log>0){
         or_lwr[idx] = std::exp(theta - zq*se_log);
         or_upr[idx] = std::exp(theta + zq*se_log);
         double z    = theta / se_log;
         or_z[idx]   = z;
         or_p[idx]   = R::pchisq(z*z, 1.0, false, false);
       } else {
         or_lwr[idx]=or_upr[idx]=or_z[idx]=or_p[idx]=NA_REAL;
       }
     }
     ++pos;
   }

    NumericMatrix cov_or(NR_rr, NR_rr);
   {
     for (int gg=0; gg<K; ++gg){
       NumericMatrix IFg = if_aj1[gg];
       const int ng = IFg.nrow();
       for (int i=0;i<ng;++i){
         std::vector<double> v(NR_rr, 0.0);
         int pos = 0;
         for (int g=0; g<K; ++g){
           if (g==refg) continue;
           for (int r=0; r<Rw; ++r){
             int idx = pos*Rw + r;
             double prj=0.0;
             for (int j=1;j<=m;++j){
               double u = clamp01(aj1(g,j),   EPS);
               double vref = clamp01(aj1(refg,j), EPS);
               if (gg == g)      prj += W[r][j]*dt[j]*( IFg(i,j)/(u*(1.0-u)) );
               else if (gg==refg) prj += W[r][j]*dt[j]*( -IFg(i,j)/(vref*(1.0-vref)) );
             }
             v[idx] = prj;
           }
           ++pos;
         }
         for (int a=0;a<NR_rr;++a){
           if (v[a]==0.0) continue;
           for (int b=a;b<NR_rr;++b){
             double s = v[a]*v[b];
             if (s!=0.0){
               cov_or(a,b)+=s;
               if (b!=a) cov_or(b,a)+=s;
             }
           }
         }
       }
     }
   }

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
     _["group"]     = or_group,
     _["weight_id"] = or_weight_id,
     _["est"] = or_est,
     _["se"]  = or_se,
     _["lwr"] = or_lwr,
     _["upr"] = or_upr,
     _["z"]   = or_z,
     _["p"]   = or_p
   );

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

namespace {
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
 //' @param n_simulation integer, Monte Carlo draws for MVN-based p-value (default 20000).
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
 List calculateMaxCombo(NumericVector z, NumericMatrix cov, int n_simulation = 20000) {
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
   for (int i=0;i<R;++i) corr(i,i) += 1e-8;
   bool ok = chol_lower(corr, L);

   double p_mc = NA_REAL;
   if (ok) {
     int exceed = 0;
     for (int b=0;b<n_simulation;++b){
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
     p_mc = (exceed + 1.0) / (n_simulation + 1.0);
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


//' Convenience wrapper: build maxCombo directly from calculateWeightedAverage() result
 //'
 //' @param out_cwa      List output from calculateWeightedAverage()
 //' @param measure "diff", "rr", or "or"  （diff_aj1 / rr_aj1 / or_aj1 に対応）
 //' @param n_simulation   Monte Carlo draws (default 20000)
 //'
 //' @return list(per_weight = data.frame, maxcombo = list(stat,p,R,corr))
 //'
 //' @details
 //' - それぞれの measure で、per-weight（× per-contrast or per-group）の Z を連結して maxCombo を計算します。
 //' - 共分散行列は out_cwa$cov_diff / cov_rr / cov_or を使用します。
 //' - NAのZは自動で落とします（共分散も対応成分で抜粋）。
 // [[Rcpp::export]]
 List callMaxCombo(List out_cwa, std::string measure = "diff", int n_simulation = 20000) {
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
   List mc = core(z, cov, n_simulation);

   return List::create(
     _["per_weight"] = df,
     _["maxcombo"]   = mc
   );
 }


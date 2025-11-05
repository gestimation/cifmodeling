#include <Rcpp.h>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

namespace {
const double EPS = 1e-9;

inline bool feq(double a, double b, double eps=EPS){
  double s = std::max({1.0, std::fabs(a), std::fabs(b)});
  return std::fabs(a-b) <= eps * s;
}
inline double clamp01(double x, double eps=1e-12){
  if (!R_finite(x)) return eps;
  if (x < eps) return eps;
  if (x > 1.0 - eps) return 1.0 - eps;
  return x;
}
inline double qnorm_from_conf(double conf) {
  return R::qnorm(0.5 + conf/2.0, 0.0, 1.0, /*lower*/1, /*log*/0);
}
inline double sdiv(double num, double den){
  if (!R_finite(num) || !R_finite(den) || std::fabs(den) < 1e-12) return 0.0;
  return num/den;
}

struct RecW {
  double t;
  int    eps;   // 0=censor, 1=cause1, 2+=competing
  int    g;     // stratum label (>=1)
  int    id;    // original index
  double w;     // frequency weight (>=0)
};
} // namespace

// [[Rcpp::export]]
Rcpp::List calculateAJ_Rcpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::Nullable<Rcpp::NumericVector> w       = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> strata  = R_NilValue,
    std::string error = "greenwood",
    std::string conf_type = "arcsin",
    bool return_if = true,
    double conf_int = 0.95
){
  const int N = t.size();

  if (epsilon.size() != N) stop("Length mismatch: epsilon.");
  if (N == 0) stop("Empty input.");

  // weights
  NumericVector W(N, 1.0);
  if (!w.isNull()) {
    NumericVector wv = as<NumericVector>(w);
    if (wv.size() == N) W = wv;
    else if (wv.size() != 0) stop("w must have length 0 or N.");
  }

  // helpers (local)
  auto to_lower = [](std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
  };
  auto canon_conf = [&](std::string s){
    s = to_lower(s);
    if (s=="n" || s=="none") return std::string("none");
    if (s=="p" || s=="plain" || s=="linear") return std::string("plain");
    if (s=="log") return std::string("log");
    if (s=="logit") return std::string("logit");
    if (s=="log-log") return std::string("log-log");
    if (s=="a" || s=="arcsin" || s=="arcsine-square root") return std::string("arcsine-square root");
    return std::string("arcsine-square root");
  };
  const std::string CONF = canon_conf(conf_type);

  // --- NEW: normalize error & flags ---
  auto to_lower2 = [](std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
  };
  const std::string ERROR = to_lower2(error);

  bool error_tsiatis = (ERROR == "tsiatis");
  bool error_aalen   = (ERROR == "aalen");
  bool error_delta   = (ERROR == "delta");
  bool error_if      = (ERROR == "if") || return_if;
  bool error_cif     = error_aalen || error_delta || error_if;

  auto any_of_int = [](const IntegerVector& v, std::function<bool(int)> pred){
    for (int x : v) if (pred(x)) return true;
    return false;
  };
  const bool has_cause1    = any_of_int(epsilon, [](int e){ return e==1; });
  const bool has_competing = any_of_int(epsilon, [](int e){ return e>=2; });
  const bool return_aj     = (has_cause1 && (has_competing || error_cif));
  const bool need_any_cif  = return_aj;

  std::vector<double> combined_aj;  combined_aj.reserve(N);
  std::vector<double> combined_km;  combined_km.reserve(N);

  // strata
  IntegerVector G(N, 1);
  Rcpp::CharacterVector levs;
  if (!strata.isNull()){
    IntegerVector gs = as<IntegerVector>(strata);
    if (gs.size() == N) {
      G = gs;
      if (gs.hasAttribute("levels")) levs = gs.attr("levels");
    } else if (gs.size() != 0) {
      stop("strata must have length 0 or N.");
    }
  }

  // basic checks
  for (int i = 0; i < N; ++i){
    if (!R_finite(t[i]) || t[i] < 0.0) stop("t must be finite and nonnegative.");
    if (epsilon[i] < 0) stop("epsilon must be >=0.");
    if (G[i] < 1) stop("strata labels must be >=1.");
    if (!R_finite(W[i]) || W[i] < 0.0) stop("weights must be nonnegative and finite.");
  }

  // unique strata codes (respect factor levels if present)
  IntegerVector guniq;
  if (levs.size() > 0) {
    int gmax = 0; for (int i=0;i<N;++i) if (G[i] > gmax) gmax = G[i];
    std::vector<char> present(gmax+1, 0);
    for (int i=0;i<N;++i) if (G[i] >= 1) present[G[i]] = 1;

    std::vector<int> tmp; tmp.reserve(levs.size());
    for (int k=0; k<levs.size(); ++k) {
      int code = k + 1;
      if (code < (int)present.size() && present[code]) tmp.push_back(code);
    }
    guniq = wrap(tmp);
  } else {
    guniq = sort_unique(G);
  }
  const int K = guniq.size();

  if (has_competing) {
    if (error_tsiatis || (!error_delta && !error_if && !error_aalen)) {
      error_tsiatis = false;
      error_aalen   = false;
      error_delta   = true;
      error_if      = error_if;
      error_cif     = true;
    }
  }

  const std::string CONF_CIF = (
    has_competing ?
  ( (CONF=="plain" || CONF=="logit" || CONF=="arcsine-square root") ?
  CONF : std::string("arcsine-square root") )
      : CONF );

  int gmax2 = 0; for (int i=0;i<N;++i) if (G[i] > gmax2) gmax2 = G[i];
  std::vector<int> lab2k(gmax2 + 1, -1);
  for (int k=0;k<K;++k) lab2k[guniq[k]] = k;

  std::vector< std::vector<int> > idx_of_strata(K);
  for (int i=0;i<N;i++){
    int kg = lab2k[G[i]];
    if (kg<0) stop("Internal strata mapping failed.");
    idx_of_strata[kg].push_back(i);
  }

  // combined outputs
  std::vector<double> combined_times, combined_surv, combined_n_risk, combined_n_event, combined_n_censor;
  std::vector<double> combined_std_err, combined_std_err_km, combined_std_err_aj, combined_high, combined_low;
  std::vector<int>    combined_n_stratum, combined_u_stratum;
  Rcpp::List          IF_AJ_output(K);
  IntegerVector       strata_levels_out = guniq;

  // one-shot quantile
  const double z = qnorm_from_conf(conf_int);
  const double PI = 3.14159265358979323846;

  // per-stratum work
  for (int kg=0; kg<K; ++kg){
    const std::vector<int>& ids = idx_of_strata[kg];
    const int n_g = (int)ids.size();

    if (n_g==0){
      combined_n_stratum.push_back(0);
      combined_u_stratum.push_back(0);
      IF_AJ_output[kg] = NumericMatrix(0,0);
      continue;
    }

    // make and sort records within stratum
    std::vector<RecW> v; v.reserve(n_g);
    for (int ii=0; ii<n_g; ++ii){
      int i = ids[ii];
      v.push_back({ (double)t[i], (int)epsilon[i], (int)G[i], i, (double)W[i] });
    }
    std::sort(v.begin(), v.end(), [](const RecW& a, const RecW& b){
      if (!feq(a.t, b.t)) return a.t < b.t;
      return (a.eps >= 1) > (b.eps >= 1);
    });

    // unique time grid and weighted counts
    std::vector<double> t_all, w_event_all, w_censor_all, w_total_all;
    t_all.reserve(n_g); w_event_all.reserve(n_g); w_censor_all.reserve(n_g); w_total_all.reserve(n_g);

    for (size_t p=0; p<v.size(); ){
      double curT = v[p].t;
      double sumE = 0.0, sumC = 0.0, sumTot = 0.0;
      size_t q=p;
      while (q<v.size() && feq(v[q].t, curT)){
        if (v[q].eps >= 1) sumE += v[q].w; else sumC += v[q].w;
        sumTot += v[q].w;
        ++q;
      }
      t_all.push_back(curT);
      w_event_all.push_back(sumE);
      w_censor_all.push_back(sumC);
      w_total_all.push_back(sumTot);
      p=q;
    }
    const int Uall = (int)t_all.size();

    // unweighted counts by unique time
    std::vector<double> cnt_event_all(Uall, 0.0),
                        cnt_censor_all(Uall, 0.0);
    {
      size_t p2 = 0;
      for (int j = 0; j < Uall; ++j) {
        double curT = t_all[j];
        double e = 0.0, c = 0.0, tot = 0.0;
        while (p2 < v.size() && feq(v[p2].t, curT)) {
          if (v[p2].eps >= 1) e += 1.0; else c += 1.0;
          ++p2;
        }
        cnt_event_all[j]  = e;
        cnt_censor_all[j] = c;
      }
    }

    std::vector<double> Y_all(Uall, 0.0);
    {
      double suf = 0.0;
      for (int j = Uall - 1; j >= 0; --j) { suf += w_total_all[j]; Y_all[j] = suf; }
    }

    // cause-1 weighted counts on the same grid (only if any CIF is needed)
    std::vector<double> w_event1_all(Uall, 0.0);
    if (need_any_cif) {
      for (size_t p=0; p<v.size(); ){
        double curT = v[p].t;
        double sum1 = 0.0;
        size_t q=p;
        while (q<v.size() && feq(v[q].t, curT)){
          if (v[q].eps == 1) sum1 += v[q].w;
          ++q;
        }
        int idx = (int)(std::lower_bound(t_all.begin(), t_all.end(), curT) - t_all.begin());
        if (idx>=0 && idx<Uall && feq(t_all[idx], curT)) w_event1_all[idx] = sum1;
        p=q;
      }
    }

    // event indices (ANY / cause-1)
    std::vector<int> ev_any_idx, ev_c1_idx;
    ev_any_idx.reserve(Uall); ev_c1_idx.reserve(Uall);
    for (int j=0;j<Uall;++j){
      if (w_event_all[j]  > 0.0) ev_any_idx.push_back(j);
      if (need_any_cif && w_event1_all[j] > 0.0) ev_c1_idx.push_back(j);
    }
    const int M_any = (int)ev_any_idx.size();
    const int M_c1  = (need_any_cif ? (int)ev_c1_idx.size() : 0);

    // AJ for cause-1 on ANY grid (only if CIF is needed)
    std::vector<double> dl_any(M_any,0.0), dl_c1, aj1_any;
    if (need_any_cif){
      dl_c1.resize(M_any, 0.0);
      aj1_any.resize(M_any, 0.0);

      for (int m=0;m<M_any;++m){
        int j = ev_any_idx[m];
        double Yj = Y_all[j], dj = w_event_all[j];
        dl_any[m] = (Yj>0.0) ? (dj/Yj) : 0.0;
      }
      std::vector<char> is_c1_at_all(Uall, 0);
      for (int u=0; u<M_c1; ++u) is_c1_at_all[ ev_c1_idx[u] ] = 1;
      for (int m=0;m<M_any;++m){
        int j = ev_any_idx[m];
        double Yj = Y_all[j];
        dl_c1[m] = (is_c1_at_all[j] && Yj>0.0) ? (w_event1_all[j] / Yj) : 0.0;
      }
      double Sprev = 1.0, Fprev = 0.0;
      for (int m=0;m<M_any;++m){
        double dF1 = Sprev * dl_c1[m];
        Fprev += dF1;
        aj1_any[m] = Fprev;
        Sprev *= (1.0 - dl_any[m]);
      }
    } else {
      // 使わない場合でも dl_any は上で必要（KM のため）
      for (int m=0;m<M_any;++m){
        int j = ev_any_idx[m];
        double Yj = Y_all[j], dj = w_event_all[j];
        dl_any[m] = (Yj>0.0) ? (dj/Yj) : 0.0;
      }
    }

    std::vector<double> S_any(M_any, 1.0), varKM_any(M_any, 0.0);
    double accKM = 0.0;
    for (int m=0; m<M_any; ++m){
      int j = ev_any_idx[m];
      double Yj = Y_all[j], dj = w_event_all[j];
      double Sprev = (m==0? 1.0 : S_any[m-1]);
      double fac   = (Yj>0.0) ? (1.0 - dj/Yj) : 1.0;
      S_any[m]     = Sprev * fac;

      if (dj>0.0 && Yj>0.0){
        if (error_tsiatis) accKM += dj / (Yj*Yj);
        else {
          if (Yj > dj) accKM += dj / (Yj*(Yj-dj));
          else         accKM  = std::numeric_limits<double>::infinity();
        }
      }
      double S2 = S_any[m] * S_any[m];
      varKM_any[m] = S2 * accKM;
    }

    std::vector<int> all2any(Uall, 0);
    {
      int c_any = 0;
      for (int j=0; j<Uall; ++j){
        if (w_event_all[j] > 0.0) ++c_any;
        all2any[j] = c_any;
      }
    }
    std::vector<double> S_all(Uall, 1.0), SE_all(Uall, 0.0);
    for (int j=0;j<Uall;++j){
      int c = all2any[j];
      if (c>0){
        S_all[j] = S_any[c-1];
        double v = varKM_any[c-1];
        SE_all[j] = (R_finite(v) && v>=0.0) ? std::sqrt(v) : R_PosInf;
      } else { S_all[j]=1.0; SE_all[j]=0.0; }
    }

    std::vector<double> F1_all(Uall, 0.0);
    if (need_any_cif){
      for (int j = 0; j < Uall; ++j){
        int c = all2any[j];
        F1_all[j] = (c>0 ? aj1_any[c-1] : 0.0);
      }
    }

    std::vector<double> se_cif_if(Uall, 0.0);
    NumericMatrix IF_AJ_all;
    if (error_if) {
      NumericMatrix IF_AJ_any(n_g, M_any);
      for (int r=0; r<n_g; ++r){
        int i = ids[r];
        double Ti = t[i], wi = W[i];
        int Ei = epsilon[i];
        double xs_prev = 0.0, xf_prev = 0.0;
        for (int m=0; m<M_any; ++m){
          int j_all = ev_any_idx[m];
          double tj = t_all[j_all];
          double Yj = Y_all[j_all];
          double dlam_any = dl_any[m];
          double dlam_c1  = (need_any_cif ? dl_c1[m] : 0.0);

          double dm_any = 0.0, dm_c1 = 0.0;
          if (Ti + EPS >= tj){
            dm_any = wi * ((feq(Ti,tj) && (Ei >= 1)) ? (1.0 - dlam_any) : (0.0 - dlam_any));
            dm_c1  = wi * ((feq(Ti,tj) && (Ei == 1)) ? (1.0 - dlam_c1 ) : (0.0 - dlam_c1 ));
          }
          double xs_now = xs_prev + sdiv(dm_any, Yj);
          double Sprev  = (m==0? 1.0 : S_any[m-1]);
          double term   = sdiv(dm_c1, Yj) - xs_prev * dlam_c1;
          double xf_now = xf_prev + Sprev * term;

          xs_prev = xs_now; xf_prev = xf_now;
          IF_AJ_any(r, m) = xf_now;
        }
      }

      double sum_w = 0.0; for (int ii = 0; ii < n_g; ++ii) sum_w += W[ ids[ii] ];
      const double denom = (sum_w > 0.0 ? sum_w : (double)n_g);

      std::vector<long double> ss(Uall, 0.0L);
      if (return_if) IF_AJ_all = NumericMatrix(n_g, Uall);

      for (int j = 0; j < Uall; ++j) {
        int c = all2any[j];
        for (int r = 0; r < n_g; ++r) {
          double v = (c==0 ? 0.0 : IF_AJ_any(r, c-1));
          ss[j] += (long double)v * v;
          if (return_if) IF_AJ_all(r,j) = v;
        }
      }
      for (int j=0; j<Uall; ++j)
        se_cif_if[j] = std::sqrt((double)(ss[j] / ((long double)denom * denom)));
    }
    IF_AJ_output[kg] = (error_if && return_if) ? IF_AJ_all : Rcpp::NumericMatrix(0,0);

    std::vector<double> se_cif_aalen(Uall, 0.0), se_cif_delta(Uall, 0.0);

    if (error_aalen || error_delta) {
      std::vector<double> km_time(M_any);
      for (int m=0;m<M_any;++m) km_time[m] = t_all[ ev_any_idx[m] ];

      std::vector<double> var_aalen, var_delta;
      if (error_aalen) var_aalen.assign(M_c1, 0.0);
      if (error_delta) var_delta.assign(M_c1, 0.0);

      std::vector<double> c1_time(M_c1), c1_n1(M_c1), c1_n2(M_c1), c1_Y(M_c1), c1_F(M_c1);
      for (int u=0; u<M_c1; ++u){
        int j = ev_c1_idx[u];
        c1_time[u] = t_all[j];
        c1_n1[u]   = w_event1_all[j];
        c1_n2[u]   = w_event_all[j] - w_event1_all[j];
        c1_Y[u]    = Y_all[j];
        int count_any_before = 0;
        for (int m=0;m<M_any;++m){ if (ev_any_idx[m]==j){ count_any_before = m+1; break; } }
        c1_F[u] = (count_any_before>0 ? aj1_any[count_any_before-1] : 0.0);
      }

      if (error_aalen) {
        double first_cum=0.0, second_cum=0.0, third_cum=0.0;
        for (int i=0;i<M_c1;++i){
          int idx_count = 0;
          for (int m=0;m<M_any;++m){ if (km_time[m] < c1_time[i] - EPS) ++idx_count; }
          int idx_cap = std::min(M_c1-1, idx_count);
          if (idx_cap != 0){
            double Spre = (idx_count==0? 1.0 : S_any[idx_count-1]);
            double F_idx = c1_F[idx_cap];
            double Fi    = c1_F[i];

            double n1 = c1_n1[i], n2 = c1_n2[i], Y  = c1_Y[i];
            double Ym1   = std::max(Y-1.0, 1e-12);
            double Ymn12 = std::max(Y - n1 - n2, 1e-12);
            double Y2    = std::max(Y*Y, 1e-12);

            first_cum  += sdiv( (F_idx-Fi)*(F_idx-Fi) * (n1+n2),  Ym1*Ymn12 );
            second_cum += sdiv( (Spre*Spre) * n1 * (Y-n1),        (Y2*Ym1)   );
            third_cum  += sdiv( (F_idx-Fi)*Spre * n1*(Y-n1),      (Y*(Ymn12)*Ym1) );

            var_aalen[i] = first_cum + second_cum - 2.0*third_cum;
          } else var_aalen[i] = 0.0;
        }
      }

      if (error_delta) {
        double first_cum=0.0, second_cum=0.0, third_cum=0.0;
        for (int i=0;i<M_c1;++i){
          int idx_count = 0;
          for (int m=0;m<M_any;++m){ if (km_time[m] < c1_time[i] - EPS) ++idx_count; }
          int idx_cap = std::min(M_c1-1, idx_count);
          if (idx_cap != 0){
            double Spre = (idx_count==0? 1.0 : S_any[idx_count-1]);
            double F_idx = c1_F[idx_cap];
            double Fi    = c1_F[i];

            double n1 = c1_n1[i], n2 = c1_n2[i], Y  = c1_Y[i];
            double Ymn12 = std::max(Y - n1 - n2, 1e-12);
            double Y2    = std::max(Y*Y, 1e-12);
            double Y3    = std::max(Y*Y*Y, 1e-12);

            first_cum  += sdiv( (F_idx-Fi)*(F_idx-Fi) * (n1+n2),  (Y*Ymn12) );
            second_cum += sdiv( (Spre*Spre) * n1 * (Y-n1),        Y3        );
            third_cum  += sdiv( (F_idx-Fi)*Spre * n1,             Y2        );

            var_delta[i] = first_cum + second_cum - 2.0*third_cum;
          } else var_delta[i] = 0.0;
        }
      }

      // carry to ALL grid
      std::vector<int> all2c1(Uall, 0);
      {
        int c_c1 = 0;
        for (int j = 0; j < Uall; ++j) {
          if (w_event1_all[j] > 0.0) ++c_c1;
          all2c1[j] = c_c1;
        }
      }

      for (int j = 0; j < Uall; ++j) {
        int c = all2c1[j];          // 0..M_c1
        if (c == 0) {
          if (error_aalen && M_c1 > 0) se_cif_aalen[j] = std::sqrt(std::max(0.0, var_aalen[0]));
          if (error_delta && M_c1 > 0) se_cif_delta[j] = std::sqrt(std::max(0.0, var_delta[0]));
        } else {
          int use = std::min(c, M_c1);  // 末尾越え防止
          if (error_aalen) se_cif_aalen[j] = std::sqrt(std::max(0.0, var_aalen[use - 1]));
          if (error_delta) se_cif_delta[j] = std::sqrt(std::max(0.0, var_delta[use - 1]));
        }
      }
    }

    if (Uall >= 2) {
      if (error_aalen) {
        std::vector<double> tmp = se_cif_aalen;          // copy
        for (int j = 0; j < Uall - 1; ++j) se_cif_aalen[j] = tmp[j + 1];
        se_cif_aalen[Uall - 1] = tmp[Uall - 1];          // keep last (no t_max+1)
      }
      if (error_delta) {
        std::vector<double> tmp = se_cif_delta;          // copy
        for (int j = 0; j < Uall - 1; ++j) se_cif_delta[j] = tmp[j + 1];
        se_cif_delta[Uall - 1] = tmp[Uall - 1];          // keep last
      }
    }

    std::vector<double> SE_aj(Uall, 0.0);
    if      (error_aalen) SE_aj = se_cif_aalen;
    else if (error_delta) SE_aj = se_cif_delta;
    else if (error_if)    SE_aj = se_cif_if;
    else                  SE_aj = SE_all;

    std::vector<double> high_all(Uall, NA_REAL), low_all(Uall, NA_REAL);
    if (CONF_CIF != "none" || CONF != "none"){
      for (int j=0;j<Uall;++j){
        if (return_aj){
          double Fj  = clamp01(F1_all[j], 1e-15);
          double SEj = SE_aj[j];
          if (!R_finite(SEj)){ high_all[j]=low_all[j]=NA_REAL; continue; }
          if (CONF_CIF=="plain"){
            double lo = std::max(0.0, Fj - z*SEj);
            double hi = std::min(1.0, Fj + z*SEj);
            low_all[j]  = 1.0 - hi;
            high_all[j] = 1.0 - lo;
          } else if (CONF_CIF=="logit"){
            double se = SEj / std::max(1e-15, Fj*(1.0 - Fj));
            double lo = Fj / ( Fj + (1.0 - Fj)*std::exp( z*se) );
            double hi = Fj / ( Fj + (1.0 - Fj)*std::exp(-z*se) );
            lo = std::max(0.0, std::min(1.0, lo));
            hi = std::max(0.0, std::min(1.0, hi));
            low_all[j]  = 1.0 - hi;
            high_all[j] = 1.0 - lo;
          } else { // if (CONF=="arcsine-square root" || CONF=="a" || CONF=="arcsin"){
            double denom = 2.0 * std::sqrt( std::max(1e-15, Fj*(1.0-Fj)) );
            double se    = (denom>0.0 ? SEj/denom : R_PosInf);
            double ang   = std::asin( std::sqrt(Fj) );
            double lo    = std::sin( std::max(0.0, ang - z*se) );
            double hi    = std::sin( std::min(PI/2, ang + z*se) );
            lo = lo*lo; hi = hi*hi;
            low_all[j]  = 1.0 - hi;
            high_all[j] = 1.0 - lo;
          }
        } else {
          double Sj  = clamp01(S_all[j], 1e-15);
          double SEj = SE_all[j];
          if (!R_finite(SEj) || Sj<=0.0 || Sj>=1.0){ high_all[j]=low_all[j]=NA_REAL; continue; }
          if (CONF=="plain"){
            double lo = std::max(0.0, Sj - z*SEj);
            double hi = std::min(1.0, Sj + z*SEj);
            low_all[j]=lo; high_all[j]=hi;
          } else if (CONF=="log"){
            double se = SEj / Sj;
            double lo = Sj * std::exp(-z*se);
            double hi = Sj * std::exp( z*se);
            low_all[j]=std::max(0.0, lo);
            high_all[j]=std::min(1.0, hi);
          } else if (CONF=="log-log"){
            double se = SEj / ( Sj * std::fabs(std::log(Sj)) );
            double lo = std::exp( - std::exp( std::log(-std::log(Sj)) + z*se ) );
            double hi = std::exp( - std::exp( std::log(-std::log(Sj)) - z*se ) );
            low_all[j]=std::max(0.0, std::min(1.0, lo));
            high_all[j]=std::max(0.0, std::min(1.0, hi));
          } else if (CONF=="logit"){
            double se = SEj / ( Sj * (1.0 - Sj) );
            double lo = Sj / ( Sj + (1.0 - Sj)*std::exp( z*se) );
            double hi = Sj / ( Sj + (1.0 - Sj)*std::exp(-z*se) );
            low_all[j]=std::max(0.0, std::min(1.0, lo));
            high_all[j]=std::max(0.0, std::min(1.0, hi));
          } else { // arcsine-square root
            double denom = 2.0 * std::sqrt( std::max(1e-15, Sj*(1.0-Sj)) );
            double se    = (denom>0.0 ? SEj/denom : R_PosInf);
            double ang   = std::asin( std::sqrt(Sj) );
            double lo    = std::sin( std::max(0.0, ang - z*se) );
            double hi    = std::sin( std::min(PI/2, ang + z*se) );
            low_all[j]   = lo*lo;
            high_all[j]  = hi*hi;
          }
        }
      }
    }

    combined_times.reserve(combined_times.size() + Uall);
    combined_surv .reserve(combined_surv .size() + Uall);
    combined_n_risk.reserve(combined_n_risk.size() + Uall);
    combined_n_event.reserve(combined_n_event.size() + Uall);
    combined_n_censor.reserve(combined_n_censor.size() + Uall);
    combined_std_err.reserve(combined_std_err.size() + Uall);
    combined_std_err_km.reserve(combined_std_err_km.size() + Uall);
    combined_std_err_aj.reserve(combined_std_err_aj.size() + Uall);
    combined_low.reserve(combined_low.size() + Uall);
    combined_high.reserve(combined_high.size() + Uall);

    for (int j=0;j<Uall;++j){
      combined_times.push_back(t_all[j]);
      combined_n_risk.push_back(Y_all[j]);
      combined_n_event.push_back( cnt_event_all[j] );
      combined_n_censor.push_back(cnt_censor_all[j]);
      if (return_aj) combined_surv.push_back(1.0 - F1_all[j]);
      else           combined_surv.push_back(S_all[j]);
      if (return_aj) combined_std_err.push_back(SE_aj[j]);
      else           combined_std_err.push_back(SE_all[j]);
      if (return_aj) combined_aj.push_back(F1_all[j]);
      combined_std_err_km.push_back(SE_all[j]);
      if (return_aj) combined_std_err_aj.push_back(SE_aj[j]);
      combined_low.push_back(low_all[j]);
      combined_high.push_back(high_all[j]);
    }
    combined_n_stratum.push_back(n_g);
    combined_u_stratum.push_back(Uall);
  }

  std::string out_type, out_method;
  if (return_aj){
    out_type   = "aalen-johansen";
    out_method = "Aalen-Johansen";
  } else {
    out_type   = "kaplan-meier";
    out_method = "Kaplan-Meier";
  }

  SEXP strata_out, strata_levels_out_sexp;
  if (K==1){
    strata_out = Rcpp::IntegerVector(0);
    strata_levels_out_sexp = Rcpp::IntegerVector(0);
  } else {
    strata_out = wrap(combined_u_stratum);
    strata_levels_out_sexp = wrap(strata_levels_out);
  }

  Rcpp::NumericVector aj_out, std_err_aj_out, km_out, std_err_km_out;
  if (has_competing) {
    aj_out        = wrap(combined_aj);
    std_err_aj_out= wrap(combined_std_err_aj);
  } else {
    aj_out        = Rcpp::NumericVector(0);
    std_err_aj_out= Rcpp::NumericVector(0);
  }
  if (has_competing) {
    std_err_km_out = Rcpp::NumericVector(0);
    km_out         = Rcpp::NumericVector(0);
  } else {
    std_err_km_out = wrap(combined_std_err);
    km_out = Rcpp::NumericVector(0);
  }

  if (has_competing) {
    std_err_aj_out = wrap(
      [&](){
        return Rcpp::wrap(combined_std_err_aj); }() );
  }

  Rcpp::List out = Rcpp::List::create(
    _["time"]               = wrap(combined_times),
    _["surv"]               = wrap(combined_surv),
    _["aj"]                 = aj_out,
    _["km"]                 = km_out,
    _["n.risk"]             = wrap(combined_n_risk),
    _["n"]                  = wrap(combined_n_stratum),
    _["n.event"]            = wrap(combined_n_event),
    _["n.censor"]           = wrap(combined_n_censor),
    _["unweighted.n"]       = wrap(combined_n_stratum),
    _["std.err"]            = wrap(combined_std_err),
    _["std.err.aj"]         = std_err_aj_out,
    _["std.err.km"]         = std_err_km_out,
    _["low"]                = wrap(combined_low),
    _["high"]               = wrap(combined_high),
    _["conf.type"]          = (return_aj ? CONF_CIF : CONF),
    _["strata"]             = strata_out,
    _["strata.levels"]      = strata_levels_out_sexp,
    _["type"]               = out_type,
    _["method"]             = out_method,
    _["influence.function"] = IF_AJ_output
  );
  return out;
}

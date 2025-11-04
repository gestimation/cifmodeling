#include <Rcpp.h>
using namespace Rcpp;
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List calculateKM(Rcpp::NumericVector t, Rcpp::IntegerVector d,
                       Rcpp::NumericVector w = Rcpp::NumericVector::create(),
                       Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                       std::string error = "greenwood") {
  Rcpp::List survfit_list;

  std::vector<double> combined_times;
  std::vector<double> combined_surv;
  std::vector<int> combined_n_risk;
  std::vector<int> combined_n_event;
  std::vector<int> combined_n_censor;
  std::vector<double> combined_std_err;
  std::vector<int> combined_n_stratum;
  std::vector<int> combined_u_stratum;
  Rcpp::IntegerVector strata_levels_out;

  if ((strata.size() == 0 || Rcpp::unique(strata).size() == 1) && (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1))) {

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());

    int u = unique_times.size();
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector weighted_n_event(u);
    Rcpp::IntegerVector weighted_n_censor(u);
    Rcpp::NumericVector std_err(u);

    int n_stratum = t.size();
    int n = t.size();
    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];
      double weighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_i++;
        }
      }
      weighted_n_risk[i] = weighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events++;
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event[i] += weighted_events;
      } else {
        km_i[i] = 1;
      }

      if (i > 0) {
        weighted_n_censor[i-1] = weighted_n_risk[i-1] - weighted_n_risk[i] - weighted_n_event[i-1];
      }
      if (i == u-1) {
        weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i ++;
          }
        }
        if (n_i > d_i) {
          if (error == "tsiatis") {
            sum_se += (d_i / (n_i * n_i));
          } else if (error == "greenwood") {
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        } else {
          sum_se = std::numeric_limits<double>::infinity();
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
    combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);

  } else if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {

    Rcpp::NumericVector unique_times = Rcpp::unique(t);
    std::sort(unique_times.begin(), unique_times.end());

    int u = unique_times.size();
    Rcpp::NumericVector km(u);
    Rcpp::NumericVector km_i(u);
    Rcpp::IntegerVector weighted_n_risk(u);
    Rcpp::IntegerVector weighted_n_event(u);
    Rcpp::IntegerVector weighted_n_censor(u);
    Rcpp::NumericVector std_err(u);

    int n_stratum = t.size();
    int n = t.size();
    for (int i = 0; i < u; ++i) {
      double time = unique_times[i];
      double weighted_n_i = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] >= time) {
          weighted_n_i += w[j];
        }
      }
      weighted_n_risk[i] = weighted_n_i;

      double weighted_events = 0;
      for (int j = 0; j < n; ++j) {
        if (t[j] == time && d[j] == 1) {
          weighted_events += w[j];
        }
      }

      if (weighted_n_risk[i] > 0) {
        km_i[i] = 1 - weighted_events / weighted_n_risk[i];
        weighted_n_event[i] += weighted_events;
      } else {
        km_i[i] = 1;
      }

      if (i > 0) {
        weighted_n_censor[i-1] = weighted_n_risk[i-1] - weighted_n_risk[i] - weighted_n_event[i-1];
      }
      if (i == u-1) {
        weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
      }

      if (i == 0) {
        km[i] = km_i[i];
      } else {
        km[i] = km[i - 1] * km_i[i];
      }

      double sum_se = 0;
      for (int j = 0; j <= i; ++j) {
        double n_i = weighted_n_risk[j];
        double d_i = 0;
        for (int k = 0; k < n; ++k) {
          if (t[k] == unique_times[j] && d[k] == 1) {
            d_i += w[k];
          }
        }
        if (n_i > d_i) {
          if (error == "tsiatis") {
            sum_se += (d_i / (n_i * n_i));
          } else if (error == "greenwood") {
            sum_se += (d_i / (n_i * (n_i - d_i)));
          }
        } else {
          sum_se = std::numeric_limits<double>::infinity();
        }
      }
      std_err[i] = sqrt(sum_se);
    }

    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
    combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);

  } else if (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1)) {

    Rcpp::IntegerVector strata_vec = strata;

    Rcpp::IntegerVector levels = Rcpp::sort_unique(strata_vec);

    for (int li = 0; li < levels.size(); ++li) {
      int level = levels[li];
      Rcpp::LogicalVector strata_condition = (strata_vec == level);

      Rcpp::NumericVector t_selected = t[strata_condition];
      Rcpp::IntegerVector d_selected = d[strata_condition];

      int n_stratum = t_selected.size();

      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u_stratum = unique_times.size();

      Rcpp::NumericVector km(u_stratum), km_i(u_stratum), std_err(u_stratum);
      Rcpp::IntegerVector weighted_n_event(u_stratum), weighted_n_censor(u_stratum);
      Rcpp::IntegerVector weighted_n_risk(u_stratum);

      for (int j = 0; j < u_stratum; ++j) {
        double time = unique_times[j];
        int n_at_risk = 0, n_events = 0;

        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) n_at_risk++;
          if (t_selected[k] == time && d_selected[k] == 1) n_events++;
        }
        weighted_n_risk[j] = n_at_risk;
        weighted_n_event[j] = n_events;

        km_i[j] = (n_at_risk > 0) ? (1.0 - (double)n_events / (double)n_at_risk) : 1.0;
        km[j]   = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        if (j > 0) weighted_n_censor[j-1] = weighted_n_risk[j-1] - weighted_n_risk[j] - weighted_n_event[j-1];
        if (j == u_stratum - 1) weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];

        double sum_se = 0.0;
        for (int k = 0; k <= j; ++k) {
          double n_i = (double)weighted_n_risk[k];
          double d_i = (double)weighted_n_event[k];
          if (n_i > d_i) {
            if (error == "tsiatis")      sum_se += d_i / (n_i * n_i);
            else /* greenwood */         sum_se += d_i / (n_i * (n_i - d_i));
          } else {
            sum_se = std::numeric_limits<double>::infinity();
          }
        }
        std_err[j] = std::sqrt(sum_se);
      }

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
      combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
      combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
      combined_n_stratum.push_back(n_stratum);
      combined_u_stratum.push_back(u_stratum);
      strata_levels_out = levels;
    }
  } else {

    Rcpp::IntegerVector strata_vec = strata;
    Rcpp::IntegerVector levels = Rcpp::sort_unique(strata_vec);

    for (int li = 0; li < levels.size(); ++li) {
      int level = levels[li];
      Rcpp::LogicalVector strata_condition = (strata_vec == level);

      Rcpp::NumericVector t_selected = t[strata_condition];
      Rcpp::IntegerVector d_selected = d[strata_condition];
      Rcpp::NumericVector w_selected = w[strata_condition];

      int n_stratum = t_selected.size();

      Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
      std::sort(unique_times.begin(), unique_times.end());
      int u_stratum = unique_times.size();

      Rcpp::NumericVector km(u_stratum), km_i(u_stratum), std_err(u_stratum);
      Rcpp::NumericVector weighted_n_risk(u_stratum), weighted_n_event(u_stratum), weighted_n_censor(u_stratum);

      for (int j = 0; j < u_stratum; ++j) {
        double time = unique_times[j];
        double n_at_risk = 0.0, w_events = 0.0;

        for (int k = 0; k < t_selected.size(); ++k) {
          if (t_selected[k] >= time) n_at_risk += w_selected[k];
          if (t_selected[k] == time && d_selected[k] == 1) w_events += w_selected[k];
        }
        weighted_n_risk[j]  = n_at_risk;
        weighted_n_event[j] = w_events;

        km_i[j] = (n_at_risk > 0.0) ? (1.0 - w_events / n_at_risk) : 1.0;
        km[j]   = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

        if (j > 0) weighted_n_censor[j-1] = weighted_n_risk[j-1] - weighted_n_risk[j] - weighted_n_event[j-1];
        if (j == u_stratum - 1) weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];

        double sum_se = 0.0;
        for (int k = 0; k <= j; ++k) {
          double n_i = weighted_n_risk[k];
          double d_i = 0.0;
          for (int m = 0; m < t_selected.size(); ++m) {
            if (t_selected[m] == unique_times[k] && d_selected[m] == 1) d_i += w_selected[m];
          }
          if (n_i > d_i) {
            if (error == "tsiatis")      sum_se += d_i / (n_i * n_i);
            else /* greenwood */         sum_se += d_i / (n_i * (n_i - d_i));
          } else {
            sum_se = std::numeric_limits<double>::infinity();
          }
        }
        std_err[j] = std::sqrt(sum_se);
      }

      combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
      combined_surv.insert(combined_surv.end(), km.begin(), km.end());
      combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
      combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
      combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
      combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
      combined_n_stratum.push_back(n_stratum);
      combined_u_stratum.push_back(u_stratum);
      strata_levels_out = levels;
    }
  }

  Rcpp::NumericVector all_times      = Rcpp::wrap(combined_times);
  Rcpp::NumericVector all_surv       = Rcpp::wrap(combined_surv);
  Rcpp::IntegerVector all_n_risk     = Rcpp::wrap(combined_n_risk);
  Rcpp::IntegerVector all_n_event    = Rcpp::wrap(combined_n_event);
  Rcpp::IntegerVector all_n_censor   = Rcpp::wrap(combined_n_censor);
  Rcpp::NumericVector all_std_err    = Rcpp::wrap(combined_std_err);
  Rcpp::IntegerVector all_n_stratum  = Rcpp::wrap(combined_n_stratum);
  Rcpp::IntegerVector all_u_stratum  = Rcpp::wrap(combined_u_stratum);

  survfit_list = Rcpp::List::create(
    Rcpp::_["time"]            = all_times,
    Rcpp::_["surv"]            = all_surv,
    Rcpp::_["n.risk"]          = all_n_risk,
    Rcpp::_["n"]               = combined_n_stratum,
    Rcpp::_["n.event"]         = all_n_event,
    Rcpp::_["n.censor"]        = all_n_censor,
    Rcpp::_["unweighted.n"]    = all_n_stratum,
    Rcpp::_["std.err"]         = all_std_err,
    Rcpp::_["high"]            = R_NilValue,
    Rcpp::_["low"]             = R_NilValue,
    Rcpp::_["conf.type"]       = "log-log",
    Rcpp::_["strata"]          = all_u_stratum,
    Rcpp::_["strata.levels"]   = strata_levels_out,
    Rcpp::_["type"]            = "kaplan-meier",
    Rcpp::_["method"]          = "Kaplan-Meier"
  );
  return survfit_list;
}



// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {
constexpr double EPS = 1e-12;

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

   // KM(生存)の分散推定に Tsiatis を使うか（それ以外は Greenwood）
   const bool use_tsiatis_surv = (ERROR == "tsiatis");

   // CIF の SE 計算モード
   const bool need_cif_aalen = (ERROR == "aalen");
   const bool need_cif_delta = (ERROR == "delta");
   const bool need_cif_if    = (ERROR == "if") || return_if; // 返すだけでも IF は計算する
   const bool need_any_cif   = need_cif_aalen || need_cif_delta || need_cif_if;

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

   // strata -> block index
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
   std::vector<double> combined_std_err, combined_std_err_cif, combined_high, combined_low;
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
     cnt_censor_all(Uall, 0.0),
     cnt_total_all(Uall, 0.0);
     {
       size_t p2 = 0;
       for (int j = 0; j < Uall; ++j) {
         double curT = t_all[j];
         double e = 0.0, c = 0.0, tot = 0.0;
         while (p2 < v.size() && feq(v[p2].t, curT)) {
           if (v[p2].eps >= 1) e += 1.0; else c += 1.0;
           tot += 1.0;
           ++p2;
         }
         cnt_event_all[j]  = e;
         cnt_censor_all[j] = c;
         cnt_total_all[j]  = tot;
       }
     }

     std::vector<double> Y_all(Uall, 0.0);
     {
       double suf = 0.0;
       for (int j = Uall - 1; j >= 0; --j) { suf += w_total_all[j]; Y_all[j] = suf; }
     }
     std::vector<double> Y_unw_all(Uall, 0.0);
     {
       double suf_unw = 0.0;
       for (int j = Uall - 1; j >= 0; --j) { suf_unw += cnt_total_all[j]; Y_unw_all[j] = suf_unw; }
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

     // KM S and Greenwood/Tsiatis var on ANY grid
     // KM S and Greenwood/Tsiatis var on ANY grid (compute only the one you need)
     std::vector<double> S_any(M_any, 1.0), varKM_any(M_any, 0.0);
     double accKM = 0.0;
     for (int m=0; m<M_any; ++m){
       int j = ev_any_idx[m];
       double Yj = Y_all[j], dj = w_event_all[j];
       double Sprev = (m==0? 1.0 : S_any[m-1]);
       double fac   = (Yj>0.0) ? (1.0 - dj/Yj) : 1.0;
       S_any[m]     = Sprev * fac;

       if (dj>0.0 && Yj>0.0){
         if (use_tsiatis_surv) accKM += dj / (Yj*Yj);
         else {
           if (Yj > dj) accKM += dj / (Yj*(Yj-dj));
           else         accKM  = std::numeric_limits<double>::infinity();
         }
       }
       double S2 = S_any[m] * S_any[m];
       varKM_any[m] = S2 * accKM;
     }

     // carry S/SE to ALL times (carry forward)
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

     // SE(log S) on ALL times (unweighted Greenwood/Tsiatis accumulator)
     std::vector<double> SE_log_surv(Uall, 0.0), SE_surv(Uall, 0.0);
     {
       double acc = 0.0;
       for (int j = 0; j < Uall; ++j) {
         double n_i = Y_unw_all[j];
         double d_i = cnt_event_all[j];
         if (d_i > 0.0) {
           if (n_i > d_i) {
             if (error == "tsiatis") acc += (d_i / (n_i * n_i));
             else                    acc += (d_i / (n_i * (n_i - d_i)));
           } else {
             acc = std::numeric_limits<double>::infinity();
           }
         }
         SE_log_surv[j] = (R_finite(acc) ? std::sqrt(acc) : R_PosInf);
       }
       for (int j = 0; j < Uall; ++j) SE_surv[j] = SE_all[j];
     }

      std::vector<double> se_cif_if(Uall, 0.0);
      NumericMatrix IF_AJ_all;
      if (need_cif_if) {
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
      IF_AJ_output[kg] = (need_cif_if && return_if) ? IF_AJ_all : Rcpp::NumericMatrix(0,0);


     // Aalen/Delta variance (cause-1 grid), then carry to ALL
     std::vector<double> se_cif_aalen(Uall, 0.0), se_cif_delta(Uall, 0.0);

     if (need_cif_aalen || need_cif_delta) {
       // cause-1 系の下ごしらえ
       std::vector<double> km_time(M_any);
       for (int m=0;m<M_any;++m) km_time[m] = t_all[ ev_any_idx[m] ];

       std::vector<double> var_aalen, var_delta;
       if (need_cif_aalen) var_aalen.assign(M_c1, 0.0);
       if (need_cif_delta) var_delta.assign(M_c1, 0.0);

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

       if (need_cif_aalen) {
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

       if (need_cif_delta) {
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
           if (need_cif_aalen && M_c1 > 0) se_cif_aalen[j] = std::sqrt(std::max(0.0, var_aalen[0]));
           if (need_cif_delta && M_c1 > 0) se_cif_delta[j] = std::sqrt(std::max(0.0, var_delta[0]));
         } else {
           int use = std::min(c, M_c1);  // 末尾越え防止
           if (need_cif_aalen) se_cif_aalen[j] = std::sqrt(std::max(0.0, var_aalen[use - 1]));
           if (need_cif_delta) se_cif_delta[j] = std::sqrt(std::max(0.0, var_delta[use - 1]));
         }
       }
     }

     // --- shift se_cif_* one step forward on the ALL grid (t -> t+1) ---
     if (Uall >= 2) {
       if (need_cif_aalen) {
         std::vector<double> tmp = se_cif_aalen;          // copy
         for (int j = 0; j < Uall - 1; ++j) se_cif_aalen[j] = tmp[j + 1];
         se_cif_aalen[Uall - 1] = tmp[Uall - 1];          // keep last (no t_max+1)
       }
       if (need_cif_delta) {
         std::vector<double> tmp = se_cif_delta;          // copy
         for (int j = 0; j < Uall - 1; ++j) se_cif_delta[j] = tmp[j + 1];
         se_cif_delta[Uall - 1] = tmp[Uall - 1];          // keep last
       }
     }

     // choose SE for CI on S(t) / CIF
     std::vector<double> SE_cif(Uall, 0.0);
     if      (ERROR=="aalen") SE_cif = se_cif_aalen;
     else if (ERROR=="delta") SE_cif = se_cif_delta;
     else if (ERROR=="if")    SE_cif = se_cif_if;
     else                     SE_cif = SE_all;  // ← ここは SE_surv ではなく SE_all が中身

     // CI bands (per CONF)
     std::vector<double> high_all(Uall, NA_REAL), low_all(Uall, NA_REAL);
     if (CONF != "none"){
       for (int j=0;j<Uall;++j){
         double Sj  = clamp01(S_all[j], 1e-15);
         double SEj = SE_cif[j];
         if (!R_finite(SEj) || Sj<=0.0 || Sj>=1.0){ high_all[j]=NA_REAL; low_all[j]=NA_REAL; continue; }

         if (CONF=="plain"){
           double lo   = std::max(0.0, Sj - z*SEj);
           double hi   = std::min(1.0, Sj + z*SEj);
           low_all[j]  = lo; high_all[j]=hi;
         } else if (CONF=="log"){
           double se   = SEj / Sj;
           double lo   = Sj * std::exp(-z*se);
           double hi   = Sj * std::exp( z*se);
           low_all[j]  = std::max(0.0, lo);
           high_all[j] = std::min(1.0, hi);
         } else if (CONF=="log-log"){
           double se   = SEj / ( Sj * std::fabs(std::log(Sj)) );
           double lo   = std::exp( - std::exp( std::log(-std::log(Sj)) + z*se ) );
           double hi   = std::exp( - std::exp( std::log(-std::log(Sj)) - z*se ) );
           low_all[j]  = std::max(0.0, std::min(1.0, lo));
           high_all[j] = std::max(0.0, std::min(1.0, hi));
         } else if (CONF=="logit"){
           double se   = SEj / ( Sj * (1.0 - Sj) );
           double lo   = Sj / ( Sj + (1.0 - Sj)*std::exp( z*se) );
           double hi   = Sj / ( Sj + (1.0 - Sj)*std::exp(-z*se) );
           low_all[j]  = std::max(0.0, std::min(1.0, lo));
           high_all[j] = std::max(0.0, std::min(1.0, hi));
         } else if (CONF=="arcsine-square root"){
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

     // append combined outputs
     combined_times.reserve(combined_times.size() + Uall);
     combined_surv .reserve(combined_surv .size() + Uall);
     combined_n_risk.reserve(combined_n_risk.size() + Uall);
     combined_n_event.reserve(combined_n_event.size() + Uall);
     combined_n_censor.reserve(combined_n_censor.size() + Uall);
     combined_std_err.reserve(combined_std_err.size() + Uall);
     combined_std_err_cif.reserve(combined_std_err_cif.size() + Uall);
     combined_low.reserve(combined_low.size() + Uall);
     combined_high.reserve(combined_high.size() + Uall);

     for (int j=0;j<Uall;++j){
       combined_times.push_back(t_all[j]);
       combined_surv.push_back(S_all[j]);
       combined_n_risk.push_back(Y_all[j]);
       combined_n_event.push_back( cnt_event_all[j] );
       combined_n_censor.push_back(cnt_censor_all[j]);
       combined_std_err.push_back(SE_log_surv[j]);
       combined_std_err_cif.push_back(SE_cif[j]);
       combined_low.push_back(low_all[j]);
       combined_high.push_back(high_all[j]);
     }
     combined_n_stratum.push_back(n_g);
     combined_u_stratum.push_back(Uall);
   }

   // output meta
   std::string out_type   = "kaplan-meier";
   std::string out_method = "Kaplan-Meier";
   if (error=="aalen" || error=="delta" || error=="if"){
     out_type = "aalen-johansen";
     out_method = "Aalen-Johansen";
   }

   SEXP strata_out, strata_levels_out_sexp;
   if (K==1){
     strata_out = Rcpp::IntegerVector(0);
     strata_levels_out_sexp = Rcpp::IntegerVector(0);
   } else {
     strata_out = wrap(combined_u_stratum);
     strata_levels_out_sexp = wrap(strata_levels_out);
   }

   Rcpp::List out = Rcpp::List::create(
     _["time"]               = wrap(combined_times),
     _["surv"]               = wrap(combined_surv),
     _["n.risk"]             = wrap(combined_n_risk),
     _["n"]                  = wrap(combined_n_stratum),
     _["n.event"]            = wrap(combined_n_event),
     _["n.censor"]           = wrap(combined_n_censor),
     _["unweighted.n"]       = wrap(combined_n_stratum),
     _["std.err"]            = wrap(combined_std_err),
     _["std.err.cif"]        = wrap(combined_std_err_cif),
     _["low"]                = wrap(combined_low),
     _["high"]               = wrap(combined_high),
     _["conf.type"]          = CONF,
     _["strata"]             = strata_out,
     _["strata.levels"]      = strata_levels_out_sexp,
     _["type"]               = out_type,
     _["method"]             = out_method,
     _["influence.function"] = IF_AJ_output
   );
   return out;
 }



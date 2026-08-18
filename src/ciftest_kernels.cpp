#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List ciftest_fg_iid_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector strata_id,
    int code_event1,
    int code_event2,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericVector fh_weight,
    Rcpp::NumericVector g_at_competing,
    Rcpp::NumericMatrix g_event_stratum,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    double prob_bound
) {
  const int n = t.size();
  const int p = x.ncol();
  const int m = event_times.size();
  const int h_count = hazard_time.size();
  if (epsilon.size() != n || x.nrow() != n || weights.size() != n ||
      strata_id.size() != n || g_at_competing.size() != n ||
      fh_weight.size() != m || g_event_stratum.nrow() != m ||
      hazard_stratum.size() != h_count || hazard.size() != h_count ||
      hazard_n_risk.size() != h_count) {
    stop("Fine-Gray Rcpp kernel received incompatible dimensions.");
  }

  NumericMatrix base_iid(n, p);
  NumericMatrix censor_iid(n, p);
  NumericMatrix censor_derivative(h_count, p);
  NumericMatrix xbar(m, p);
  std::fill(xbar.begin(), xbar.end(), NA_REAL);
  NumericVector event_score(p);
  NumericVector risk_total(m);
  NumericVector event_count(m);
  NumericVector censor_weight(n);
  NumericVector mean_x(p);

  for (int j = 0; j < m; ++j) {
    const double current = event_times[j];
    double s0 = 0.0;
    double d1 = 0.0;
    bool any_competing_before = false;
    std::fill(mean_x.begin(), mean_x.end(), 0.0);

    for (int i = 0; i < n; ++i) {
      const bool competing_before =
        epsilon[i] == code_event2 && t[i] < current;
      const bool subrisk = t[i] >= current || competing_before;
      double cw = 1.0;
      if (competing_before) {
        any_competing_before = true;
        const int stratum = strata_id[i] - 1;
        if (stratum < 0 || stratum >= g_event_stratum.ncol()) {
          stop("Fine-Gray Rcpp kernel received an invalid stratum code.");
        }
        const double denominator = g_at_competing[i];
        const double numerator = g_event_stratum(j, stratum);
        if (denominator <= prob_bound || numerator <= prob_bound) {
          stop("Censoring positivity is violated in the Fine-Gray risk set.");
        }
        cw = numerator / denominator;
      }
      censor_weight[i] = cw;
      const double weighted_risk = weights[i] * (subrisk ? cw : 0.0);
      s0 += weighted_risk;
      for (int q = 0; q < p; ++q) {
        mean_x[q] += x(i, q) * weighted_risk;
      }
      if (epsilon[i] == code_event1 && t[i] == current) {
        d1 += weights[i];
      }
    }
    if (!R_finite(s0) || s0 <= prob_bound) {
      stop("Fine-Gray subdistribution risk set is empty.");
    }
    for (int q = 0; q < p; ++q) {
      mean_x[q] /= s0;
      xbar(j, q) = mean_x[q];
    }
    if (d1 <= 0.0) continue;

    const double hazard1 = d1 / s0;
    const double a = fh_weight[j];
    for (int i = 0; i < n; ++i) {
      const bool is_event = epsilon[i] == code_event1 && t[i] == current;
      const bool competing_before =
        epsilon[i] == code_event2 && t[i] < current;
      const bool subrisk = t[i] >= current || competing_before;
      const double coefficient = weights[i] *
        ((is_event ? 1.0 : 0.0) -
         (subrisk ? censor_weight[i] * hazard1 : 0.0));
      for (int q = 0; q < p; ++q) {
        base_iid(i, q) += (x(i, q) - mean_x[q]) * a * coefficient;
        if (is_event) event_score[q] += a * x(i, q) * weights[i];
      }
    }
    for (int q = 0; q < p; ++q) {
      event_score[q] -= a * d1 * mean_x[q];
    }
    risk_total[j] = s0;
    event_count[j] = d1;

    if (h_count > 0 && any_competing_before) {
      for (int h = 0; h < h_count; ++h) {
        if (!(hazard_time[h] < current)) continue;
        const double one_minus_hazard = 1.0 - hazard[h];
        bool any_affected = false;
        for (int i = 0; i < n; ++i) {
          if (epsilon[i] == code_event2 && t[i] < current &&
              strata_id[i] == hazard_stratum[h] &&
              t[i] <= hazard_time[h]) {
            any_affected = true;
            break;
          }
        }
        if (!any_affected) continue;
        if (one_minus_hazard <= prob_bound) {
          stop("Censoring positivity is violated after a competing event.");
        }
        for (int q = 0; q < p; ++q) {
          double derivative_sum = 0.0;
          for (int i = 0; i < n; ++i) {
            if (epsilon[i] == code_event2 && t[i] < current &&
                strata_id[i] == hazard_stratum[h] &&
                t[i] <= hazard_time[h]) {
              const double derivative_weight =
                -censor_weight[i] / one_minus_hazard;
              derivative_sum += (x(i, q) - mean_x[q]) *
                weights[i] * derivative_weight;
            }
          }
          censor_derivative(h, q) -=
            a * d1 * derivative_sum / s0;
        }
      }
    }
  }

  if (h_count > 0) {
    for (int h = 0; h < h_count; ++h) {
      if (hazard_n_risk[h] <= prob_bound) {
        stop("Censoring hazard risk set is empty in the Fine-Gray kernel.");
      }
      for (int i = 0; i < n; ++i) {
        const bool in_stratum = strata_id[i] == hazard_stratum[h];
        const bool at_risk = in_stratum && t[i] >= hazard_time[h];
        const bool censored = in_stratum && t[i] == hazard_time[h] &&
          epsilon[i] == code_censoring;
        const double martingale = weights[i] *
          ((censored ? 1.0 : 0.0) - hazard[h] * (at_risk ? 1.0 : 0.0));
        const double multiplier = martingale / hazard_n_risk[h];
        for (int q = 0; q < p; ++q) {
          censor_iid(i, q) += multiplier * censor_derivative(h, q);
        }
      }
    }
  }

  NumericVector score(p);
  for (int q = 0; q < p; ++q) {
    for (int i = 0; i < n; ++i) score[q] += base_iid(i, q);
  }
  return List::create(
    _["score"] = score,
    _["event.score"] = event_score,
    _["score.iid.base"] = base_iid,
    _["score.iid.censor"] = censor_iid,
    _["xbar"] = xbar,
    _["risk.total"] = risk_total,
    _["event.count"] = event_count,
    _["censor.derivative"] = censor_derivative
  );
}

// Batch-only variant that evaluates several fixed FH weight processes while
// sharing the Fine-Gray risk-set and censoring-hazard traversals.
// [[Rcpp::export]]
Rcpp::List ciftest_fg_iid_multi_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector strata_id,
    int code_event1,
    int code_event2,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericMatrix fh_weight,
    Rcpp::NumericVector g_at_competing,
    Rcpp::NumericMatrix g_event_stratum,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    double prob_bound
) {
  const int n = t.size();
  const int p = x.ncol();
  const int m = event_times.size();
  const int weight_count = fh_weight.ncol();
  const int h_count = hazard_time.size();
  if (epsilon.size() != n || x.nrow() != n || weights.size() != n ||
      strata_id.size() != n || g_at_competing.size() != n ||
      fh_weight.nrow() != m || weight_count < 1 ||
      g_event_stratum.nrow() != m ||
      hazard_stratum.size() != h_count || hazard.size() != h_count ||
      hazard_n_risk.size() != h_count) {
    stop("Fine-Gray multiweight kernel received incompatible dimensions.");
  }

  NumericVector base_iid(n * p * weight_count);
  NumericVector censor_iid(n * p * weight_count);
  NumericVector censor_derivative(h_count * p * weight_count);
  base_iid.attr("dim") = IntegerVector::create(n, p, weight_count);
  censor_iid.attr("dim") = IntegerVector::create(n, p, weight_count);
  censor_derivative.attr("dim") =
    IntegerVector::create(h_count, p, weight_count);
  NumericMatrix event_score(p, weight_count);
  NumericMatrix score(p, weight_count);
  NumericMatrix xbar(m, p);
  std::fill(xbar.begin(), xbar.end(), NA_REAL);
  NumericVector risk_total(m);
  NumericVector event_count(m);
  NumericVector censor_weight(n);
  NumericVector mean_x(p);

  auto score_index = [=](int i, int q, int w) {
    return i + n * (q + p * w);
  };
  auto derivative_index = [=](int h, int q, int w) {
    return h + h_count * (q + p * w);
  };

  for (int j = 0; j < m; ++j) {
    const double current = event_times[j];
    double s0 = 0.0;
    double d1 = 0.0;
    bool any_competing_before = false;
    std::fill(mean_x.begin(), mean_x.end(), 0.0);

    for (int i = 0; i < n; ++i) {
      const bool competing_before =
        epsilon[i] == code_event2 && t[i] < current;
      const bool subrisk = t[i] >= current || competing_before;
      double cw = 1.0;
      if (competing_before) {
        any_competing_before = true;
        const int stratum = strata_id[i] - 1;
        if (stratum < 0 || stratum >= g_event_stratum.ncol()) {
          stop("Fine-Gray multiweight kernel received an invalid stratum code.");
        }
        const double denominator = g_at_competing[i];
        const double numerator = g_event_stratum(j, stratum);
        if (denominator <= prob_bound || numerator <= prob_bound) {
          stop("Censoring positivity is violated in the Fine-Gray risk set.");
        }
        cw = numerator / denominator;
      }
      censor_weight[i] = cw;
      const double weighted_risk = weights[i] * (subrisk ? cw : 0.0);
      s0 += weighted_risk;
      for (int q = 0; q < p; ++q) {
        mean_x[q] += x(i, q) * weighted_risk;
      }
      if (epsilon[i] == code_event1 && t[i] == current) d1 += weights[i];
    }
    if (!R_finite(s0) || s0 <= prob_bound) {
      stop("Fine-Gray subdistribution risk set is empty.");
    }
    for (int q = 0; q < p; ++q) {
      mean_x[q] /= s0;
      xbar(j, q) = mean_x[q];
    }
    if (d1 <= 0.0) continue;

    const double hazard1 = d1 / s0;
    for (int i = 0; i < n; ++i) {
      const bool is_event = epsilon[i] == code_event1 && t[i] == current;
      const bool competing_before =
        epsilon[i] == code_event2 && t[i] < current;
      const bool subrisk = t[i] >= current || competing_before;
      const double coefficient = weights[i] *
        ((is_event ? 1.0 : 0.0) -
         (subrisk ? censor_weight[i] * hazard1 : 0.0));
      for (int w = 0; w < weight_count; ++w) {
        const double a = fh_weight(j, w);
        for (int q = 0; q < p; ++q) {
          base_iid[score_index(i, q, w)] +=
            (x(i, q) - mean_x[q]) * a * coefficient;
          if (is_event) event_score(q, w) += a * x(i, q) * weights[i];
        }
      }
    }
    for (int w = 0; w < weight_count; ++w) {
      const double a = fh_weight(j, w);
      for (int q = 0; q < p; ++q) {
        event_score(q, w) -= a * d1 * mean_x[q];
      }
    }
    risk_total[j] = s0;
    event_count[j] = d1;

    if (h_count > 0 && any_competing_before) {
      for (int h = 0; h < h_count; ++h) {
        if (!(hazard_time[h] < current)) continue;
        const double one_minus_hazard = 1.0 - hazard[h];
        bool any_affected = false;
        for (int i = 0; i < n; ++i) {
          if (epsilon[i] == code_event2 && t[i] < current &&
              strata_id[i] == hazard_stratum[h] &&
              t[i] <= hazard_time[h]) {
            any_affected = true;
            break;
          }
        }
        if (!any_affected) continue;
        if (one_minus_hazard <= prob_bound) {
          stop("Censoring positivity is violated after a competing event.");
        }
        for (int q = 0; q < p; ++q) {
          double derivative_sum = 0.0;
          for (int i = 0; i < n; ++i) {
            if (epsilon[i] == code_event2 && t[i] < current &&
                strata_id[i] == hazard_stratum[h] &&
                t[i] <= hazard_time[h]) {
              derivative_sum += (x(i, q) - mean_x[q]) * weights[i] *
                (-censor_weight[i] / one_minus_hazard);
            }
          }
          for (int w = 0; w < weight_count; ++w) {
            censor_derivative[derivative_index(h, q, w)] -=
              fh_weight(j, w) * d1 * derivative_sum / s0;
          }
        }
      }
    }
  }

  for (int h = 0; h < h_count; ++h) {
    if (hazard_n_risk[h] <= prob_bound) {
      stop("Censoring hazard risk set is empty in the Fine-Gray kernel.");
    }
    for (int i = 0; i < n; ++i) {
      const bool in_stratum = strata_id[i] == hazard_stratum[h];
      const bool at_risk = in_stratum && t[i] >= hazard_time[h];
      const bool censored = in_stratum && t[i] == hazard_time[h] &&
        epsilon[i] == code_censoring;
      const double martingale = weights[i] *
        ((censored ? 1.0 : 0.0) - hazard[h] * (at_risk ? 1.0 : 0.0));
      const double multiplier = martingale / hazard_n_risk[h];
      for (int w = 0; w < weight_count; ++w) {
        for (int q = 0; q < p; ++q) {
          censor_iid[score_index(i, q, w)] += multiplier *
            censor_derivative[derivative_index(h, q, w)];
        }
      }
    }
  }

  for (int w = 0; w < weight_count; ++w) {
    for (int q = 0; q < p; ++q) {
      for (int i = 0; i < n; ++i) {
        score(q, w) += base_iid[score_index(i, q, w)];
      }
    }
  }
  return List::create(
    _["score"] = score,
    _["event.score"] = event_score,
    _["score.iid.base"] = base_iid,
    _["score.iid.censor"] = censor_iid,
    _["xbar"] = xbar,
    _["risk.total"] = risk_total,
    _["event.count"] = event_count,
    _["censor.derivative"] = censor_derivative
  );
}

// Prefix/suffix implementation of the Fine-Gray risk moments, base IID,
// censoring derivative, and censoring IID.  It remains a separate entry point
// so that the optimization can be enabled and rolled back independently. For
// censoring stratum s and censoring time c, the derivative evaluates
//
//   C1_s(c) R0_s(c) - C0_s(c) R1_s(c)
//
// from a competing-event prefix and an event-time suffix, avoiding the legacy
// event x censoring-time x subject traversal.
static Rcpp::List ciftest_fg_iid_prefix_core(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector strata_id,
    int code_event1,
    int code_event2,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericMatrix fh_weight,
    Rcpp::NumericVector g_at_competing,
    Rcpp::NumericMatrix g_event_stratum,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    double prob_bound
) {
  const int n = t.size();
  const int p = x.ncol();
  const int m = event_times.size();
  const int weight_count = fh_weight.ncol();
  const int h_count = hazard_time.size();
  const int stratum_count = g_event_stratum.ncol();
  if (epsilon.size() != n || x.nrow() != n || weights.size() != n ||
      strata_id.size() != n || g_at_competing.size() != n ||
      fh_weight.nrow() != m || weight_count < 1 ||
      g_event_stratum.nrow() != m || stratum_count < 1 ||
      hazard_stratum.size() != h_count || hazard.size() != h_count ||
      hazard_n_risk.size() != h_count) {
    stop("Fine-Gray prefix kernel received incompatible dimensions.");
  }

  NumericVector base_iid(n * p * weight_count);
  NumericVector censor_iid(n * p * weight_count);
  NumericVector censor_derivative(h_count * p * weight_count);
  base_iid.attr("dim") = IntegerVector::create(n, p, weight_count);
  censor_iid.attr("dim") = IntegerVector::create(n, p, weight_count);
  censor_derivative.attr("dim") =
    IntegerVector::create(h_count, p, weight_count);
  NumericMatrix event_score(p, weight_count);
  NumericMatrix score(p, weight_count);
  NumericMatrix xbar(m, p);
  NumericVector risk_total(m);
  NumericVector event_count(m);
  NumericMatrix event_x(m, p);

  auto score_index = [=](int i, int q, int w) {
    return i + n * (q + p * w);
  };
  auto derivative_index = [=](int h, int q, int w) {
    return h + h_count * (q + p * w);
  };

  for (int i = 0; i < n; ++i) {
    if (strata_id[i] < 1 || strata_id[i] > stratum_count) {
      stop("Fine-Gray prefix kernel received an invalid stratum code.");
    }
  }
  for (int h = 0; h < h_count; ++h) {
    if (hazard_stratum[h] < 1 || hazard_stratum[h] > stratum_count) {
      stop("Fine-Gray prefix kernel received an invalid hazard stratum code.");
    }
  }
  std::vector<std::vector<int> > competing_by_stratum(stratum_count);
  std::vector<std::vector<int> > hazard_by_stratum(stratum_count);
  for (int i = 0; i < n; ++i) {
    if (epsilon[i] == code_event2) {
      competing_by_stratum[strata_id[i] - 1].push_back(i);
    }
  }
  for (int h = 0; h < h_count; ++h) {
    hazard_by_stratum[hazard_stratum[h] - 1].push_back(h);
  }
  for (int s = 0; s < stratum_count; ++s) {
    std::sort(
      competing_by_stratum[s].begin(), competing_by_stratum[s].end(),
      [&](int left, int right) {
        if (t[left] != t[right]) return t[left] < t[right];
        return left < right;
      }
    );
    std::sort(
      hazard_by_stratum[s].begin(), hazard_by_stratum[s].end(),
      [&](int left, int right) {
        if (hazard_time[left] != hazard_time[right]) {
          return hazard_time[left] < hazard_time[right];
        }
        return left < right;
      }
    );
  }

  // Fine-Gray risk moments are an ordinary at-risk suffix plus a
  // stratum-specific competing-event prefix. This replaces the event x
  // subject risk traversal.
  std::vector<int> subject_order(n);
  std::iota(subject_order.begin(), subject_order.end(), 0);
  std::sort(
    subject_order.begin(), subject_order.end(),
    [&](int left, int right) {
      if (t[left] != t[right]) return t[left] < t[right];
      return left < right;
    }
  );
  int subject_pos = n - 1;
  double ordinary0 = 0.0;
  std::vector<double> ordinary1(p, 0.0);
  for (int j = m - 1; j >= 0; --j) {
    while (subject_pos >= 0 &&
           t[subject_order[subject_pos]] >= event_times[j]) {
      const int i = subject_order[subject_pos];
      ordinary0 += weights[i];
      for (int q = 0; q < p; ++q) {
        ordinary1[q] += weights[i] * x(i, q);
      }
      --subject_pos;
    }
    risk_total[j] = ordinary0;
    for (int q = 0; q < p; ++q) xbar(j, q) = ordinary1[q];
  }
  for (int s = 0; s < stratum_count; ++s) {
    const std::vector<int>& comp = competing_by_stratum[s];
    int comp_pos = 0;
    int comp_count = 0;
    double competing0 = 0.0;
    std::vector<double> competing1(p, 0.0);
    for (int j = 0; j < m; ++j) {
      while (comp_pos < static_cast<int>(comp.size()) &&
             t[comp[comp_pos]] < event_times[j]) {
        const int i = comp[comp_pos];
        const double denominator = g_at_competing[i];
        if (denominator <= prob_bound) {
          stop("Censoring positivity is violated in the Fine-Gray risk set.");
        }
        const double value = weights[i] / denominator;
        competing0 += value;
        for (int q = 0; q < p; ++q) {
          competing1[q] += value * x(i, q);
        }
        ++comp_pos;
        ++comp_count;
      }
      if (comp_count > 0) {
        const double g = g_event_stratum(j, s);
        if (g <= prob_bound) {
          stop("Censoring positivity is violated in the Fine-Gray risk set.");
        }
        risk_total[j] += g * competing0;
        for (int q = 0; q < p; ++q) {
          xbar(j, q) += g * competing1[q];
        }
      }
    }
  }

  // Event totals retain original subject order, while their risk-set moments
  // come from the sorted prefix/suffix calculation above.
  for (int i = 0; i < n; ++i) {
    if (epsilon[i] != code_event1 || weights[i] <= 0.0) continue;
    const int j = std::lower_bound(
      event_times.begin(), event_times.end(), t[i]
    ) - event_times.begin();
    if (j < m && event_times[j] == t[i]) {
      event_count[j] += weights[i];
      for (int q = 0; q < p; ++q) {
        event_x(j, q) += weights[i] * x(i, q);
      }
    }
  }
  for (int j = 0; j < m; ++j) {
    if (!R_finite(risk_total[j]) || risk_total[j] <= prob_bound) {
      stop("Fine-Gray subdistribution risk set is empty.");
    }
    for (int q = 0; q < p; ++q) xbar(j, q) /= risk_total[j];
    for (int w = 0; w < weight_count; ++w) {
      const double a = fh_weight(j, w);
      for (int q = 0; q < p; ++q) {
        event_score(q, w) += a *
          (event_x(j, q) - event_count[j] * xbar(j, q));
      }
    }
  }

  // Known-censoring IID terms use an event-time prefix for ordinary risk and
  // a G-weighted reverse suffix after a competing event.
  std::vector<double> prefix0(m * weight_count, 0.0);
  std::vector<double> prefix1(m * p * weight_count, 0.0);
  std::vector<double> suffix0((m + 1) * weight_count * stratum_count, 0.0);
  std::vector<double> suffix1(
    (m + 1) * p * weight_count * stratum_count, 0.0
  );
  auto prefix0_index = [=](int j, int w) { return j + m * w; };
  auto prefix1_index = [=](int j, int q, int w) {
    return j + m * (q + p * w);
  };
  auto suffix0_index = [=](int j, int w, int s) {
    return j + (m + 1) * (w + weight_count * s);
  };
  auto suffix1_index = [=](int j, int q, int w, int s) {
    return j + (m + 1) * (q + p * (w + weight_count * s));
  };
  for (int w = 0; w < weight_count; ++w) {
    double running0 = 0.0;
    std::vector<double> running1(p, 0.0);
    for (int j = 0; j < m; ++j) {
      const double value = fh_weight(j, w) * event_count[j] /
        risk_total[j];
      running0 += value;
      prefix0[prefix0_index(j, w)] = running0;
      for (int q = 0; q < p; ++q) {
        running1[q] += value * xbar(j, q);
        prefix1[prefix1_index(j, q, w)] = running1[q];
      }
    }
    for (int s = 0; s < stratum_count; ++s) {
      running0 = 0.0;
      std::fill(running1.begin(), running1.end(), 0.0);
      for (int j = m - 1; j >= 0; --j) {
        const double value = fh_weight(j, w) * event_count[j] /
          risk_total[j] * g_event_stratum(j, s);
        running0 += value;
        suffix0[suffix0_index(j, w, s)] = running0;
        for (int q = 0; q < p; ++q) {
          running1[q] += value * xbar(j, q);
          suffix1[suffix1_index(j, q, w, s)] = running1[q];
        }
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    const int last = std::upper_bound(
      event_times.begin(), event_times.end(), t[i]
    ) - event_times.begin() - 1;
    const int first_after = last + 1;
    int event_index = -1;
    if (epsilon[i] == code_event1) {
      event_index = std::lower_bound(
        event_times.begin(), event_times.end(), t[i]
      ) - event_times.begin();
      if (event_index >= m || event_times[event_index] != t[i]) {
        event_index = -1;
      }
    }
    for (int w = 0; w < weight_count; ++w) {
      for (int q = 0; q < p; ++q) {
        double value = 0.0;
        if (last >= 0) {
          value -= weights[i] *
            (x(i, q) * prefix0[prefix0_index(last, w)] -
             prefix1[prefix1_index(last, q, w)]);
        }
        if (event_index >= 0) {
          value += weights[i] * fh_weight(event_index, w) *
            (x(i, q) - xbar(event_index, q));
        }
        if (epsilon[i] == code_event2 && first_after < m) {
          const double denominator = g_at_competing[i];
          if (denominator <= prob_bound) {
            stop("Censoring positivity is violated in the Fine-Gray risk set.");
          }
          const int s = strata_id[i] - 1;
          value -= weights[i] / denominator *
            (x(i, q) * suffix0[suffix0_index(first_after, w, s)] -
             suffix1[suffix1_index(first_after, q, w, s)]);
        }
        base_iid[score_index(i, q, w)] = value;
      }
    }
  }

  // C0/C1 are competing-event prefixes at each censoring time. R0/R1 are
  // event-time suffixes strictly after that censoring time. Ties therefore
  // retain the same <= and < conventions as the legacy implementation.
  for (int s = 0; s < stratum_count; ++s) {
    const std::vector<int>& comp = competing_by_stratum[s];
    const std::vector<int>& haz = hazard_by_stratum[s];
    const int hs = haz.size();
    if (hs == 0) continue;
    std::vector<double> c0(hs, 0.0);
    std::vector<double> c1(hs * p, 0.0);
    std::vector<int> affected(hs, 0);
    double running_c0 = 0.0;
    std::vector<double> running_c1(p, 0.0);
    int comp_pos = 0;
    const double last_event = m > 0 ? event_times[m - 1] :
      -std::numeric_limits<double>::infinity();
    for (int local_h = 0; local_h < hs; ++local_h) {
      const int h = haz[local_h];
      if (!(hazard_time[h] < last_event)) break;
      while (comp_pos < static_cast<int>(comp.size()) &&
             t[comp[comp_pos]] <= hazard_time[h]) {
        const int i = comp[comp_pos];
        const double denominator = g_at_competing[i];
        if (denominator <= prob_bound) {
          stop("Censoring positivity is violated in the Fine-Gray risk set.");
        }
        const double value = weights[i] / denominator;
        running_c0 += value;
        for (int q = 0; q < p; ++q) {
          running_c1[q] += value * x(i, q);
        }
        ++comp_pos;
      }
      c0[local_h] = running_c0;
      affected[local_h] = comp_pos;
      for (int q = 0; q < p; ++q) {
        c1[local_h * p + q] = running_c1[q];
      }
    }

    std::vector<double> r0(weight_count, 0.0);
    std::vector<double> r1(p * weight_count, 0.0);
    int event_pos = m - 1;
    int future_event_count = 0;
    for (int local_h = hs - 1; local_h >= 0; --local_h) {
      const int h = haz[local_h];
      while (event_pos >= 0 && event_times[event_pos] > hazard_time[h]) {
        if (event_count[event_pos] > 0.0) {
          const double common = event_count[event_pos] /
            risk_total[event_pos] * g_event_stratum(event_pos, s);
          ++future_event_count;
          for (int w = 0; w < weight_count; ++w) {
            const double value = fh_weight(event_pos, w) * common;
            r0[w] += value;
            for (int q = 0; q < p; ++q) {
              r1[q + p * w] += value * xbar(event_pos, q);
            }
          }
        }
        --event_pos;
      }
      if (affected[local_h] == 0 || future_event_count == 0) continue;
      const double one_minus_hazard = 1.0 - hazard[h];
      if (one_minus_hazard <= prob_bound) {
        stop("Censoring positivity is violated after a competing event.");
      }
      for (int w = 0; w < weight_count; ++w) {
        for (int q = 0; q < p; ++q) {
          censor_derivative[derivative_index(h, q, w)] =
            (c1[local_h * p + q] * r0[w] -
             c0[local_h] * r1[q + p * w]) / one_minus_hazard;
        }
      }
    }
  }

  // The censoring martingale consists of a hazard prefix through each
  // subject's follow-up time and, for observed censoring, its jump term.
  std::vector<double> hazard_prefix(h_count * p * weight_count, 0.0);
  for (int h = 0; h < h_count; ++h) {
    if (hazard_n_risk[h] <= prob_bound) {
      stop("Censoring hazard risk set is empty in the Fine-Gray kernel.");
    }
  }
  for (int s = 0; s < stratum_count; ++s) {
    const std::vector<int>& haz = hazard_by_stratum[s];
    std::vector<double> running(p * weight_count, 0.0);
    for (int local_h = 0; local_h < static_cast<int>(haz.size()); ++local_h) {
      const int h = haz[local_h];
      for (int w = 0; w < weight_count; ++w) {
        for (int q = 0; q < p; ++q) {
          running[q + p * w] += hazard[h] / hazard_n_risk[h] *
            censor_derivative[derivative_index(h, q, w)];
          hazard_prefix[derivative_index(h, q, w)] = running[q + p * w];
        }
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    const std::vector<int>& haz = hazard_by_stratum[strata_id[i] - 1];
    const int through = std::upper_bound(
      haz.begin(), haz.end(), t[i],
      [&](double value, int h) { return value < hazard_time[h]; }
    ) - haz.begin();
    for (int w = 0; w < weight_count; ++w) {
      for (int q = 0; q < p; ++q) {
        double value = 0.0;
        if (through > 0) {
          value -= weights[i] *
            hazard_prefix[derivative_index(haz[through - 1], q, w)];
        }
        if (epsilon[i] == code_censoring) {
          int equal = through - 1;
          while (equal >= 0 && hazard_time[haz[equal]] == t[i]) {
            const int h = haz[equal];
            value += weights[i] / hazard_n_risk[h] *
              censor_derivative[derivative_index(h, q, w)];
            --equal;
          }
        }
        censor_iid[score_index(i, q, w)] = value;
      }
    }
  }

  for (int w = 0; w < weight_count; ++w) {
    for (int q = 0; q < p; ++q) {
      for (int i = 0; i < n; ++i) {
        score(q, w) += base_iid[score_index(i, q, w)];
      }
    }
  }
  return List::create(
    _["score"] = score,
    _["event.score"] = event_score,
    _["score.iid.base"] = base_iid,
    _["score.iid.censor"] = censor_iid,
    _["xbar"] = xbar,
    _["risk.total"] = risk_total,
    _["event.count"] = event_count,
    _["censor.derivative"] = censor_derivative
  );
}

// [[Rcpp::export]]
Rcpp::List ciftest_fg_iid_prefix_multi_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector strata_id,
    int code_event1,
    int code_event2,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericMatrix fh_weight,
    Rcpp::NumericVector g_at_competing,
    Rcpp::NumericMatrix g_event_stratum,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    double prob_bound
) {
  return ciftest_fg_iid_prefix_core(
    t, epsilon, x, weights, strata_id, code_event1, code_event2,
    code_censoring, event_times, fh_weight, g_at_competing,
    g_event_stratum, hazard_time, hazard_stratum, hazard, hazard_n_risk,
    prob_bound
  );
}

// [[Rcpp::export]]
Rcpp::List ciftest_fg_iid_prefix_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector strata_id,
    int code_event1,
    int code_event2,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericVector fh_weight,
    Rcpp::NumericVector g_at_competing,
    Rcpp::NumericMatrix g_event_stratum,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    double prob_bound
) {
  NumericMatrix weight_matrix(event_times.size(), 1);
  for (int j = 0; j < event_times.size(); ++j) {
    weight_matrix(j, 0) = fh_weight[j];
  }
  List result = ciftest_fg_iid_prefix_core(
    t, epsilon, x, weights, strata_id, code_event1, code_event2,
    code_censoring, event_times, weight_matrix, g_at_competing,
    g_event_stratum, hazard_time, hazard_stratum, hazard, hazard_n_risk,
    prob_bound
  );
  const int n = t.size();
  const int p = x.ncol();
  const int h_count = hazard_time.size();
  NumericVector base_array = result["score.iid.base"];
  NumericVector censor_array = result["score.iid.censor"];
  NumericVector derivative_array = result["censor.derivative"];
  NumericMatrix score_matrix = result["score"];
  NumericMatrix event_score_matrix = result["event.score"];
  NumericMatrix base_iid(n, p);
  NumericMatrix censor_iid(n, p);
  NumericMatrix censor_derivative(h_count, p);
  NumericVector score(p);
  NumericVector event_score(p);
  for (int q = 0; q < p; ++q) {
    score[q] = score_matrix(q, 0);
    event_score[q] = event_score_matrix(q, 0);
    for (int i = 0; i < n; ++i) {
      base_iid(i, q) = base_array[i + n * q];
      censor_iid(i, q) = censor_array[i + n * q];
    }
    for (int h = 0; h < h_count; ++h) {
      censor_derivative(h, q) = derivative_array[h + h_count * q];
    }
  }
  result["score"] = score;
  result["event.score"] = event_score;
  result["score.iid.base"] = base_iid;
  result["score.iid.censor"] = censor_iid;
  result["censor.derivative"] = censor_derivative;
  return result;
}

// [[Rcpp::export]]
Rcpp::List ciftest_augmentation_iid_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector censor_stratum_id,
    Rcpp::IntegerVector augmentation_cell_id,
    Rcpp::IntegerVector working_cell_id,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    Rcpp::NumericVector hazard_g_left,
    Rcpp::NumericVector working_survival,
    Rcpp::NumericVector working_cif2,
    Rcpp::NumericVector h_process,
    double prob_bound
) {
  const int n = t.size();
  const int h_count = hazard_time.size();
  IntegerVector h_dim = h_process.attr("dim");
  IntegerVector working_dim = working_survival.attr("dim");
  if (h_dim.size() != 3 || working_dim.size() != 2) {
    stop("Augmentation Rcpp kernel received arrays without expected dimensions.");
  }
  const int cell_count = h_dim[0];
  const int event_count = h_dim[1];
  const int p = h_dim[2];
  const int working_h = working_dim[0];
  const int working_count = working_dim[1];
  if (epsilon.size() != n || weights.size() != n ||
      censor_stratum_id.size() != n || augmentation_cell_id.size() != n ||
      working_cell_id.size() != n || event_times.size() != event_count ||
      hazard_stratum.size() != h_count || hazard.size() != h_count ||
      hazard_n_risk.size() != h_count || hazard_g_left.size() != h_count ||
      working_h != h_count || working_cif2.size() != working_survival.size()) {
    stop("Augmentation Rcpp kernel received incompatible dimensions.");
  }

  NumericMatrix augment_iid(n, p);
  NumericMatrix hstar(n, p);
  NumericVector hbar(p);
  double max_raw_center_error = 0.0;
  double min_working_survival = 1.0;
  double min_censoring_survival = 1.0;

  auto h_value = [&](int cell, int event, int q) {
    return h_process[cell + cell_count * (event + event_count * q)];
  };
  auto working_value = [&](const NumericVector& value, int h, int cell) {
    return value[h + h_count * cell];
  };

  for (int h = 0; h < h_count; ++h) {
    const double current = hazard_time[h];
    const int stratum = hazard_stratum[h];
    const double g_left = hazard_g_left[h];
    min_censoring_survival = std::min(min_censoring_survival, g_left);
    if (hazard_n_risk[h] <= prob_bound) continue;
    int event_index = std::lower_bound(
      event_times.begin(), event_times.end(), current
    ) - event_times.begin();
    std::fill(hstar.begin(), hstar.end(), 0.0);
    std::fill(hbar.begin(), hbar.end(), 0.0);
    double risk_weight = 0.0;

    for (int i = 0; i < n; ++i) {
      if (censor_stratum_id[i] != stratum) continue;
      const bool at_risk = t[i] >= current;
      if (at_risk) risk_weight += weights[i];
      const int working_cell = working_cell_id[i] - 1;
      const int augmentation_cell = augmentation_cell_id[i] - 1;
      if (working_cell < 0 || working_cell >= working_count ||
          augmentation_cell < 0 || augmentation_cell >= cell_count) {
        stop("Augmentation Rcpp kernel received an invalid cell code.");
      }
      const double survival = working_value(
        working_survival, h, working_cell
      );
      const double cif2 = working_value(working_cif2, h, working_cell);
      min_working_survival = std::min(min_working_survival, survival);
      if (event_index < event_count && cif2 > 0.0) {
        bool active = false;
        for (int q = 0; q < p; ++q) {
          if (std::abs(h_value(augmentation_cell, event_index, q)) >
              std::sqrt(std::numeric_limits<double>::epsilon())) {
            active = true;
            break;
          }
        }
        if (active && (survival <= prob_bound || g_left <= prob_bound)) {
          stop("Working-model positivity is violated in the closed-form augmentation.");
        }
        if (active) {
          const double ratio = cif2 / (survival * g_left);
          for (int q = 0; q < p; ++q) {
            hstar(i, q) = ratio *
              h_value(augmentation_cell, event_index, q);
          }
        }
      }
      if (at_risk) {
        for (int q = 0; q < p; ++q) {
          hbar[q] += hstar(i, q) * weights[i];
        }
      }
    }
    if (risk_weight <= prob_bound) continue;
    for (int q = 0; q < p; ++q) hbar[q] /= risk_weight;

    NumericVector raw_increment(p);
    NumericVector centered_increment(p);
    for (int i = 0; i < n; ++i) {
      const bool in_stratum = censor_stratum_id[i] == stratum;
      const bool at_risk = in_stratum && t[i] >= current;
      const bool censored = in_stratum && t[i] == current &&
        epsilon[i] == code_censoring;
      const double martingale = weights[i] *
        ((censored ? 1.0 : 0.0) - hazard[h] * (at_risk ? 1.0 : 0.0));
      for (int q = 0; q < p; ++q) {
        const double centered = hstar(i, q) - hbar[q];
        raw_increment[q] += hstar(i, q) * martingale;
        centered_increment[q] += centered * martingale;
        augment_iid(i, q) += centered * martingale;
      }
    }
    for (int q = 0; q < p; ++q) {
      max_raw_center_error = std::max(
        max_raw_center_error,
        std::abs(raw_increment[q] - centered_increment[q])
      );
    }
  }
  return List::create(
    _["score.iid.augment"] = augment_iid,
    _["augmentation.centering.error"] = max_raw_center_error,
    _["minimum.working.survival"] = min_working_survival,
    _["minimum.censoring.survival"] = min_censoring_survival
  );
}

// Batch-only multiweight augmentation. The working probabilities, censoring
// martingales, and risk-set means are traversed once for all FH weights.
// [[Rcpp::export]]
Rcpp::List ciftest_augmentation_iid_multi_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector censor_stratum_id,
    Rcpp::IntegerVector augmentation_cell_id,
    Rcpp::IntegerVector working_cell_id,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    Rcpp::NumericVector hazard_g_left,
    Rcpp::NumericVector working_survival,
    Rcpp::NumericVector working_cif2,
    Rcpp::NumericVector h_process,
    double prob_bound
) {
  const int n = t.size();
  const int h_count = hazard_time.size();
  IntegerVector h_dim = h_process.attr("dim");
  IntegerVector working_dim = working_survival.attr("dim");
  if (h_dim.size() != 4 || working_dim.size() != 2) {
    stop("Multiweight augmentation kernel received unexpected array dimensions.");
  }
  const int cell_count = h_dim[0];
  const int event_count = h_dim[1];
  const int p = h_dim[2];
  const int weight_count = h_dim[3];
  const int working_h = working_dim[0];
  const int working_count = working_dim[1];
  if (epsilon.size() != n || weights.size() != n ||
      censor_stratum_id.size() != n || augmentation_cell_id.size() != n ||
      working_cell_id.size() != n || event_times.size() != event_count ||
      hazard_stratum.size() != h_count || hazard.size() != h_count ||
      hazard_n_risk.size() != h_count || hazard_g_left.size() != h_count ||
      working_h != h_count || working_cif2.size() != working_survival.size() ||
      weight_count < 1) {
    stop("Multiweight augmentation kernel received incompatible dimensions.");
  }

  NumericVector augment_iid(n * p * weight_count);
  augment_iid.attr("dim") = IntegerVector::create(n, p, weight_count);
  NumericVector hstar(n * p * weight_count);
  NumericMatrix hbar(p, weight_count);
  NumericVector max_raw_center_error(weight_count);
  double min_working_survival = 1.0;
  double min_censoring_survival = 1.0;

  auto score_index = [=](int i, int q, int w) {
    return i + n * (q + p * w);
  };
  auto h_value = [=](int cell, int event, int q, int w) {
    return h_process[cell + cell_count *
      (event + event_count * (q + p * w))];
  };
  auto working_value = [=](const NumericVector& value, int h, int cell) {
    return value[h + h_count * cell];
  };

  for (int h = 0; h < h_count; ++h) {
    const double current = hazard_time[h];
    const int stratum = hazard_stratum[h];
    const double g_left = hazard_g_left[h];
    min_censoring_survival = std::min(min_censoring_survival, g_left);
    if (hazard_n_risk[h] <= prob_bound) continue;
    const int event_index = std::lower_bound(
      event_times.begin(), event_times.end(), current
    ) - event_times.begin();
    std::fill(hstar.begin(), hstar.end(), 0.0);
    std::fill(hbar.begin(), hbar.end(), 0.0);
    double risk_weight = 0.0;

    for (int i = 0; i < n; ++i) {
      if (censor_stratum_id[i] != stratum) continue;
      const bool at_risk = t[i] >= current;
      if (at_risk) risk_weight += weights[i];
      const int working_cell = working_cell_id[i] - 1;
      const int augmentation_cell = augmentation_cell_id[i] - 1;
      if (working_cell < 0 || working_cell >= working_count ||
          augmentation_cell < 0 || augmentation_cell >= cell_count) {
        stop("Multiweight augmentation kernel received an invalid cell code.");
      }
      const double survival = working_value(
        working_survival, h, working_cell
      );
      const double cif2 = working_value(working_cif2, h, working_cell);
      min_working_survival = std::min(min_working_survival, survival);
      if (event_index < event_count && cif2 > 0.0) {
        for (int w = 0; w < weight_count; ++w) {
          bool active = false;
          for (int q = 0; q < p; ++q) {
            if (std::abs(h_value(augmentation_cell, event_index, q, w)) >
                std::sqrt(std::numeric_limits<double>::epsilon())) {
              active = true;
              break;
            }
          }
          if (active && (survival <= prob_bound || g_left <= prob_bound)) {
            stop("Working-model positivity is violated in the multiweight augmentation.");
          }
          if (active) {
            const double ratio = cif2 / (survival * g_left);
            for (int q = 0; q < p; ++q) {
              hstar[score_index(i, q, w)] = ratio *
                h_value(augmentation_cell, event_index, q, w);
            }
          }
        }
      }
      if (at_risk) {
        for (int w = 0; w < weight_count; ++w) {
          for (int q = 0; q < p; ++q) {
            hbar(q, w) += hstar[score_index(i, q, w)] * weights[i];
          }
        }
      }
    }
    if (risk_weight <= prob_bound) continue;
    for (int w = 0; w < weight_count; ++w) {
      for (int q = 0; q < p; ++q) hbar(q, w) /= risk_weight;
    }

    NumericMatrix raw_increment(p, weight_count);
    NumericMatrix centered_increment(p, weight_count);
    for (int i = 0; i < n; ++i) {
      const bool in_stratum = censor_stratum_id[i] == stratum;
      const bool at_risk = in_stratum && t[i] >= current;
      const bool censored = in_stratum && t[i] == current &&
        epsilon[i] == code_censoring;
      const double martingale = weights[i] *
        ((censored ? 1.0 : 0.0) - hazard[h] * (at_risk ? 1.0 : 0.0));
      for (int w = 0; w < weight_count; ++w) {
        for (int q = 0; q < p; ++q) {
          const double raw = hstar[score_index(i, q, w)];
          const double centered = raw - hbar(q, w);
          raw_increment(q, w) += raw * martingale;
          centered_increment(q, w) += centered * martingale;
          augment_iid[score_index(i, q, w)] += centered * martingale;
        }
      }
    }
    for (int w = 0; w < weight_count; ++w) {
      for (int q = 0; q < p; ++q) {
        max_raw_center_error[w] = std::max(
          max_raw_center_error[w],
          std::abs(raw_increment(q, w) - centered_increment(q, w))
        );
      }
    }
  }
  return List::create(
    _["score.iid.augment"] = augment_iid,
    _["augmentation.centering.error"] = max_raw_center_error,
    _["minimum.working.survival"] = min_working_survival,
    _["minimum.censoring.survival"] = min_censoring_survival
  );
}

// Cell-aggregated prefix implementation of the closed-form augmentation.
// H*(c) is constant within an augmentation cell, so risk-set centering can be
// performed from cell-specific reverse risk sums. Subject-level martingale
// contributions then reduce to a cell-specific censoring-hazard prefix plus
// the observed censoring jump.
static Rcpp::List ciftest_augmentation_iid_prefix_core(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector censor_stratum_id,
    Rcpp::IntegerVector augmentation_cell_id,
    Rcpp::IntegerVector working_cell_id,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    Rcpp::NumericVector hazard_g_left,
    Rcpp::NumericVector working_survival,
    Rcpp::NumericVector working_cif2,
    Rcpp::NumericVector h_process,
    double prob_bound
) {
  const int n = t.size();
  const int h_count = hazard_time.size();
  IntegerVector h_dim = h_process.attr("dim");
  IntegerVector working_dim = working_survival.attr("dim");
  if (h_dim.size() != 4 || working_dim.size() != 2) {
    stop("Prefix augmentation kernel received unexpected array dimensions.");
  }
  const int cell_count = h_dim[0];
  const int event_count = h_dim[1];
  const int p = h_dim[2];
  const int weight_count = h_dim[3];
  const int working_h = working_dim[0];
  const int working_count = working_dim[1];
  if (epsilon.size() != n || weights.size() != n ||
      censor_stratum_id.size() != n || augmentation_cell_id.size() != n ||
      working_cell_id.size() != n || event_times.size() != event_count ||
      hazard_stratum.size() != h_count || hazard.size() != h_count ||
      hazard_n_risk.size() != h_count || hazard_g_left.size() != h_count ||
      working_h != h_count || working_cif2.size() != working_survival.size() ||
      cell_count < 1 || working_count < 1 || weight_count < 1) {
    stop("Prefix augmentation kernel received incompatible dimensions.");
  }

  int stratum_count = 0;
  for (int i = 0; i < n; ++i) {
    if (censor_stratum_id[i] < 1 || augmentation_cell_id[i] < 1 ||
        augmentation_cell_id[i] > cell_count || working_cell_id[i] < 1 ||
        working_cell_id[i] > working_count) {
      stop("Prefix augmentation kernel received an invalid cell code.");
    }
    stratum_count = std::max(stratum_count, censor_stratum_id[i]);
  }
  for (int h = 0; h < h_count; ++h) {
    if (hazard_stratum[h] < 1) {
      stop("Prefix augmentation kernel received an invalid hazard stratum.");
    }
    stratum_count = std::max(stratum_count, hazard_stratum[h]);
  }

  std::vector<int> cell_stratum(cell_count, -1);
  std::vector<int> cell_working(cell_count, -1);
  std::vector<std::vector<int> > subjects_by_cell(cell_count);
  std::vector<std::vector<int> > cells_by_stratum(stratum_count);
  std::vector<std::vector<int> > hazards_by_stratum(stratum_count);
  for (int i = 0; i < n; ++i) {
    const int cell = augmentation_cell_id[i] - 1;
    const int stratum = censor_stratum_id[i] - 1;
    const int working = working_cell_id[i] - 1;
    if (cell_stratum[cell] < 0) {
      cell_stratum[cell] = stratum;
      cell_working[cell] = working;
    } else if (cell_stratum[cell] != stratum ||
               cell_working[cell] != working) {
      stop("An augmentation cell maps to inconsistent nuisance cells.");
    }
    subjects_by_cell[cell].push_back(i);
  }
  for (int cell = 0; cell < cell_count; ++cell) {
    if (cell_stratum[cell] < 0) {
      stop("Prefix augmentation kernel received an empty augmentation cell.");
    }
    cells_by_stratum[cell_stratum[cell]].push_back(cell);
    std::sort(
      subjects_by_cell[cell].begin(), subjects_by_cell[cell].end(),
      [&](int left, int right) {
        if (t[left] != t[right]) return t[left] < t[right];
        return left < right;
      }
    );
  }
  for (int h = 0; h < h_count; ++h) {
    hazards_by_stratum[hazard_stratum[h] - 1].push_back(h);
  }
  for (int s = 0; s < stratum_count; ++s) {
    std::sort(
      hazards_by_stratum[s].begin(), hazards_by_stratum[s].end(),
      [&](int left, int right) {
        if (hazard_time[left] != hazard_time[right]) {
          return hazard_time[left] < hazard_time[right];
        }
        return left < right;
      }
    );
  }

  auto cell_h_index = [=](int h, int cell) {
    return h + h_count * cell;
  };
  auto process_index = [=](int cell, int event, int q, int w) {
    return cell + cell_count *
      (event + event_count * (q + p * w));
  };
  auto centered_index = [=](int cell, int h, int q, int w) {
    return cell + cell_count * (h + h_count * (q + p * w));
  };
  auto score_index = [=](int i, int q, int w) {
    return i + n * (q + p * w);
  };
  auto working_value = [=](const NumericVector& value, int h, int cell) {
    return value[h + h_count * cell];
  };

  std::vector<double> cell_risk(h_count * cell_count, 0.0);
  std::vector<double> cell_censor(h_count * cell_count, 0.0);
  for (int cell = 0; cell < cell_count; ++cell) {
    const std::vector<int>& subjects = subjects_by_cell[cell];
    const std::vector<int>& hazards = hazards_by_stratum[cell_stratum[cell]];
    int subject_pos = static_cast<int>(subjects.size()) - 1;
    double running_risk = 0.0;
    for (int local_h = static_cast<int>(hazards.size()) - 1;
         local_h >= 0; --local_h) {
      const int h = hazards[local_h];
      while (subject_pos >= 0 &&
             t[subjects[subject_pos]] >= hazard_time[h]) {
        running_risk += weights[subjects[subject_pos]];
        --subject_pos;
      }
      cell_risk[cell_h_index(h, cell)] = running_risk;
    }
    for (int subject : subjects) {
      if (epsilon[subject] != code_censoring) continue;
      auto first = std::lower_bound(
        hazards.begin(), hazards.end(), t[subject],
        [&](int h, double value) { return hazard_time[h] < value; }
      );
      while (first != hazards.end() && hazard_time[*first] == t[subject]) {
        cell_censor[cell_h_index(*first, cell)] += weights[subject];
        ++first;
      }
    }
  }

  std::vector<double> raw(
    cell_count * h_count * p * weight_count, 0.0
  );
  std::vector<double> centered(
    cell_count * h_count * p * weight_count, 0.0
  );
  NumericVector max_raw_center_error(weight_count);
  double min_working_survival = 1.0;
  double min_censoring_survival = 1.0;
  const double active_bound = std::sqrt(
    std::numeric_limits<double>::epsilon()
  );

  for (int h = 0; h < h_count; ++h) {
    const int stratum = hazard_stratum[h] - 1;
    const double g_left = hazard_g_left[h];
    min_censoring_survival = std::min(min_censoring_survival, g_left);
    if (hazard_n_risk[h] <= prob_bound) continue;
    const int event_index = std::lower_bound(
      event_times.begin(), event_times.end(), hazard_time[h]
    ) - event_times.begin();
    double risk_weight = 0.0;
    std::vector<double> hbar(p * weight_count, 0.0);

    for (int cell : cells_by_stratum[stratum]) {
      const int working = cell_working[cell];
      const double survival = working_value(
        working_survival, h, working
      );
      const double cif2 = working_value(working_cif2, h, working);
      min_working_survival = std::min(min_working_survival, survival);
      const double risk = cell_risk[cell_h_index(h, cell)];
      risk_weight += risk;
      if (event_index >= event_count || cif2 <= 0.0) continue;
      for (int w = 0; w < weight_count; ++w) {
        bool active = false;
        for (int q = 0; q < p; ++q) {
          if (std::abs(h_process[
                process_index(cell, event_index, q, w)
              ]) > active_bound) {
            active = true;
            break;
          }
        }
        if (active && (survival <= prob_bound || g_left <= prob_bound)) {
          if (weight_count == 1) {
            stop("Working-model positivity is violated in the closed-form augmentation.");
          }
          stop("Working-model positivity is violated in the multiweight augmentation.");
        }
        if (!active) continue;
        const double ratio = cif2 / (survival * g_left);
        for (int q = 0; q < p; ++q) {
          const double value = ratio * h_process[
            process_index(cell, event_index, q, w)
          ];
          raw[centered_index(cell, h, q, w)] = value;
          hbar[q + p * w] += risk * value;
        }
      }
    }
    if (risk_weight <= prob_bound) continue;
    for (double& value : hbar) value /= risk_weight;

    std::vector<double> raw_increment(p * weight_count, 0.0);
    std::vector<double> centered_increment(p * weight_count, 0.0);
    for (int cell : cells_by_stratum[stratum]) {
      const double martingale =
        cell_censor[cell_h_index(h, cell)] -
        hazard[h] * cell_risk[cell_h_index(h, cell)];
      for (int w = 0; w < weight_count; ++w) {
        for (int q = 0; q < p; ++q) {
          const int index = centered_index(cell, h, q, w);
          const int local = q + p * w;
          centered[index] = raw[index] - hbar[local];
          raw_increment[local] += raw[index] * martingale;
          centered_increment[local] += centered[index] * martingale;
        }
      }
    }
    for (int w = 0; w < weight_count; ++w) {
      for (int q = 0; q < p; ++q) {
        const int local = q + p * w;
        max_raw_center_error[w] = std::max(
          max_raw_center_error[w],
          std::abs(raw_increment[local] - centered_increment[local])
        );
      }
    }
  }

  std::vector<double> prefix(
    cell_count * h_count * p * weight_count, 0.0
  );
  for (int cell = 0; cell < cell_count; ++cell) {
    const std::vector<int>& hazards = hazards_by_stratum[cell_stratum[cell]];
    std::vector<double> running(p * weight_count, 0.0);
    for (int h : hazards) {
      for (int w = 0; w < weight_count; ++w) {
        for (int q = 0; q < p; ++q) {
          const int local = q + p * w;
          const int index = centered_index(cell, h, q, w);
          running[local] += hazard[h] * centered[index];
          prefix[index] = running[local];
        }
      }
    }
  }

  NumericVector augment_iid(n * p * weight_count);
  augment_iid.attr("dim") = IntegerVector::create(n, p, weight_count);
  for (int i = 0; i < n; ++i) {
    const int cell = augmentation_cell_id[i] - 1;
    const std::vector<int>& hazards = hazards_by_stratum[cell_stratum[cell]];
    const int through = std::upper_bound(
      hazards.begin(), hazards.end(), t[i],
      [&](double value, int h) { return value < hazard_time[h]; }
    ) - hazards.begin();
    for (int w = 0; w < weight_count; ++w) {
      for (int q = 0; q < p; ++q) {
        double value = 0.0;
        if (through > 0) {
          value -= weights[i] * prefix[
            centered_index(cell, hazards[through - 1], q, w)
          ];
        }
        if (epsilon[i] == code_censoring) {
          int equal = through - 1;
          while (equal >= 0 && hazard_time[hazards[equal]] == t[i]) {
            value += weights[i] * centered[
              centered_index(cell, hazards[equal], q, w)
            ];
            --equal;
          }
        }
        augment_iid[score_index(i, q, w)] = value;
      }
    }
  }

  return List::create(
    _["score.iid.augment"] = augment_iid,
    _["augmentation.centering.error"] = max_raw_center_error,
    _["minimum.working.survival"] = min_working_survival,
    _["minimum.censoring.survival"] = min_censoring_survival
  );
}

// [[Rcpp::export]]
Rcpp::List ciftest_augmentation_iid_prefix_multi_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector censor_stratum_id,
    Rcpp::IntegerVector augmentation_cell_id,
    Rcpp::IntegerVector working_cell_id,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    Rcpp::NumericVector hazard_g_left,
    Rcpp::NumericVector working_survival,
    Rcpp::NumericVector working_cif2,
    Rcpp::NumericVector h_process,
    double prob_bound
) {
  return ciftest_augmentation_iid_prefix_core(
    t, epsilon, weights, censor_stratum_id, augmentation_cell_id,
    working_cell_id, code_censoring, event_times, hazard_time,
    hazard_stratum, hazard, hazard_n_risk, hazard_g_left,
    working_survival, working_cif2, h_process, prob_bound
  );
}

// [[Rcpp::export]]
Rcpp::List ciftest_augmentation_iid_prefix_kernel_cpp(
    Rcpp::NumericVector t,
    Rcpp::IntegerVector epsilon,
    Rcpp::NumericVector weights,
    Rcpp::IntegerVector censor_stratum_id,
    Rcpp::IntegerVector augmentation_cell_id,
    Rcpp::IntegerVector working_cell_id,
    int code_censoring,
    Rcpp::NumericVector event_times,
    Rcpp::NumericVector hazard_time,
    Rcpp::IntegerVector hazard_stratum,
    Rcpp::NumericVector hazard,
    Rcpp::NumericVector hazard_n_risk,
    Rcpp::NumericVector hazard_g_left,
    Rcpp::NumericVector working_survival,
    Rcpp::NumericVector working_cif2,
    Rcpp::NumericVector h_process,
    double prob_bound
) {
  IntegerVector h_dim = h_process.attr("dim");
  if (h_dim.size() != 3) {
    stop("Single-weight prefix augmentation requires a three-dimensional H process.");
  }
  NumericVector h_array(h_process.size());
  std::copy(h_process.begin(), h_process.end(), h_array.begin());
  h_array.attr("dim") = IntegerVector::create(
    h_dim[0], h_dim[1], h_dim[2], 1
  );
  List result = ciftest_augmentation_iid_prefix_core(
    t, epsilon, weights, censor_stratum_id, augmentation_cell_id,
    working_cell_id, code_censoring, event_times, hazard_time,
    hazard_stratum, hazard, hazard_n_risk, hazard_g_left,
    working_survival, working_cif2, h_array, prob_bound
  );
  const int n = t.size();
  const int p = h_dim[2];
  NumericVector iid_array = result["score.iid.augment"];
  NumericMatrix iid(n, p);
  for (int q = 0; q < p; ++q) {
    for (int i = 0; i < n; ++i) iid(i, q) = iid_array[i + n * q];
  }
  NumericVector error = result["augmentation.centering.error"];
  result["score.iid.augment"] = iid;
  result["augmentation.centering.error"] = error[0];
  return result;
}

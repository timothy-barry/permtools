#include <Rcpp.h>
#include <boost/random.hpp>
#include <boost/math/distributions/binomial.hpp>
using boost::math::binomial;
#include <stdexcept>
#include "utilities.h"
#include <algorithm>
using namespace Rcpp;

double compute_treated_mean(int n_trt, double n_trt_double, const std::vector<int>& trt_idxs, const NumericVector& y) {
  double sum = 0;
  for (int i = 0; i < n_trt; i ++) sum += y[trt_idxs[i]];
  return sum/n_trt_double;
}


// [[Rcpp::export]]
List run_bc_multiple_test(IntegerVector x, List Y, int h, double alpha) {
  // define variables and objects
  int n = x.length(), m = Y.length(), max_n_losses_active_set;
  std::vector<bool> active_set(m, true), futility_set(m, false), rejected_set(m, false);
  std::vector<double> stop_times(m), original_statistics(m), p_values(m);
  std::vector<int> n_losses(m, 0);
  double n_trt_double, curr_test_stat, t = 0, h_doub = static_cast<double>(h), m_doub = static_cast<double>(m), threshold, n_in_active_set;

  // populate the trt_idxs vector
  std::vector<int> trt_idxs;
  for (int i = 0; i < x.length(); i++) if (x[i] == 1) trt_idxs.push_back(i);
  int n_trt = trt_idxs.size();
  n_trt_double = static_cast<double>(n_trt);
  if (n_trt == 0) throw std::invalid_argument("Zero treatment units.");

  // compute the original test statistics
  for (int i = 0; i < m; i++) {
    original_statistics[i] = compute_treated_mean(n_trt, n_trt_double, trt_idxs, Y(i));
  }

  // define objects related to random permutations
  std::vector<int> random_samp(n);
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);
  std::vector<double> i_doub_array(n_trt);
  for (int i = 0; i < n_trt; i ++) i_doub_array[i] = static_cast<double>(i);

  // iterate through time
  for (int k = 0; k < 50000; k++) {
    // increment time
    t++;
    // generate a random permutation
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n, n_trt);
    // iterate over hypotheses
    for (int i = 0; i < m; i ++) {
      // if in the active set, compute the test statistic and update gamma
      if (active_set[i]) {
        // compute the test statistic
        curr_test_stat = compute_treated_mean(n_trt, n_trt_double, random_samp, Y(i));
        // determine whether we have a loss; if so, increment n_losses
        n_losses[i] += (curr_test_stat >= original_statistics[i] ? 1.0 : 0.0);
      }
    }

    // move elements from active set into futility set
    max_n_losses_active_set = 0;
    n_in_active_set = 0;
    for (int i = 0; i < m; i++) {
      if (active_set[i]) {
        if (n_losses[i] == h) { // hit loss limit
          active_set[i] = false;
          futility_set[i] = true;
          stop_times[i] = t;
        } else {
          n_in_active_set ++;
          if (n_losses[i] > max_n_losses_active_set) max_n_losses_active_set = n_losses[i];
        }
      }
    }

    // check for rejection of active set hypotheses
    threshold = t + h_doub - (h_doub * m_doub)/(n_in_active_set * alpha);
    if (static_cast<double>(max_n_losses_active_set) <= threshold) {
      for (int i = 0; i < m; i ++) {
        if (active_set[i]) {
          active_set[i] = false;
          rejected_set[i] = true;
          stop_times[i] = t;
        }
      }
      break;
    }
  }
  // compute the p-values
  for (int i = 0; i < m; i ++) {
    p_values[i] = h_doub/(stop_times[i] - static_cast<double>(n_losses[i]) + h);
  }

  return List::create(Named("p_values") = p_values, Named("rejected") = rejected_set);
}

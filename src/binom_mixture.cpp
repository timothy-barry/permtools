#include <Rcpp.h>
#include <boost/random.hpp>
#include <boost/math/distributions/binomial.hpp>
using boost::math::binomial;
#include <stdexcept>
#include "utilities.h"
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List run_binom_mixture_perm_test(IntegerVector x, NumericVector y, double alpha) {
  // define variables
  int n = x.length(), n_trt = 0, n_cntrl = 0, min_permutations_futility = 50,
    max_permutations = 10000, rank = 1, n_permutations = 0;
  double total_sum = 0, trt_sum = 0, cntrl_sum = 0, orig_test_stat = 0, curr_test_stat = 0,
    n_trt_double = 0, n_cntrl_double = 0, c = 0.9 * alpha, min_wealth_threshold = alpha,
    curr_wealth = 0, max_wealth_threshold = 1/alpha;
  bool stop = false;
  IntegerVector trt_idxs;

  // populate the trt_idxs vector
  for (int i = 0; i < x.length(); i++) if (x[i] == 1) trt_idxs.push_back(i);
  n_trt = trt_idxs.length(); n_cntrl = n - n_trt;
  if (n_trt == 0) throw std::invalid_argument("Zero treatment units.");
  n_trt_double = static_cast<double>(n_trt); n_cntrl_double = static_cast<double>(n_cntrl);

  // define objects related to random permutations
  std::vector<int> random_samp(n);
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);
  std::vector<double> i_doub_array(n_trt);
  for (int i = 0; i < n_trt; i ++) i_doub_array[i] = static_cast<double>(i);

  // compute original test statistic
  for (int i = 0; i < n_trt; i++) trt_sum += y[trt_idxs[i]];
  for (int i = 0; i < n; i++) total_sum += y[i];
  cntrl_sum = total_sum - trt_sum;
  orig_test_stat = trt_sum/n_trt_double - cntrl_sum/n_cntrl_double;

  // iterate
  stop = false;
  while (!stop) {
    // increment n_permutations
    n_permutations ++;

    // sample permutation
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n, n_trt);

    // compute current test statistic
    trt_sum = 0;
    for (int i = 0; i < n_trt; i++) trt_sum += y[random_samp[i]];
    cntrl_sum = total_sum - trt_sum;
    curr_test_stat = trt_sum/n_trt_double - cntrl_sum/n_cntrl_double;

    // check for a loss
    if (curr_test_stat >= orig_test_stat) rank ++;

    // compute and evaluate wealth
    curr_wealth = cdf(complement(binomial(n_permutations + 1, c), rank - 1))/c;

    if (curr_wealth >= max_wealth_threshold) {
      stop = true;
    } else if (curr_wealth <= min_wealth_threshold && n_permutations >= min_permutations_futility) {
      stop = true;
    } else if (n_permutations == max_permutations) {
      stop = true;
      if (curr_wealth >= distribution(generator)/alpha) curr_wealth = 1/alpha;
    }
  }
  double p = std::min(1/curr_wealth, 1.0);
  return List::create(Named("p") = p, Named("n_permutations") = n_permutations);
}

// [[Rcpp::export]]
double compute_binomial_right_tail_prob(int n, double p, int k) {
  return cdf(complement(binomial(n, p), k));
}

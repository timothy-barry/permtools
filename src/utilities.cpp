#include <boost/random.hpp>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


void draw_wor_sample(boost::random::mt19937& generator, boost::random::uniform_real_distribution<double>& distribution, const std::vector<double>& i_doub_array, std::vector<int>& random_samp, int n, int n_trt) {
  for (int i = 0; i < n; i ++) random_samp[i] = i;
  double u, temp, n_doub = static_cast<double>(n);
  int pos;
  // perform the swap
  for (int i = 0; i < n_trt; i ++) {
    u = distribution(generator);
    pos = floor((n_doub - i_doub_array[i]) * u);
    temp = random_samp[pos];
    random_samp[pos] = random_samp[n - i - 1];
    random_samp[n - i - 1] = temp;
  }
  return;
}

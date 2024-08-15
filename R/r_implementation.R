# this script contains an R implementation of the binomial mixture and Besag-Clifford adaptive permutation tests.
# for the time being, we assume a binary treatment variable and a difference in means test statistic.

#' Run binomial mixture permutation test
#'
#' @param x a binary, integer treatment vector of 1s and 0s
#' @param y a numeric response vector
#'
#' @return a list containing (i) the p-value and (ii) the numer of permutations
#' @export
#'
#' @examples
#' set.seed(4)
#' n <- 100
#' y <- rnorm(100)
#' x <- y > median(y)
#' swap_idxs <- sample(x = seq_along(x), size = 35)
#' x[swap_idxs] <- !x[swap_idxs]
#' target_alpha <- 0.05
#' res <- run_binom_mixture_perm_test(x, y, target_alpha)
run_binom_mixture_perm_test_r <- function(x, y, target_alpha = 0.05, b = 0.9) {
  set.seed(4)
  # compute n nonzero entries
  n <- length(x)
  n_trt <- sum(x)
  n_cntrl <- n - n_trt
  sum_y <- sum(y)

  # compute original test statistic
  trt_sum <- sum(y[x])
  cntrl_sum <- sum_y - trt_sum
  orig_test_stat <- trt_sum/n_trt - cntrl_sum/n_cntrl

  # initialize tuning parameters
  c <- b * target_alpha
  min_permutations_futility <- 50L
  max_permutations <- 10000L
  min_wealth_threshold <- target_alpha

  # initialize vars related to wealth
  max_wealth <- 0
  max_wealth_threshold <- 1/target_alpha
  rank <- 1L
  n_permutations <- 0L
  stop <- FALSE

  # iterate
  while (!stop) {
    # increment n_permutations
    n_permutations <- n_permutations + 1L
    # permute x
    x_perm <- sample(x = seq_len(n), size = n_trt, replace = FALSE)
    # compute the current test statistic
    trt_sum <- sum(y[x_perm])
    cntrl_sum <- sum_y - trt_sum
    curr_test_stat <- trt_sum/n_trt - cntrl_sum/n_cntrl

    # check for a loss
    if (curr_test_stat >= orig_test_stat) rank <- rank + 1L

    # if n_permutations is multiple of 5, compute wealth and update stopping rule
    curr_wealth <- (1 - pbinom(rank - 1L, n_permutations + 1L, c))/c
    if (curr_wealth > max_wealth) max_wealth <- curr_wealth
    # stop if (i) rejection, (ii) futility, (iii) max permutations
    if (curr_wealth >= max_wealth_threshold) {
      stop <- TRUE
    } else if (curr_wealth <= min_wealth_threshold && n_permutations >= min_permutations_futility) {
      stop <- TRUE
    } else if (n_permutations == max_permutations) {
      stop <- TRUE
      if (curr_wealth >= runif(1)/target_alpha) max_wealth <- 1/target_alpha
    }
  }
  out <- list(p = min(1, 1/max_wealth), n_permutations = n_permutations)
  return(out)
}


#' Run binomial mixture permutation test
#'
#' @param x a binary, integer treatment vector of 1s and 0s
#' @param y a numeric response vector
#' @param h tuning parameter; the number of losses to incur until stopping the algorithm
#'
#' @return a list containing (i) the p-value and (ii) the numer of permutations
#' @export
#'
#' @examples
#' set.seed(4)
#' n <- 100
#' y <- rnorm(100)
#' x <- y > median(y)
#' swap_idxs <- sample(x = seq_along(x), size = 35)
#' x[swap_idxs] <- !x[swap_idxs]
#' res <- run_bc_perm_test_r(x, y)
run_bc_perm_test_r <- function(x, y, h = 15L) {
  set.seed(4)
  # compute n nonzero entries
  n <- length(x)
  n_trt <- sum(x)
  n_cntrl <- n - n_trt
  sum_y <- sum(y)

  # compute original test statistic
  trt_sum <- sum(y[x])
  cntrl_sum <- sum_y - trt_sum
  orig_test_stat <- trt_sum/n_trt - cntrl_sum/n_cntrl

  # set vars
  n_permutations <- 0L
  n_losses <- 0L
  max_permutations <- 10000L

  # iterate
  stop <- FALSE
  while (!stop) {
    # increment n_permutations
    n_permutations <- n_permutations + 1L
    # permute x
    x_perm <- sample(x = seq_len(n), size = n_trt, replace = FALSE)
    # compute the current test statistic
    trt_sum <- sum(y[x_perm])
    cntrl_sum <- sum_y - trt_sum
    curr_test_stat <- trt_sum/n_trt - cntrl_sum/n_cntrl
    # check for a loss
    if (curr_test_stat >= orig_test_stat) n_losses <- n_losses + 1L
    # determine stop
    if (n_losses == h) {
      p <- n_losses/n_permutations
      stop <- TRUE
    } else if (n_permutations == max_permutations) {
      p <- (n_losses + 1)/(max_permutations + 1)
      stop <- TRUE
    }
  }
  return(list(p = p, n_permutations = n_permutations))
}

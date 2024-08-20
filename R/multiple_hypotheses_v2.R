#' Run BC multiple testing procedure
#'
#' Runs the permutation test based on the anytime-valid BC p-value and BH correction.
#'
#' @param x an n by m matrix of treatments
#' @param Y an n by m matrix of responses
#' @param h the tuning parameter anytime-balid BC p-value
#'
#' @return a list containing (1) the 1-based indices of the hypotheses that are significant, (2) the number of permutations for each hypothesis, and (3) the p-value for each hypothesis
#' @export
#'
#' @examples
#' n <- 100L
#' m <- 1000L
#' x <- sample(x = c(0L, 1L), size = n, replace = TRUE)
#' under_null <- sample(c(TRUE, FALSE), size = m, prob = c(0.95, 0.05), replace = TRUE)
#' Y <- sapply(X = seq_len(m), FUN = function(i) {
#'   rnorm(n = n, mean = 2 + (if (under_null[i]) 0 else 2 * x), sd = 1)
#' })
#' h <- 15L; alpha <- 0.1
#' res_r <- run_bc_multiple_test_r_v2(x, Y, alpha = alpha)
#' mean(under_null[res_r |> filter(significant) |> pull(hypothesis_idx)])
#'
#' # c++ version
#' Y_list <- apply(X = Y, MARGIN = 2, FUN = function(col) col, simplify = FALSE)
#' system.time(res <- run_bc_multiple_test(x, Y_list, h, alpha))
#' mean(under_null[which(res)])
run_bc_multiple_test_r_v2 <- function(x, Y, h = 15L, alpha = 0.1) {
  set.seed(4)
  # initialize other variables
  n <- nrow(Y) # sample size
  m <- ncol(Y) # number of hypotheses
  t <- 0L # current time

  # initialize active, futility, and rejected sets
  active_set <- !logical(m)
  futility_set <- logical(m)
  rejected_set <- logical(m)
  stop_times <- integer(m)
  original_test_statistics <- integer(length = m)
  n_losses <- integer(m)

  # function to compute the test statistic
  compute_test_stat <- function(x, y) mean(y[x == 1L]) - mean(y[x == 0L])

  # compute the original test statistics
  for (i in seq_len(m)) original_test_statistics[i] <- compute_test_stat(x, Y[,i])

  # iterate through time
  for (k in seq_len(50000)) {
    # increment t
    t <- t + 1L
    # generate a random permutation of 1...n
    perm_sample <- sample(x = seq(1L, n))
    # iterate over the hypotheses
    for (i in seq_len(m)) {
      # if in the active set, compute the test statistic and update gamma
      if (active_set[i]) {
        # compute the test statistic
        curr_test_stat <- compute_test_stat(x = x[perm_sample], y = Y[,i])
        # determine whether we have a loss; if so, increment n_losses
        n_losses[i] <- n_losses[i] + as.integer(curr_test_stat >= original_test_statistics[i])
      }
    }

    # a. move hypotheses that have failed due to futility from active set to futility set
    # additionally,
    max_n_losses_active_set <- 0L
    n_in_active_set <- 0L
    for (i in seq_len(m)) {
      if (active_set[i]) {
        if (n_losses[i] == h) { # hit loss limit
          active_set[i] <- FALSE
          futility_set[i] <- TRUE
          stop_times[i] <- t
        } else {
          n_in_active_set <- n_in_active_set + 1L
          if (n_losses[i] > max_n_losses_active_set) max_n_losses_active_set <- n_losses[i]
        }
      }
    }

    # b. check for rejection of active set hypotheses
    threshold <- t + h - (h * m)/(n_in_active_set * alpha)
    if (max_n_losses_active_set <= threshold) {
      for (i in seq_len(m)) {
        if (active_set[i]) {
          active_set[i] <- FALSE
          rejected_set[i] <- TRUE
          stop_times[i] <- t
        }
      }
      break
    }
  }
  p <- h/(stop_times - n_losses + h)
  data.frame(hypothesis_idx = seq_len(m),
             p = p,
             significant = rejected_set,
             stop_time = stop_times,
             n_losses = n_losses) |> dplyr::arrange(p)
}

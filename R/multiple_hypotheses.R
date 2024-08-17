#' Run BC multiple testing procedure
#'
#' Runs the permutation test based on the anytime-valid BC p-value and BH correction.
#'
#' @param X an n by m matrix of treatments
#' @param Y an n by m matrix of responses
#' @param h the tuning parameter anytime-balid BC p-value
#'
#' @return a list containing (1) the 1-based indices of the hypotheses that are significant, (2) the number of permutations for each hypothesis, and (3) the p-value for each hypothesis
#' @export
#'
#' @examples
#' n <- 100L
#' m <- 5000L
#' X <- matrix(sample(x = c(0L, 1L), size = n * m, replace = TRUE), nrow = n, ncol = m)
#' under_null <- sample(c(TRUE, FALSE), size = m, prob = c(0.95, 0.05), replace = TRUE)
#' Y <- sapply(X = seq_len(m), FUN = function(i) {
#'  col <- X[,i]
#'  Y <- rnorm(n = n, mean = 2 + (if (under_null[i]) 0 else 2 * col), sd = 1)
#' })
#' h <- 15L; alpha <- 0.1
#' r <- sapply(seq(1, m), function(i) cor(X[,i], Y[,i]))
#' res <- run_bc_multiple_test(X, Y, alpha = alpha)
run_bc_multiple_test <- function(X, Y, h = 15L, alpha = 0.1) {
  set.seed(4)
  # initialize other variables
  n <- nrow(X) # sample size
  m <- ncol(X) # number of hypotheses
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
  for (i in seq_len(m)) original_test_statistics[i] <- compute_test_stat(X[,i], Y[,i])

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
        curr_test_stat <- compute_test_stat(x = X[perm_sample, i], y = Y[,i])
        # determine whether we have a loss; if so, increment n_losses
        loss_indicator <- as.integer(curr_test_stat >= original_test_statistics[i])
        n_losses[i] <- n_losses[i] + loss_indicator
      }
    }
    # a. check for futility
    for (i in seq_len(m)) {
      if (active_set[i]) {
        if (n_losses[i] == h) { # hit loss limit
          active_set[i] <- FALSE
          futility_set[i] <- TRUE
          stop_times[i] <- t
        }
      }
    }

    # b. apply bh
    p <- h/(t - n_losses + h)
    rejections <- p.adjust(p, "BH") < alpha
    for (i in seq_len(m)) {
      if (rejections[i]) {
        active_set[i] <- FALSE
        rejected_set[i] <- TRUE
        stop_times[i] <- t
      }
    }

    # c. determine whether to continue
    if (!any(active_set)) break
  }
  data.frame(hypothesis_idx = seq_len(m), p = p,
             significant = rejected_set, stop_time = stop_times,
             n_losses = n_losses) |>
    dplyr::arrange(p)
}

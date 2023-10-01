#' Michelle Cleary, s1979093, michellecleary
#' Add your own function definitions at the top of this file.

#' Poisson model confidence interval construction
#' (non-approximated ratio)
#'
#' @param y vector of observations
#' @param alpha percentage level of significance of interval
#'
#' @return (a, b), the confidence interval for lambda

conf_int <- function(y, alpha) {
  # Compute number of observations
  n <- length(y)
  # Compute ybar
  ybar <- mean(y)
  # Compute lambda^hat
  lambda_hat <- ybar
  # Compute the quantiles a and b
  a <- qnorm(alpha / 2)
  b <- qnorm(1 - alpha / 2)
  # Compute the confidence interval
  lower <- max(0, 1 / 2 * (a^2 / n + 2 * lambda_hat + a / n * sqrt(a^2 + 4 * n * lambda_hat)))
  upper <- 1 / 2 * (b^2 / n + 2 * lambda_hat + b / n * sqrt(b^2 + 4 * n * lambda_hat))
  interval <- c(lower, upper)
  interval
}


#' Poisson model confidence interval construction
#' (approximated ratio)
#'
#' @param y vector of observations
#' @param alpha percentage level of significance of interval
#'
#' @return (a, b), the confidence interval for lambda

approx_conf_int <- function(y, alpha) {
  # Compute number of observations
  n <- length(y)
  # Compute yabr
  ybar <- mean(y)
  # Compute lambda^hat
  lambda_hat <- ybar
  # Compute the quantiles a and b
  a <- qnorm(alpha / 2)
  b <- qnorm(1 - alpha / 2)
  # Compute the confidence interval
  lower <- max(0, lambda_hat - b * sqrt(lambda_hat / n))
  upper <- lambda_hat - a * sqrt(lambda_hat / n)
  interval <- c(lower, upper)
  interval
}


#' Log prior density
#'
#' Evaluate the logarithm of the joint prior density p(<U+03B8>)
#' for each parameter in <U+03B8>
#'
#' @param theta theta parameter vector
#' @param params vector of gamma parameters
#'
#'
log_prior_density <- function(theta, params) {
  # Compute the univariate log prior densities p(theta_i), i=1,....,4
  p_theta_1 <- dnorm(theta[1], mean = 0, sd = sqrt(params[1]), log = TRUE)
  p_theta_2 <- dnorm(theta[2], mean = 1, sd = sqrt(params[2]), log = TRUE)
  p_theta_3 <- dlogexp(theta[3], rate = params[3], log = TRUE)
  p_theta_4 <- dlogexp(theta[4], rate = params[4], log = TRUE)
  # Compute the joint prior densities p(theta)
  p_theta <- p_theta_1 + p_theta_2 + p_theta_3 + p_theta_4
  p_theta
}

#' Log-likelihood
#'
#' Evaluate the observation log-likelihood p(y|theta) for the
#' model yi~Normal[beta_1+beta_2*x_i, beta_3+beta_4*x_i^2)]
#'
#' @param theta theta parameter vector
#' @param x CAD weight for observation i
#' @param y actual weight for observation i
#'
#'

log_like <- function(theta, x, y) {
  # Compute the number of observations
  n <- length(x)
  # Initialise the log-likelihood p(y|theta) as 0
  p_y_theta <- 0
  # Iterate over each observation
  for (i in seq_len(n)) {
    # Compute the log-likelihood for observation i
    p_y_i_theta <- dnorm(y[i],
      mean = theta[1] + theta[2] * x[i],
      sd = sqrt(exp(theta[3]) + exp(theta[4]) * x[i]^2), log = TRUE
    )
    # Add this value to the log-likelihood p(y|theta)
    p_y_theta <- p_y_theta + p_y_i_theta
  }
  p_y_theta
}

#' Log posterior density
#'
#' Evaluate the logarithm of the posterior density p(<U+03B8>|y),
#' apart from some unevaluated normalisation constant
#'
#' @param theta theta parameter vector
#' @param x CAD weight for observation i
#' @param y actual weight for observation i
#' @param params vector of gamma parameters
#'
#'
log_posterior_density <- function(theta, x, y, params) {
  # Compute the log of the prior density p(theta)
  p_theta <- log_prior_density(theta, params)
  # Compute the log-likelihood p(y|theta)
  p_y_theta <- log_like(theta, x, y)
  # Compute the log of the posterior density p(theta|y)
  p_theta_y <- p_theta + p_y_theta
  p_theta_y
}

#' Importance sampling function
#'
#' @param N number of samples to generate
#' @param mu mean vector for the importance distribution
#' @param S covariance matrix
#'
#' @return A data frame with 5 columns containing the beta_i
#' samples, for i = 1,...,4, and the normalised
#' log-importance-weights
#'
do_importance <- function(N, mu, S, x, y, params) {
  # Simulate an importance sample for theta
  sample_theta <- rmvnorm(N, mean = mu, sigma = S)
  # Compute the log-weights
  log_weights <- log_prior_density(sample_theta, params) +
    log_like(sample_theta, x, y) -
    dmvnorm(sample_theta, mean = mu, sigma = S, log = TRUE)
  # Compute the normalisation constant
  normalisation_const <- log_sum_exp(log_weights)
  # Subtract the normalisation constant from the log-weights such that sum(exp(log_weights)) = 1
  log_weights <- log_weights - normalisation_const
  # Take the exponential of the sample values for theta_3 and theta_4 as we want an importance sample for beta
  sample_theta[, 3] <- exp(sample_theta[, 3])
  sample_theta[, 4] <- exp(sample_theta[, 4])
  # Create a data frame containing the beta samples and the normalised log-weights
  df <- data.frame(sample_theta, log_weights)
  # Set the column names
  colnames(df) <- c(paste("beta", seq_len(ncol(df) - 1), sep = ""), "log_weights")
  df
}

#' Credible interval construction
#'
#' @param x vector of importance samples
#' @param weights importance weights
#' @param prob intended coverage probability

make_CI <- function(x, prob, weights) {
  # Compute the significance level alpha
  alpha <- 1 - prob
  # Compute the credible interval
  CI <- wquantile(x, probs = c(alpha / 2, 1 - (alpha / 2)), weights = weights)
  data.frame(CI)
}

#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}


#' Interval coverage estimation
#'
#' @param approx_method a function object taking arguments x and alpha,
#'   returning a 2-element vector.
#' @param non_approx_method a function object taking arguments x and alpha,
#'   returning a 2-element vector.
#' @param N The number of simulation replications to use for the coverage estimate
#' @param alpha 1-alpha is the intended coverage probability
#' @param n The sample size
#' @param lambda The true lambda values for the Poisson(lambda) model
#'
#' @return A 2-element vector, each element a scalar between 0 and 1,
#' estimating the coverage probabilities for the `approx_method` and
#' `non_approx_method` intervals for the given parameters.


estimate_coverage <- function(approx_method,
                              non_approx_method,
                              N = 10000,
                              alpha = 0.1,
                              n = 2,
                              lambda = 3) {
  cover_approx <- 0
  cover_non_approx <- 0
  for (loop in seq_len(N)) {
    y <- rpois(n, lambda)
    ci_approx <- approx_method(y, alpha)
    ci_non_approx <- non_approx_method(y, alpha)
    cover_approx <- cover_approx +
      ((ci_approx[1] <= lambda) && (lambda <= ci_approx[2]))
    cover_non_approx <- cover_non_approx +
      ((ci_non_approx[1] <= lambda) && (lambda <= ci_non_approx[2]))
  }
  c(cover_approx / N, cover_non_approx / N)
}

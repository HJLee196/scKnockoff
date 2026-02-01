#' Estimate the number of spiked eigenvalues.
#'
#' @param data A numeric matrix of dimension n x p.
#' @param alpha A numeric scalar in (0, 0.5) specifying the trimming proportion
#'   for the eigenvalue index set used in fitting.
#' @param beta_q A numeric scalar in (0, 0.5) specifying the tail probability
#'   used to construct the empirical upper envelope of Monte Carlo null eigenvalues.
#' @param n_mc A positive integer specifying the number of Monte Carlo simulations.
#' @param verbose Logical scalar; if \code{TRUE}, print progress and timing
#'   messages during computation.
#'
#' @return A list with components:
#' \describe{
#'   \item{L_mat}{Monte Carlo eigenvalue matrix (n_mc x p), scaled by \code{1/n}.}
#'   \item{l_true}{Eigenvalues of \code{t(data) \%*\% data / n}.}
#'   \item{sigma2_est}{Estimated scaling factor (slope) used to match spectra.}
#'   \item{theta_est}{Estimated Gamma shape parameter.}
#'   \item{result_mat}{Matrix of estimated ranks evaluated at quantile
#'     levels from 0.05 to 0.95.}
#'   \item{alpha}{The input \code{alpha}.}
#'   \item{beta_q}{The input \code{beta_q}.}
#'   \item{estimated_rank}{Estimated rank at the quantile level \code{beta_q}.}
#' }
#'
#' @references
#' This method is based on Algorithm 2 and GetQT1 in:
#'
#' Ke, Z. T., Ma, Y., & Lin, X. (2023).
#' \emph{Estimation of the Number of Spiked Eigenvalues in a Covariance Matrix
#' by Bulk Eigenvalue Matching Analysis}.
#' Journal of the American Statistical Association, 118(541), 374--392.
#' \doi{10.1080/01621459.2021.1933497}
#'
#' The implementation is adapted from the BEMA reference code
#' (\url{https://github.com/ZhengTracyKe/BEMA})
#' and is otherwise identical, with only minor implementation-level
#' modifications for computational efficiency.
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' r <- 5
#' p <- 30
#'
#' X <- matrix(rnorm(n * r), n, r)
#' B <- matrix(rnorm(p * r), p, r)
#' err <- matrix(rnorm(n * p), n, p)
#'
#' XBt <- X %*% t(B) + err
#'
#' out <- SQM_Gamma(data = XBt, alpha = 0.2, beta_q = 0.1, n_mc = 30, verbose = TRUE)
#'
#' @export
SQM_Gamma = function(data, alpha = 0.2, beta_q = 0.1, n_mc = 500, verbose = FALSE) {
  # data
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("`data` must be a numeric matrix or data.frame.")
  }
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    stop("`data` must contain only numeric values.")
  }
  if (anyNA(data)) {
    stop("`data` must not contain NA values.")
  }
  # alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a numeric scalar in (0, 0.5).")
  }
  # beta_q
  if (!is.numeric(beta_q) || length(beta_q) != 1 || beta_q <= 0 || beta_q >= 0.5) {
    stop("`beta_q` must be a numeric scalar in (0, 0.5).")
  }
  # n_mc
  if (!is.numeric(n_mc) || length(n_mc) != 1 || is.na(n_mc) ||
      n_mc <= 0 || n_mc != as.integer(n_mc)) {
    stop("`n_mc` must be a positive integer.")
  }

  data <- as.matrix(data)

  n <- nrow(data)
  p <- ncol(data)

  S <- crossprod(data) / n


  l <- eigen(S, symmetric = TRUE, only.values = TRUE)$values

  # Computing Time: optimize()
  start_time_opt <- Sys.time()
  o <- stats::optimize(
    f = loss_internal,
    interval = c(0.1, 50),
    l = l, alpha = alpha, n = n, p = p
  )
  a <- o$minimum
  end_time_opt <- Sys.time()

  opt_total_seconds <- as.numeric(difftime(end_time_opt, start_time_opt, units = "secs"))
  opt_minutes <- floor(opt_total_seconds / 60)
  opt_seconds <- round(opt_total_seconds %% 60)
  if (verbose){
    message(sprintf("Spent Computation Time (Optimize): %02d:%02d", opt_minutes, opt_seconds))
  }

  L <- matrix(0, nrow = n_mc, ncol = p)

  # Computing Time: Monte Carlo
  start_time_nsim <- Sys.time()

  for (i in seq_len(n_mc)) {

    sd_candidate <- matrix(sqrt(stats::rgamma(p, a, a)), n, p, byrow = TRUE)
    x1 <- matrix(stats::rnorm(n = n * p), n, p) * sd_candidate

    if (p <= n) {
      l1 <- eigen(crossprod(x1), symmetric = TRUE, only.values = TRUE)$values
      L[i, ] = l1 / n
    } else {
      l1 <- eigen(tcrossprod(x1), symmetric = TRUE, only.values = TRUE)$values
      L[i, 1:n] = l1 / n
    }

    if (verbose && ((i %% 10) == 0)) {
      end_time_nsim <- Sys.time()
      nsim_total_seconds <- as.numeric(difftime(end_time_nsim, start_time_nsim, units = "secs"))
      nsim_minutes <- floor(nsim_total_seconds / 60)
      nsim_seconds <- round(nsim_total_seconds %% 60)
      message(sprintf("Iteration %d", i))
      message(sprintf("Spent Computation Time: %02d:%02d", nsim_minutes, nsim_seconds))
      start_time_nsim <- Sys.time()
    }
  }

  l_mean <- numeric(p)
  l_q_hi <- numeric(p)
  l_med  <- numeric(p)

  for (j in seq_len(p)) {
    l_mean[j] <- mean(L[, j])
    l_q_hi[j] <- as.numeric(stats::quantile(L[, j], 1 - beta_q))
    l_med[j]  <- as.numeric(stats::quantile(L[, j], 0.5))
  }

  k <- floor(min(p, n) * alpha) : floor(min(p, n) * (1 - alpha))
  s1 <- as.numeric(stats::lm(l[k] ~ l_mean[k] - 1)$coef[[1]])
  estimated_rank <- sum(l > max(l_q_hi * s1))

  if (verbose) {
    message(sprintf("Beta: %s", beta_q))
    message(sprintf("Estimated Rank: %d", estimated_rank))
  }

  L_max <- apply(L, 1, max)
  quantile_seq <- seq(from = 0.05, to = 0.95, by = 0.05)
  Ls_quantile <- stats::quantile(L_max, probs = quantile_seq) * s1

  result_vec <- sapply(
    Ls_quantile,
    FUN = function(Ls_temp, l_temp) sum(l_temp > Ls_temp),
    l_temp = l
  )
  result_mat <- cbind(Rank = as.integer(result_vec), Quantile = names(result_vec))

  result_list <-
    list(
      L_mat = L,
      l_true = l,
      sigma2_est = s1,
      theta_est = a,
      result_mat = result_mat,
      alpha = alpha,
      beta_q = beta_q,
      estimated_rank = estimated_rank
    )

  return(result_list)
}

# ---- Internal helper (not exported) -----------------------------------------
#' @noRd
loss_internal = function(proposal, l, alpha, n, p) {

  L <- matrix(0, nrow = 10, ncol = p)

  for (i in 1:10) {

    sd_candidate <- matrix(sqrt(stats::rgamma(p, proposal, proposal)), n, p, byrow = TRUE)
    x1 <- matrix(stats::rnorm(n = n * p), n, p) * sd_candidate

    if (p <= n) {
      l1 <- eigen(crossprod(x1), symmetric = TRUE, only.values = TRUE)$values
      L[i, ] <- l1 / n
    } else {
      l1 <- eigen(tcrossprod(x1), symmetric = TRUE, only.values = TRUE)$values
      L[i, 1:n] <- l1 / n
    }
  }

  l1 <- colMeans(L)

  k <- floor(min(p, n) * alpha) : floor(min(p, n) * (1 - alpha))
  s1 <- as.numeric(stats::lm(l[k] ~ l1[k] - 1)$coef[[1]])
  l1 <- s1 * l1

  return(sum((l1 - l)[k]^2))
}

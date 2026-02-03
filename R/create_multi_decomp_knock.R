
# Input:
# np_data.matrix.imp: The imputed data.
# Xl: Random matrix of latent factors.
# Bl: Matrix of deterministic factor loadings.
# err: Estimated variances (of the error matrix).
# PC: Number of latent factors, excluding additional covariates.
# m: Number of knockoff copies.

# Output:
# Multi-decomp knockoffs, which are not rescaled yet.

#' Multiple knockoff construction using matrix decomposition
#'
#' @param np_data.matrix.imp A numeric matrix of dimension n x p containing the
#'   imputed data.
#' @param Xl A numeric matrix of dimension n x k containing latent factors, where
#'   k = PC plus the number of additional covariates.
#' @param Bl A numeric matrix of dimension k x p containing factor loadings.
#' @param err A numeric vector of length p containing estimated error variances.
#' @param PC A positive integer specifying the number of latent factors,
#'   excluding additional covariates.
#' @param m A positive integer specifying the number of knockoffs to be constructed.
#'
#' @return
#' A list of length \code{m}. Each element is a numeric matrix of dimension
#' n x p containing one knockoff copy constructed using the matrix
#' decomposition-based covariance estimate and the multiple knockoff procedure.
#' The returned knockoffs are not rescaled; rescaling should be performed
#' separately (e.g., using \code{rescale_knockoff}).
#'
#' @family create
#'
#' @details A multiple knockoff constructor that takes advantage of a matrix decomposition,
#' creating more efficient knockoffs for high-dimensional datasets.
#'
#' @references
#' Roquero Gimenez, J., and Zou, J. (2019).
#' \emph{Improving the Stability of the Knockoff Procedure: Multiple Simultaneous
#' Knockoffs and Entropy Maximization}.
#' In \emph{Proceedings of the 22nd International Conference on Artificial Intelligence
#' and Statistics}, pp. 2184--2192.
#'
#' He, Z., Chu, B., Yang, J., Gu, J., Chen, Z., Liu, L., Morrison, T.,
#' Belloy, M. E., Qi, X., Hejazi, N., Mathur, M., Le Guen, Y., Tang, H.,
#' Hastie, T., Ionita-Laza, I., Sabatti, C., and Candes, E. (2024).
#' \emph{Beyond Guilty by Association at Scale: Searching for Causal Variants
#' on the Basis of Genome-Wide Summary Statistics}.
#' bioRxiv preprint.
#' \doi{10.1101/2024.02.28.582621}
#'
#' @examples
#' set.seed(2024)
#' # Create dataset
#' Xl <- matrix(rnorm(1000),nrow = 100)
#' Bl <- matrix(rnorm(800),ncol=80)
#' err <- abs(rnorm(80,sd=3))
#'
#' E <- matrix(NA,nrow = 100,ncol = 80)
#' for(i in 1:80){
#'   E[,i] <- rnorm(100, sd = err[i])
#' }
#'
#' np_data <- Xl%*%Bl + E
#'
#' rownames(np_data) <- paste0("r", 1:100)
#' colnames(np_data) <- paste0("c", 1:80)
#'
#' # Create missing data
#' miss <- sample(1:prod(dim(np_data)),floor(prod(dim(np_data))*0.6))
#' np_data[miss] <- 0
#'
#' np_data_exp <- matrix(1,nrow=100,ncol=80)
#' np_data_exp[miss] <- 0
#'
#' np_data_exp.count <- apply(np_data_exp,2,sum)
#'
#' # Center target matrix, but only center the expressed parts.
#' np_data.avg <- (colSums(np_data, na.rm = TRUE))/(np_data_exp.count-1)
#' np_data_centered <- sweep(np_data,2,np_data.avg,FUN = "-")*(np_data_exp)
#'
#' # Impute missing values
#' np_data_imp <- sc_softImpute(np_data_centered,np_data_exp,Xl[,1:3],PC=min(np_data_exp.count)-5)
#'
#' # multidecomp knockoff construction
#' m <- 5
#' multi_decomp_knock <- create_multi_decomp_knock(np_data_imp$X_imp,
#'                                                 np_data_imp$Xl,
#'                                                 np_data_imp$Bl,
#'                                                 np_data_imp$err,
#'                                                 np_data_imp$PC,
#'                                                 m)
#'
#' # Rescale each of the knockoff copies
#' multi_decomp_knock <- rescale_knockoff(multi_decomp_knock, np_data_exp, np_data.avg)
#'
#' @export
create_multi_decomp_knock <- function(np_data.matrix.imp,Xl,Bl,err,PC,m){
  # np_data.matrix.imp
  if (!is.matrix(np_data.matrix.imp) && !is.data.frame(np_data.matrix.imp)) {
    stop("`np_data.matrix.imp` must be a numeric matrix or data.frame.", call. = FALSE)
  }

  np_data.matrix.imp <- as.matrix(np_data.matrix.imp)
  if (!is.numeric(np_data.matrix.imp)) {
    stop("`np_data.matrix.imp` must contain only numeric values.", call. = FALSE)
  }
  if (anyNA(np_data.matrix.imp)) {
    stop("`np_data.matrix.imp` must not contain NA values.", call. = FALSE)
  }

  n <- nrow(np_data.matrix.imp)
  p <- ncol(np_data.matrix.imp)

  # Xl
  if (!is.matrix(Xl) && !is.data.frame(Xl)) {
    stop("`Xl` must be a numeric matrix or data.frame.", call. = FALSE)
  }
  Xl <- as.matrix(Xl)
  if (!is.numeric(Xl)) {
    stop("`Xl` must contain only numeric values.", call. = FALSE)
  }
  if (nrow(Xl) != n) {
    stop("`Xl` must have the same number of rows as `np_data.matrix.imp`.", call. = FALSE)
  }

  # Bl
  if (!is.matrix(Bl) && !is.data.frame(Bl)) {
    stop("`Bl` must be a numeric matrix or data.frame.", call. = FALSE)
  }
  Bl <- as.matrix(Bl)
  if (!is.numeric(Bl)) {
    stop("`Bl` must contain only numeric values.", call. = FALSE)
  }
  if (ncol(Bl) != p) {
    stop("`Bl` must have ncol(Bl) = p = ncol(np_data.matrix.imp).", call. = FALSE)
  }

  # err
  if (!is.numeric(err) || length(err) != p || anyNA(err)) {
    stop("`err` must be a numeric vector of length p = ncol(np_data.matrix.imp) with no NA.", call. = FALSE)
  }
  if (any(err < 0)) {
    stop("`err` must be non-negative (variances).", call. = FALSE)
  }

  # PC
  if (!is.numeric(PC) || length(PC) != 1 || is.na(PC) || PC <= 0 || PC != as.integer(PC)) {
    stop("`PC` must be a positive integer.", call. = FALSE)
  }
  if (ncol(Xl) < PC) {
    stop("`Xl` must have at least `PC` columns.", call. = FALSE)
  }
  if (nrow(Bl) != ncol(Xl)) {
    stop("`nrow(Bl)` must equal `ncol(Xl)`.", call. = FALSE)
  }

  # m
  if (!is.numeric(m) || length(m) != 1 || is.na(m) || m <= 0 || m != as.integer(m)) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }
  if (m < 2) {
    stop(
      "`m` must be at least 2. For a single knockoff (m = 1), use `create_decomp_knock()` instead.",
      call. = FALSE
    )
  }

  is_no_additional_covariates <- FALSE
  if (dim(Xl)[2] == PC){ ## Without conditioning on X ## "dim(Xl)[2] == PC" indicates that there is no additional covariate.
    is_no_additional_covariates <- TRUE
  }

  W <- diag(err)
  if (is_no_additional_covariates){
    sigma.Xl <- stats::cov(Xl)
    Bl <- t(Bl)
    sigma.Y <- W + (Bl)%*%sigma.Xl%*%t(Bl)
    np_data.matrix.imp0 <- np_data.matrix.imp
  } else { ## Conditioning on X
    Q_index <- 1:(ncol(Xl)-PC)
    A_index <- (max(Q_index)+1):(max(Q_index)+PC)

    X_mat <- Xl[,Q_index]
    A_hat <- Xl[,A_index]

    B0_hat <- t(Bl)[,Q_index]
    B1_hat <- t(Bl)[,A_index]

    Gamma_hat <- t(solve(t(X_mat) %*% X_mat) %*% t(X_mat) %*% A_hat)
    U_hat <- A_hat - (X_mat %*% t(Gamma_hat))
    SigmaA_hat <- (t(U_hat) %*% U_hat)/nrow(U_hat)

    np_data.matrix.imp0 <- np_data.matrix.imp - (X_mat %*% t(B0_hat)) - (X_mat %*% t(Gamma_hat) %*% t(B1_hat))
    sigma.Y <- W + (B1_hat %*% SigmaA_hat %*% t(B1_hat))
  }

  # the number depends on m, (m+1)/m. E.g. when m = 5, it's 1.2
  # (-0.05 to be conservative) taking a number slightly smaller than
  # (m+1)/m in order to avoid possible computational issues.
  S_multi_decomp <- diag(((m+1)/m-0.05)*(err))
  sigma.Y_inv <- solve(sigma.Y)

  C_Cov <- (2*S_multi_decomp) - (S_multi_decomp %*% sigma.Y_inv %*% S_multi_decomp)
  V1_Cov <- C_Cov - ((m - 1)/m)*S_multi_decomp
  V1 <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = V1_Cov)

  V2_array <- array(dim = c(n, p, m))
  for (i_m in 1:m)
  {
    V2_array[,,i_m] <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = S_multi_decomp)
  }

  V2bar <- apply(X = V2_array, MARGIN = c(1,2), FUN = mean)
  V2_array <- sweep(V2_array, c(1, 2), V2bar, "-")
  V2_array <- sweep(V2_array, c(1, 2), V1, "+")

  np_data.knockoff.sav0_list <- list()
  for (i_m in 1:m)
  {
    np_data.knockoff <-(np_data.matrix.imp0 %*% t(diag(1,p,p) - (S_multi_decomp %*% sigma.Y_inv)) ) + V2_array[,,i_m]

    colnames(np_data.knockoff) <- colnames(np_data.matrix.imp)
    rownames(np_data.knockoff) <- rownames(np_data.matrix.imp)

    ## scale back the knockoff, but not fully rescaled
    if (is_no_additional_covariates){
      np_data.knockoff.sav0 <- np_data.knockoff
    } else {
      np_data.knockoff.sav0 <- np_data.knockoff + (X_mat %*% t(B0_hat)) + (X_mat %*% t(Gamma_hat) %*% t(B1_hat))
    }

    np_data.knockoff.sav0_list[[i_m]] <- np_data.knockoff.sav0
  }

  return(np_data.knockoff.sav0_list)
}

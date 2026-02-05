
# Input:
#  np_data.matrix.imp: The imputed data.
#  Xl: Random matrix of latent factors.
#  Bl: Matrix of deterministic factor loadings.
#  err: Estimated variances (of the error matrix).
#  m: Number of knockoff copies.

# Output:
# Multi-lr knockoffs, which are not rescaled yet.

#' Multiple knockoff construction using low rank structure
#'
#' @param np_data.matrix.imp A numeric matrix of dimension n x p containing the
#'   imputed data.
#' @param Xl A numeric matrix of dimension n x k containing latent factors.
#' @param Bl A numeric matrix of dimension k x p containing factor loadings.
#' @param err A numeric vector of length p containing estimated error variances.
#' @param m A positive integer specifying the number of knockoffs to be constructed.
#'
#' @return
#' A list of length \code{m}. Each element is a numeric matrix of dimension
#' n x p containing one knockoff copy constructed using the low-rank structure.
#' The knockoffs are constructed as \eqn{X_l B_l + E}, where \eqn{E} is a noise
#' matrix with independent columns and variances specified by \code{err}.
#' The returned matrix is not rescaled; rescaling should be performed
#' separately (e.g., using \code{rescale_knockoff}).
#'
#' @family create
#'
#' @details Multi-LR knockoffs that can be constructed based on the results from \code{sc_softImpute}.
#'
#' @references
#' Fan, Y., Lv, J., Sharifvaghefi, M., and Uematsu, Y. (2020).
#' \emph{IPAD: Stable Interpretable Forecasting with Knockoffs Inference}.
#' \emph{Journal of the American Statistical Association}, 115(532), 1822--1834.
#'
#' @examples
#' set.seed(2024)
#' # Create dataset
#' X <- matrix(rnorm(1000),nrow = 100)
#' Bl <- matrix(rnorm(800),ncol=80)
#' err <- abs(rnorm(80,sd=3))
#'
#' E <- matrix(NA,nrow = 100,ncol = 80)
#' for(i in 1:80){
#'   E[,i] = rnorm(100, sd = err[i])
#' }
#'
#' np_data <- X%*%Bl + E
#'
#' rownames(np_data) <- paste0("r", 1:100)
#' colnames(np_data) <- paste0("c", 1:80)
#'
#' # Create missing data
#' miss <- sample(1:prod(dim(np_data)), floor(prod(dim(np_data))*0.6))
#' np_data[miss] <- 0
#'
#' np_data_exp <- matrix(1, nrow=100, ncol=80)
#' np_data_exp[miss] <- 0
#'
#' np_data_exp.count <- apply(np_data_exp,2,sum)
#'
#' # Center target matrix, but only center the expressed parts.
#' np_data.avg <- (colSums(np_data, na.rm = TRUE))/(np_data_exp.count-1)
#' np_data_centered <- sweep(np_data,2,np_data.avg,FUN = "-")*(np_data_exp)
#'
#' # Impute missing values
#' np_data_imp <- sc_softImpute(np_data_centered, np_data_exp, X[,1:3], PC = min(np_data_exp.count)-5)
#'
#' # Multi-lr Knockoff construction
#' m <- 5
#' multi_lr_knock <- create_multi_lr_knock(np_data_imp$X_imp,
#'                                         np_data_imp$Xl,
#'                                         np_data_imp$Bl,
#'                                         np_data_imp$err,
#'                                         m)
#'
#' multi_lr_knock <- rescale_knockoff(multi_lr_knock, np_data_exp, np_data.avg)
#'
#' @export
create_multi_lr_knock <- function(np_data.matrix.imp, Xl, Bl, err, m){
  Xl <- as.matrix(Xl)
  Bl <- as.matrix(Bl)

  if (!is.numeric(Xl) || !is.numeric(Bl)) {
    stop("`Xl` and `Bl` must be numeric matrices.", call. = FALSE)
  }
  if (!is.numeric(err) || anyNA(err) || any(err < 0)) {
    stop("`err` must be a non-negative numeric vector with no NA.", call. = FALSE)
  }
  if (ncol(Xl) != nrow(Bl)) {
    stop("`ncol(Xl)` must equal `nrow(Bl)`.", call. = FALSE)
  }
  if (length(err) != ncol(Bl)) {
    stop("`length(err)` must equal `ncol(Bl)`.", call. = FALSE)
  }

  # m
  if (!is.numeric(m) || length(m) != 1 || is.na(m) || m <= 0 || m != as.integer(m)) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }
  if (m < 2) {
    stop(
      "`m` must be at least 2. For a single knockoff (m = 1), use `create_lr_knock()` instead.",
      call. = FALSE
    )
  }

  np_data.knockoff.sav0_list <- list()

  for (i_m in 1:m)
  {
    E <- matrix(NA, nrow = dim(Xl)[1], ncol = dim(Bl)[2])

    for(i in 1:dim(Bl)[2]){
      E[,i] = stats::rnorm(dim(Xl)[1], mean = 0,sd = sqrt(err[i]))
    }

    np_data.knockoff <- Xl%*%(Bl) + E

    colnames(np_data.knockoff) <- colnames(np_data.matrix.imp)
    rownames(np_data.knockoff) <- rownames(np_data.matrix.imp)

    np_data.knockoff.sav0_list[[i_m]] = np_data.knockoff
  }

  return(np_data.knockoff.sav0_list)
}

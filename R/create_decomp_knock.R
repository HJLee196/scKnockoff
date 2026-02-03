
# Input:
# np_data.matrix.imp: The imputed data.
# Xl: Random matrix of latent factors.
# Bl: Matrix of deterministic factor loadings.
# err: Estimated variances (of the error matrix).
# PC: Number of latent factors, excluding additional covariates.

# Output:
# Decomp knockoffs, which are not rescaled yet.

#' Knockoff construction using matrix decomposition
#'
#' @param np_data.matrix.imp A numeric matrix of dimension n x p containing the
#'   imputed data.
#' @param Xl A numeric matrix of dimension n x k containing latent factors, where
#'   k = PC plus the number of additional covariates.
#' @param Bl A numeric matrix of dimension k x p containing factor loadings.
#' @param err A numeric vector of length p containing estimated error variances.
#' @param PC A positive integer specifying the number of latent factors,
#'   excluding additional covariates.
#' @param decomp Logical scalar; if \code{TRUE}, generate knockoffs without solving an approximate
#'  semidefinite program (ASDP). If \code{FALSE}, generate knockoffs with ASDP.
#'
#' @return
#' If successful, a numeric matrix of dimension \eqn{n \times p} containing
#' one set of knockoff variables constructed using either a matrix
#' decomposition-based covariance estimate or an ASDP-based approach (see
#' Details). The returned matrix is not rescaled; rescaling should be
#' performed separately (e.g., using \code{rescale_knockoff}). If knockoff
#' construction fails, the function returns \code{NULL} and issues a
#' warning.
#'
#' @family create
#'
#' @details
#' Constructs knockoff variables using second-order Gaussian knockoff
#' generation. When \code{decomp = TRUE}, the method uses a matrix
#' decomposition-based covariance estimator to reduce computational and
#' memory costs in high-dimensional settings compared to full
#' covariance-based approaches. When \code{decomp = FALSE}, the method constructs knockoffs using an
#' approximate semidefinite programming (ASDP) approach.
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
#' # Decomp Knockoff construction
#' decomp_knock <- create_decomp_knock(np_data_imp$X_imp,
#'                                     np_data_imp$Xl,
#'                                     np_data_imp$Bl,
#'                                     np_data_imp$err,
#'                                     np_data_imp$PC)
#'
#' decomp_knock <- rescale_knockoff(decomp_knock, np_data_exp, np_data.avg)
#'
#' @export
create_decomp_knock <- function(np_data.matrix.imp,Xl,Bl,err,PC,decomp=TRUE){
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

  # decomp
  if (!is.logical(decomp) || length(decomp) != 1 || is.na(decomp)) {
    stop("`decomp` must be a single logical value.", call. = FALSE)
  }

  W <- diag(err)
  if (dim(Xl)[2] == PC) { ## Without conditioning on X ## "dim(Xl)[2] == PC" indicates that there is no additional covariate.
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

  ##generate second-order gaussian knockoff
  skip_to_next <- FALSE

  if(decomp){
    # diag_s can be as large as 2*W. It is set at 1.95*W to avoid
    # potential computational errors.
    tryCatch(np_data.knockoff0 <- knockoff::create.gaussian(np_data.matrix.imp0,
                                                            mu = rep(0, length.out = p),
                                                            Sigma = sigma.Y,
                                                            method = "asdp",
                                                            diag_s = 1.95*W),
             error = function(e) {skip_to_next <<- TRUE})}
  else{
    tryCatch(np_data.knockoff0 <- knockoff::create.gaussian(np_data.matrix.imp0,
                                                            mu = rep(0, length.out = p),
                                                            Sigma = sigma.Y,
                                                            method = "asdp"),
             error = function(e) {skip_to_next <<- TRUE})}

  if(skip_to_next) {
    warning("Knockoff construction failed.", call. = FALSE)
    return(NULL) }

  if (dim(Xl)[2] != PC) {
    np_data.knockoff0 <-
      np_data.knockoff0 + (X_mat %*% t(B0_hat)) + (X_mat %*% t(Gamma_hat) %*% t(B1_hat))
  }

  return(np_data.knockoff0)
}

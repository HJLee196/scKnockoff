
#' softImpute with additional covariates
#'
#' @param np_data.matrix A centered numeric matrix of dimension n x p containing
#'   the original variables.
#' @param np_data.matrix.exp A numeric matrix of dimension n x p indicating
#'   observed entries, with 1 for observed values and 0 for missing values.
#' @param X A numeric matrix of additional covariates. If \code{NULL}, the
#'   function reduces to a standard \code{softImpute}-based imputation.
#' @param PC A positive integer specifying the number of latent factors.
#' @param lambda A penalty parameter used for imputation; by default,
#'   \code{softImpute::lambda0()} is used.
#' @param lam_ratio A ratio applied to \code{lambda0}; only used when
#'   \code{lambda = NULL} and is set to 0.2 by default.
#' @param max_it A positive integer specifying the maximum number of iterations;
#'   default is 100.
#'
#' @return
#' A list containing the imputed data and estimated model components obtained
#' from the adjusted softImpute algorithm with additional covariates.
#' The returned list contains:
#' \describe{
#'  \item{X_imp}{A numeric matrix of dimension n x p containing the imputed data
#'     with error added.}
#'  \item{Xl}{A numeric matrix of dimension n x k containing estimated latent
#'     factors, where \code{k = ncol(X) + PC}.}
#'  \item{Bl}{A numeric matrix of dimension k x p containing estimated factor
#'     loadings, where \code{k = ncol(X) + PC}.}
#'  \item{err}{A numeric vector of length p containing estimated error variances.}
#'  \item{PC}{The number of latent factors used in the algorithm.}
#' }
#'
#' @details
#' This function implements an adjusted softImpute algorithm that allows the
#' inclusion of additional covariates. Each variable must have more than
#' \code{ncol(X) + PC + 1} non-missing observations.
#' \code{ncol(X)} is treated as 0 when \code{X} is \code{NULL}.
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
#' np_data_imp <- sc_softImpute(np_data_centered,np_data_exp,X[,1:3],PC=min(np_data_exp.count)-5)
#'
#' @export
sc_softImpute <- function(np_data.matrix,
                          np_data.matrix.exp,
                          X = NULL,
                          PC,
                          lambda = NULL,
                          lam_ratio = 0.2,
                          max_it = 100){

  if (!is.matrix(np_data.matrix) &&
      !is.data.frame(np_data.matrix) &&
      !inherits(np_data.matrix, "Matrix")) {
    stop("`np_data.matrix` must be a matrix-like object.", call. = FALSE)
  }
  np_data.matrix <- as.matrix(np_data.matrix)
  if (!is.numeric(np_data.matrix)) {
    stop("`np_data.matrix` must contain only numeric values.", call. = FALSE)
  }
  if (anyNA(np_data.matrix)) {
    stop("`np_data.matrix` must not contain NA values.", call. = FALSE)
  }

  n <- nrow(np_data.matrix)
  p <- ncol(np_data.matrix)

  # np_data.matrix.exp
  if (!is.matrix(np_data.matrix.exp) &&
      !is.data.frame(np_data.matrix.exp) &&
      !inherits(np_data.matrix.exp, "Matrix")) {
    stop("`np_data.matrix.exp` must be a matrix-like object.",
         call. = FALSE)
  }
  np_data.matrix.exp <- as.matrix(np_data.matrix.exp)
  if (!is.numeric(np_data.matrix.exp)) {
    stop("`np_data.matrix.exp` must contain only numeric values.", call. = FALSE)
  }
  if (anyNA(np_data.matrix.exp)) {
    stop("`np_data.matrix.exp` must not contain NA values.", call. = FALSE)
  }
  if (!all(dim(np_data.matrix.exp) == c(n, p))) {
    stop("`np_data.matrix.exp` must have the same dimensions as `np_data.matrix`.", call. = FALSE)
  }
  if (!all(np_data.matrix.exp %in% c(0, 1))) {
    stop("`np_data.matrix.exp` must be a 0/1 matrix indicating observed entries.", call. = FALSE)
  }

  # PC
  if (!is.numeric(PC) || length(PC) != 1L || is.na(PC) ||
      PC <= 0 || PC != as.integer(PC)) {
    stop("`PC` must be a positive integer.", call. = FALSE)
  }
  if (PC > min(n, p)) {
    stop("`PC` must be less than or equal to min(nrow(np_data.matrix), ncol(np_data.matrix)).",
         call. = FALSE)
  }

  # X (additional covariates)

  np_data.matrix.exp.count <- Matrix::colSums(np_data.matrix.exp)

  if (!is.null(X)) {
    if (!is.matrix(X) && !is.data.frame(X)) {
      stop("`X` must be NULL or a numeric matrix/data.frame.", call. = FALSE)
    }
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("`X` must contain only numeric values.", call. = FALSE)
    }
    if (anyNA(X)) {
      stop("`X` must not contain NA values.", call. = FALSE)
    }
    if (nrow(X) != n) {
      stop("`X` must have the same number of rows as `np_data.matrix`.", call. = FALSE)
    }
    k_num <- ncol(X) # number of additional covariates
  } else {
    k_num <- 0
  }

  if (min(np_data.matrix.exp.count) - PC - k_num - 1 <= 0) {
    stop(
      sprintf("Each variable must have more than PC + %d + 1 non-missing observations.", k_num),
      call. = FALSE
    )
  }

  # lam_ratio
  if (!is.numeric(lam_ratio) || length(lam_ratio) != 1L ||
      is.na(lam_ratio) || lam_ratio <= 0) {
    stop("`lam_ratio` must be a positive numeric scalar.", call. = FALSE)
  }

  # lambda
  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda <= 0) {
      stop("`lambda` must be NULL or a positive numeric scalar.", call. = FALSE)
    }
  } else {
    lambda <- softImpute::lambda0(np_data.matrix) * lam_ratio
  }

  # max_it
  if (!is.numeric(max_it) || length(max_it) != 1L ||
      is.na(max_it) || max_it <= 0 || max_it != as.integer(max_it)) {
    stop("`max_it` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(X)){
    # if (requireNamespace('Seurat', quietly=T)){
    #   X_0 <- Seurat::RunPCA(t(np_data.matrix), assay='RNA', slot = "data", npcs = PC, verbose = F)@cell.embeddings
    # } else {
    X_0 <- stats::prcomp(np_data.matrix, center = FALSE, rank. = PC)$x
    # }

    X_0 <- cbind(X,X_0)

    # Initial Values
    X_0 <- sweep(X_0,2,colMeans(X_0),FUN = "-")
    B_0 <- solve(t(X_0)%*%X_0 + lambda*diag(ncol = ncol(X_0), nrow = ncol(X_0)))%*%t(X_0)%*%np_data.matrix

    np_data.matrix_l <- lm_impute(Yl_imp = np_data.matrix,
                                  Yl = np_data.matrix,
                                  Yl_exp = np_data.matrix.exp,
                                  Bl = B_0,
                                  Xl = X_0,
                                  PC = PC,
                                  q_ncol = k_num,
                                  lambda = lambda)

    np_data.matrix_l <- lm_impute(Yl_imp = np_data.matrix_l[[1]],
                                  Yl = np_data.matrix,
                                  Yl_exp = np_data.matrix.exp,
                                  Bl = np_data.matrix_l[[3]],
                                  Xl = np_data.matrix_l[[2]],
                                  PC = PC,
                                  q_ncol = k_num,
                                  lambda = lambda)

    # stop_crit_compare <- (sum((np_data.matrix_l[[1]] - np_data.matrix_l[[2]] %*% t(np_data.matrix_l[[3]]))^2*
    #                     np_data.matrix.exp))
    stop_crit_compare <- sum((np_data.matrix_l[[1]] - np_data.matrix)^2)

    stop_crit <- TRUE
    i <- 1
    while (stop_crit) {
      np_data.matrix_l <- lm_impute(Yl_imp = np_data.matrix_l[[1]],
                                    Yl = np_data.matrix,
                                    Yl_exp = np_data.matrix.exp,
                                    Bl = np_data.matrix_l[[3]],
                                    Xl = np_data.matrix_l[[2]],
                                    PC = PC,
                                    q_ncol = k_num,
                                    lambda = lambda)

      # loop_compare <- sum((np_data.matrix_l[[1]] - np_data.matrix_l[[2]] %*% (np_data.matrix_l[[3]]))^2*
      #                       np_data.matrix.exp)
      loop_compare <- sum((np_data.matrix_l[[1]] - np_data.matrix)^2)

      # Iterate until the reconstruction error drops below 10% of the initial error
      stop_crit <- (loop_compare > stop_crit_compare/10) && (i < max_it)

      i <- i + 1
    }

    np_data.matrix <- np_data.matrix_l[[1]]
  } else {
    warning("X is NULL; using softImpute instead.", call. = FALSE)
    sI_np_data.matrix <- np_data.matrix
    sI_np_data.matrix[np_data.matrix.exp == 0] <- NA_real_

    result_softImpute <- softImpute::softImpute(x = sI_np_data.matrix,
                                                rank.max = PC,
                                                lambda = lambda,
                                                maxit = max_it)

    np_data.matrix <- softImpute::complete(sI_np_data.matrix, object = result_softImpute)
    np_data.matrix_l <- list(np_data.matrix,
                             result_softImpute$u %*% diag(sqrt(result_softImpute$d)),
                             diag(sqrt(result_softImpute$d)) %*% t(result_softImpute$v))

    # np_data.matrix.imp$X_imp[np_data.matrix.exp] %>% head()
    # np_data.matrix[np_data.matrix.exp] %>% head()
    # dim( result_softImpute$u)
    # dim( result_softImpute$d)
    # complete_mat = np_data.matrix.imp$Xl %*% np_data.matrix.imp$Bl
    # np_data.matrix.imp$X_imp[!np_data.matrix.exp] %>% head()
    # complete_mat[!np_data.matrix.exp] %>% head()
  }

  # Add noise to the imputed values
  err <- Matrix::colSums((np_data.matrix - np_data.matrix_l[[2]] %*% (np_data.matrix_l[[3]]))^2)
  err <- err/(np_data.matrix.exp.count-PC-k_num-1)

  # n <- dim(np_data.matrix)[1] # defined above
  for(i in 1:length(err)){
    np_data.matrix[,i] <-  np_data.matrix[,i] + stats::rnorm(n,0,sqrt(err[i]))*(1-np_data.matrix.exp[,i])
  }

  return(list(
    X_imp = np_data.matrix,
    Xl = np_data.matrix_l[[2]],
    Bl = np_data.matrix_l[[3]],
    err = err,
    PC = PC))
}

# ---- Internal helper (not exported) -----------------------------------------
# Internal function used in sc_softImpute.
#' @noRd
lm_impute <- function(Yl_imp,
                      Xl,
                      Bl,
                      Yl,
                      Yl_exp,
                      PC,
                      q_ncol,
                      lambda){

  ## Fit linear model
  n <- nrow(Yl_imp)
  p <- ncol(Yl_imp)
  Q_index = 1:q_ncol
  A_index = (1:PC) + q_ncol

  # Fix Bl
  B_q <- matrix(Bl[Q_index,], nrow = q_ncol)
  B_A <- matrix(Bl[A_index,], nrow = PC)
  Q <- matrix(Xl[,Q_index], ncol = q_ncol)

  Y_hat <- Xl%*%Bl
  Yl_imp <- (Yl - Y_hat)*Yl_exp + Y_hat - Q%*%(B_q)
  I <- diag(1,PC,PC)
  Al <- solve((B_A)%*%t(B_A) + lambda*I)%*%(B_A)%*%t(Yl_imp)
  Al <- t(Al)

  # Fix Al
  Xl[,A_index] <- Al

  Y_hat <- Xl%*%Bl
  Yl_imp <- (Yl - Y_hat)*Yl_exp + Y_hat
  I <- diag(1,PC+q_ncol,PC+q_ncol)
  Bl <- solve(t(Xl)%*%Xl + lambda*I)%*%t(Xl)%*%Yl_imp

  Y_hat <- Xl%*%Bl
  Yl_imp <- (Yl - Y_hat)*Yl_exp + Y_hat

  return(list(Yl_imp,Xl,Bl))
}

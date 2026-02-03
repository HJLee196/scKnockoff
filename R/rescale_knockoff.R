#' Rescale knockoff variables
#'
#' @param X_knk A numeric matrix of dimension n x p containing knockoff
#'   variables, or a list of such n x p knockoff matrices.
#' @param np_data.matrix.exp A numeric matrix of dimension n x p indicating
#'   observed entries, with 1 for observed values and 0 for missing values.
#' @param center A numeric vector of length p giving the column-wise means
#'   of the observed elements in the original data.
#'
#' @return
#' A numeric matrix or a list of numeric matrices containing knockoff
#' variables re-centered to match the original data.
#'
#' @details
#' This function re-centers knockoff variables using the column-wise means
#' computed from the observed elements of the original data. No rescaling
#' of variability is performed.
#' Entries corresponding to missing values are set to zero.
#'
#' @examples
#' ## Simple matrix example
#' set.seed(1)
#' n <- 4
#' p <- 3
#'
#' X_knk <- matrix(rnorm(n * p), n, p)
#' np_data.matrix.exp <- matrix(c(1, 1, 0,
#'                                1, 0, 1,
#'                                1, 1, 1,
#'                                0, 1, 1),
#'                              n, p, byrow = TRUE)
#' center <- colMeans(matrix(rnorm(n * p), n, p))
#'
#' rescale_knockoff(X_knk, np_data.matrix.exp, center)
#'
#' ## List input example
#' X_knk_list <- list(X_knk, X_knk + 1)
#' rescale_knockoff(X_knk_list, np_data.matrix.exp, center)
#'
#' @export
rescale_knockoff <- function(X_knk, np_data.matrix.exp, center) {

  if (is.list(X_knk)) {
    return(lapply(X_knk, rescale_knockoff,
                  np_data.matrix.exp = np_data.matrix.exp,
                  center = center))
  }

  # Basic input checks (lightweight but helpful)
  if (!is.matrix(X_knk) || !is.numeric(X_knk)) {
    stop("`X_knk` must be a numeric matrix or a list of numeric matrices.")
  }
  np_data.matrix.exp = as.matrix(np_data.matrix.exp)
  if (!is.matrix(np_data.matrix.exp) || !all(dim(np_data.matrix.exp) == dim(X_knk))) {
    stop("`np_data.matrix.exp` must be a matrix with the same dimensions as `X_knk`.")
  }
  if (length(center) != ncol(X_knk)) {
    stop("`center` must have length equal to ncol(X_knk).")
  }

  # Re-center (column-wise)
  X_knk <- sweep(X_knk, 2, center, FUN = "+")

  # Keep only observed entries; set missing entries to 0 (as encoded by np_data.matrix.exp)
  X_knk <- X_knk * np_data.matrix.exp

  return(X_knk)
}

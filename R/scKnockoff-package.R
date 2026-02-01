#' scKnockoff: Knockoff-based differential analysis for scRNA-seq data
#'
#' Provides methods for differential analysis of single-cell RNA sequencing
#' (scRNA-seq) data under the knockoff framework. The package implements
#' large-scale conditional independence screening to control false discovery
#' rates in high-dimensional settings, with support for imputing missing data
#' via sc-softImpute incorporating additional latent variables and a matrix
#' decomposition-based knockoff construction.
#'
#' @details
#' The scKnockoff package is designed for high-dimensional single-cell data,
#' where the number of features can be much larger than the number of samples.
#' It supports missing data imputation using sc-softImpute with additional
#' latent variables to improve stability, and employs a matrix
#' decomposition-based knockoff construction that enables efficient inference
#' in high-dimensional settings.
#'
#' @name scKnockoff-package
#' @aliases scKnockoff
"_PACKAGE"

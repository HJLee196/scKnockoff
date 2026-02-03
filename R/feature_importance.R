
# Calculate the importance statistics using the Seurat package.
# return W_imp.

# Input:
# A single knockoff or a list of knockoffs (for eBH or multiple knockoffs)
#   The knockoffs should be rescaled!
# bonf: Whether we want to generate bonferroni corrected importance statistics

# Output:
# Either a vector of importance statistics or a list of importance statistics.

#' Calculate feature importance based on p-values.
#'
#' @param np_data.sub A \code{Seurat} object containing the original data.
#' @param np_data.knockoff A numeric matrix representing a single knockoff copy,
#'   or a list of such matrices for eBH or multiple knockoff procedures.
#'   The knockoffs must be rescaled prior to calling this function.
#' @param bonf Logical scalar; if \code{TRUE}, use Bonferroni-adjusted p-values
#'   to compute feature importance statistics. Default is \code{TRUE}.
#' @param test.use A character string specifying the test used to compute
#'   feature importance statistics. This argument corresponds to
#'   \code{test.use} in \code{Seurat::FindMarkers}. In addition, the lasso
#'   coefficient difference statistic is supported by setting
#'   \code{test.use = "LCD"}.
#' @param latent.vars A character vector specifying covariates to be included
#'   when computing feature importance. These variables must be present in
#'   the \code{Seurat} object.
#' @param ident.1 Same as \code{ident.1} in \code{Seurat::FindMarkers}. If
#'   \code{NULL}, the first level of \code{Idents(np_data.sub)} is used.
#' @param slot_name A character string specifying which assay slot (or layer) to use.
#'   This argument is passed to \code{Seurat::FindMarkers} as \code{slot}.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers}.
#'
#' @return
#' A numeric matrix of dimension \eqn{(m+1) \times p} when \code{np_data.knockoff}
#' is a list of \code{m} knockoff copies, or \eqn{2 \times p} when a single
#' knockoff matrix is provided. The first row corresponds to the original data.
#'
#' @family select
#'
#' @details Uses \code{Seurat::FindMarkers} to compute p-value-based feature
#'   importance statistics. The \code{"MAST"} option requires the \pkg{MAST}
#'   package (\url{https://github.com/RGLab/MAST}).
#'
#' @examples
#' # Read a dataset from the Seurat package
#' library(Seurat)
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_small <- Seurat::CreateSeuratObject(counts = pbmc_raw)
#'
#' # Prepare dataset
#' np_data <- t(as.matrix(pbmc_small[['RNA']]$counts))
#' np_data <- ifelse(np_data>0,log(np_data),0)
#' pbmc_small[['RNA']]$data = t(np_data)
#'
#' np_data_exp <- t(as.matrix(pbmc_small[['RNA']]$counts))
#' np_data_exp <- ifelse(np_data_exp>0,1,0)
#'
#' np_data_exp.count <- apply(np_data_exp,2,sum)
#'
#' # Trim down data
#' discard <- which(np_data_exp.count<20)
#'
#' np_data <- np_data[,-discard]
#' np_data_exp <- np_data_exp[,-discard]
#' np_data_exp.count <- np_data_exp.count[-discard]
#' pbmc_small <- subset(pbmc_small,features = colnames(np_data))
#'
#'
#' # Center target matrix, but only center the expressed parts.
#' np_data.avg <- (colSums(np_data, na.rm = TRUE))/(np_data_exp.count-1)
#' np_data_centered <- sweep(np_data,2,np_data.avg,FUN = "-")*(np_data_exp)
#'
#' # CDR as additional covariate
#' X <- matrix(apply(np_data_exp,1,sum),ncol=1)
#'
#' # Impute missing values
#' min(np_data_exp.count)
#' np_data_imp <- sc_softImpute(np_data_centered,np_data_exp,X,PC=15)
#'
#' # Decomp Knockoff construction
#' set.seed(1010)
#' decomp_knk <- create_decomp_knock(np_data_imp$X_imp,
#'                                   np_data_imp$Xl,
#'                                   np_data_imp$Bl,
#'                                   np_data_imp$err,
#'                                   np_data_imp$PC)
#'
#' decomp_knk <- rescale_knockoff(decomp_knk,np_data_exp,np_data.avg)
#'
#'
#' # Generate synthetic signal
#' np_data.scale <- apply(np_data, 2, scale)
#'
#' n <- dim(np_data.scale)[1]
#' p <- dim(np_data.scale)[2]
#'
#' norm_coef <- rep(0, p)
#' tmp_coef <- rep(0,20)
#'
#' sign_strength <- 3
#' for(i in 1:20){
#'   new_coef <- 0
#'   # truncate the coefficients which are too small
#'   while(isTRUE((new_coef)^2 < (sign_strength^2)*2*log(p)/n)){
#'     new_coef <- rnorm(1,0,sign_strength*sqrt(2*log(p)/n))
#'   }
#'
#'   tmp_coef[i] <- new_coef
#' }
#'
#' sig_coef <- sample(1:p,20)
#' norm_coef[sig_coef] <- tmp_coef
#'
#' label_p <- 1/(1 + exp(- np_data.scale %*% norm_coef))
#'
#' # Create random labels
#' label_binom <- rbinom(n = length(label_p), size = 1, prob = label_p)
#' label_ad <- ifelse((label_binom>=0.5),"y1", "y2")
#'
#' pbmc_small <- Seurat::SetIdent(object = pbmc_small, value = label_ad)
#'
#'
#' W_imp <- feature_importance(pbmc_small,decomp_knk,bonf=TRUE,
#'                            ident.1 = 'y1', ident.2 = 'y2',
#'                            test.use ="MAST")
#'
#' # Won't be able to select anything due to the small sample size
#' selected <- knockfilter_select(W_imp)
#'
#' @export
feature_importance <- function(np_data.sub,
                               np_data.knockoff,
                               bonf=TRUE,
                               test.use,
                               latent.vars = NULL,
                               ident.1 = NULL,
                               slot_name = "data", ...){
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' must be installed.", call. = FALSE)
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' must be installed.", call. = FALSE)
  }

  # np_data.sub
  if (!inherits(np_data.sub, "Seurat")) {
    stop("`np_data.sub` must be a Seurat object.", call. = FALSE)
  }

  # test.use
  if (missing(test.use) || !is.character(test.use) || length(test.use) != 1L || is.na(test.use)) {
    stop("`test.use` must be a single non-missing character string.", call. = FALSE)
  }

  # bonf
  if (!is.logical(bonf) || length(bonf) != 1L || is.na(bonf)) {
    stop("`bonf` must be a single logical value.", call. = FALSE)
  }

  # slot_name
  if (!is.character(slot_name) || length(slot_name) != 1L || is.na(slot_name) || !nzchar(slot_name)) {
    stop("`slot_name` must be a single non-empty character string.", call. = FALSE)
  }

  # latent.vars
  if (!is.null(latent.vars)) {
    if (!is.character(latent.vars) || anyNA(latent.vars)) {
      stop("`latent.vars` must be NULL or a character vector with no NA.", call. = FALSE)
    }
  }

  # ident.1 (optional)
  if (!is.null(ident.1)) {
    if (!is.character(ident.1) || length(ident.1) != 1L || is.na(ident.1) || !nzchar(ident.1)) {
      stop("`ident.1` must be NULL or a single non-empty character string.", call. = FALSE)
    }
  }

  # To match the p-values according to the variables
  # since FindMarkers() will return sorted p-values
  feature.names <- rownames(Seurat::GetAssayData(np_data.sub, layer = slot_name))
  label_Idents <- Seurat::Idents(np_data.sub)

  if (is.null(ident.1)) {
    ident.1 <- levels(label_Idents)[1]

    warning(
      sprintf(
        "`ident.1` was not specified; using '%s' as the reference identity.",
        ident.1
      ),
      call. = FALSE
    )

    ident.else <- paste0(levels(label_Idents)[-1], collapse = "_")

    label_Idents <- as.character(label_Idents)
    label_Idents[label_Idents != ident.1] <- ident.else
    label_Idents <- factor(label_Idents, levels = c(ident.1, ident.else))
  }

  # For multiple knockoffs, could be multiple independent
  # (single) knockoffs or multiple knockoffs.
  if(is.list(np_data.knockoff)){

    W_imp <- matrix(0,nrow = 1+length(np_data.knockoff), ncol = dim(np_data.knockoff[[1]])[2])

    if(test.use != "LCD"){
      # Find p-values for original data
      result <- Seurat::FindMarkers(np_data.sub, test.use = test.use,
                                    ident.1 = ident.1, latent.vars = latent.vars,
                                    slot = slot_name,
                                    min.pct = 0, logfc.threshold = 0, verbose = FALSE, ...)

      if(bonf){
        W_imp[1,] <- -log(result$p_val_adj[match(feature.names, rownames(result))]) # 0<= -log(p-value)

      } else {
        W_imp[1,] <- -log(result$p_val[match(feature.names, rownames(result))]) # 0<= -log(p-value)
      }


      # Importance statistics
      for(i in 1:length(np_data.knockoff)){
        # Row names and column names in new.data have to match those of np_data.sub!
        np_data.sub <-
          Seurat::SetAssayData(object = np_data.sub,
                               layer = slot_name,
                               new.data = Matrix::t(np_data.knockoff[[i]]))

        # Find p-values for knockoffs
        result <- Seurat::FindMarkers(np_data.sub, test.use = test.use,
                                      ident.1 = ident.1, latent.vars = latent.vars,
                                      slot = slot_name,
                                      min.pct = 0, logfc.threshold = 0, verbose = FALSE, ...)

        if(bonf){
          W_imp[1+i,] <- -log(result$p_val_adj[match(feature.names, rownames(result))]) # 0<= -log(p-value)
        } else {
          W_imp[1+i,] <- -log(result$p_val[match(feature.names, rownames(result))]) # 0<= -log(p-value)
        }
      }
    } else if (test.use == "LCD") {
      np_data.matrix <-
        Matrix::t(Seurat::GetAssayData(np_data.sub, layer = slot_name))

      m_kos <- length(np_data.knockoff)
      n <- nrow(np_data.matrix)

      n_total = n*(m_kos + 1) #kos + original

      np_data.matrix.sum <- Matrix::colSums(np_data.matrix)
      np_data.knockoff.sum <-
        Reduce(`+`, lapply(np_data.knockoff, \(x) Matrix::colSums(x)))

      np_data.all.mean <- (np_data.matrix.sum + np_data.knockoff.sum)/n_total

      np_data.matrix.sq_sum <- Matrix::colSums(np_data.matrix*np_data.matrix)
      np_data.knockoff.sq_sum <-
        Reduce(`+`, lapply(np_data.knockoff, \(x) Matrix::colSums(x*x)))

      np_data.all.sq_mean <- (np_data.matrix.sq_sum + np_data.knockoff.sq_sum)/n_total

      np_data.all.var <- np_data.all.sq_mean - np_data.all.mean^2
      np_data.matrix_sd <- sqrt((n_total/(n_total - 1)) * pmax(np_data.all.var, 0))

      np_data.full <-
        sweep(x = np_data.matrix, MARGIN = 2, STATS = np_data.matrix_sd, FUN = "/")

      for (m_ko in 1:m_kos)
      {
        np_data.knockoff.scale <-
          sweep(x = np_data.knockoff[[m_ko]], MARGIN = 2, STATS = np_data.matrix_sd, FUN = "/")
        colnames(np_data.knockoff.scale) <- paste0(colnames(np_data.matrix), "_", m_ko)
        np_data.full <- cbind(np_data.full, np_data.knockoff.scale)
      }

      # generate X
      X <- NULL
      if (!is.null(latent.vars) && length(latent.vars) > 0) {
        for (imp_vars_temp in latent.vars) {

          # meta.data column as vector
          v <- np_data.sub[[imp_vars_temp]][, 1, drop = TRUE]

          # if non-numeric: one-hot encode
          if (!is.numeric(v)) {
            mm <- stats::model.matrix(~ x, data = data.frame(x = v))  # includes intercept
            mm <- mm[, -1, drop = FALSE]  # drop intercept -> k-1 dummies
            # optional: prefix colnames to avoid collisions
            colnames(mm) <- paste0(imp_vars_temp, ":", colnames(mm))
            vmat <- mm
          } else {
            vmat <- matrix(v, ncol = 1)
            colnames(vmat) <- imp_vars_temp
          }

          X <- if (is.null(X)) vmat else cbind(X, vmat)
        }
      }

      if(is.null(X)){
        X_center <- NULL
      } else {
        X_center <- scale(X)
      }

      np_data.full <-
        cbind(X_center, np_data.full)

      cv_lasso_fit <-
        glmnet::cv.glmnet(y = label_Idents, x = np_data.full,
                          nfolds = 10,
                          alpha = 1, # alpha = 1 means lasso
                          family = "binomial",
                          standardize = FALSE)

      cv_lasso_fit_best_lambda <- cv_lasso_fit[["lambda.1se"]] #lambda.1se #lambda.min
      cv_lasso_fit_coef <- Matrix::t(stats::coef(cv_lasso_fit, s = cv_lasso_fit_best_lambda))
      cv_lasso_fit_coef <- cv_lasso_fit_coef[,colnames(cv_lasso_fit_coef) != "(Intercept)", drop = F]

      cv_lasso_fit_coef_orig <- cv_lasso_fit_coef[,match(colnames(np_data.matrix), colnames(cv_lasso_fit_coef))]

      W_imp[1, ] <- abs(cv_lasso_fit_coef_orig)
      for (m_ko in 1:m_kos)
      {
        colnames_mth_ko <- paste0(colnames(np_data.matrix), "_", m_ko)
        cv_lasso_fit_coef_mth_ko <- cv_lasso_fit_coef[,match(colnames_mth_ko, colnames(cv_lasso_fit_coef))]
        W_imp[1 + m_ko, ] <- abs(cv_lasso_fit_coef_mth_ko)
      }
    }

    # For single knockoffs
  } else {
    W_imp <- matrix(0,nrow = 2, ncol = dim(np_data.knockoff)[2])

    if (test.use != "LCD"){
      # Find p-values for original data
      result <- Seurat::FindMarkers(np_data.sub, test.use = test.use,
                                    ident.1 = ident.1, latent.vars = latent.vars,
                                    slot = slot_name, min.pct = 0, logfc.threshold=0, verbose = FALSE, ...)

      if(bonf){
        W_imp[1,] <- -log(result$p_val_adj[match(feature.names, rownames(result))]) # 0<= -log(p-value)
      } else {
        W_imp[1,] <- -log(result$p_val[match(feature.names, rownames(result))]) # 0<= -log(p-value)
      }

      # Row names and column names in new.data have to match those of np_data.sub!
      np_data.sub <-
        Seurat::SetAssayData(object = np_data.sub,
                             layer = slot_name,
                             new.data = Matrix::t(np_data.knockoff))

      # Find p-values for knockoffs
      result <- Seurat::FindMarkers(np_data.sub, test.use = test.use,
                                    ident.1 = ident.1, latent.vars = latent.vars,
                                    slot = slot_name, min.pct = 0,logfc.threshold=0, verbose = FALSE, ...)

      if(bonf){
        W_imp[2,] <- -log(result$p_val_adj[match(feature.names, rownames(result))]) # 0<= -log(p-value)

      } else {
        W_imp[2,] <- -log(result$p_val[match(feature.names, rownames(result))]) # 0<= -log(p-value)
      }
    } else if (test.use == "LCD"){
      np_data.matrix <-
        Matrix::t(Seurat::GetAssayData(np_data.sub, layer = slot_name))

      m_kos <- 1
      n <- nrow(np_data.matrix)

      n_total = n*(m_kos + 1) #kos + original

      np_data.matrix.sum <- Matrix::colSums(np_data.matrix)
      np_data.knockoff.sum <- Matrix::colSums(np_data.knockoff)

      np_data.all.mean <- (np_data.matrix.sum + np_data.knockoff.sum)/n_total

      np_data.matrix.sq_sum <- Matrix::colSums(np_data.matrix * np_data.matrix)
      np_data.knockoff.sq_sum <- Matrix::colSums(np_data.knockoff * np_data.knockoff)

      np_data.all.sq_mean <- (np_data.matrix.sq_sum + np_data.knockoff.sq_sum)/n_total

      np_data.all.var <- np_data.all.sq_mean - np_data.all.mean^2
      np_data.matrix_sd <- sqrt((n_total/(n_total - 1)) * pmax(np_data.all.var, 0))


      np_data.matrix_scale <-
        sweep(x = np_data.matrix, 2, np_data.matrix_sd, "/")

      np_data.knockoff.scale <-
        sweep(x = np_data.knockoff, 2, np_data.matrix_sd, "/")

      colnames(np_data.knockoff.scale) <- paste0(colnames(np_data.matrix_scale), "_1")

      # generate X
      X <- NULL
      if (!is.null(latent.vars) && length(latent.vars) > 0) {
        for (imp_vars_temp in latent.vars) {

          # meta.data column as vector
          v <- np_data.sub[[imp_vars_temp]][, 1, drop = TRUE]

          # if non-numeric: one-hot encode
          if (!is.numeric(v)) {
            mm <- stats::model.matrix(~ x, data = data.frame(x = v))  # includes intercept
            mm <- mm[, -1, drop = FALSE]  # drop intercept -> k-1 dummies
            # optional: prefix colnames to avoid collisions
            colnames(mm) <- paste0(imp_vars_temp, ":", colnames(mm))
            vmat <- mm
          } else {
            vmat <- matrix(v, ncol = 1)
            colnames(vmat) <- imp_vars_temp
          }

          X <- if (is.null(X)) vmat else cbind(X, vmat)
        }
      }

      if(is.null(X)){
        X_center <- NULL
      } else {
        X_center <- scale(X)
      }

      np_data.full <- cbind(X_center,
                           np_data.matrix_scale,
                           np_data.knockoff.scale)

      cv_lasso_fit <-
        glmnet::cv.glmnet(y = label_Idents, x = np_data.full,
                          nfolds = 10,
                          alpha = 1, # alpha = 1 means lasso
                          family = "binomial",
                          standardize = FALSE)

      cv_lasso_fit_best_lambda <- cv_lasso_fit[["lambda.1se"]] #lambda.1se #lambda.min
      cv_lasso_fit_coef <- Matrix::t(stats::coef(cv_lasso_fit, s = cv_lasso_fit_best_lambda))
      cv_lasso_fit_coef <- cv_lasso_fit_coef[,colnames(cv_lasso_fit_coef) != "(Intercept)", drop = F]

      cv_lasso_fit_coef_orig <- cv_lasso_fit_coef[,match(colnames(np_data.matrix.scale), colnames(cv_lasso_fit_coef))]
      cv_lasso_fit_coef_ko <- cv_lasso_fit_coef[,match(colnames(np_data.knockoff.scale), colnames(cv_lasso_fit_coef))]

      W_imp[1,] <- abs(cv_lasso_fit_coef_orig)
      W_imp[2,] <- abs(cv_lasso_fit_coef_ko)
    }

  }

  if (sum(is.na(W_imp)) != 0) {
    warning(
      "NA values were detected during feature importance computation. ",
      "Please check the p-values returned by `Seurat::FindMarkers()`.",
      call. = FALSE
    )
  }

  return(W_imp)
}

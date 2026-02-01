#' Variable selection using decomp knockoffs with a Seurat dataset
#'
#' @param Seurat_data A Seurat object
#' @param PC A positive integer specifying the number of latent factors.
#' @param latent.vars_imp A character vector specifying covariates to be included
#'   when performing \code{sc_softImpute()}. These variables must be present in
#'   the \code{Seurat} object.
#' @param latent.vars_comp A character vector specifying covariates to be included
#'   when computing feature importance. These variables must be present in
#'   the \code{Seurat} object.
#' @param ident.1 Same as \code{ident.1} in \code{Seurat::FindMarkers}. If
#'   \code{NULL}, the first level of \code{Idents(np_data.sub)} is used.
#' @param ident.2 Same as \code{ident.2} in \code{Seurat::FindMarkers}. If
#'   \code{NULL}, \code{feature_importance()} uses the default behavior of
#'   \code{Seurat::FindMarkers}.
#' @param slot_name A character string specifying which assay slot (or layer) to use.
#'   This argument is passed to \code{Seurat::FindMarkers} as \code{slot}.
#' @param test.use A character string specifying the test used to compute
#'   feature importance statistics. This argument corresponds to
#'   \code{test.use} in \code{Seurat::FindMarkers}. In addition, the lasso
#'   coefficient difference statistic is supported by setting
#'   \code{test.use = "LCD"}.
#' @param q A numeric scalar in (0, 1) specifying the target FDR level.
#' @param m A positive integer specifying the number of knockoffs to be constructed;
#'   default is 1.
#'   If \code{m = 1}, a single knockoff is constructed using
#'   \code{create_decomp_knock()}.
#'   If \code{m >= 2}, multiple knockoffs are constructed using
#'   \code{create_multi_decomp_knock()}.
#' @param max_it A positive integer specifying the maximum number of iterations;
#'   default is 100.
#'
#' @details
#' This function provides a minimal, end-to-end example of variable selection
#' using matrix decomposition-based knockoffs on a \code{Seurat} dataset.
#' The normalized expression matrix (specified by \code{slot_name}) is centered
#' gene-wise using means computed over observed entries (as indicated by the
#' count-based expression pattern). The centered matrix is then imputed via
#' \code{sc_softImpute()}, followed by knockoff construction and feature
#' importance calculation using \code{feature_importance()}.
#' Variables are selected by \code{knockfilter_select()} to control the FDR at
#' level \code{q}.
#'
#' @return
#' A list containing the indices of selected variables and the corresponding
#' variable names.
#'
#' @examples
#' # Read a dataset from the Seurat package
#' library(Seurat)
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_raw <- as(as.matrix(pbmc_raw), "dgCMatrix")
#' pbmc_small <- Seurat::CreateSeuratObject(counts = pbmc_raw)
#'
#' pbmc_small.count <- Seurat::GetAssayData(object = pbmc_small, layer = "count")
#' pbmc_small.data <-
#'   Seurat::NormalizeData(object = pbmc_small.count,
#'                         normalization.method = "LogNormalize",
#'                         scale.factor = 10000, #' Default = 10000,
#'                         margin = 1,
#'                         verbose = FALSE)
#'
#' pbmc_small = Seurat::SetAssayData(object = pbmc_small, layer = "data", new.data = pbmc_small.data)
#'
#' pbmc_small.exp <- t(pbmc_small.count > 0)
#' CDR <- rowMeans(pbmc_small.exp)
#'
#' pbmc_small$CDR = CDR
#'
#' Seurat::Idents(pbmc_small) = as.factor(rep(c("A", "B"), each = length(CDR)/2))
#'
#' set.seed(1010)
#' full_process(Seurat_data = pbmc_small,
#'              PC = 3,
#'              ident.1 = "A",
#'              ident.2 = "B",
#'              latent.vars_comp = "CDR",
#'              latent.vars_imp = "CDR",
#'              test.use = "LCD",
#'              m = 1)
#'
#'
#' @export
full_process <- function(Seurat_data,
                         latent.vars_imp = NULL,
                         latent.vars_comp = NULL,
                         PC,
                         ident.1 = NULL,
                         ident.2 = NULL,
                         slot_name = "data",
                         test.use,
                         q = 0.1,
                         m = 1,
                         max_it = 100){

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' must be installed.", call. = FALSE)
  }

  # Input validation ----------------------------------------------------

  # Seurat_data
  if (!inherits(Seurat_data, "Seurat")) {
    stop("`Seurat_data` must be a Seurat object.", call. = FALSE)
  }

  # PC
  if (!is.numeric(PC) || length(PC) != 1L || is.na(PC) || PC <= 0 || PC != as.integer(PC)) {
    stop("`PC` must be a positive integer.", call. = FALSE)
  }

  # m
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m <= 0 || m != as.integer(m)) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }
  if (m < 1) {
    stop("`m` must be at least 1.", call. = FALSE)
  }

  # q
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0 || q >= 1) {
    stop("`q` must be a numeric scalar in (0, 1).", call. = FALSE)
  }

  # max_it
  if (!is.numeric(max_it) || length(max_it) != 1L || is.na(max_it) ||
      max_it <= 0 || max_it != as.integer(max_it)) {
    stop("`max_it` must be a positive integer.", call. = FALSE)
  }

  # slot_name
  if (!is.character(slot_name) || length(slot_name) != 1L || is.na(slot_name) || !nzchar(slot_name)) {
    stop("`slot_name` must be a single non-empty character string.", call. = FALSE)
  }

  # test.use
  if (missing(test.use) || !is.character(test.use) || length(test.use) != 1L || is.na(test.use) || !nzchar(test.use)) {
    stop("`test.use` must be a single non-empty character string.", call. = FALSE)
  }

  # latent.vars_imp / latent.vars_comp
  if (!is.null(latent.vars_imp)) {
    if (!is.character(latent.vars_imp) || anyNA(latent.vars_imp)) {
      stop("`latent.vars_imp` must be NULL or a character vector with no NA.", call. = FALSE)
    }
  }
  if (!is.null(latent.vars_comp)) {
    if (!is.character(latent.vars_comp) || anyNA(latent.vars_comp)) {
      stop("`latent.vars_comp` must be NULL or a character vector with no NA.", call. = FALSE)
    }
  }

  # ident.1 / ident.2
  if (!is.null(ident.1)) {
    if (!is.character(ident.1) || length(ident.1) != 1L || is.na(ident.1) || !nzchar(ident.1)) {
      stop("`ident.1` must be NULL or a single non-empty character string.", call. = FALSE)
    }
  }
  if (!is.null(ident.2)) {
    if (!is.character(ident.2) || length(ident.2) != 1L || is.na(ident.2) || !nzchar(ident.2)) {
      stop("`ident.2` must be NULL or a single non-empty character string.", call. = FALSE)
    }
  }

  # Check that requested metadata vars exist
  all_meta <- colnames(Seurat_data[[]])
  if (!is.null(latent.vars_imp)) {
    miss_imp <- setdiff(latent.vars_imp, all_meta)
    if (length(miss_imp) > 0) {
      stop(
        paste0("The following `latent.vars_imp` are not in `Seurat_data@meta.data`: ",
               paste(miss_imp, collapse = ", ")),
        call. = FALSE
      )
    }
  }
  if (!is.null(latent.vars_comp)) {
    miss_comp <- setdiff(latent.vars_comp, all_meta)
    if (length(miss_comp) > 0) {
      stop(
        paste0("The following `latent.vars_comp` are not in `Seurat_data@meta.data`: ",
               paste(miss_comp, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  # Check ident levels
  id_levels <- levels(as.factor(Seurat::Idents(Seurat_data)))
  if (length(id_levels) < 2) {
    stop("`Seurat_data` must have at least two identity groups in `Idents()`.", call. = FALSE)
  }
  if (!is.null(ident.1) && !(ident.1 %in% id_levels)) {
    stop("`ident.1` is not a level in `Idents(Seurat_data)`.", call. = FALSE)
  }
  if (!is.null(ident.2) && !(ident.2 %in% id_levels)) {
    stop("`ident.2` is not a level in `Idents(Seurat_data)`.", call. = FALSE)
  }
  if (!is.null(ident.1) && !is.null(ident.2) && identical(ident.1, ident.2)) {
    stop("`ident.1` and `ident.2` must be different.", call. = FALSE)
  }

  if (Seurat::DefaultAssay(Seurat_data) != "RNA"){
    warning(
      "This method is designed for scRNA-seq data. ",
      "Using a default assay other than 'RNA' may lead to unexpected results.",
      call. = FALSE
    )
  }
  # Take a subset of genes
  feature.names <- rownames(Seurat_data)

  np_data.count <- Seurat::GetAssayData(object = Seurat_data, layer = "count")

  np_data.matrix.exp <- Matrix::t(np_data.count > 0) # indicator matrix for expressed genes,
  np_data.matrix.exp.count <- Matrix::colSums(np_data.matrix.exp)

  # # Calculate cellular detection rate (CDR)
  # np_data.matrix <- np_data.count # do not transpose here. Keep the genes as rows.
  #
  # CDR <- rowMeans(np_data.matrix.exp) # expressed genes in each CELL
  # Seurat_data$CDR <- CDR

  # generate X
  X <- NULL
  if (!is.null(latent.vars_imp) && length(latent.vars_imp) > 0) {
    for (imp_vars_temp in latent.vars_imp) {

      # meta.data column as vector
      v <- Seurat_data[[imp_vars_temp]][, 1, drop = TRUE]

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

  # Exclude empty variables
  if (is.null(X)){
    X_p <- 0
  } else {
    X_p <- ncol(X)
  }
  feature.exp.ind <- which(np_data.matrix.exp.count > (PC + X_p + 1))

  if (length(feature.exp.ind) == 0) {
    stop(
      "No features remain after filtering by non-missing counts. ",
      "Try reducing `PC` or using fewer covariates in `latent.vars_imp`.",
      call. = FALSE
    )
  }

  feature.names.new <- feature.names[feature.exp.ind]  # length(feature.names.new)

  np_data.matrix.exp <- np_data.matrix.exp[,feature.exp.ind]
  np_data.matrix.exp.count <- np_data.matrix.exp.count[feature.exp.ind]

  # Refresh subset
  np_data.sub <- subset(x = Seurat_data, features = feature.names.new)

  np_data.matrix <- Seurat::GetAssayData(object = np_data.sub, layer = slot_name) # Normalized data matrix
  np_data.matrix <- Matrix::t(np_data.matrix)

  # Center target matrix, but only center the expressed parts.
  np_data.matrix.fix <- (Matrix::colSums(np_data.matrix*np_data.matrix.exp))/np_data.matrix.exp.count # na.rm = T is removed
  np_data.matrix <- sweep(np_data.matrix,2,np_data.matrix.fix,FUN = "-")*(np_data.matrix.exp)

  #### Imputation ####
  np_data.matrix.exp <- np_data.matrix.exp*1
  np_data.matrix.imp <- sc_softImpute(np_data.matrix = np_data.matrix,
                                      np_data.matrix.exp = np_data.matrix.exp,
                                      X = X,
                                      PC = PC,
                                      max_it = max_it)


  #### Decomp Knockoff construction ####
  if (m == 1){
    decomp_knk <-
      create_decomp_knock(np_data.matrix.imp = np_data.matrix.imp$X_imp,
                          Xl = np_data.matrix.imp$Xl,
                          Bl = np_data.matrix.imp$Bl,
                          err = np_data.matrix.imp$err,
                          decomp = TRUE,
                          PC = np_data.matrix.imp$PC)
  } else {
    decomp_knk <-
      create_multi_decomp_knock(np_data.matrix.imp = np_data.matrix.imp$X_imp,
                                Xl = np_data.matrix.imp$Xl,
                                Bl = np_data.matrix.imp$Bl,
                                err = np_data.matrix.imp$err,
                                PC = np_data.matrix.imp$PC,
                                m = m)
  }


  decomp_knk <- rescale_knockoff(decomp_knk,
                                 np_data.matrix.exp,
                                 np_data.matrix.fix)

  #### Feature calculation and variable selection ####

  W_imp <- feature_importance(np_data.sub = np_data.sub,
                              np_data.knockoff = decomp_knk,
                              bonf = TRUE,
                              ident.1 = ident.1, ident.2 = ident.2, slot_name = slot_name,
                              test.use = test.use, latent.vars = latent.vars_comp)


  selected <- knockfilter_select(W_imp = W_imp, q = q, m = m)

  genes <- if (is.null(selected)) NULL else feature.names.new[selected]
  return(list(selected = selected, genes = genes))
}

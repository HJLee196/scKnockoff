#' Generate a toy Seurat object for differential analysis
#'
#' This function generates a toy single-cell RNA-seq dataset with disease
#' labels, observed covariates, latent factors, and ground-truth signal genes.
#' The simulated latent expression layer follows
#' \eqn{G = X B0 + A B1 + E}, where \eqn{G} is interpreted as a Gaussian
#' log-mean expression matrix. Observed UMI-like counts are generated from
#' \eqn{Y_{ij} \sim Poisson(\mu_{ij})}, where
#' \eqn{\mu_{ij} = \exp(G_{ij})}, up to cell-specific library size scaling.
#'
#' @param n_genes Number of genes.
#' @param donors_per_group Number of donors in each disease group.
#' @param cell_types Character vector of cell types.
#' @param cells_per_donor_per_type Number of cells per donor and cell type.
#' @param n_batch Number of batches.
#' @param n_signals Number of ground-truth disease-associated genes.
#' @param signal_strength Fold change for disease-associated genes on the
#'   expression scale.
#' @param r Number of latent factors.
#' @param sigma Standard deviation of Gaussian noise in the log-expression layer.
#' @param mean_umi Target mean library size for the generated count matrix.
#' @param lib_size_sd Standard deviation of the log-library-size factor.
#' @param age_effect_sd Standard deviation of age effects across genes.
#' @param batch_effect_sd Standard deviation of batch effects across genes.
#' @param cell_type_effect_sd Standard deviation of cell-type effects across genes.
#' @param loading_prop Proportion of genes with nonzero loadings for each latent
#'   factor.
#' @param loading_sd Standard deviation of nonzero latent-factor loadings.
#' @param normalize Logical; whether to run \code{Seurat::NormalizeData()}.
#' @param scale_data Logical; whether to run \code{Seurat::ScaleData()}.
#' @param seed Random seed. If \code{NULL}, the random seed is not set.
#'
#' @return
#' A Seurat object containing simulated UMI-like counts and cell-level metadata.
#' The ground-truth signal indices are stored in \code{object@misc$true_signal},
#' and the names of the ground-truth signal genes are stored in
#' \code{object@misc$true_signal_genes}. The object also stores simulation
#' components such as \code{B0}, \code{B1}, \code{A}, \code{X}, \code{G},
#' \code{mu}, and the true library sizes in \code{object@misc}.
#'
#' @export
make_toy_seurat_ad <- function(
    n_genes = 1000,
    donors_per_group = 6,
    cell_types = c("Astro", "Micro", "Neuron"),
    cells_per_donor_per_type = 50,
    n_batch = 2,
    n_signals = 50,
    signal_strength = 2,
    r = 5,
    sigma = 0.5,
    mean_umi = 5000,
    lib_size_sd = 0.4,
    age_effect_sd = 0.05,
    batch_effect_sd = 0.15,
    cell_type_effect_sd = 0.25,
    loading_prop = 0.2,
    loading_sd = 0.7,
    normalize = TRUE,
    scale_data = TRUE,
    seed = 1
) {
  # Required packages
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(
      "Seurat is required to build the toy Seurat object. Install Seurat or skip this example.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Matrix is required to create a sparse count matrix.",
      call. = FALSE
    )
  }
  
  # Input checks
  .check_positive_integer(n_genes, "n_genes")
  .check_positive_integer(donors_per_group, "donors_per_group")
  .check_positive_integer(cells_per_donor_per_type, "cells_per_donor_per_type")
  .check_positive_integer(n_batch, "n_batch")
  .check_positive_integer(r, "r")
  .check_nonnegative_integer(n_signals, "n_signals")
  
  if (n_signals > n_genes) {
    stop("`n_signals` cannot be larger than `n_genes`.", call. = FALSE)
  }
  
  if (!is.character(cell_types) || length(cell_types) == 0L ||
      any(is.na(cell_types)) || any(cell_types == "")) {
    stop(
      "`cell_types` must be a non-empty character vector without missing or empty values.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(cell_types)) {
    stop("`cell_types` must not contain duplicated values.", call. = FALSE)
  }
  
  if (!is.numeric(signal_strength) || length(signal_strength) != 1L ||
      is.na(signal_strength) || signal_strength <= 1) {
    stop("`signal_strength` must be a number greater than 1.", call. = FALSE)
  }
  
  .check_nonnegative_number(sigma, "sigma")
  .check_positive_number(mean_umi, "mean_umi")
  .check_nonnegative_number(lib_size_sd, "lib_size_sd")
  .check_nonnegative_number(age_effect_sd, "age_effect_sd")
  .check_nonnegative_number(batch_effect_sd, "batch_effect_sd")
  .check_nonnegative_number(cell_type_effect_sd, "cell_type_effect_sd")
  .check_nonnegative_number(loading_sd, "loading_sd")
  
  if (!is.numeric(loading_prop) || length(loading_prop) != 1L ||
      is.na(loading_prop) || loading_prop <= 0 || loading_prop > 1) {
    stop(
      "`loading_prop` must be a number greater than 0 and at most 1.",
      call. = FALSE
    )
  }
  
  .check_logical_scalar(normalize, "normalize")
  .check_logical_scalar(scale_data, "scale_data")
  
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L ||
        is.na(seed) || seed != as.integer(seed)) {
      stop("`seed` must be an integer or NULL.", call. = FALSE)
    }
    set.seed(seed)
  }
  
  # 1. Metadata
  n_types <- length(cell_types)
  n_donors <- 2L * donors_per_group
  donors <- paste0("D", seq_len(n_donors))
  
  donor_disease <- rep(c("Control", "AD"), each = donors_per_group)
  names(donor_disease) <- donors
  
  batch_levels <- paste0("B", seq_len(n_batch))
  
  # Balanced batch assignment across donors when possible.
  batch_by_donor <- sample(rep(batch_levels, length.out = n_donors))
  names(batch_by_donor) <- donors
  
  donor_age <- round(stats::runif(n_donors, 65, 90))
  names(donor_age) <- donors
  
  n_cells <- n_donors * n_types * cells_per_donor_per_type
  
  meta <- data.frame(
    cell_id = paste0("cell", seq_len(n_cells)),
    donor = rep(donors, each = n_types * cells_per_donor_per_type),
    cell_type = rep(
      rep(cell_types, each = cells_per_donor_per_type),
      times = n_donors
    ),
    stringsAsFactors = FALSE
  )
  
  meta$disease <- unname(donor_disease[meta$donor])
  meta$batch <- unname(batch_by_donor[meta$donor])
  meta$age <- unname(donor_age[meta$donor])
  
  meta$disease <- factor(meta$disease, levels = c("Control", "AD"))
  meta$batch <- factor(meta$batch, levels = batch_levels)
  meta$cell_type <- factor(meta$cell_type, levels = cell_types)
  
  genes <- paste0("gene", seq_len(n_genes))
  
  # 2. Observed covariates X
  #
  # This is the X in G = X B0 + A B1 + E.
  # The diseaseAD coefficient defines the ground-truth disease-associated genes.
  X <- stats::model.matrix(
    ~ disease + batch + age + cell_type,
    data = meta
  )
  
  # Standardize age column only, if present.
  if ("age" %in% colnames(X)) {
    X[, "age"] <- as.numeric(scale(X[, "age"]))
  }
  
  q <- ncol(X)
  
  # 3. Construct B0
  B0 <- matrix(0, nrow = q, ncol = n_genes)
  rownames(B0) <- colnames(X)
  colnames(B0) <- genes
  
  # Baseline log-expression.
  #
  # Since G is interpreted as log-expression, this intercept should be
  # sufficiently positive to avoid generating too many near-zero means after
  # the inverse log transform.
  baseline_log_expression <- stats::rnorm(n_genes, mean = 1.5, sd = 1)
  B0["(Intercept)", ] <- baseline_log_expression
  
  # Small age effects.
  if ("age" %in% rownames(B0)) {
    B0["age", ] <- stats::rnorm(n_genes, mean = 0, sd = age_effect_sd)
  }
  
  # Batch effects.
  batch_rows <- grep("^batch", rownames(B0), value = TRUE)
  if (length(batch_rows) > 0L) {
    for (rr in batch_rows) {
      B0[rr, ] <- stats::rnorm(n_genes, mean = 0, sd = batch_effect_sd)
    }
  }
  
  # Cell-type effects.
  cell_type_rows <- grep("^cell_type", rownames(B0), value = TRUE)
  if (length(cell_type_rows) > 0L) {
    for (rr in cell_type_rows) {
      B0[rr, ] <- stats::rnorm(n_genes, mean = 0, sd = cell_type_effect_sd)
    }
  }
  
  # Disease signals.
  signal_genes <- integer(0)
  signal_gene_names <- character(0)
  
  if (n_signals > 0L) {
    signal_genes <- sort(sample(seq_len(n_genes), n_signals))
    signal_gene_names <- genes[signal_genes]
    
    disease_row <- "diseaseAD"
    
    if (!disease_row %in% rownames(B0)) {
      stop("Could not find `diseaseAD` in the design matrix.", call. = FALSE)
    }
    
    n_up <- ceiling(n_signals / 2)
    n_down <- n_signals - n_up
    
    # Since G is log-expression, these are log fold changes.
    disease_effect <- c(
      rep(log(signal_strength), n_up),
      rep(-log(signal_strength), n_down)
    )
    
    disease_effect <- sample(disease_effect)
    
    B0[disease_row, signal_genes] <- disease_effect
  }
  
  # 4. Latent factors A and loadings B1
  A <- matrix(
    stats::rnorm(n_cells * r),
    nrow = n_cells,
    ncol = r
  )
  rownames(A) <- meta$cell_id
  colnames(A) <- paste0("LF", seq_len(r))
  
  B1 <- matrix(0, nrow = r, ncol = n_genes)
  rownames(B1) <- colnames(A)
  colnames(B1) <- genes
  
  n_loading_genes <- floor(loading_prop * n_genes)
  
  if (n_loading_genes < 1L) {
    stop(
      "`loading_prop` is too small for the given `n_genes`; at least one loading gene is required.",
      call. = FALSE
    )
  }
  
  for (k in seq_len(r)) {
    idx <- sample(seq_len(n_genes), n_loading_genes)
    B1[k, idx] <- stats::rnorm(length(idx), mean = 0, sd = loading_sd)
  }
  
  # 5. Gaussian log-expression layer
  #
  # G = X B0 + A B1 + E
  E <- matrix(
    stats::rnorm(n_cells * n_genes, mean = 0, sd = sigma),
    nrow = n_cells,
    ncol = n_genes
  )
  
  G <- X %*% B0 + A %*% B1 + E
  rownames(G) <- meta$cell_id
  colnames(G) <- genes
  
  # 6. Generate observed counts from the Gaussian log-expression layer
  #
  # Here G is the underlying Gaussian log-mean expression matrix.
  # We convert it to a positive mean count matrix using the exponential
  # transform and then generate UMI-like counts from a Poisson model.
  mu <- exp(G)
  
  # Add cell-specific library size factors.
  library_size_factor <- exp(
    stats::rnorm(n_cells, mean = 0, sd = lib_size_sd)
  )
  
  mu <- mu * library_size_factor
  
  # Rescale so that the expected mean library size is approximately mean_umi.
  current_mean_library <- mean(rowSums(mu))
  
  if (!is.finite(current_mean_library) || current_mean_library <= 0) {
    stop(
      "The generated mean expression matrix is degenerate. Try increasing the baseline expression or decreasing noise.",
      call. = FALSE
    )
  }
  
  mu <- mu * mean_umi / current_mean_library
  
  counts_dense <- matrix(
    stats::rpois(n_cells * n_genes, lambda = as.vector(mu)),
    nrow = n_cells,
    ncol = n_genes
  )
  
  rownames(counts_dense) <- meta$cell_id
  colnames(counts_dense) <- genes
  
  counts <- t(counts_dense)
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  
  library_size <- Matrix::colSums(counts)
  
  # 7. Build Seurat object
  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "toyAD"
  )
  
  seu$donor <- meta$donor
  seu$disease <- meta$disease
  seu$cell_type <- meta$cell_type
  seu$batch <- meta$batch
  seu$age <- meta$age
  seu$library_size_true <- as.numeric(library_size)
  
  Seurat::Idents(seu) <- "disease"
  
  if (isTRUE(normalize)) {
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  }
  
  if (isTRUE(scale_data)) {
    seu <- Seurat::ScaleData(seu, features = rownames(seu), verbose = FALSE)
  }
  
  # 8. Store ground truth and simulation components
  beta_disease <- rep(0, n_genes)
  names(beta_disease) <- genes
  
  if (n_signals > 0L) {
    beta_disease[signal_genes] <- B0["diseaseAD", signal_genes]
  }
  
  seu@misc$true_signal <- signal_genes
  seu@misc$true_signal_genes <- signal_gene_names
  seu@misc$beta_disease <- beta_disease
  seu@misc$B0 <- B0
  seu@misc$B1 <- B1
  seu@misc$A <- A
  seu@misc$X <- X
  seu@misc$G <- G
  seu@misc$mu <- mu
  seu@misc$library_size_factor <- library_size_factor
  seu@misc$library_size_true <- as.numeric(library_size)
  
  seu@misc$simulation_model <- paste(
    "G = X B0 + A B1 + E;",
    "G is the Gaussian log-mean expression layer;",
    "counts_ij ~ Poisson(mu_ij), where mu_ij = exp(G_ij)",
    "with cell-specific library size scaling."
  )
  
  seu
}
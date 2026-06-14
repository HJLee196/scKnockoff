Introduction to scKnockoff - Easy
================
Hyunjae Lee

# 1. Load scKnockoff

``` r
# if (!requireNamespace("scKnockoff", quietly = TRUE)) {
#   devtools::install_github("HJLee196/scKnockoff")
# }

library(scKnockoff)
```

# 2. Create dataset

We generate a synthetic single-cell RNA-seq dataset from the latent
factor model $G = X B_0 + A B_1 + E$, where $X$ contains observed
covariates such as disease status, batch, age, and cell type. Disease
status is assigned at the donor level, and a subset of genes is assigned
nonzero disease effects in $B_0$, serving as the ground-truth signals.
The latent expression matrix $G$ is transformed into cell-specific
relative expression probabilities, from which observed UMI counts are
sampled using a multinomial model.

``` r
make_toy_seurat_ad <- function(
    n_genes = 1000,
    donors_per_group = 6,
    cell_types = c("Astro", "Micro", "Neuron"),
    cells_per_donor_per_type = 50,
    n_batch = 2,
    # disease-associated genes
    n_signals = 50,
    # fold-change on relative abundance scale
    signal_strength = 2,
    # latent dimension
    r = 5,
    # noise level in Gaussian latent layer
    sigma = 0.5,
    # library size distribution
    mean_umi = 5000,
    lib_size_sd = 0.4,
    # effect sizes for observed covariates
    age_effect_sd = 0.05,
    batch_effect_sd = 0.15,
    cell_type_effect_sd = 0.25,
    # loading sparsity
    loading_prop = 0.2,
    loading_sd = 0.7,
    normalize = TRUE,
    scale_data = TRUE,
    seed = 1
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required to build the toy Seurat object. Install Seurat or skip this example.")
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required to create a sparse count matrix.")
  }

  set.seed(seed)

  # 1. Metadata

  n_types <- length(cell_types)
  n_donors <- 2 * donors_per_group
  donors <- paste0("D", seq_len(n_donors))

  donor_disease <- rep(c("Control", "AD"), each = donors_per_group)
  names(donor_disease) <- donors

  batch_levels <- paste0("B", seq_len(n_batch))
  batch_by_donor <- sample(batch_levels, size = n_donors, replace = TRUE)
  names(batch_by_donor) <- donors

  donor_age <- round(runif(n_donors, 65, 90))
  names(donor_age) <- donors

  n_cells <- n_donors * n_types * cells_per_donor_per_type

  meta <- data.frame(
    cell_id = paste0("cell", seq_len(n_cells)),
    donor = rep(donors, each = n_types * cells_per_donor_per_type),
    cell_type = rep(rep(cell_types, each = cells_per_donor_per_type), times = n_donors),
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

  # This is the X in G = X B0 + A B1 + E.
  # diseaseAD coefficient defines true disease-associated genes.

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

  # Baseline gene abundance through intercept.
  # This makes genes have heterogeneous baseline expression.
  baseline_log_abundance <- rnorm(n_genes, mean = 0, sd = 1)
  B0["(Intercept)", ] <- baseline_log_abundance

  # Small age effects.
  if ("age" %in% rownames(B0)) {
    B0["age", ] <- rnorm(n_genes, mean = 0, sd = age_effect_sd)
  }

  # Batch effects.
  batch_rows <- grep("^batch", rownames(B0), value = TRUE)
  if (length(batch_rows) > 0L) {
    for (rr in batch_rows) {
      B0[rr, ] <- rnorm(n_genes, mean = 0, sd = batch_effect_sd)
    }
  }

  # Cell-type effects.
  cell_type_rows <- grep("^cell_type", rownames(B0), value = TRUE)
  if (length(cell_type_rows) > 0L) {
    for (rr in cell_type_rows) {
      B0[rr, ] <- rnorm(n_genes, mean = 0, sd = cell_type_effect_sd)
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

  n_loading_genes <- max(1L, floor(loading_prop * n_genes))

  for (k in seq_len(r)) {
    idx <- sample(seq_len(n_genes), n_loading_genes)
    B1[k, idx] <- stats::rnorm(length(idx), mean = 0, sd = loading_sd)
  }

  # 5. Gaussian latent expression layer

  # G = X B0 + A B1 + E

  E <- matrix(
    stats::rnorm(n_cells * n_genes, mean = 0, sd = sigma),
    nrow = n_cells,
    ncol = n_genes
  )

  G <- X %*% B0 + A %*% B1 + E
  rownames(G) <- meta$cell_id
  colnames(G) <- genes

  # 6. Library sizes and multinomial counts

  library_size <- round(
    exp(stats::rnorm(n_cells, mean = log(mean_umi), sd = lib_size_sd))
  )

  library_size[library_size < 1L] <- 1L

  counts <- matrix(
    0L,
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(genes, meta$cell_id)
  )

  for (i in seq_len(n_cells)) {
    # Stable softmax.
    eta <- G[i, ]
    eta <- eta - max(eta)

    prob <- exp(eta)
    prob <- prob / sum(prob)

    counts[, i] <- as.integer(
      stats::rmultinom(
        n = 1,
        size = library_size[i],
        prob = prob
      )
    )
  }

  counts <- Matrix::Matrix(counts, sparse = TRUE)

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
  seu$library_size_true <- library_size

  Seurat::Idents(seu) <- "disease"

  if (isTRUE(normalize)) {
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  }

  if (isTRUE(scale_data)) {
    seu <- Seurat::ScaleData(seu, features = rownames(seu), verbose = FALSE)
  }

  # 8. Store ground truth

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
  seu@misc$library_size_true <- library_size
  seu@misc$simulation_model <- "G = X B0 + A B1 + E; counts_i ~ Multinomial(N_i, softmax(G_i))"

  seu
}

Seurat_toy = make_toy_seurat_ad(n_genes = 100, 
                                donors_per_group = 6,
                                n_signals = 50,
                                # fold-change on relative abundance scale
                                signal_strength = 2,
                                # latent dimension
                                r = 5,
                                seed = 1)
```

# 3. Add the cellular detection rate (CDR) as a latent variable in the toy dataset.

In addition to `batch`, `cell_type`, and `age`, we compute and include
the cellular detection rate (CDR) as a latent variable. Here, `CDR`
summarizes the overall detection level of each cell and is commonly used
to adjust for cell-level technical variation in single-cell RNA-seq
data.

``` r
feature.names <- rownames(Seurat_toy)

np_data.count <- Seurat::GetAssayData(object = Seurat_toy, layer = "count")

np_data.matrix.exp <- Matrix::t(np_data.count > 0) # indicator matrix for expressed genes,
np_data.matrix.exp.count <- Matrix::colSums(np_data.matrix.exp)

# Calculate cellular detection rate (CDR)
np_data.matrix <- np_data.count # do not transpose here. Keep the genes as rows.

CDR <- Matrix::rowMeans(np_data.matrix.exp) # expressed genes in each CELL
Seurat_toy$CDR <- CDR
```

# 4. Select the significant genes using `full_process`

We then apply `full_process()` to identify genes associated with disease
status. In this example, we compare cells from the AD group against
cells from the Control group by setting `ident.1 = "AD"` and
`ident.2 = "Control"`. Setting `PC = NULL` allows the number of
principal components to be estimated automatically using the bulk
eigenvalue matching analysis (BEMA) method (Ke et al. 2023). The
covariates `batch`, `cell_type`, `age`, and `CDR` are included in both
`latent.vars_imp` and `latent.vars_comp` to adjust for potential
confounding effects during the imputation and comparison steps.

We use `test.use = "LCD"` to construct feature-importance statistics
based on the Lasso coefficient difference statistic, and `m = 1`
specifies that a single knockoff copy is generated for each gene.

``` r
result_LCD = 
  full_process(Seurat_data = Seurat_toy,
               PC = NULL,
               ident.1 = "AD",
               ident.2 = "Control",
               latent.vars_imp = c("batch", "cell_type", "age", "CDR"),
               latent.vars_comp = c("batch", "cell_type", "age", "CDR"),
               test.use = "LCD",
               m = 1)
```

We then evaluate the false discovery proportion (FDP) and power of the
selected genes.

``` r
true_signal <- Seurat_toy@misc$true_signal

print("FDP and Power")
#> [1] "FDP and Power"
print((length(result_LCD$selected) - sum(result_LCD$selected %in% true_signal))/length(result_LCD$selected)) # FDR
#> [1] 0
print(sum(result_LCD$selected %in% true_signal)/length(true_signal)) # Power
#> [1] 0.94
```

The feature-importance statistics can also be constructed using other
testing methods supported by `Seurat::FindMarkers()`. In addition,
setting `m > 1` enables the use of multiple knockoffs. In the following
example, we use `test.use = "MAST"` and generate five knockoff copies
for each gene by setting `m = 5`. Here, we set `PC = 5`, corresponding
to the true number of latent factors used to generate the toy dataset.

``` r
result_MAST = 
  full_process(Seurat_data = Seurat_toy,
               PC = 5,
               ident.1 = "AD",
               ident.2 = "Control",
               latent.vars_imp = c("batch", "cell_type", "age", "CDR"),
               latent.vars_comp = c("batch", "cell_type", "age", "CDR"),
               test.use = "MAST",
               m = 5)
```

We then evaluate the false discovery proportion (FDP) and power of the
selected genes.

``` r
print("FDP and Power")
#> [1] "FDP and Power"
print((length(result_MAST$selected) - sum(result_MAST$selected %in% true_signal))/length(result_MAST$selected)) # FDR
#> [1] 0
print(sum(result_MAST$selected %in% true_signal)/length(true_signal)) # Power
#> [1] 0.94
```

The ground-truth signal genes and the genes selected by LCD and MAST are
shown below:

``` r
cat("Ground-truth signal genes:\n")
#> Ground-truth signal genes:
cat(Seurat_toy@misc$true_signal_genes, sep = ", ")
#> gene2, gene4, gene5, gene6, gene7, gene8, gene9, gene10, gene11, gene12, gene13, gene16, gene17, gene19, gene21, gene22, gene23, gene24, gene29, gene33, gene34, gene39, gene41, gene42, gene43, gene46, gene58, gene60, gene63, gene65, gene67, gene68, gene69, gene70, gene72, gene74, gene75, gene78, gene79, gene81, gene82, gene85, gene86, gene88, gene90, gene92, gene93, gene95, gene97, gene100
cat("\n\n")

cat("Genes selected by LCD (m=1):\n")
#> Genes selected by LCD (m=1):
cat(result_LCD$variables_name, sep = ", ")
#> gene2, gene4, gene5, gene6, gene8, gene9, gene10, gene11, gene12, gene13, gene16, gene17, gene19, gene21, gene22, gene23, gene24, gene29, gene33, gene39, gene41, gene42, gene43, gene46, gene58, gene60, gene63, gene65, gene67, gene69, gene70, gene72, gene74, gene75, gene78, gene79, gene81, gene82, gene85, gene86, gene88, gene90, gene92, gene93, gene95, gene97, gene100
cat("\n\n")

cat("Genes selected by MAST (m=5):\n")
#> Genes selected by MAST (m=5):
cat(result_MAST$variables_name, sep = ", ")
#> gene2, gene4, gene5, gene6, gene8, gene9, gene10, gene11, gene12, gene13, gene16, gene17, gene19, gene21, gene22, gene23, gene24, gene29, gene33, gene34, gene39, gene42, gene43, gene46, gene58, gene60, gene63, gene65, gene67, gene69, gene70, gene72, gene74, gene75, gene78, gene79, gene81, gene82, gene85, gene86, gene88, gene90, gene92, gene93, gene95, gene97, gene100
cat("\n")
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ke2023" class="csl-entry">

Ke, Zheng Tracy, Yucong Ma, and Xihong Lin. 2023. “Estimation of the
Number of Spiked Eigenvalues in a Covariance Matrix by Bulk Eigenvalue
Matching Analysis.” *Journal of the American Statistical Association*
118 (541): 374–92.

</div>

</div>

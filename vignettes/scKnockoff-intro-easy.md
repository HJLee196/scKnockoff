Introduction to scKnockoff: Integrated Workflow
================

This vignette demonstrates a simplified workflow for applying
`scKnockoff` to a toy single-cell RNA-seq dataset. The main analysis
steps, including imputation, knockoff generation, feature-importance
calculation, and gene selection, are combined into a single function
call using `full_process()`. In addition, when `PC = NULL`,
`full_process()` automatically estimates the number of latent factors
using the Bulk Eigenvalue Matching Analysis (BEMA) method proposed by
(Ke et al. 2023).

# 1. Load scKnockoff

``` r
# if (!requireNamespace("scKnockoff", quietly = TRUE)) {
#   devtools::install_github("HJLee196/scKnockoff")
# }

library(scKnockoff)
```

# 2. Create Dataset

For convenience, the package provides `make_toy_seurat_ad()`, which
generates a toy Seurat object under a Gaussian latent-factor model. The
synthetic single-cell RNA-seq dataset is generated from

$$G = X B_0 + A B_1 + E,$$

where $G$ represents the latent Gaussian log-expression matrix, $X$
contains observed covariates such as disease status, batch, age, and
cell type, $A$ contains Gaussian latent factors, and $E$ is Gaussian
noise. Disease status is assigned at the donor level, and nonzero
disease effects are assigned to a subset of genes in $B_0$, defining the
ground-truth signals. To construct a Seurat object, the latent log-mean
expression matrix $G$ is converted to a positive mean expression matrix
by setting $\mu_{ij} = \exp(G_{ij})$, and observed UMI-like counts are
generated from a Poisson model.

``` r

Seurat_toy = make_toy_seurat_ad(n_genes = 200, 
                                donors_per_group = 6,
                                n_signals = 50,
                                # fold-change on relative abundance scale
                                signal_strength = 2,
                                # latent dimension
                                r = 5,
                                # standard deviation of noise
                                sigma = 2,
                                seed = 2)
```

# 3. Add the Cellular Detection Rate (CDR)

In addition to `batch`, `cell_type`, and `age`, we compute and include
the cellular detection rate (CDR) as an observed covariate. Here, `CDR`
summarizes the overall detection level of each cell and is commonly used
to adjust for cell-level technical variation in single-cell RNA-seq
data.

``` r
feature.names <- rownames(Seurat_toy)

np_data.count <- Seurat::GetAssayData(object = Seurat_toy, layer = "count")

# Cell-by-gene indicator matrix: 1 if a gene is expressed in a cell.
np_data.matrix.exp <- Matrix::t(np_data.count > 0)
# Number of cells in which each gene is expressed.
np_data.matrix.exp.count <- Matrix::colSums(np_data.matrix.exp)

# Calculate cellular detection rate (CDR)
# Fraction of expressed genes in each cell.
CDR <- Matrix::rowMeans(np_data.matrix.exp)

Seurat_toy$CDR <- CDR
```

# 4. Select Significant Genes Using `full_process`

We then apply `full_process()` to identify genes associated with disease
status. In this example, we compare cells from the AD group against
cells from the Control group by setting `ident.1 = "AD"` and
`ident.2 = "Control"`. Setting `PC = NULL` allows the number of
principal components to be estimated automatically using the BEMA method
(Ke et al. 2023). The covariates `batch`, `cell_type`, `age`, and `CDR`
are included in both `latent.vars_imp` and `latent.vars_comp` to adjust
for potential confounding effects during the imputation and comparison
steps.

We use `test.use = "LCD"` to construct feature-importance statistics
based on the lasso coefficient difference statistic, and `m = 1`
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

compute_fdp <- function(selected, true_signal) {
  if (length(selected) == 0) return(0)
  (length(selected) - sum(selected %in% true_signal)) / length(selected)
}

compute_power <- function(selected, true_signal) {
  sum(selected %in% true_signal) / length(true_signal)
}

lcd_results <- data.frame(
  Method = "LCD (m = 1)",
  FDP = compute_fdp(result_LCD$selected, true_signal),
  Power = compute_power(result_LCD$selected, true_signal)
)

lcd_results
#>        Method       FDP Power
#> 1 LCD (m = 1) 0.1403509  0.98
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

We again evaluate the false discovery proportion (FDP) and power of the
selected genes.

``` r
results <- data.frame(
  Method = c("LCD (m = 1)",
             "MAST (m = 5)"),
  FDP = c(compute_fdp(result_LCD$selected, true_signal),
          compute_fdp(result_MAST$selected, true_signal)),
  Power = c(compute_power(result_LCD$selected, true_signal),
            compute_power(result_MAST$selected, true_signal))
)

results
#>         Method       FDP Power
#> 1  LCD (m = 1) 0.1403509  0.98
#> 2 MAST (m = 5) 0.0000000  0.80
```

The ground-truth signal genes and the genes selected by LCD and MAST are
shown below:

``` r
cat("Ground-truth signal genes:\n")
#> Ground-truth signal genes:
cat(Seurat_toy@misc$true_signal_genes, sep = ", ")
#> gene3, gene8, gene17, gene18, gene22, gene28, gene29, gene31, gene34, gene39, gene41, gene44, gene46, gene48, gene61, gene67, gene68, gene71, gene78, gene84, gene90, gene91, gene92, gene96, gene98, gene99, gene101, gene104, gene106, gene112, gene117, gene119, gene126, gene128, gene130, gene138, gene143, gene150, gene154, gene159, gene162, gene163, gene164, gene173, gene180, gene183, gene185, gene188, gene191, gene200
cat("\n\n")

cat("Genes selected by LCD (m = 1):\n")
#> Genes selected by LCD (m = 1):
cat(result_LCD$variables_name, sep = ", ")
#> gene2, gene3, gene6, gene7, gene8, gene17, gene18, gene22, gene28, gene29, gene31, gene34, gene39, gene41, gene44, gene46, gene48, gene61, gene67, gene68, gene71, gene78, gene84, gene89, gene90, gene91, gene92, gene96, gene98, gene99, gene101, gene104, gene106, gene109, gene111, gene117, gene118, gene119, gene126, gene128, gene130, gene138, gene143, gene150, gene152, gene154, gene159, gene162, gene163, gene164, gene173, gene180, gene183, gene185, gene188, gene191, gene200
cat("\n\n")

cat("Genes selected by MAST (m = 5):\n")
#> Genes selected by MAST (m = 5):
cat(result_MAST$variables_name, sep = ", ")
#> gene8, gene17, gene18, gene28, gene31, gene34, gene39, gene44, gene46, gene61, gene67, gene68, gene78, gene84, gene90, gene91, gene98, gene99, gene101, gene104, gene106, gene117, gene119, gene126, gene128, gene130, gene138, gene143, gene150, gene154, gene159, gene162, gene163, gene164, gene173, gene180, gene183, gene185, gene188, gene200
cat("\n")
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ke2023" class="csl-entry">

Ke, Zheng Tracy, Yucong Ma, and Xihong Lin. 2023. “Estimation of the
Number of Spiked Eigenvalues in a Covariance Matrix by Bulk Eigenvalue
Matching Analysis.” *Journal of the American Statistical Association*
118 (541): 374–92.

</div>

</div>

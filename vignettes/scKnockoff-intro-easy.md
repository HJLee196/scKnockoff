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
#> [1] 0.84
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
#> [1] 0.86
```

The ground-truth signal genes and the genes selected by LCD and MAST are
shown below:

``` r
cat("Ground-truth signal genes:\n")
#> Ground-truth signal genes:
cat(Seurat_toy@misc$true_signal_genes, sep = ", ")
#> gene2, gene4, gene5, gene6, gene7, gene8, gene9, gene10, gene11, gene12, gene13, gene16, gene17, gene19, gene21, gene22, gene24, gene29, gene33, gene34, gene39, gene41, gene42, gene43, gene44, gene46, gene57, gene58, gene60, gene65, gene66, gene69, gene70, gene72, gene74, gene75, gene76, gene77, gene79, gene81, gene84, gene86, gene88, gene89, gene90, gene93, gene94, gene95, gene97, gene100
cat("\n\n")

cat("Genes selected by LCD (m=1):\n")
#> Genes selected by LCD (m=1):
cat(result_LCD$variables_name, sep = ", ")
#> gene2, gene4, gene5, gene6, gene8, gene9, gene10, gene11, gene12, gene16, gene17, gene19, gene21, gene22, gene24, gene29, gene33, gene34, gene39, gene42, gene43, gene44, gene46, gene57, gene58, gene60, gene65, gene66, gene69, gene70, gene72, gene76, gene79, gene81, gene84, gene86, gene89, gene93, gene94, gene95, gene97, gene100
cat("\n\n")

cat("Genes selected by MAST (m=5):\n")
#> Genes selected by MAST (m=5):
cat(result_MAST$variables_name, sep = ", ")
#> gene2, gene4, gene5, gene6, gene8, gene9, gene10, gene11, gene12, gene16, gene17, gene19, gene21, gene22, gene24, gene29, gene33, gene34, gene39, gene41, gene42, gene46, gene57, gene58, gene60, gene65, gene66, gene69, gene70, gene72, gene74, gene75, gene76, gene79, gene81, gene84, gene86, gene90, gene93, gene94, gene95, gene97, gene100
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

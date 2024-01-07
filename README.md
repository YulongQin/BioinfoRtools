
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BioinfoRtools

<!-- badges: start -->
<!-- badges: end -->

The BioinfoRtools is A comprehensive bioinformatics R package toolbox

## Installation

You can install the development version of BioinfoRtools like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
devtools::install_github("YulongQin/BioinfoRtools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BioinfoRtools)
#> 
#> 载入需要的程辑包：Biobase
#> 载入需要的程辑包：BiocGenerics
#> 
#> 载入程辑包：'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
#>     table, tapply, union, unique, unsplit, which.max, which.min
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> 载入程辑包：'DynDoc'
#> The following object is masked from 'package:BiocGenerics':
#> 
#>     path
#> 
## basic example code
# ?DEGs_DESeq2
# ?DEGs_Wilcoxon
# data(countdata)
# data(coldata)
res <- DEGs_DESeq2(
  countdata=countdata,
  group=coldata$condition,
  case = "treated",
  ctrl = "untreated",
  covariate = NULL,
  updown_thres = 1,
  rate_sap = 0.1,
  OUT_df_index = FALSE,
  updown_index = FALSE,
  p_value = 0.05,
  padj_value = 0.05,
  grp_nm = "DEGs_DESeq2",
  dir_nm = "DEGs_DESeq2"
)
#> Warning in dir.create(output_dir, recursive = T):
#> '.\outputdata\DEGs_DESeq2\DEGs_DESeq2'已存在
#> converting counts to integer mode
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
```

Here’s how to use some of the graphing functions`PCA` and `HTmap`:

``` r
PCA(countdata=countdata,
    group=coldata$condition, 
    grp_nm = "PCA", 
    scale = T)
#> Too few points to calculate an ellipse
```

<img src="man/figures/README-PCA-1.png" width="100%" />

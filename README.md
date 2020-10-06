
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DELocal

<!-- badges: start -->

<!-- badges: end -->

The goal of DELocal is to …

![neighbor](Fig_2_solid.svg.png)

## Installation

You can install the released version of DELocal from
[CRAN](https://CRAN.R-project.org) with:

``` r
devtools::install_github("dasroy/delocal")
```

## Example

This is a basic example which shows you how to solve a common problem:

    #> class: DESeqDataSet 
    #> dim: 52183 14 
    #> metadata(1): version
    #> assays(1): counts
    #> rownames(52183): ENSMUSG00000000001 ENSMUSG00000000003 ...
    #>   ENSMUSG00000114967 ENSMUSG00000114968
    #> rowData names(0):
    #> colnames(14): ME14.E1M1R ME14.E2M1R ... ME13.E9M1R ME13.EXM1L
    #> colData names(1): condition

## Getting gene chromosomal location

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

    #> Loading required package: dplyr
    #> 
    #> Attaching package: 'dplyr'
    #> The following object is masked from 'package:matrixStats':
    #> 
    #>     count
    #> The following object is masked from 'package:Biobase':
    #> 
    #>     combine
    #> The following objects are masked from 'package:GenomicRanges':
    #> 
    #>     intersect, setdiff, union
    #> The following object is masked from 'package:GenomeInfoDb':
    #> 
    #>     intersect
    #> The following objects are masked from 'package:IRanges':
    #> 
    #>     collapse, desc, intersect, setdiff, slice, union
    #> The following objects are masked from 'package:S4Vectors':
    #> 
    #>     first, intersect, rename, setdiff, setequal, union
    #> The following objects are masked from 'package:BiocGenerics':
    #> 
    #>     combine, intersect, setdiff, union
    #> The following objects are masked from 'package:stats':
    #> 
    #>     filter, lag
    #> The following objects are masked from 'package:base':
    #> 
    #>     intersect, setdiff, setequal, union
    #> Loading required package: gtools
    #> Loading required package: limma
    #> 
    #> Attaching package: 'limma'
    #> The following object is masked from 'package:DESeq2':
    #> 
    #>     plotMA
    #> The following object is masked from 'package:BiocGenerics':
    #> 
    #>     plotMA
    #>                        logFC       AveExpr         t      P.Value    adj.P.Val
    #> ENSMUSG00000032334 1891.3317  1.299278e-13  22.33965 2.656896e-11 9.598419e-07
    #> ENSMUSG00000037217  509.6374  4.060244e-15  21.73943 3.681434e-11 9.598419e-07
    #> ENSMUSG00000015501 1256.9704 -4.060244e-15  19.61168 1.258837e-10 1.828783e-06
    #> ENSMUSG00000029468  131.1903 -3.552714e-15  19.43399 1.402845e-10 1.828783e-06
    #> ENSMUSG00000021196  451.0717  1.624098e-14  17.43213 5.090738e-10 4.545027e-06
    #> ENSMUSG00000049420 -428.4620 -1.421085e-14 -17.39244 5.229679e-10 4.545027e-06
    #>                           B
    #> ENSMUSG00000032334 14.80633
    #> ENSMUSG00000037217 14.59157
    #> ENSMUSG00000015501 13.74295
    #> ENSMUSG00000029468 13.66530
    #> ENSMUSG00000021196 12.70721
    #> ENSMUSG00000049420 12.68654

## Trying with summerizedexperiment

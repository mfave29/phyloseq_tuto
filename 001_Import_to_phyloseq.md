R Notebook
================

  - [Imports de donnée de Dada2 vers
    PhyloSeq](#imports-de-donnée-de-dada2-vers-phyloseq)
      - [Chargement du “Dada2\_tutorial”](#chargement-du-dada2_tutorial)
      - [Installation de BiocManager et de
        phyloseq](#installation-de-biocmanager-et-de-phyloseq)
      - [Vérification de la version de
        phyloseq](#vérification-de-la-version-de-phyloseq)
      - [Vérification de la version de
        Biostrings](#vérification-de-la-version-de-biostrings)
      - [Vérification et importation de
        ggplot2](#vérification-et-importation-de-ggplot2)
  - [Vizualize alpha-diversity](#vizualize-alpha-diversity)
  - [Ordinate](#ordinate)

# Imports de donnée de Dada2 vers PhyloSeq

## Chargement du “Dada2\_tutorial”

``` r
load("Dada2_tutorial_FinalEnv")
```

## Installation de BiocManager et de phyloseq

``` r
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
```

    ## Loading required namespace: BiocManager

``` r
BiocManager::install("phyloseq")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phyloseq'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'isoband', 'lme4'

## Vérification de la version de phyloseq

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.34.0'

## Vérification de la version de Biostrings

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## [1] '2.58.0'

## Vérification et importation de ggplot2

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.2'

``` r
theme_set(theme_bw())
```

``` r
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

``` r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
```

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 232 taxa by 6 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 232 reference sequences ]

# Vizualize alpha-diversity

``` r
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](001_Import_to_phyloseq_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
No difference

# Ordinate

``` r
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.08574537 
    ## Run 1 stress 0.0894286 
    ## Run 2 stress 0.08942862 
    ## Run 3 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04283619  max resid 0.1433491 
    ## Run 4 stress 0.09421602 
    ## Run 5 stress 0.0894286 
    ## Run 6 stress 0.08002299 
    ## ... Procrustes: rmse 4.589663e-06  max resid 1.563771e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 1.683942e-06  max resid 4.40338e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.08002299 
    ## ... Procrustes: rmse 3.523666e-06  max resid 9.695695e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.08942864 
    ## Run 10 stress 0.08002299 
    ## ... Procrustes: rmse 3.515738e-06  max resid 1.078929e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.08942864 
    ## Run 12 stress 0.08942866 
    ## Run 13 stress 0.1216669 
    ## Run 14 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 1.417483e-06  max resid 4.263064e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.09421602 
    ## Run 16 stress 0.1323637 
    ## Run 17 stress 0.08942869 
    ## Run 18 stress 0.1233096 
    ## Run 19 stress 0.09421601 
    ## Run 20 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 3.055629e-06  max resid 9.255715e-06 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```

![](001_Import_to_phyloseq_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Clear separation between the early and late samples

Bar plot

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```

![](001_Import_to_phyloseq_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

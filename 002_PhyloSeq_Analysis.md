Analyse - Phyloseq
================

  - [PhyloSeq](#phyloseq)
      - [Chargement des données à partir d’un jeu de
        données](#chargement-des-données-à-partir-dun-jeu-de-données)
      - [Filtration taxonomique](#filtration-taxonomique)
      - [Création de table, traits de chaque
        phylum](#création-de-table-traits-de-chaque-phylum)
      - [Exploration des prévalences des caractéristiques, soit le
        nombre d’échantillons dans lesquels un taxon apparaît au moins
        une
        fois.](#exploration-des-prévalences-des-caractéristiques-soit-le-nombre-déchantillons-dans-lesquels-un-taxon-apparaît-au-moins-une-fois.)
      - [Filtrage](#filtrage)
      - [Filtrage de prévalence](#filtrage-de-prévalence)
      - [Définition du seuil de prévalence à 5 % des échantillons
        totaux.](#définition-du-seuil-de-prévalence-à-5-des-échantillons-totaux.)
      - [Agglomération](#agglomération)
      - [Combine all features that descend from the same
        genus](#combine-all-features-that-descend-from-the-same-genus)
      - [Installation du package gridExtra pour avoir grid.arrange
        permettant le
        regroupement](#installation-du-package-gridextra-pour-avoir-grid.arrange-permettant-le-regroupement)
      - [Transformation de la valeur
        d’abondance](#transformation-de-la-valeur-dabondance)
      - [Transformation de l’abondance
        relative](#transformation-de-labondance-relative)
  - [Sous-ensemble par la taxonomie](#sous-ensemble-par-la-taxonomie)
      - [Installation des packages pour une analyse statistique plus
        approfondie](#installation-des-packages-pour-une-analyse-statistique-plus-approfondie)
      - [Prétraitement - Âge des
        souris](#prétraitement---âge-des-souris)
      - [Analyse en PCoA avec
        Bray-Curtis](#analyse-en-pcoa-avec-bray-curtis)
  - [Projections d’ordination](#projections-dordination)
      - [Les positions d’échantillon produites par un PCoA en utilisant
        Unifrac
        pondéré](#les-positions-déchantillon-produites-par-un-pcoa-en-utilisant-unifrac-pondéré)
  - [PCA sur les rangs](#pca-sur-les-rangs)
      - [On peut maintenant faire une
        PCA.](#on-peut-maintenant-faire-une-pca.)
  - [Supervised learning](#supervised-learning)
  - [Analyses basées sur des
    graphiques](#analyses-basées-sur-des-graphiques)
      - [Création et traçage des
        graphiques](#création-et-traçage-des-graphiques)
  - [Tests à deux échantillons basés sur des
    graphiques](#tests-à-deux-échantillons-basés-sur-des-graphiques)
      - [Minimum Spanning Tree (MST)](#minimum-spanning-tree-mst)
      - [Voisins les plus proches](#voisins-les-plus-proches)
  - [Modélisation linéaire](#modélisation-linéaire)
  - [Tests multiples hiérarchiques](#tests-multiples-hiérarchiques)
  - [Techniques polyvalentes](#techniques-polyvalentes)
      - [Récupération et filtrage des
        données](#récupération-et-filtrage-des-données)
  - [Conclusion](#conclusion)

# PhyloSeq

Amplicon bioinformatics : from raw reads to tables

## Chargement des données à partir d’un jeu de données

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

    ## Loading required package: phyloseq

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

Il y a 389 taxons et 360 échantillons, se regroupant en 6 rangs
taxonomiques.

## Filtration taxonomique

``` r
# Show available ranks in the dataset
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

On voit les noms du classement des espèces.

## Création de table, traits de chaque phylum

``` r
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

On voit le nombre des taxons. Les non-identifiés doivent être éliminés.
On a 327 taxons différents chez les Firmicutes.

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

## Exploration des prévalences des caractéristiques, soit le nombre d’échantillons dans lesquels un taxon apparaît au moins une fois.

``` r
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

Nous voyons la prévalence des phyla. Deinococcus Thermus est apparu dans
un peu plus de un pour cent des échantillons et Fusobactéries dans
seulement 2 échantillons. Il faut donc les filtrer. Pour donner un
exemple, on a 1562 séquences d’Actinobacteries et environ 120 types
d’Actinobacteries différentes.

## Filtrage

``` r
# Define phyla to filter
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

Il s’agit d’une étape de filtrage, où l’on a supprimé de l’analyse
Deinococcus thermus et les Fusobacteria.

## Filtrage de prévalence

Non supervisée. Cette étape de filtrage peut être appliquée même dans
des contextes où l’annotation taxinomique n’est pas disponible ou n’est
pas fiable.

``` r
library(ggplot2)
```

``` r
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Figure : prévalence des taxons en rapport avec le total. Nous voyons que
la prévalence des Firmicutes est très importante. Aucune séparation
naturelle n’est immédiatement évidente. Ce serait bien de définir un
seuil de prévalence dans un intervalle entre 0 et 10%. C’est pour cela
qu’on va définir le seuil de prévalence à 5% afin de ne garder que ceux
ayant une prévalence supérieure.

## Définition du seuil de prévalence à 5 % des échantillons totaux.

Nous fixons le seuil de prévalence à 5%

``` r
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

Nous gardons les taxa dont la prévalence est supérieure au seuil de
prévalence de 0.05

## Agglomération

Il est utile d’agglomérer les caractéristiques correspondant aux données
étroitement en lien car on sait qu’il y a beaucoup d’espèces ou de
sous-espèces redondantes d’un point de vue fonctionnel dans une
communauté. L’agglomération taxonomique a l’avantage d’être plus facile
pour définir en avance.

## Combine all features that descend from the same genus

Combien de genres seront présents après une étape de filtration ?

``` r
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

Il y a 49 genres qui seront présents après le filtrage.

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

Définition d’une distance phylogénétique de 0.4 pour analyse.

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

plot\_tree de phyloseq permet de comparer les données originales non
filtrées, l’arbre avant et l’arbre après l’agglomération.

``` r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```

## Installation du package gridExtra pour avoir grid.arrange permettant le regroupement

``` r
install.packages("gridExtra")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
library(gridExtra)
```

``` r
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Sur la gauche, on voit l’arbre original, au milieu celui des rangs des
genres et à droite l’agglomération phylogénétique avec la distance
phylogénétique de 0.4. On peut constater que l’arbre original contenait
plus de branches que les deux autres.

## Transformation de la valeur d’abondance

Transformation des données du microbiome pour tenir compte de différents
paramètres (taille de la library, la variance, l’échelle…). On veut un
graphique d’abondance relative.

``` r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

## Transformation de l’abondance relative

La transformation dans ce cas convertit les dénombrements de chaque
échantillon en leurs fréquences, souvent appelées proportions ou
abondances relatives.

``` r
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

``` r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

On obtient les valeurs d’abondance avant et après la transformation. Les
abondance originales sont en haut et les abondances relatives en bas. On
voit les abondances relatives transformées selon le sexe.

# Sous-ensemble par la taxonomie

Nous spécifions un rang taxonomique plus précis afin d’expliquer
l’abondance bimodale précédente.

``` r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Nous pouvons voir les abondances relatives des Lactobacillales selon le
sexe et le genre d’hôte. Nous voyons une abondance relative des
Lactobacillus plus importante que celle des Streptococcus. La
distribution bimodale précédente était en fait du à un mélange de deux
genres dans les Lactobacillales avec les Lactobacillus et les
Streptococcus.

## Installation des packages pour une analyse statistique plus approfondie

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}
```

    ## Skipping install of 'phyloseqGraphTest' from a github remote, the SHA1 (3fb6c274) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
.inst <- .bioc_packages %in% installed.packages()
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

## Prétraitement - Âge des souris

``` r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Histogramme des groupes d’âge, 3 parties distinctes

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Histogramme comparant les profondeurs de lecture brutes et transformées.
La transformation pourrait être suffisante par normalisation des données
d’abondance. En fait, elle ne l’est pas et donc, nous faisons une
transformation de stabilisation de la variance.

## Analyse en PCoA avec Bray-Curtis

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGCTTTGCAAGTCAGAAGTGAAATCCATGGGCTTAACCCATGAACTGCTTTTGAAACTGCAGAGCTTGAGTGGAGTAGAGGTAGGCGGAATTCCCGGTGTAGCGGTGAAATGCGTAGAGATCGGGAGGAACACCAGTGGCGAAGGCGGCCTGCTGGGCTCTAACTGACGCTGAGGCACGAAAGCGTGGGTAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

Analyse d’ordination avec le log d’abondance. On retrouve des
aberrances. Il s’agit des échantillons des femelles 5 et 6 au jour 165
et des échantillons des mâles 3, 4, 5 et 6 au jour 175. Nous les
retirerons, car nous nous intéressons principalement aux relations entre
les non-aberrants.

Vérification des valeur aberrantes des femelles.

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

Les échantillons aberrants sont dominés par un seul ASV (Amplicon
Sequence Variant).

# Projections d’ordination

Calcul des ordinations avec les valeurs aberrantes supprimées.

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

Suppression des échantillons avec moins de 1000 lectures

``` r
which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

Nous avons donc les différents échantillons contnenant moins de 1000
reads.

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

PCoA en utilisant la dissemblance de Bray-Curtis

``` r
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

PcoA - Bray-Curtis On voit qu’il y a un effet d’âge important qui est
cohérent avec le reste et que les litters ne sont à priori pas
discriminatoires, du moins il y a un effet de l’âge plus important.

``` r
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Figure : Un graphique DPCoA incorpore des informations phylogénétiques,
mais est dominé par le premier axe. La variabilité en plus de l’âge est
la composition phylogénétique de chaque échantillon. Le premier axe
explique 75% de la variabilité, environ 9 fois celle du deuxième axe;
cela se traduit par la forme allongée du tracé d’ordination.

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Graphique : taxons responsables des axes 1 et 2 C’est une analyse en
fonction des phylum.

Ainsi, la différence de similarité est due à l’âge et aux phylums.

## Les positions d’échantillon produites par un PCoA en utilisant Unifrac pondéré

``` r
library(grid)
```

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTGTCCGGATTCACTGGGCGTAAAGAGCACGTAGGCGGTTATTTAAGTCAGGTGTGAAAGTTTTCGGCTCAACCGGAAAAGTGCACTTGAAACTGGATAACTTGAGAATCGGAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAAGAACACCGGTGGCGAAGGCGGCTTACTGGACGATTACTGACGCTGAGGTGCGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

Utilisé UniFrac pondéré permet de mesurer quantitativement la
dissimilitude communautaire qui intègre les relations phylogénétiques
entre les échantillons.Il y a réellement un effet d’âge.

# PCA sur les rangs

On ignore les abondances brutes et on remplace avec des rangs car il
peut être difficile d’identifier une transformation qui ramène les
données à la normalité. On utilise une version transformée par rang de
données pour effectuer l’ACP. Création matrice représentant les
abondances par rang où le microbe avec le plus petit dans un échantillon
est mappé au rang 1, au deuxième rang le plus petit 2, etc.

``` r
library(phyloseq)
```

``` r
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

``` r
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

Ici, on a fait le choix où les microbes dont le rang est inférieur à un
certain seuil doivent être liés à 1.

``` r
library(ggplot2)
```

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

Graphique : transformation du seuil de rang. On visualise les abondances
selon le rang du seuil.

## On peut maintenant faire une PCA.

``` r
library(ade4)
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
col_scores <- col_scores %>%
  left_join(tax)
```

    ## Joining, by = "seq"

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

Les résultats sont similaires aux analyses PCoA sans transformation,
ainsi nous avons une bonne confiance de l’analyse de nos données. On
voit donc que les triangles (mélange de données) sont rassemblées et que
les ronds sont dispersés. La vraie variance ici est donc l’ordre plutôt
que la litter par exemple. \# Correspondance canonique CCpnA :
ordination d’une espèce par table d’échantillons. Le but de la création
de biplots est de déterminer quels types de communautés bactériennes
sont les plus importants dans différents types d’échantillons de souris.
Il peut être plus facile d’interpréter ces biplots lorsque l’ordre entre
les échantillons reflète les caractéristiques de l’échantillon - les
variations d’âge ou de statut de portée dans les données de souris, par
exemple - et cela est au cœur de la conception de CCpnA.

``` r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

``` r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

# Supervised learning

Nous appliquerons des techniques supervisées pour essayer de prédire
l’âge à partir de la composition du microbiome. Diviser les données
en ensembles d’apprentissage et de test, les affectations étant
effectuées à la souris plutôt que par échantillon, pour garantir que
l’ensemble de test simule de manière réaliste la collecte de nouvelles
données. Une fois que nous avons divisé les données, nous pouvons
utiliser la fonction train pour ajuster le modèle PLS.

``` r
library(caret)
```

    ## Loading required package: lattice

``` r
sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
# take 8 mice at random to be the training set, and the remaining 4 the test set
trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
plsFit <- train(age ~ ., data = training,
                method = "pls", preProc = "center")
```

Nous pouvons prédire les étiquettes de classe sur l’ensemble de test.

``` r
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)
```

    ##            
    ## plsClasses  (0,100] (100,400]
    ##   (0,100]        67         1
    ##   (100,400]       5        45

Prédiction

``` r
library(randomForest)
```

    ## randomForest 4.6-14

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

``` r
rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)
```

    ##            
    ## rfClasses   (0,100] (100,400]
    ##   (0,100]        70         2
    ##   (100,400]       2        44

``` r
library(vegan)
```

    ## Loading required package: permute

    ## This is vegan 2.5-7

    ## 
    ## Attaching package: 'vegan'

    ## The following object is masked from 'package:caret':
    ## 
    ##     tolerance

``` r
pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"

pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],
                                pls_biplot$scores)

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)
```

``` r
ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ age2) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

Figure : PLS produit une représentation biplot conçue pour séparer les
échantillons par une variable de réponse

La projection est choisie avec une référence explicite à la variable
d’âge regroupée. Plus précisément, PLS identifie un sous-espace pour
maximiser la discrimination entre les classes, et le biplot affiche des
projections d’échantillons et des coefficients ASV par rapport à ce
sous-espace.

``` r
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])
ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 1, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

Graphique : distance entre les échantillons. Une distance a été calculée
entre les échantillons en fonction de la fréquence à laquelle
l’échantillon se produit dans la même partition d’arbre dans un
random forest’s bootstrapping. Si une paire d’échantillons se produit
fréquemment dans la même cloison, la paire se voit attribuer une faible
distance.

``` r
as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
```

    ## [1] "Lachnospiraceae" "Roseburia"

``` r
impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund)) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

L’espèce appartenant au genre Rosbeuria est le microorganisme qui a le
plus d’influence dans la prédiction aléatoire de forest. Son abondance
est basse dans les âges de 0 à 100 jours et beaucoup plus élevée dans
les jours entre 100 et 400.

# Analyses basées sur des graphiques

## Création et traçage des graphiques

Création d’un réseau avec la dissimilarité de Jaccard à .35. On ajoute
un attribut aux sommets en indiquant de quelle souris provient
l’échantillon et de quelle portée se trouvait la souris.

``` r
library("phyloseqGraphTest")
library("igraph")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:vegan':
    ## 
    ##     diversity

    ## The following object is masked from 'package:permute':
    ## 
    ##     permute

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library("ggnetwork")
net <- make_network(ps, max.dist=0.35)
sampledata <- data.frame(sample_data(ps))
V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$litter <- sampledata[names(V(net)), "family_relationship"]

net_graph <- ggnetwork(net)
ggplot(net_graph, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter),  size = 3 ) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(0.5,"line")) +
  guides(col = guide_legend(override.aes = list(size = .5)))
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

Graphique : réseau réalisé avec la matrice de dissimilarité de Jaccard.
Il y a un regroupement des échantillons par souris et par portée.

# Tests à deux échantillons basés sur des graphiques

On a proposé l’utilisation d’un arbre MST basé sur les distances entre
les échantillons, puis en comptant le nombre d’arêtes sur l’arbre qui se
trouvaient entre les échantillons dans différents groupes.

## Minimum Spanning Tree (MST)

Avec une dissemblance de Jaccard. Nous voulons savoir si les deux
portées viennent de la même distribution.

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "mst")
gt$pval
```

    ## [1] 0.002

``` r
plotNet1=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet1, plotPerm1)
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

Graphique et l’histogramme de permutation obtenus à partir de l’arbre
couvrant minimal avec similitude Jaccard

Il y a une petite p-value (0.004) et nous rejetons l’hypothèse nulle que
les deux échantillons proviennent de la même distribution. A travers
l’histogramme, on voit que les échantillons sont groupés par portée
plus que ce à quoi nous nous attendrions par hasard.

## Voisins les plus proches

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
```

``` r
plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

Figure : k=1 nearest-neighbor network and permutation histogram

Si une paire d’échantillons a un lien entre eux dans le graphe du plus
proche voisin, ils sont extrêmement susceptibles d’être dans la même.

# Modélisation linéaire

Etude de la relation entre la diversité de la communauté microbienne de
la souris et la les variables d’âge et de portée.

Calcul de la diversité de Shannon associé à chaque échantillon.

``` r
library("nlme")
```

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
library("reshape2")
ps_alpha_div <- estimate_richness(ps, split = TRUE, measure = "Shannon")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>%
  as.factor()
ps_samp <- sample_data(ps) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ps_alpha_div, by = "SampleID") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")

# reorder's facet from lowest to highest diversity
diversity_means <- ps_samp %>%
  group_by(host_subject_id) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
ps_samp$host_subject_id <- factor(ps_samp$host_subject_id)
#                                  diversity_means$host_subject_id)
```

``` r
alpha_div_model <- lme(fixed = alpha_diversity ~ age_binned, data = ps_samp,
                       random = ~ 1 | host_subject_id)
```

``` r
new_data <- expand.grid(host_subject_id = levels(ps_samp$host_subject_id),
                        age_binned = levels(ps_samp$age_binned))
new_data$pred <- predict(alpha_div_model, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2
# fitted values, with error bars
ggplot(ps_samp %>% left_join(new_data)) +
  geom_errorbar(aes(x = age_binned, ymin = pred - 2 * sqrt(pred_var),
                    ymax = pred + 2 * sqrt(pred_var)),
                col = "#858585", size = .1) +
  geom_point(aes(x = age_binned, y = alpha_diversity,
                 col = family_relationship), size = 0.8) +
  facet_wrap(~host_subject_id) +
  scale_y_continuous(limits = c(2.4, 4.6), breaks = seq(0, 5, .5)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Binned Age", y = "Shannon Diversity", color = "Litter") +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)),
        axis.text.x = element_text(angle = -90, size = 6),
        axis.text.y = element_text(size = 6))
```

    ## Joining, by = c("host_subject_id", "age_binned")

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

Graphique : diversité de Shannon selon l’échantillon, l’âge sur deux
portées différentes. On voit que l’âge joue un rôle sur cette différence
de diversité. Il y a plus de diversité pour les souris âgée de 100 à
200.

# Tests multiples hiérarchiques

On teste l’association entre l’abondance microbienne et l’âge.

``` r
library("BiocGenerics")
```

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:igraph':
    ## 
    ##     normalize, path, union

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     combine

    ## The following object is masked from 'package:ade4':
    ## 
    ##     score

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

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

``` r
library("reshape2")
library("DESeq2")
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     collapse

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

``` r
#New version of DESeq2 needs special levels
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship = gsub(" ", "", sample_data(ps)$family_relationship)
ps_dds <- phyloseq_to_deseq2(ps, design = ~ age_binned + family_relationship)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
# geometric mean, set to zero when all coordinates are zero
geo_mean_protected <- function(x) {
  if (all(x == 0)) {
    return (0)
  }
  exp(mean(log(x[x != 0])))
}

geoMeans <- apply(counts(ps_dds), 1, geo_mean_protected)
ps_dds <- estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds <- estimateDispersions(ps_dds)
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
abund <- getVarianceStabilizedData(ps_dds)
```

Raccourcissement des noms de chaque taxon

``` r
short_names <- substr(rownames(abund), 1, 5)%>%
  make.names(unique = TRUE)
rownames(abund) <- short_names
```

``` r
abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(pslog)),
                               sample = rownames(otu_table(pslog)),
                               type = "log(1 + x)"))

ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

Abondance totale transformée DEseq2 dans chaque échantillon

``` r
library("structSSI")
el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$age_binned)
```

La fonction treePvalues permet de tester si ces abondances agrégées sont
significativement différentes entre les environnements

``` r
hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)
```

    ## Number of hypotheses: 764 
    ## Number of tree discoveries: 579 
    ## Estimated tree FDR: 1 
    ## Number of tip discoveries: 280 
    ## Estimated tips FDR: 1 
    ## 
    ##  hFDR adjusted p-values: 
    ##                 unadjp         adjp adj.significance
    ## GCAAG.95  1.861873e-82 3.723745e-82              ***
    ## GCAAG.70  1.131975e-75 2.263950e-75              ***
    ## GCAAG.187 5.148758e-59 1.029752e-58              ***
    ## GCAAG.251 3.519276e-50 7.038553e-50              ***
    ## GCAAG.148 1.274481e-49 2.548962e-49              ***
    ## GCAAG.30  9.925218e-49 1.985044e-48              ***
    ## GCGAG.76  1.722591e-46 3.445183e-46              ***
    ## GCAAG.167 6.249050e-43 1.249810e-42              ***
    ## 255       8.785479e-40 1.757096e-39              ***
    ## GCAAG.64  2.727610e-36 5.455219e-36              ***
    ## [only 10 most significant hypotheses shown] 
    ## --- 
    ## Signif. codes:  0 '***' 0.015 '**' 0.15 '*' 0.75 '.' 1.5 '-' 1

``` r
#interactive part: not run
plot(hfdr_res, height = 5000) # opens in a browser
```

La figure qui s’affiche sur une autre page est un sous-arbre avec de
nombreuses bactéries différentiellement abondantes. Les noeuds sont
colorés selon la p-value en représentant les associations les plus
fortes aux plus faibles. L’association entre le groupe d’âge et
l’abondance bactérienne n’est présente que dans quelques groupes
taxonomiques isolés, mais elle est assez forte dans ces groupes.

``` r
tax <- tax_table(pslog)[, c("Family", "Genus")] %>%
  data.frame()
tax$seq <- short_names
```

``` r
options(digits=3)
hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>% head(10)
```

    ## Joining, by = "seq"

    ##             Family            Genus       seq   unadjp     adjp
    ## 1  Lachnospiraceae             <NA>  GCAAG.95 1.86e-82 3.72e-82
    ## 2  Lachnospiraceae        Roseburia  GCAAG.70 1.13e-75 2.26e-75
    ## 3  Lachnospiraceae Clostridium_XlVa GCAAG.187 5.15e-59 1.03e-58
    ## 4  Lachnospiraceae             <NA> GCAAG.251 3.52e-50 7.04e-50
    ## 5  Lachnospiraceae Clostridium_XlVa GCAAG.148 1.27e-49 2.55e-49
    ## 6  Lachnospiraceae             <NA>  GCAAG.30 9.93e-49 1.99e-48
    ## 7  Ruminococcaceae     Ruminococcus  GCGAG.76 1.72e-46 3.45e-46
    ## 8  Lachnospiraceae Clostridium_XlVa GCAAG.167 6.25e-43 1.25e-42
    ## 9  Lachnospiraceae        Roseburia  GCAAG.64 2.73e-36 5.46e-36
    ## 10            <NA>             <NA>   GCAAG.1 5.22e-35 1.04e-34
    ##    adj.significance
    ## 1               ***
    ## 2               ***
    ## 3               ***
    ## 4               ***
    ## 5               ***
    ## 6               ***
    ## 7               ***
    ## 8               ***
    ## 9               ***
    ## 10              ***

les bactéries les plus fortement associées appartiennent à la famille
des Lachnospiraceae

# Techniques polyvalentes

## Récupération et filtrage des données

On fait une analyse adaptée aussi bien aux comparaisons entre
échantillons qu’à l’identification de traits présentant des variations
intéressantes. 12 échantillons ont été obtenus et 20609 OTU mais 96% des
entrées du tableau d’abondance microbienne sont nulles. Ainsi, le code
ci-dessous récupère et filtre ces données.

``` r
metab <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/metabolites.csv",row.names = 1)
microbe_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/microbe.rda")
load(microbe_connect)
microbe
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 20609 taxa and 12 samples ]
    ## tax_table()   Taxonomy Table:    [ 20609 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 20609 tips and 20607 internal nodes ]

On voit que microbe est un objet phyloseq. On filtre les microbes et les
métabolites d’intérêt. Nous les transformons ensuite.

``` r
library("genefilter")
```

    ## 
    ## Attaching package: 'genefilter'

    ## The following objects are masked from 'package:MatrixGenerics':
    ## 
    ##     rowSds, rowVars

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     rowSds, rowVars

``` r
keep_ix <- rowSums(metab == 0) <= 3
metab <- metab[keep_ix, ]
microbe <- prune_taxa(taxa_sums(microbe) > 4, microbe)
microbe <- filter_taxa(microbe, filterfun(kOverA(3, 2)), TRUE)
metab <- log(1 + metab, base = 10)
X <- otu_table(microbe)
X[X > 50] <- 50
dim(X)
```

    ## [1] 174  12

``` r
dim(metab)
```

    ## [1] 405  12

X et metab ont 12 colonnes, ce sont les échantillons. On fait ensuite
une CCA où il peut y avoir des ensembles d’entités dans des tables de
données de grande dimension, où il peut y avoir plus d’entités mesurées
que d’échantillons.

``` r
library(PMA)
cca_res <- CCA(t(X),  t(metab), penaltyx = .15, penaltyz = .15)
```

    ## 123456789101112131415

``` r
cca_res
```

    ## Call: CCA(x = t(X), z = t(metab), penaltyx = 0.15, penaltyz = 0.15)
    ## 
    ## 
    ## Num non-zeros u's:  5 
    ## Num non-zeros v's:  15 
    ## Type of x:  standard 
    ## Type of z:  standard 
    ## Penalty for x: L1 bound is  0.15 
    ## Penalty for z: L1 bound is  0.15 
    ## Cor(Xu,Zv):  0.974

5 microorganismes et 15 métabolites ont été sélectionnés, selon leur
capacité à expliquer les covariations. Il y a une corrélation de 0,974
entre les deux tableaux. Les données microbiennes et métaboliques
reflètent des signaux sous-jacents similaires.

``` r
combined <- cbind(t(X[cca_res$u != 0, ]),
                  t(metab[cca_res$v != 0, ]))
pca_res <- dudi.pca(combined, scannf = F, nf = 3)
```

``` r
genotype <- substr(rownames(pca_res$li), 1, 2)
sample_type <- substr(rownames(pca_res$l1), 3, 4)
feature_type <- grepl("\\.", colnames(combined))
feature_type <- ifelse(feature_type, "Metabolite", "OTU")
sample_info <- data.frame(pca_res$li, genotype, sample_type)
feature_info <- data.frame(pca_res$c1,
                           feature = substr(colnames(combined), 1, 6))
```

``` r
ggplot() +  geom_point(data = sample_info,
            aes(x = Axis1, y = Axis2, col = sample_type, shape = genotype), size = 3) + 
  geom_label_repel(data = feature_info,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = feature_type),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = feature_info,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = feature_type),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed(sqrt(pca_res$eig[2] / pca_res$eig[2])) +
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pca_res$eig[1] / sum(pca_res$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pca_res$eig[2] / sum(pca_res$eig), 2)),
       fill = "Feature Type", col = "Sample Type")
```

![](002_PhyloSeq_Analysis_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

Graphique : c’est un triplot PCA où il y a les différents types
d’échantillons et les caractéristiques multidomaines (métabolites et
OTU). On peut donc comparer les échantillons et on peut caractériser
l’influence des caractéristiques.

# Conclusion

Le langage R permet de débruiter, d’identifier et de normaliser des
données de séquençage avec des modèles probabilistes avec des
paramètres que l’on a choisi. On a pu visualiser les différentes
abondances, prévalences mais surtout quels étaient les paramètres qui
permettaient les différentes variabilités. On a pu constater que l’âge
était un critère important pour la variabilité, ce qui a été confirmé
avec différentes expériences probabilistes.

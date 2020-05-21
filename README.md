# Summary

You can find a preprint version of the manuscript  
Ines Assum & Julia Krause et al., Tissue-specific multiOMICs analysis of atrial fibrillation, *bioRxiv* (2020)  
here:
https://doi.org/10.1101/2020.04.06.021527

**Abstract:**  
Genome-wide association studies (GWAS) for atrial fibrillation (AF) have uncovered numerous disease-associated variants. Their underlying molecular mechanisms, especially consequences for mRNA and protein expression remain largely elusive. Thus, novel multiOMICs approaches are needed for deciphering the underlying molecular networks. Here, we integrated genomics, transcriptomics, and proteomics of human atrial tissue which allowed for identifying widespread effects of genetic variants on both transcript (cis eQTL) and protein (cis pQTL) abundance. We further established a novel targeted trans QTL approach based on polygenic risk scores to identify candidates for AF core genes. Using this approach, we identified two trans eQTLs and four trans pQTLs for AF GWAS hits, and elucidated the role of the transcription factor NKX2-5 as a link between the GWAS SNP rs9481842 and AF. Altogether, we present an integrative multiOMICs method to uncover trans-acting networks in small datasets and provide a rich resource of atrial tissue-specific regulatory variants for transcript and protein levels for cardiovascular disease gene prioritization.



We address a key hypothesis about the existence of core genes as postulated in the omnigenic model by [Liu et al.](https://doi.org/10.1016/j.cell.2019.04.014), *Cell* (2019). Core genes are central genes with trans-associations to GWAS loci, whose expression levels directly affect a disease phenotype. Here we sought to identify candidate core genes for AF to understand the contribution of trans-genetic effects in the pathology of AF. To prioritize genes satisfying the properties predicted by the omnigenic model, we evaluated the accumulation of trans-effects, their relevance in gene regulatory networks, and the disease association by the following strategy:


i) Evaluate cumulated trans-effects of disease-associated variants on expression by ranking genes based on their correlation of mRNA and protein abundance with a PRS (eQTS) as proposed by [Võsa et al.](https://doi.org/10.1101/447367), *bioRxiv* (2018)

ii) Identify genes sharing molecular function and representing biological networks that propagate genetic trans-effects to core genes by pathway enrichment analysis (GSEA) on the eQTS rankings. Genes driving the enrichment of multiple gene sets were selected as core gene candidates.

iii) Establish the link between the core gene candidates and the disease based on a significant trans eQTL GWAS hit.



# Run your own analysis

Want to try our approach? Let's get started with a short [tutorial](https://colab.research.google.com/drive/1ZSQF1Lh86tVgIlfrK10EB3VUlNDoG9fs?usp=sharing) on google colab (15min installation and <5min runtime)!

## Installation instructions
Analysis was done using R 3.4.1.

In general, four packages need to be installed
* [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
* [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
* [MatrixEQTL](https://github.com/andreyshabalin/MatrixEQTL)
* [eQTLpipeline](https://github.com/matthiasheinig/eQTLpipeline) (optional)

```
install.packages("devtools")
BiocManager::install("fgsea", dependencies = T, clean = T)
library(devtools)
devtools::install_github("andreyshabalin/MatrixEQTL", force=T)
devtools::install_github("matthiasheinig/eQTLpipeline", force=T)
```

Additional packages are required to run the cis QTL pipeline with PEER analysis. For this, we provide a conda environment running on fedora 25, as PEER frequently causes problems when being installed as a R package. A docker will be supplied soon.   
For installation, please
* create a new [conda]() environment, using the [yml-file](https://raw.githubusercontent.com/heiniglab/symatrial/master/envs/r341peer.yml) `conda env create -f r341peer.yml`
* activate the environment `conda activate r341peer`
* start `R` and install packages MatrixEQTL and eQTLpipeline from github by running
```
library(devtools)
devtools::install_github("andreyshabalin/MatrixEQTL", force=T)
devtools::install_github("matthiasheinig/eQTLpipeline", force=T)
```



## Tutorial

A short example of how to run our pre-selection approach with a short example dataset can be found as a [R markdown](https://github.com/heiniglab/symatrial/blob/master/example_data/PRSenrichQTL_tutorial.Rmd), [html document](https://github.com/heiniglab/symatrial/blob/master/example_data/PRSenrichQTL_tutorial.html) or as a google [colab notebook](https://colab.research.google.com/drive/1ZSQF1Lh86tVgIlfrK10EB3VUlNDoG9fs?usp=sharing).

If you want to run the analysis, in general the following R packages are required:
* [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
* [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
* [MatrixEQTL](https://github.com/andreyshabalin/MatrixEQTL)
* [eQTLpipeline](https://github.com/matthiasheinig/eQTLpipeline)

Please also download gene set annotations here:
[GO biological processes](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.1/c5.bp.v7.1.symbols.gmt)


## Cis QTL analysis

Cis QTL analyses were performed 

You can find analysis code for the cis QTL analysis in the folder
[qtl_pipeline](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline)

containing numerous scripts for
* preprocessing
* running PEER and then QTL analysis
* functional annotations
* [comparison](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline/analysis/comparison/) comparisons to other datasets (e.g. GTEx)

### Functional annotations
* [build snp-gene pair annotations](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline/functional_analysis/)
* [functions to evaluate enrichment of functional elements](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline/functional_analysis/enrichment_analysis/)
* [derivation of RBP binding sites](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline/functional_analysis/eclip_preprocessing.R)
* [NKX2-5 binding sites and TF activity](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline/analysis/TF_activity.Rmd)

### GWAS annotations
* [calculate GWAS overlaps and enrichments](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline/analysis/gwas_imputed/script.R)


## Trans analysis

### Genome-wide polygenic scores for AF and CAD

Code for the computation of the polygenic risk score for AF on both our cohort and the 1000 genomes individuals can be found  [here](https://github.com/heiniglab/symatrial/blob/master/scripts/PRS_trans_analyses/polygenic_risk_scores.R).

The final list of 109 tested SNPs derived by [pruning](https://github.com/heiniglab/symatrial/blob/master/scripts/qtl_pipeline/analysis/gwas_imputed/AF_snp_pruning.R) all SNPs annotated with AF in the GWAS catalog are supplied [here](https://github.com/heiniglab/symatrial/blob/master/example_data/AF_SNPs_pruned_hg19.txt).

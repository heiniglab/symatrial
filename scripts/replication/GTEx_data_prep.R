# Setup ----

get.duplicates <- function(df, cols){
  df2 <- data.frame(table(df[, cols]))
  df2 <- df2[df2$Freq>1, ]
  colnames(df2) <- c(cols, "Freq")
  df <- merge(df,
              df2[cols],
              all.y=T)
}

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(readxl)

HMGU.blue <- "#003E6E"
mygray <- "#C6DDEA"
col.paired <- brewer.pal(n = 11, "Paired")
col.set <- col.paired[c(2,9,8)]

## GTEx data directory ----

# project based settings for data paths (incl. get.path function)
source("../helper/helper.R")
local=F
gtex <- paste0(get.path("replication", local), "GTEx_v8/")

# Load data ----

# Histology info: https://gtexportal.org/home/histologyPage
# (manual selection + csv download)

# Sources:
# GTEx v8: https://gtexportal.org/home/datasets
# Proteomics: PXD016999, supplement https://doi.org/10.1016/j.cell.2020.08.036

## Variable dictionaries:
# https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx
# https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx

## Sample annotations:
# https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
# https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

## RNAseq
# Gene level TPM:
# https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
# Transcript TPM:
# https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
# Normalized expression matrices for QTL analysis (archive with matrices per tissue):
# https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar

## Proteomics
# Annotations: https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc2.xlsx
# Data: https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc3.xlsx


## GTEx sample annotations ----

# GTEx sample names:
# GTEX-12WSD-1226    -SM-9KMWI
# GTEX-donor-tissueID-SM-measurementID
# donor 12WSD has
# e.g. RNAseq sample     GTEX-12WSD-1226-SM-5HL9Q and
#      proteomics sample GTEX-12WSD-1226-SM-9KNJG
sample.anno <- read.csv(paste0(gtex,
                               "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
                        stringsAsFactors = F, sep = "\t", h=T)
sample.anno$SAMPID2 <- gsub("\\-", "\\.", sample.anno$SAMPID)
sample.anno$sample <- gsub("(^[^\\-]+\\-[^\\-]+).*$", "\\1", sample.anno$SAMPID)
sample.anno$sample2 <- gsub("\\-", "\\.", sample.anno$sample)

sample.anno2 <- read.csv(paste0(gtex,
                                "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"),
                         stringsAsFactors = F, sep = "\t", h=T)
sample.anno <- merge(sample.anno, sample.anno2,
                     by.x = "sample", by.y = "SUBJID",
                     all = T)

# see variable dictionary
colnames(sample.anno)

# tissue types
unique(sample.anno$SMTSD)

# select other columns of interest or tissue types
samples.haa <- sample.anno[sample.anno$SMTSD=="Heart - Atrial Appendage",
                           c("SAMPID", "SAMPID2", "sample", "sample2",
                             "SEX", "AGE", "DTHHRDY",
                             "SMGEBTCH", "SMCENTER", "SMPTHNTS",
                             "SMRIN", "SMTS", "SMTSD", "SMUBRID", "SMTSISCH",
                             "SMTSPAX")]


## GTEx RNAseq data ----
expr.norm <- read.csv(paste0(get.path("replication", local),
                             "GTEx_v8/GTEx_Analysis_v8_eQTL_expression_matrices/",
                             "Heart_Atrial_Appendage.v8.normalized_expression.bed"),
                      stringsAsFactors = F, h = T, sep = "\t")
rownames(expr.norm) <- expr.norm$gene_id
anno.gene <- expr.norm[, 1:4]
anno.gene$gene_id2 <- gsub("\\..*", "", anno.gene$gene_id)
expr.norm[, 1:4] <- NULL

# GTEx proteomics ----

gtex.prot <- data.frame(read_excel(paste0(gtex,
                                          "/PXD016999/1-s2.0-S0092867420310783-mmc3.xlsx"),
                                   sheet = "C protein normalized abundance",
                                   na = "NA"))

gtex.prot.samples <- data.frame(t(gtex.prot[3:1, -c(1:2)]))
colnames(gtex.prot.samples) <- c("ID", "tissue", "tag")
gtex.prot.samples$sample <- gsub("(^[^\\-]+\\-[^\\-]+).*$", "\\1", gtex.prot.samples$ID)
rownames(gtex.prot.samples) <- NULL
gtex.prot.samples <- gtex.prot.samples[gtex.prot.samples$tissue!="reference", ]
gtex.prot.samples$ID2 <- gsub("\\-", "\\.", gtex.prot.samples$ID)
table(gtex.prot.samples$tissue)

colnames(gtex.prot) <- gtex.prot[3, ]
gtex.prot <- gtex.prot[-c(1:3), grep("gene.id.full|gene.id|GTEX", colnames(gtex.prot))]

gtex.prot.anno <- data.frame(read_excel(paste0(gtex,
                                               "/PXD016999/1-s2.0-S0092867420310783-mmc3.xlsx"),
                                        sheet = "G protein TS score",
                                        na = "NA", skip = 2))
gtex.prot.anno <- merge(gtex.prot.anno[, c("ensembl_id", "entrez_id", "hgnc_name", "hgnc_symbol")],
                        gtex.prot[, c("gene.id", "gene.id.full")],
                        by.x = "ensembl_id", by.y = "gene.id", 
                        all.y = T, sort = F)
rownames(gtex.prot.anno) <- gtex.prot.anno$ensembl_id
rownames(gtex.prot) <- gtex.prot$gene.id
gtex.prot[, c("gene.id", "gene.id.full")] <- NULL
gtex.prot.anno <- gtex.prot.anno[rownames(gtex.prot), ]

gtex.prot.aa <- gtex.prot[rownames(gtex.prot) %in% gtex.prot.anno[gtex.prot.anno$hgnc_symbol %in% genes, "ensembl_id"],
                          gtex.prot.samples[gtex.prot.samples$tissue=="Heart - Atrial Appendage", "ID"]]
unique(gtex.prot.samples[gtex.prot.samples$tissue=="Heart - Atrial Appendage", "sample"])
sum(unique(gtex.prot.samples[gtex.prot.samples$tissue=="Heart - Atrial Appendage", "sample"]) %in% df.norm$sample)
colnames(gtex.prot.aa)
rownames(gtex.prot.aa) <- gtex.prot.anno[rownames(gtex.prot.aa), "hgnc_symbol"]
gtex.prot.aa[genes, ]

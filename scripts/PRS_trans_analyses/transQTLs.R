## calculate transQTLs
setwd("~/work/symAtrial_QTL/scripts/analysis/gwas_imputed")
source("../../helper/helper.R")
source("../../preprocessing/correction/normalization/norm.R")
local=F

library(ggplot2)
library(ggpubr)
library(VariantAnnotation)
library(GenomicRanges)
library(eQTLpipeline)

# AF GWAS hits (pruned) trans eQTLs for core gene candidates -------------------
snp.list <- paste0(get.path("gwas2", local),
                   "AF_snps_pruned.txt")
geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_t_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalog_snps_pruned_mRNA.txt")
system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
             ">", geno2, sep=" "))
snplocs <- read.csv(get.path("snplocs imputed", local),
                    sep = "\t", h=T, stringsAsFactors = F)
colnames(snplocs) <- c("snps", "chr", "pos")

locs.gene <- as.data.frame(readRDS(paste0(get.path("locations", local),
                                          "gencode.v31lift37.gene.level.locations.RDS")))
colnames(locs.gene) <- c("chr.gene", "gene.start", "gene.end", "gene.length", "gene.strand", "gene")

fgsea.trans <- readRDS(paste0(get.path("dataset", local),
                              "risk_score/correlations/",
                              "fgsea_assoc_transcriptomics_GPS_covs_RIN.RDS"))
df1 <- fgsea.trans[["GObp_high"]]
lead.trans <- as.data.frame(table(unlist(df1[df1$padj<0.05, "leadingEdge"])),
                            stringsAsFactors=F)
dim(lead.trans)
hist(lead.trans$Freq, breaks = 100)
lead.trans <- lead.trans[lead.trans$Freq>=16, ]
dim(lead.trans)
expr <- readRDS(paste0(get.path("dataset", local),
                       "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expr) <- sub(".RAA", "", colnames(expr))
expr <- expr[rownames(expr) %in% lead.trans$Var1, geno.names]
expr <- correct.linewise(expr)
identical(colnames(expr), geno.names)
dim(expr)

tpheno <- readRDS(paste0(get.path("dataset", local),
                         "AFHRI_B_transcriptomics_phenos_symbol.RDS"))
tpheno$externID2 <- sub(".RAA", "", tpheno$externID)
rownames(tpheno) <- tpheno$externID2
covs <- data.frame(t(tpheno[geno.names,
                            c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "RIN", "fibro.score")]))
identical(colnames(covs), geno.names)

trans.core <- trans.qtl(prefix = paste0(get.path("results", local),
                                        "imputed/trans/PRS_core_mRNA_RIN_AF_GWAS_snps_trans_eQTL_25"),
                        genotype_file_name = geno2,
                        expression_file_name = expr,
                        covariates_file_name = covs,
                        snp.pos = snplocs,
                        threshold = 1,
                        compute.all = T,
                        min.pv.by.genesnp = T,
                        save.memory = F,
                        load.qtls = T)

trans.core$all$eqtls$chr <- factor(trans.core$all$eqtls$chr,
                                   levels = c(1:22),
                                   ordered = T)
trans.core$all$eqtls <- merge(trans.core$all$eqtls,
                              locs.gene,
                              all.x = T)
saveRDS(trans.core,
        file = paste0(get.path("results", local),
                      "imputed/trans/PRS_core_mRNA_RIN_AF_GWAS_snps_trans_eQTL_results_25.RDS"))

trans.core <- readRDS(paste0(get.path("results", local),
                             "imputed/trans/PRS_core_mRNA_RIN_AF_GWAS_snps_trans_eQTL_results_25.RDS"))
manhattan.qtl(trans.core) +
  ggtitle("Trans eQTL associations between AF GWAS SNPs (109) and leading edge transcripts (25)")
table <- trans.core$all$eqtls
print.data.frame(table[order(table$pvalue)[1:10], ],
                 row.names = F, digits = 4)

# GWAS SNPs (pruned), proteins for trans leading edge --------------------------
snp.list <- paste0(get.path("gwas2", local),
                   "AF_snps_pruned_prot.txt")
geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_p_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalogue_snps_pruned_prot.txt")
system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
             ">", geno2, sep=" "))
snplocs <- read.csv(get.path("snplocs imputed", local),
                    sep = "\t", h=T, stringsAsFactors = F)
colnames(snplocs) <- c("snps", "chr", "pos")

locs.gene <- as.data.frame(readRDS(paste0(get.path("locations", local),
                                          "gencode.v31lift37.gene.level.locations.RDS")))
colnames(locs.gene) <- c("chr.gene", "gene.start", "gene.end", "gene.length", "gene.strand", "gene")

fgsea.trans <- readRDS(paste0(get.path("dataset", local),
                              "risk_score/correlations/",
                              "fgsea_assoc_transcriptomics_GPS_covs_RIN.RDS"))
df1 <- fgsea.trans[["GObp_high"]]
lead.trans <- as.data.frame(table(unlist(df1[df1$padj<0.05, "leadingEdge"])),
                            stringsAsFactors=F)
dim(lead.trans)
lead.trans <- lead.trans[lead.trans$Freq>=16, ]
dim(lead.trans)
expr <- readRDS(paste0(get.path("dataset", local),
                       "AFHRI_B_proteomics_QC_symbol.RDS"))
colnames(expr) <- sub(".RAA", "", colnames(expr))
expr <- expr[rownames(expr) %in% lead.trans$Var1, geno.names]
expr <- correct.linewise(expr)
identical(colnames(expr), geno.names)
dim(expr)

ppheno <- readRDS(paste0(get.path("dataset", local),
                         "AFHRI_B_proteomics_phenos_symbol.RDS"))
ppheno$externID2 <- sub(".RAA", "", ppheno$externID)
rownames(ppheno) <- ppheno$externID2
covs <- data.frame(t(ppheno[geno.names,
                            c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "Protein.c..ug.ul.", "fibro.score")]))
identical(colnames(covs), geno.names)

trans.core.prot <- trans.qtl(prefix = paste0(get.path("results", local),
                                             "imputed/trans/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_pQTL_25"),
                             genotype_file_name = geno2,
                             expression_file_name = expr,
                             covariates_file_name = covs,
                             snp.pos = snplocs,
                             threshold = 1,
                             compute.all = T,
                             min.pv.by.genesnp = T,
                             save.memory = F,
                             load.qtls = T)

trans.core.prot$all$eqtls$chr <- factor(trans.core.prot$all$eqtls$chr,
                                        levels = c(1:22),
                                        ordered = T)
trans.core.prot$all$eqtls <- merge(trans.core.prot$all$eqtls,
                                   locs.gene,
                                   all.x = T)
saveRDS(trans.core.prot,
        file = paste0(get.path("results", local),
                      "imputed/trans/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_pQTL_results_25.RDS"))

trans.core.prot <- readRDS(paste0(get.path("results", local),
                                  "imputed/trans/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_pQTL_results_25.RDS"))
pqtl <- manhattan.qtl(trans.core.prot) +
  ggtitle("Proteins for transcript leading edge")
table <- trans.core.prot$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[which(!duplicated(table$gene))[1:10], ],
                 row.names = F, digits = 4)


# GWAS SNPs (pruned), proteins for protein leading edge ------------------------
snp.list <- paste0(get.path("gwas2", local),
                   "AF_snps_pruned.txt")
geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_p_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalogue_snps_pruned_prot.txt")
system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
             ">", geno2, sep=" "))
snplocs <- read.csv(get.path("snplocs imputed", local),
                    sep = "\t", h=T, stringsAsFactors = F)
colnames(snplocs) <- c("snps", "chr", "pos")

locs.gene <- as.data.frame(readRDS(paste0(get.path("locations", local),
                                          "gencode.v31lift37.gene.level.locations.RDS")))
colnames(locs.gene) <- c("chr.gene", "gene.start", "gene.end", "gene.length", "gene.strand", "gene")

fgsea.prot <- readRDS(paste0(get.path("dataset", local),
                             "risk_score/correlations/",
                             "fgsea_assoc_proteomics_GPS_covs_conc.RDS"))
df1 <- fgsea.prot[["GObp_all"]]
lead.prot <- as.data.frame(table(unlist(df1[df1$padj<0.05, "leadingEdge"])),
                           stringsAsFactors=F)
dim(lead.prot)

expr <- readRDS(paste0(get.path("dataset", local),
                       "AFHRI_B_proteomics_QC_symbol.RDS"))
colnames(expr) <- sub(".RAA", "", colnames(expr))
expr <- expr[rownames(expr) %in% lead.prot$Var1, geno.names]
expr <- correct.linewise(expr)
identical(colnames(expr), geno.names)
dim(expr)

ppheno <- readRDS(paste0(get.path("dataset", local),
                         "AFHRI_B_proteomics_phenos_symbol.RDS"))
ppheno$externID2 <- sub(".RAA", "", ppheno$externID)
rownames(ppheno) <- ppheno$externID2
covs <- data.frame(t(ppheno[geno.names,
                            c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "Protein.c..ug.ul.", "fibro.score")]))
identical(colnames(covs), geno.names)

prot.core <- trans.qtl(prefix = paste0(get.path("results", local),
                                       "imputed/trans/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL_145"),
                       genotype_file_name = geno2,
                       expression_file_name = expr,
                       covariates_file_name = covs,
                       snp.pos = snplocs,
                       threshold = 1,
                       compute.all = T,
                       min.pv.by.genesnp = T,
                       save.memory = F,
                       load.qtls = T)

prot.core$all$eqtls$chr <- factor(prot.core$all$eqtls$chr,
                                  levels = c(1:22),
                                  ordered = T)
prot.core$all$eqtls <- merge(prot.core$all$eqtls,
                             locs.gene,
                             all.x = T)

saveRDS(prot.core,
        file = paste0(get.path("results", local),
                      "imputed/trans/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL_results_145.RDS"))

# prot.core <- readRDS(paste0(get.path("results", local),
#                                   "imputed/trans/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL_results_145.RDS"))
manhattan.qtl(prot.core) +
  ggtitle("Trans pQTL associations between AF GWAS SNPs (109) and leading edge proteins (145)")
table <- prot.core$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[which(!duplicated(table$gene))[1:10], ],
                 row.names = F, digits = 4)
print.data.frame(table[table$FDR<0.2, ],
                 row.names = F, digits = 4)


# GWAS SNPs (pruned), transcripts for protein leading edge -------------
ssnp.list <- paste0(get.path("gwas2", local),
                    "AF_snps_pruned.txt")
geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_t_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalog_snps_pruned_mRNA.txt")
# system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
#              ">", geno2, sep=" "))
snplocs <- read.csv(get.path("snplocs imputed", local),
                    sep = "\t", h=T, stringsAsFactors = F)
colnames(snplocs) <- c("snps", "chr", "pos")

locs.gene <- as.data.frame(readRDS(paste0(get.path("locations", local),
                                          "gencode.v31lift37.gene.level.locations.RDS")))
colnames(locs.gene) <- c("chr.gene", "gene.start", "gene.end", "gene.length", "gene.strand", "gene")

fgsea.prot <- readRDS(paste0(get.path("dataset", local),
                             "risk_score/correlations/",
                             "fgsea_assoc_proteomics_GPS_covs_conc.RDS"))
df1 <- fgsea.prot[["GObp_all"]]
lead.prot <- as.data.frame(table(unlist(df1[df1$padj<0.05, "leadingEdge"])),
                           stringsAsFactors=F)
dim(lead.prot)

expr <- readRDS(paste0(get.path("dataset", local),
                       "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expr) <- sub(".RAA", "", colnames(expr))
lead.prot$Var1[!(lead.prot$Var1 %in% rownames(expr))]
expr <- expr[rownames(expr) %in% lead.prot$Var1, geno.names]
expr <- correct.linewise(expr)
identical(colnames(expr), geno.names)
dim(expr)

tpheno <- readRDS(paste0(get.path("dataset", local),
                         "AFHRI_B_transcriptomics_phenos_symbol.RDS"))
tpheno$externID2 <- sub(".RAA", "", tpheno$externID)
rownames(tpheno) <- tpheno$externID2
covs <- data.frame(t(tpheno[geno.names,
                            c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "RIN", "fibro.score")]))
identical(colnames(covs), geno.names)

core.prot.trans <- trans.qtl(prefix = paste0(get.path("results", local),
                                             "imputed/trans/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_eQTL_145"),
                             genotype_file_name = geno2,
                             expression_file_name = expr,
                             covariates_file_name = covs,
                             snp.pos = snplocs,
                             threshold = 1,
                             compute.all = T,
                             min.pv.by.genesnp = T,
                             save.memory = F,
                             load.qtls = T)

core.prot.trans$all$eqtls$chr <- factor(core.prot.trans$all$eqtls$chr,
                                        levels = c(1:22),
                                        ordered = T)
core.prot.trans$all$eqtls <- merge(core.prot.trans$all$eqtls,
                                   locs.gene,
                                   all.x = T)

saveRDS(core.prot.trans,
        file = paste0(get.path("results", local),
                      "imputed/trans/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_eQTL_results_145.RDS"))

core.prot.trans <- readRDS(paste0(get.path("results", local),
                                  "imputed/trans/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_eQTL_results_145.RDS"))
manhattan.qtl(core.prot.trans) +
  ggtitle("Transcripts for protein leading edge")
table <- core.prot.trans$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[which(!duplicated(table$gene))[1:10], ],
                 row.names = F, digits = 4)


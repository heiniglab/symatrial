setwd("~/work/symAtrial_QTL/scripts/polygenic_risk_scores/")
source("../helper/helper.R")
source("../preprocessing/correction/normalization/norm.R")
local=T

library(ggplot2)
library(ggpubr)
library(VariantAnnotation)
library(GenomicRanges)
library(eQTLpipeline)


# AF GWAS hits (no proxies, pruned) trans eQTLs for core genes ---------------------------------------
snp.list <- paste0(get.path("gwas2", local),
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

fgsea.eQTS <- readRDS(paste0(get.path("results", local),
                             "imputed/trans/eQTS/",
                             "fgsea_eQTS_percAF_covs_RIN.RDS"))
lead.trans <- fgsea.eQTS[["leadingEdge"]]
dim(lead.trans)
#hist(lead.trans$Freq, breaks = 100)
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
                                        "imputed/trans/final/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_eQTL_25"),
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
# saveRDS(trans.core,
#         file = paste0(get.path("results", local),
#                       "imputed/trans/final/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_eQTL_results_25.RDS"))

trans.core <- readRDS(paste0(get.path("results", local),
                             "imputed/trans/final/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_eQTL_results_25.RDS"))

manhattan.qtl(trans.core) +
  ggtitle("trans associations between AF GWAS SNPs (109) and leading edge genes (25)")
table <- trans.core$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[table$FDR<0.2, ],
                 row.names = F, digits = 4)


# only GWAS SNPs, no proxies pruned, proteins for trans leading edge ---------------------------------------
snp.list <- paste0(get.path("gwas2", local),
                   "AF_snps_pruned_prot.txt")
geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_p_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalogue_snps_pruned_prot.txt")
# system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
#              ">", geno2, sep=" "))

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
                                             "imputed/trans/final/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_pQTL_25"),
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
# saveRDS(trans.core.prot,
#         file = paste0(get.path("results", local),
#                       "imputed/trans/final/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_pQTL_results_25.RDS"))

trans.core.prot <- readRDS(paste0(get.path("results", local),
                                  "imputed/trans/final/PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_pQTL_results_25.RDS"))
pqtl <- manhattan.qtl(trans.core.prot) +
  ggtitle("...")
table <- trans.core.prot$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[which(!duplicated(table$gene))[1:10], ],
                 row.names = F, digits = 4)


# only GWAS SNPs, no proxies pruned, proteins for protein leading edge ---------------------------------------
snp.list <- paste0(get.path("gwas2", local),
                   "AF_snps_pruned.txt")
geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_p_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalogue_snps_pruned_prot.txt")
# system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
#              ">", geno2, sep=" "))
# snplocs <- read.csv(get.path("snplocs imputed", local),
#                     sep = "\t", h=T, stringsAsFactors = F)
# colnames(snplocs) <- c("snps", "chr", "pos")
# 
# locs.gene <- as.data.frame(readRDS(paste0(get.path("locations", local),
#                                           "gencode.v31lift37.gene.level.locations.RDS")))
# colnames(locs.gene) <- c("chr.gene", "gene.start", "gene.end", "gene.length", "gene.strand", "gene")

fgsea.pQTS <- readRDS(paste0(get.path("results", local),
                             "imputed/trans/pQTS/",
                             "fgsea_pQTS_percAF_covs_conc.RDS"))
lead.prot <- fgsea.pQTS[["leadingEdge"]]
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
                                       "imputed/trans/final/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL_168"),
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

# saveRDS(prot.core,
#         file = paste0(get.path("results", local),
#                       "imputed/trans/final/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL_results_168.RDS"))
prot.core <- readRDS(paste0(get.path("results", local),
                            "imputed/trans/final/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL_results_168.RDS"))

manhattan.qtl(prot.core) +
  ggtitle("trans associations between AF GWAS SNPs (109) and protein for leading edge proteins (168)")
table <- prot.core$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[table$FDR<0.2, ],
                 row.names = F, digits = 4)

# only GWAS SNPs, no proxies pruned, transcripts for protein leading edge ---------------------------------------
snp.list <- paste0(get.path("gwas2", local),
                   "AF_snps_pruned.txt")
geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_t_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalog_snps_pruned_mRNA.txt")
# system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
#              ">", geno2, sep=" "))

dim(lead.prot)

expr <- readRDS(paste0(get.path("dataset", local),
                       "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expr) <- sub(".RAA", "", colnames(expr))
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
                                             "imputed/trans/final/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_eQTL_168"),
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

# saveRDS(core.prot.trans,
#         file = paste0(get.path("results", local),
#                       "imputed/trans/final/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_eQTL_results_168.RDS"))

core.prot.trans <- readRDS(paste0(get.path("results", local),
                                  "imputed/trans/final/PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_eQTL_results_168.RDS"))
manhattan.qtl(core.prot.trans) +
  ggtitle("...")
table <- core.prot.trans$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[which(!duplicated(table$gene))[1:10], ],
                 row.names = F, digits = 4)



trans.table <- trans.core$all$eqtls
trans.table <- trans.table[order(trans.table$pvalue), ]
trans.table$snps <- gsub("\\:.*", "", trans.table$snps)

prot.table <- prot.core$all$eqtls
prot.table <- prot.table[order(prot.table$pvalue), ]
prot.table$snps <- gsub("\\:.*", "", prot.table$snps)

print.data.frame(rbind(trans.table[trans.table$FDR<0.2, ],
                       prot.table[prot.table$FDR<0.2, ]),
                 row.names = F, digits = 4)

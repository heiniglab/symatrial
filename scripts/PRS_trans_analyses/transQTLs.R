setwd("~/work/symAtrial_QTL/scripts/analysis/gwas_imputed")
pics <- "~/work/symAtrial_QTL/scripts/analysis/gwas_imputed/plots/"
source("../../helper/helper.R")
source("../../preprocessing/correction/normalization/norm.R")
source("../../preprocessing/imputation/run_imputation.R")
local=F

library(ggplot2)
library(ggpubr)
library(VariantAnnotation)
library(GenomicRanges)
library(eQTLpipeline)

snplocs <- read.csv(get.path("snplocs imputed", local),
                    sep = "\t", h=T, stringsAsFactors = F)
colnames(snplocs) <- c("snps", "chr", "pos")

locs.gene <- as.data.frame(readRDS(paste0(get.path("locations", local),
                                          "gencode.v31lift37.gene.level.locations.RDS")))
colnames(locs.gene) <- c("chr.gene", "gene.start", "gene.end", "gene.length", "gene.strand", "gene")

snp.list <- paste0(get.path("gwas2", local),
                   "AF_snps_pruned.txt")

out.dir <- paste0(get.path("results", local), "imputed/trans/new/")

# AF GWAS hits (no proxies, pruned) trans eQTLs for transcriptomics leading edge ------------------------
prefix.eqtl <- paste0(out.dir,
                      "PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_eQTL")

geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_t_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalog_snps_pruned_mRNA.txt")
# system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
#              ">", geno2, sep=" "))

fgsea.eQTS.file <- paste0(get.path("results", local),
                          "imputed/trans/eQTS_cis/",
                          "fgsea_eQTS_percAF_covs_RIN_cis_eQTLs.RDS")
fgsea.trans <- readRDS(fgsea.eQTS.file)
df1 <- fgsea.trans[["GObp_all"]]
lead.trans <- as.data.frame(table(unlist(df1[df1$padj<0.05, "leadingEdge"])),
                            stringsAsFactors=F)
dim(lead.trans)
#hist(lead.trans$Freq, breaks = 100)
lead.trans <- lead.trans[lead.trans$Freq>=14, ]
dim(lead.trans)
expr <- readRDS(paste0(get.path("dataset", local),
                       "AFHRI_B_transcriptomics_QC_symbol.RDS"))
#lead.trans$Var1[!(lead.trans$Var1 %in% rownames(expr))]
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

trans.core <- trans.qtl(prefix = prefix.eqtl,
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
#         file = paste0(prefix.eqtl, "_results.RDS"))

# trans.core <- readRDS(paste0(prefix.eqtl, "_results.RDS"))

manhattan.qtl(trans.core) +
  ggtitle("trans associations between AF GWAS SNPs (108) and leading edge genes (23)")
table <- trans.core$all$eqtls
print.data.frame(table[order(table$pvalue)[1:10], ],
                 row.names = F, digits = 4)
print.data.frame(table[table$FDR<0.2, ],
                 row.names = F, digits = 4)

# AF GWAS hits (no proxies, pruned) trans pQTLs for transcriptomics leading edge ---------------------------------------
prefix.pqtl <- paste0(out.dir,
                      "PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_pQTL")

geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_p_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalogue_snps_pruned_prot.txt")

fgsea.trans <- readRDS(fgsea.eQTS.file)
df1 <- fgsea.trans[["GObp_all"]]
lead.trans <- as.data.frame(table(unlist(df1[df1$padj<0.05, "leadingEdge"])),
                            stringsAsFactors=F)
dim(lead.trans)
lead.trans <- lead.trans[lead.trans$Freq>=14, ]

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
                            c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "Protein.c..ug.ul.", "fibro.score.imp.prot")]))
identical(colnames(covs), geno.names)

trans.core.prot <- trans.qtl(prefix = prefix.pqtl,
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
#         file = paste0(prefix.pqtl, "_results.RDS"))

# trans.core.prot <- readRDS(paste0(prefix.pqtl, "_results.RDS"))


# AF GWAS hits (no proxies, pruned) trans pQTLs for proteomics leading edge ---------------------------------------
prefix.pqtl <- paste0(out.dir,
                      "PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL")

geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_p_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalogue_snps_pruned_prot.txt")
# system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", snp.list, geno,
#              ">", geno2, sep=" "))

fgsea.pQTS.file <- paste0(get.path("results", local),
                          "imputed/trans/pQTS_cis/",
                          "fgsea_pQTS_percAF_covs_conc_cis_pQTLs.RDS")
fgsea.prot <- readRDS(fgsea.pQTS.file)
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
                            c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "Protein.c..ug.ul.", "fibro.score.imp.prot")]))
identical(colnames(covs), geno.names)

prot.core <- trans.qtl(prefix = prefix.pqtl,
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
        file = paste0(prefix.pqtl, "_results.RDS"))

# prot.core <- readRDS(paste0(prefix.pqtl, "_results.RDS"))

manhattan.qtl(prot.core) +
  ggtitle("trans associations between AF GWAS SNPs (108) and protein for leading edge proteins (152)")
table <- prot.core$all$eqtls
table <- table[order(table$pvalue), ]
print.data.frame(table[which(!duplicated(table$gene))[1:10], ],
                 row.names = F, digits = 4)
print.data.frame(table[table$FDR<0.2, ],
                 row.names = F, digits = 4)

# AF GWAS hits (no proxies, pruned) trans eQTLs for proteomics leading edge ------------------------
prefix.eqtl <- paste0(out.dir,
                      "PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_eQTL")

geno <- paste0(get.path("genotype", local),
               "genotype_imputed_common_samples_t_raa.txt")
geno.names <- colnames(read.csv(geno, sep="\t", nrows=5))
geno.names <- geno.names[!(geno.names=="X")]
geno2 <- paste0(get.path("results", local),
                "imputed/trans/AF_GWAS_catalog_snps_pruned_mRNA.txt")

fgsea.prot <- readRDS(fgsea.pQTS.file)
df1 <- fgsea.prot[["GObp_all"]]
lead.prot <- as.data.frame(table(unlist(df1[df1$padj<0.05, "leadingEdge"])),
                           stringsAsFactors=F)

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

prot.core.trans <- trans.qtl(prefix = prefix.eqtl,
                             genotype_file_name = geno2,
                             expression_file_name = expr,
                             covariates_file_name = covs,
                             snp.pos = snplocs,
                             threshold = 1,
                             compute.all = T,
                             min.pv.by.genesnp = T,
                             save.memory = F,
                             load.qtls = T)

prot.core.trans$all$eqtls$chr <- factor(prot.core.trans$all$eqtls$chr,
                                        levels = c(1:22),
                                        ordered = T)
prot.core.trans$all$eqtls <- merge(prot.core.trans$all$eqtls,
                                   locs.gene,
                                   all.x = T)
saveRDS(prot.core.trans,
        file = paste0(prefix.eqtl, "_results.RDS"))

# prot.core.trans <- readRDS(paste0(prefix.eqtl, "_results.RDS"))



out.dir <- paste0(get.path("results", local), "imputed/trans/new/")
prefix.eqtl <- paste0(out.dir,
                      "PRS_core_mRNA_RIN_AF_GWAS_snps_pruned_trans_eQTL")

if(file.exists(paste0(prefix.eqtl, "_results.RDS"))){
  trans.core <- readRDS(paste0(prefix.eqtl, "_results.RDS"))
}
manhattan.qtl(trans.core) +
  ggtitle("trans associations between AF GWAS SNPs (108) and leading edge genes (23)")
e.table <- trans.core$all$eqtls
print.data.frame(e.table[order(e.table$pvalue)[1:10], ],
                 row.names = F, digits = 4)
print.data.frame(e.table[e.table$FDR<0.2, ],
                 row.names = F, digits = 4)


prefix.pqtl <- paste0(out.dir,
                      "PRS_core_prot_conc_AF_GWAS_snps_pruned_trans_pQTL")
if(file.exists(paste0(prefix.pqtl, "_results.RDS"))){
  prot.core <- readRDS(paste0(prefix.pqtl, "_results.RDS"))
}
manhattan.qtl(prot.core) +
  ggtitle("trans associations between AF GWAS SNPs (108) and protein for leading edge proteins (152)")
p.table <- prot.core$all$eqtls
p.table <- p.table[order(p.table$pvalue), ]
print.data.frame(p.table[which(!duplicated(p.table$gene))[1:10], ],
                 row.names = F, digits = 4)
print.data.frame(p.table[p.table$FDR<0.2, ],
                 row.names = F, digits = 4)


setwd("~/work/symAtrial_QTL/scripts/analysis")

source("../helper/helper.R")
source("peer_tables/peer_result.R")
source("boxplots/boxplots.R")

library(R.utils)
library(GenomicRanges)
library(GenABEL)

local=F

map.seq.levels <- function(gr){
  newStyle <- mapSeqlevels(seqlevels(gr), "UCSC")
  newStyle <- newStyle[!is.na(newStyle)]
  gr <- renameSeqlevels(gr, newStyle)
  return(gr)
}

# Roselli GWAS -----------------------------------------------------------------
AF_GWAS <- readRDS(file=paste0(get.path("gwas2", local),
                               "Roselli2018_AF_HRC_GWAS_ALLv11.RDS"))
AF_GWAS <- AF_GWAS[AF_GWAS$P.value<5e-8, ]
rownames(AF_GWAS) <- AF_GWAS$MarkerName
AF_GWAS.gr <- GRanges(
  seqnames = Rle(AF_GWAS$chr),
  ranges = IRanges(AF_GWAS$pos, width = 1,
                   names = paste(AF_GWAS$MarkerName, AF_GWAS$Allele1, AF_GWAS$Allele2, sep = ":")),
  score = AF_GWAS$P.value,
  type = "AF GWAS",
  snp = AF_GWAS$MarkerName,
  rs_id = AF_GWAS$MarkerName,
  allele1 = AF_GWAS$Allele1,
  allele2 = AF_GWAS$Allele2,
  std.err = AF_GWAS$StdErr,
  beta = AF_GWAS$Effect,
  pvalue = AF_GWAS$P.value
)
AF_GWAS.gr <- map.seq.levels(AF_GWAS.gr)
saveRDS(AF_GWAS.gr,
        file = paste0(get.path("gwas2", local),
                      "Roselli2018_AF_HRC_GWAS_significant_v11_5e-8_GRanges.RDS"))
rtracklayer::export.bed(con = paste0(get.path("gwas2", local),
                                     "Roselli2018_AF_HRC_GWAS_significant_v11_5e-8.bed"),
                        object = AF_GWAS.gr, index = T)
rtracklayer::export.bed(con = paste0(get.path("gwas2", local),
                                     "Roselli2018_AF_HRC_GWAS_significant_v11_5e-8.bed"),
                        object = AF_GWAS.gr, index = F)

# significant QTLs -------------------------------------------------------------
## eQTLs ----
eqtls <- read.csv(paste0(get.path("results", local),
                         paste("imputed", "cis", "final",
                               "eQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                               sep="/")),
                  sep = "\t", h = T, stringsAsFactors = F)
eqtls.gr <- GRanges(
  seqnames = Rle(eqtls$chr),
  ranges = IRanges(eqtls$variant_pos,
                   end = eqtls$variant_pos,
                   names = paste0(eqtls$snpid, "_", eqtls$gene)),
  score = eqtls$FDR,
  type = "eQTL",
  snp = eqtls$snpid,
  rs_id = eqtls$rs_id,
  allele1 = eqtls$Allele1,
  allele2 = eqtls$Allele2,
  gene = eqtls$gene,
  gene_id = eqtls$gene_id,
  statistic = eqtls$statistic,
  pvalue = eqtls$pvalue,
  FDR = eqtls$FDR,
  beta = eqtls$beta
)
eqtls.gr <- map.seq.levels(eqtls.gr)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "eQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = eqtls.gr, index = T)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "eQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = eqtls.gr, index = F)
saveRDS(eqtls.gr,
        paste0(get.path("results", local),
               paste("imputed", "cis", "final",
                     "eQTL_right_atrial_appendage_allpairs.significant.FDR_GRanges.RDS",
                     sep="/")))

## pQTLs ----
pqtls <- read.csv(paste0(get.path("results", local),
                         paste("imputed", "cis", "final",
                               "pQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                               sep="/")),
                  sep = "\t", h = T, stringsAsFactors = F)
pqtls.gr <- GRanges(
  seqnames = Rle(pqtls$chr),
  ranges = IRanges(pqtls$variant_pos,
                   end = pqtls$variant_pos,
                   names = paste0(pqtls$snpid, "_", pqtls$gene)),
  score = pqtls$FDR,
  type = "pQTL",
  snp = pqtls$snpid,
  rs_id = pqtls$rs_id,
  allele1 = pqtls$Allele1,
  allele2 = pqtls$Allele2,
  gene = pqtls$gene,
  gene_id = pqtls$gene_id,
  statistic = pqtls$statistic,
  pvalue = pqtls$pvalue,
  FDR = pqtls$FDR,
  beta = pqtls$beta
)
pqtls.gr <- map.seq.levels(pqtls.gr)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "pQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = pqtls.gr, index = T)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "pQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = pqtls.gr, index = F)
saveRDS(pqtls.gr,
        paste0(get.path("results", local),
               paste("imputed", "cis", "final",
                     "pQTL_right_atrial_appendage_allpairs.significant.FDR_GRanges.RDS",
                     sep="/")))

## res eQTLs ----
reseqtls <- read.csv(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "res_eQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                                  sep="/")),
                     sep = "\t", h = T, stringsAsFactors = F)
reseqtls.gr <- GRanges(
  seqnames = Rle(reseqtls$chr),
  ranges = IRanges(reseqtls$variant_pos,
                   end = reseqtls$variant_pos,
                   names = paste0(reseqtls$snpid, "_", reseqtls$gene)),
  score = reseqtls$FDR,
  type = "res eQTL",
  snp = reseqtls$snpid,
  rs_id = reseqtls$rs_id,
  allele1 = reseqtls$Allele1,
  allele2 = reseqtls$Allele2,
  gene = reseqtls$gene,
  gene_id = reseqtls$gene_id,
  statistic = reseqtls$statistic,
  pvalue = reseqtls$pvalue,
  FDR = reseqtls$FDR,
  beta = reseqtls$beta
)
reseqtls.gr <- map.seq.levels(reseqtls.gr)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "res_eQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = reseqtls.gr, index = T)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "res_eQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = reseqtls.gr, index = F)
saveRDS(reseqtls.gr,
        paste0(get.path("results", local),
               paste("imputed", "cis", "final",
                     "res_eQTL_right_atrial_appendage_allpairs.significant.FDR_GRanges.RDS",
                     sep="/")))

## res pQTLs ----
respqtls <- read.csv(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "res_pQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                                  sep="/")),
                     sep = "\t", h = T, stringsAsFactors = F)
respqtls.gr <- GRanges(
  seqnames = Rle(respqtls$chr),
  ranges = IRanges(respqtls$variant_pos,
                   end = respqtls$variant_pos,
                   names = paste0(respqtls$snpid, "_", respqtls$gene)),
  score = respqtls$FDR,
  type = "eres pQTL",
  snp = respqtls$snpid,
  rs_id = respqtls$rs_id,
  allele1 = respqtls$Allele1,
  allele2 = respqtls$Allele2,
  gene = respqtls$gene,
  gene_id = respqtls$gene_id,
  statistic = respqtls$statistic,
  pvalue = respqtls$pvalue,
  FDR = respqtls$FDR,
  beta = respqtls$beta
)
respqtls.gr <- map.seq.levels(respqtls.gr)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "res_pQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = respqtls.gr, index = T)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "res_pQTL_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = respqtls.gr, index = F)
saveRDS(respqtls.gr,
        paste0(get.path("results", local),
               paste("imputed", "cis", "final",
                     "res_pQTL_right_atrial_appendage_allpairs.significant.FDR_GRanges.RDS",
                     sep="/")))

## ratioQTLs ----
ratioqtls <- read.csv(paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "ratios_right_atrial_appendage_allpairs.significant.FDR.txt",
                                   sep="/")),
                      sep = "\t", h = T, stringsAsFactors = F)
ratioqtls.gr <- GRanges(
  seqnames = Rle(ratioqtls$chr),
  ranges = IRanges(ratioqtls$variant_pos,
                   end = ratioqtls$variant_pos,
                   names = paste0(ratioqtls$snpid, "_", ratioqtls$gene)),
  score = ratioqtls$FDR,
  type = "ratioQTL",
  snp = ratioqtls$snpid,
  rs_id = ratioqtls$rs_id,
  allele1 = ratioqtls$Allele1,
  allele2 = ratioqtls$Allele2,
  gene = ratioqtls$gene,
  gene_id = ratioqtls$gene_id,
  statistic = ratioqtls$statistic,
  pvalue = ratioqtls$pvalue,
  FDR = ratioqtls$FDR,
  beta = ratioqtls$beta
)
ratioqtls.gr <- map.seq.levels(ratioqtls.gr)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "ratios_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = ratioqtls.gr, index = T)
rtracklayer::export.bed(con = paste0(get.path("results", local),
                                     paste("imputed", "cis", "final",
                                           "ratios_right_atrial_appendage_allpairs.significant.FDR.bed",
                                           sep="/")),
                        object = ratioqtls.gr, index = F)
saveRDS(ratioqtls.gr,
        paste0(get.path("results", local),
               paste("imputed", "cis", "final",
                     "ratios_right_atrial_appendage_allpairs.significant.FDR_GRanges.RDS",
                     sep="/")))

# trans QTLs -------------------------------------------------------------------

df.trans <- readRDS(paste0(get.path("results", local),
                           "imputed/trans/PRS_core_mRNA_RIN_AF_GWAS_snps_trans_QTL_results.RDS"))
df.prot <- readRDS(paste0(get.path("results", local),
                          "imputed/trans/PRS_core_prot_conc_AF_GWAS_snps_trans_QTL_results.RDS"))
df <- rbind(df.trans,
            df.prot[, colnames(df.trans)])
df[, c("FDR.eqtl", "FDR.pqtl")] <- NULL
dim(df.trans)
dim(df.prot)
dim(df)
df <- df[!duplicated(df), ]
dim(df)

df.trans2 <- df[which(df$pvalue.eqtl<0.01), ]
df.trans.bedpe <- data.frame(chrom1 = as.character(paste0("chr", df.trans2$chr)),
                             start1 = df.trans2$pos-1,
                             end1 = df.trans2$pos,
                             chrom2 = as.character(df.trans2$chr.gene),
                             start2 = df.trans2$gene.start-1,
                             end2 = df.trans2$gene.end,
                             name = paste0(df.trans2$snps, "_", df.trans2$gene),
                             score = df.trans2$pvalue.eqtl,
                             strand1 = ".",
                             strand2 = as.character(df.trans2$gene.strand),
                             stringsAsFactors = F)
write.table(df.trans.bedpe,
            file = paste0(get.path("results", local),
                          "imputed/trans/PRS_core_AF_GWAS_snps_trans_eQTLs.bedpe"),
            quote = F, row.names = F, col.names = F, sep = "\t")
df.trans.bed <- rbind(data.frame(chrom = as.character(paste0("chr", df.trans2$chr)),
                                 start = df.trans2$pos-1,
                                 end = df.trans2$pos,
                                 name = paste0(df.trans2$gene, "_", df.trans2$snps),
                                 score = df.trans2$pvalue.eqtl,
                                 strand = ".",
                                 stringsAsFactors = F),
                      data.frame(chrom = as.character(df.trans2$chr.gene),
                                 start = df.trans2$gene.start-1,
                                 end = df.trans2$gene.end,
                                 name = paste0(df.trans2$snps, "_", df.trans2$gene),
                                 score = df.trans2$pvalue.eqtl,
                                 strand = as.character(df.trans2$gene.strand),
                                 stringsAsFactors = F))
write.table(df.trans.bed,
            file = paste0(get.path("results", local),
                          "imputed/trans/PRS_core_AF_GWAS_snps_trans_eQTLs.bed"),
            quote = F, row.names = F, col.names = F, sep = "\t")

df.prot2 <- df[which(df$pvalue.pqtl<0.01), ]
df.prot.bedpe <- data.frame(chrom1 = as.character(paste0("chr", df.prot2$chr)),
                            start1 = df.prot2$pos-1,
                            end1 = df.prot2$pos,
                            chrom2 = as.character(df.prot2$chr.gene),
                            start2 = df.prot2$gene.start-1,
                            end2 = df.prot2$gene.end,
                            name = paste0(df.prot2$snps, "_", df.prot2$gene),
                            score = df.prot2$pvalue.pqtl,
                            strand1 = ".",
                            strand2 = as.character(df.prot2$gene.strand),
                            stringsAsFactors = F)
write.table(df.prot.bedpe,
            file = paste0(get.path("results", local),
                          "imputed/trans/PRS_core_AF_GWAS_snps_trans_pQTL.bedpe"),
            quote = F, row.names = F, col.names = F, sep = "\t")
df.prot.bed <- rbind(data.frame(chrom = as.character(paste0("chr", df.prot2$chr)),
                                start = df.prot2$pos-1,
                                end = df.prot2$pos,
                                name = paste0(df.prot2$gene, "_", df.prot2$snps),
                                score = df.prot2$pvalue.pqtl,
                                strand = ".",
                                stringsAsFactors = F),
                     data.frame(chrom = as.character(df.prot2$chr.gene),
                                start = df.prot2$gene.start-1,
                                end = df.prot2$gene.end,
                                name = paste0(df.prot2$snps, "_", df.prot2$gene),
                                score = df.prot2$pvalue.pqtl,
                                strand = as.character(df.prot2$gene.strand),
                                stringsAsFactors = F))
write.table(df.prot.bed,
            file = paste0(get.path("results", local),
                          "imputed/trans/PRS_core_AF_GWAS_snps_trans_pQTL.bed"),
            quote = F, row.names = F, col.names = F, sep = "\t")


# functional QTL categories ----------------------------------------------------

# QTL <- readRDS(paste0(get.path("results", local), "imputed/cis/",
#                             "QTL_res_all_snp_group_anno.RDS"))
# trans.qtl <- QTL[QTL$gene %in% QTL[QTL$snps %in% df[which(df$pvalue.eqtl<0.01 |
#                                                             df$pvalue.pqtl<0.01), "snps"], "gene"], ]
# saveRDS(trans.qtl,
#         paste0(get.path("results", local),
#                "imputed/cis/QTL_res_trans_QTL_genes_group_anno.RDS"))

sharedQTL <- readRDS(paste0(get.path("results", local), "imputed/cis/",
                            "QTL_res_shared_QTL_genes_group_anno.RDS"))[, c("snps", "gene", "chr", "snp.pos")]
ind.eQTL <- readRDS(paste0(get.path("results", local), "imputed/cis/",
                           "QTL_res_independent_eQTL_genes_group_anno.RDS"))[, c("snps", "gene", "chr", "snp.pos")]
ind.pQTL <- readRDS(paste0(get.path("results", local), "imputed/cis/",
                           "QTL_res_independent_pQTL_genes_group_anno.RDS"))[, c("snps", "gene", "chr", "snp.pos")]
transQTL <- readRDS(paste0(get.path("results", local),
                           "imputed/cis/QTL_res_trans_QTL_genes_group_anno.RDS"))[, c("snps", "gene", "chr", "snp.pos")]


# measured SNPs ----------------------------------------------------------------
snps <- rbind(sharedQTL, ind.eQTL, ind.pQTL, transQTL)
# snps$score <- snps$MAF
# snps <- snps[!duplicated(snps$snps), ]
# snps$name <- paste0(snps$snps, "_", snps$gene)
# rownames(snps) <- snps$name
snps <- snps[!duplicated(snps$snps), ]
rownames(snps) <- snps$snps
load(paste0(get.path("genotype", local), "AFHRI_B_imputed_gwaa_data.RData"))
info <- summary(AFHRI_B_imp)
info$SNP <- rownames(info)
snps[, c("score", "AA", "AB", "BB")] <- info[snps$snps, c("Q.2", "P.11", "P.12", "P.22")]
snps.gr <- GRanges(
  seqnames = snps$chr,
  ranges = IRanges(start = snps$snp.pos,
                   width = 1,
                   names = paste0(snps$snps, "_", snps$AA, "/", snps$AB, "/", snps$BB)),
  score = snps$score
)
#mcols(snps.gr, use.names = T) <- snps
snps.gr <- map.seq.levels(snps.gr)

rtracklayer::export.bed(con = paste0(get.path("results", local), "imputed/cis/",
                                     "measured_SNPs_functional_regions.bed"),
                        object = snps.gr, index = F)
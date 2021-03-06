# prepare textfiles with our QTL results
# first look into gtex_prep.R, plasma_pQTL_prep.R and mQTL_prep.R
# which IDs to add for further comparisons

setwd("~/work/symAtrial_QTL/scripts/analysis/comparison")

source("../../helper/helper.R")

local=F

# QTL results: TODO: add metabolomics, ratios
eqtl.table <- readRDS(paste0(get.path("results", local),
                             paste("imputed", "cis", "linear", "eqtl", "normalized", "factors_no_cov", "eqtl_nk12.RDS", sep="/")))$cis$eqtls
eqtl.table$snps <- as.character(eqtl.table$snps)
eqtl.table$gene <- as.character(eqtl.table$gene)

if(T){
  eqtl.all.table <- readRDS(paste0(get.path("results", local),
                                   paste("imputed", "cis", "linear",
                                         "eqtl", "normalized", "factors_all",
                                         "eqtl_nk06.RDS", sep="/")))$cis$eqtls
  eqtl.all.table$snps <- as.character(eqtl.all.table$snps)
  eqtl.all.table$gene <- as.character(eqtl.all.table$gene)
  
  eqtl.pop.table <- readRDS(paste0(get.path("results", local),
                                   paste("imputed", "cis", "linear",
                                         "eqtl", "normalized", "factors_pop",
                                         "eqtl_nk06.RDS", sep="/")))$cis$eqtls
  eqtl.pop.table$snps <- as.character(eqtl.pop.table$snps)
  eqtl.pop.table$gene <- as.character(eqtl.pop.table$gene)
}

pqtl.table <- readRDS(paste0(get.path("results", local),
                             paste("imputed", "cis", "linear", "pqtl", "normalized", "factors_no_cov", "eqtl_nk10.RDS", sep="/")))$cis$eqtls
pqtl.table$snps <- as.character(pqtl.table$snps)
pqtl.table$gene <- as.character(pqtl.table$gene)

pqtl.res.table <- readRDS(paste0(get.path("results", local),
                                 paste("imputed", "cis", "linear", "pqtl_res", "normalized", "factors_no_cov", "eqtl_nk12.RDS", sep="/")))$cis$eqtls
pqtl.res.table$snps <- as.character(pqtl.res.table$snps)
pqtl.res.table$gene <- as.character(pqtl.res.table$gene)

eqtl.res.table <- readRDS(paste0(get.path("results", local),
                                 paste("imputed", "cis", "linear", "eqtl_res", "normalized", "factors_no_cov", "eqtl_nk08.RDS", sep="/")))$cis$eqtls
eqtl.res.table$snps <- as.character(eqtl.res.table$snps)
eqtl.res.table$gene <- as.character(eqtl.res.table$gene)

ratio.qtl.table <- readRDS(paste0(get.path("results", local),
                                  paste("imputed", "cis", "linear", "ratios", "normalized", "factors_fibro", "eqtl_nk09.RDS", sep="/")))$cis$eqtls
ratio.qtl.table$snps <- as.character(ratio.qtl.table$snps)
ratio.qtl.table$gene <- as.character(ratio.qtl.table$gene)

snps <- read.table(paste0(get.path("genotype", local),
                          "AFHRI-B_imputed.bim"),
                   header=F, sep = "\t", stringsAsFactors = F)
snps$V3 <-NULL
colnames(snps) <- c("chr", "snpid", "variant_pos", "ref", "alt")
snps$variant_id <- paste(snps$chr, snps$variant_pos,
                         snps$ref, snps$alt, sep=":")
snps$rs_id <- NA
snps$rs_id[grep("rs", snps$snpid)] <- gsub("\\:.*", "", snps$snpid[grep("rs", snps$snpid)])

gtex.ids <- readRDS(paste0(get.path("gtex", local),
                           "tmp/match_gtex_AFHRI_B_SNPs.RDS"))
# gtex.genes <- readRDS(paste0(get.path("gtex", local),
#                              "tmp/ensembl_ids_AFHRIB_genes.RDS"))
ensembl.ids <- readRDS(paste0(get.path("locations", local),
                              "ensembl_ids_AFHRIB_QTLgenes.RDS"))
colnames(ensembl.ids)
# [1] "ensembl_gene_id"    "chromosome_name"    "start_position"    
# [4] "end_position"       "hgnc_symbol"        "external_gene_name"

# add also gtex anno
# GTEx eQTL file header:
# gene_id	variant_id	tss_distance	ma_samples	ma_count	maf	pval_nominal	slope	slope_se
# ENSG00000177757.1	1_54490_G_A_b37	-698261	68	72	0.154506	0.66026	-0.0509269	0.115697
eqtl.table.gtex <- merge(snps, eqtl.table,
                         by.x=c("snpid"), by.y=c("snps"),
                         all.y = T)
eqtl.table.gtex <- merge(eqtl.table.gtex, gtex.ids,
                         by.x=c("snpid"), by.y=c("snpid"),
                         all.x = T)
# eqtl.table.gtex <- merge(eqtl.table.gtex, gtex.genes[, c("gene_id", "symbol")],
#                          by.x=c("gene"), by.y=c("symbol"),
#                          all.x = T)
eqtl.table.gtex <- merge(eqtl.table.gtex, ensembl.ids[, c("ensembl_gene_id", "external_gene_name")],
                         by.x=c("gene"), by.y=c("external_gene_name"),
                         all.x = T)
colnames(eqtl.table.gtex) <- c("gene", "snpid", "chr", "variant_pos", "Allele1", "Allele2", "variant_id", "rs_id",
                               "statistic", "pvalue", "FDR", "beta",
                               "gtex.variant_id", "gtex.match_alleles", "gene_id")

write.table(eqtl.table.gtex,
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "eQTL_right_atrial_appendage_allpairs.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(eqtl.table.gtex[eqtl.table.gtex$pvalue < 1e-5, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "eQTL_right_atrial_appendage_allpairs.significant.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(eqtl.table.gtex[eqtl.table.gtex$FDR < 0.05, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "eQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
saveRDS(eqtl.table.gtex,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "eQTL_right_atrial_appendage_allpairs.RDS",
                          sep="/")))
head(eqtl.table.gtex)
dim(eqtl.table.gtex)
dim(eqtl.table)
rm(eqtl.table, eqtl.table.gtex)

if(T){
  eqtl.all.table.gtex <- merge(snps, eqtl.all.table,
                               by.x=c("snpid"), by.y=c("snps"),
                               all.y = T)
  eqtl.all.table.gtex <- merge(eqtl.all.table.gtex, gtex.ids,
                               by.x=c("snpid"), by.y=c("snpid"),
                               all.x = T)
  eqtl.all.table.gtex <- merge(eqtl.all.table.gtex, ensembl.ids[, c("ensembl_gene_id", "external_gene_name")],
                               by.x=c("gene"), by.y=c("external_gene_name"),
                               all.x = T)
  colnames(eqtl.all.table.gtex) <- c("gene", "snpid", "chr", "variant_pos", "Allele1", "Allele2", "variant_id", "rs_id",
                                     "statistic", "pvalue", "FDR", "beta",
                                     "gtex.variant_id", "gtex.match_alleles", "gene_id")
  
  write.table(eqtl.all.table.gtex,
              file=paste0(get.path("results", local),
                          paste("imputed", "cis", "final",
                                "eQTL_cov_all_right_atrial_appendage_allpairs.txt",
                                sep="/")),
              row.names=F, col.names=T, quote=F, sep="\t")
  write.table(eqtl.all.table.gtex[eqtl.all.table.gtex$pvalue < 1e-5, ],
              file=paste0(get.path("results", local),
                          paste("imputed", "cis", "final",
                                "eQTL_cov_all_right_atrial_appendage_allpairs.significant.txt",
                                sep="/")),
              row.names=F, col.names=T, quote=F, sep="\t")
  write.table(eqtl.all.table.gtex[eqtl.all.table.gtex$FDR < 0.05, ],
              file=paste0(get.path("results", local),
                          paste("imputed", "cis", "final",
                                "eQTL_cov_all_right_atrial_appendage_allpairs.significant.FDR.txt",
                                sep="/")),
              row.names=F, col.names=T, quote=F, sep="\t")
  saveRDS(eqtl.all.table.gtex,
          file=paste0(get.path("results", local),
                      paste("imputed", "cis", "final",
                            "eQTL_cov_all_right_atrial_appendage_allpairs.RDS",
                            sep="/")))
  head(eqtl.all.table.gtex)
  dim(eqtl.all.table.gtex)
  dim(eqtl.all.table)
  rm(eqtl.all.table, eqtl.all.table.gtex)
  
  
  eqtl.pop.table.gtex <- merge(snps, eqtl.pop.table,
                               by.x=c("snpid"), by.y=c("snps"),
                               all.y = T)
  eqtl.pop.table.gtex <- merge(eqtl.pop.table.gtex, gtex.ids,
                               by.x=c("snpid"), by.y=c("snpid"),
                               all.x = T)
  eqtl.pop.table.gtex <- merge(eqtl.pop.table.gtex, ensembl.ids[, c("ensembl_gene_id", "external_gene_name")],
                               by.x=c("gene"), by.y=c("external_gene_name"),
                               all.x = T)
  colnames(eqtl.pop.table.gtex) <- c("gene", "snpid", "chr", "variant_pos", "Allele1", "Allele2", "variant_id", "rs_id",
                                     "statistic", "pvalue", "FDR", "beta",
                                     "gtex.variant_id", "gtex.match_alleles", "gene_id")
  
  write.table(eqtl.pop.table.gtex,
              file=paste0(get.path("results", local),
                          paste("imputed", "cis", "final",
                                "eQTL_cov_pop_right_atrial_appendage_allpairs.txt",
                                sep="/")),
              row.names=F, col.names=T, quote=F, sep="\t")
  write.table(eqtl.pop.table.gtex[eqtl.pop.table.gtex$pvalue < 1e-5, ],
              file=paste0(get.path("results", local),
                          paste("imputed", "cis", "final",
                                "eQTL_cov_pop_right_atrial_appendage_allpairs.significant.txt",
                                sep="/")),
              row.names=F, col.names=T, quote=F, sep="\t")
  write.table(eqtl.pop.table.gtex[eqtl.pop.table.gtex$FDR < 0.05, ],
              file=paste0(get.path("results", local),
                          paste("imputed", "cis", "final",
                                "eQTL_cov_pop_right_atrial_appendage_allpairs.significant.FDR.txt",
                                sep="/")),
              row.names=F, col.names=T, quote=F, sep="\t")
  saveRDS(eqtl.pop.table.gtex,
          file=paste0(get.path("results", local),
                      paste("imputed", "cis", "final",
                            "eQTL_cov_pop_right_atrial_appendage_allpairs.RDS",
                            sep="/")))
  head(eqtl.pop.table.gtex)
  dim(eqtl.pop.table.gtex)
  dim(eqtl.pop.table)
  rm(eqtl.pop.table, eqtl.pop.table.gtex)
}

# add gtex and plasma pQTL anno
# plasma pQTL file header:
# SOMAMER_ID	VARIANT_ID	chromosome	position	Allele1	Allele2	Effect	StdErr	log.P.	variant_id	variant_id2	plasma.match_alleles	gene
# ACP1.3858.5.1	2_11320_23	2	11320	a	g	-0.4611	0.0268	-65.8	2:11320:G:A	2:11320:A:G	-1	ACP1
pqtl.table.gtex <- merge(snps, pqtl.table,
                         by.x=c("snpid"), by.y=c("snps"),
                         all.y = T)
pqtl.table.gtex <- merge(pqtl.table.gtex, gtex.ids,
                         by.x=c("snpid"), by.y=c("snpid"),
                         all.x = T)
pqtl.table.gtex <- merge(pqtl.table.gtex, ensembl.ids[, c("ensembl_gene_id", "external_gene_name")],
                         by.x=c("gene"), by.y=c("external_gene_name"),
                         all.x = T)
colnames(pqtl.table.gtex) <- c("gene", "snpid", "chr", "variant_pos", "Allele1", "Allele2", "variant_id", "rs_id",
                               "statistic", "pvalue", "FDR", "beta",
                               "gtex.variant_id", "gtex.match_alleles", "gene_id")

write.table(pqtl.table.gtex,
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "pQTL_right_atrial_appendage_allpairs.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(pqtl.table.gtex[pqtl.table.gtex$pvalue < 1e-5, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "pQTL_right_atrial_appendage_allpairs.significant.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(pqtl.table.gtex[pqtl.table.gtex$FDR < 0.05, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "pQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
saveRDS(pqtl.table.gtex,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "pQTL_right_atrial_appendage_allpairs.RDS",
                          sep="/")))
head(pqtl.table.gtex)
dim(pqtl.table.gtex)
dim(pqtl.table)
rm(pqtl.table, pqtl.table.gtex)


## residual eQTLs
eqtl.res.table.gtex <- merge(snps, eqtl.res.table,
                             by.x=c("snpid"), by.y=c("snps"),
                             all.y = T)
eqtl.res.table.gtex <- merge(eqtl.res.table.gtex, gtex.ids,
                             by.x=c("snpid"), by.y=c("snpid"),
                             all.x = T)
eqtl.res.table.gtex <- merge(eqtl.res.table.gtex, ensembl.ids[, c("ensembl_gene_id", "external_gene_name")],
                             by.x=c("gene"), by.y=c("external_gene_name"),
                             all.x = T)
colnames(eqtl.res.table.gtex) <- c("gene", "snpid", "chr", "variant_pos", "Allele1", "Allele2", "variant_id", "rs_id",
                                   "statistic", "pvalue", "FDR", "beta",
                                   "gtex.variant_id", "gtex.match_alleles", "gene_id")

write.table(eqtl.res.table.gtex,
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "res_eQTL_right_atrial_appendage_allpairs.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(eqtl.res.table.gtex[eqtl.res.table.gtex$pvalue < 1e-5, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "res_eQTL_right_atrial_appendage_allpairs.significant.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(eqtl.res.table.gtex[eqtl.res.table.gtex$FDR < 0.05, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "res_eQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
saveRDS(eqtl.res.table.gtex,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "res_eQTL_right_atrial_appendage_allpairs.RDS",
                          sep="/")))
head(eqtl.res.table.gtex)
dim(eqtl.res.table.gtex)
dim(eqtl.res.table)
rm(eqtl.res.table.gtex, eqtl.res.table)


## residual pQTLs
pqtl.res.table.gtex <- merge(snps, pqtl.res.table,
                             by.x=c("snpid"), by.y=c("snps"),
                             all.y = T)
pqtl.res.table.gtex <- merge(pqtl.res.table.gtex, gtex.ids,
                             by.x=c("snpid"), by.y=c("snpid"),
                             all.x = T)
pqtl.res.table.gtex <- merge(pqtl.res.table.gtex, ensembl.ids[, c("ensembl_gene_id", "external_gene_name")],
                             by.x=c("gene"), by.y=c("external_gene_name"),
                             all.x = T)
colnames(pqtl.res.table.gtex) <- c("gene", "snpid", "chr", "variant_pos", "Allele1", "Allele2", "variant_id", "rs_id",
                                   "statistic", "pvalue", "FDR", "beta",
                                   "gtex.variant_id", "gtex.match_alleles", "gene_id")

write.table(pqtl.res.table.gtex,
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "res_pQTL_right_atrial_appendage_allpairs.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(pqtl.res.table.gtex[pqtl.res.table.gtex$pvalue < 1e-5, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "res_pQTL_right_atrial_appendage_allpairs.significant.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(pqtl.res.table.gtex[pqtl.res.table.gtex$FDR < 0.05, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "res_pQTL_right_atrial_appendage_allpairs.significant.FDR.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
saveRDS(pqtl.res.table.gtex,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "res_pQTL_right_atrial_appendage_allpairs.RDS",
                          sep="/")))
head(pqtl.res.table.gtex)
dim(pqtl.res.table.gtex)
dim(pqtl.res.table)
rm(pqtl.res.table.gtex, pqtl.res.table)


## ratio QTLs
ratio.qtl.table.gtex <- merge(snps, ratio.qtl.table,
                              by.x=c("snpid"), by.y=c("snps"),
                              all.y = T)
ratio.qtl.table.gtex <- merge(ratio.qtl.table.gtex, gtex.ids,
                              by.x=c("snpid"), by.y=c("snpid"),
                              all.x = T)
ratio.qtl.table.gtex <- merge(ratio.qtl.table.gtex, ensembl.ids[, c("ensembl_gene_id", "external_gene_name")],
                              by.x=c("gene"), by.y=c("external_gene_name"),
                              all.x = T)
colnames(ratio.qtl.table.gtex) <- c("gene", "snpid", "chr", "variant_pos", "Allele1", "Allele2", "variant_id", "rs_id",
                                    "statistic", "pvalue", "FDR", "beta",
                                    "gtex.variant_id", "gtex.match_alleles", "gene_id")

write.table(ratio.qtl.table.gtex,
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "ratios_right_atrial_appendage_allpairs.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(ratio.qtl.table.gtex[ratio.qtl.table.gtex$pvalue < 1e-5, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "ratios_right_atrial_appendage_allpairs.significant.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
write.table(ratio.qtl.table.gtex[ratio.qtl.table.gtex$FDR < 0.05, ],
            file=paste0(get.path("results", local),
                        paste("imputed", "cis", "final",
                              "ratios_right_atrial_appendage_allpairs.significant.FDR.txt",
                              sep="/")),
            row.names=F, col.names=T, quote=F, sep="\t")
saveRDS(ratio.qtl.table.gtex,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "ratios_right_atrial_appendage_allpairs.RDS",
                          sep="/")))
head(ratio.qtl.table.gtex)
dim(ratio.qtl.table.gtex)
dim(ratio.qtl.table)
rm(ratio.qtl.table.gtex, ratio.qtl.table)



# Source data - QTL clumps ------

setwd("~/work/symAtrial_QTL/scripts/analysis/comparison")

source("../../helper/helper.R")

local=F

# add LD clumps eQTL ----
pairs <- readRDS(file=paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "eQTL_right_atrial_appendage_allpairs.RDS",
                                   sep="/")))
clump <- readRDS(file = paste0(get.path("results", local),
                               "imputed/cis/eQTL_clump_relaxed.RDS"))

pairs$LDclumpFDR <- NA
genes <- unique(pairs[pairs$FDR < 0.05, "gene"])
for(gene in genes){
  clump.info <- clump[["FDR"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpFDR"] <- lead
  }
}

pairs$LDclumpP <- NA
genes <- unique(pairs[pairs$pvalue < 1e-5, "gene"])
for(gene in genes){
  clump.info <- clump[["P"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpP"] <- lead
  }
}

saveRDS(pairs,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "eQTL_right_atrial_appendage_allpairs_clump.RDS",
                          sep="/")))
write.table(pairs,
            file=paste0(get.path("results", local), "/source_data/",
                        "eQTL_right_atrial_appendage_allpairs_clump.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")


# add LD clumps pQTL ----
pairs <- readRDS(file=paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "pQTL_right_atrial_appendage_allpairs.RDS",
                                   sep="/")))
clump <- readRDS(file = paste0(get.path("results", local),
                               "imputed/cis/pQTL_clump_relaxed.RDS"))

pairs$LDclumpFDR <- NA
genes <- unique(pairs[pairs$FDR < 0.05, "gene"])
for(gene in genes){
  clump.info <- clump[["FDR"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpFDR"] <- lead
  }
}

pairs$LDclumpP <- NA
genes <- unique(pairs[pairs$pvalue < 1e-5, "gene"])
for(gene in genes){
  clump.info <- clump[["P"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpP"] <- lead
  }
}

saveRDS(pairs,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "pQTL_right_atrial_appendage_allpairs_clump.RDS",
                          sep="/")))
write.table(pairs,
            file=paste0(get.path("results", local), "/source_data/",
                        "pQTL_right_atrial_appendage_allpairs_clump.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")


# add LD clumps ratioQTL ----
pairs <- readRDS(file=paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "ratios_right_atrial_appendage_allpairs.RDS",
                                   sep="/")))
clump <- readRDS(file = paste0(get.path("results", local),
                               "imputed/cis/ratioQTL_clump_relaxed.RDS"))

pairs$LDclumpFDR <- NA
genes <- unique(pairs[pairs$FDR < 0.05, "gene"])
for(gene in genes){
  clump.info <- clump[["FDR"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpFDR"] <- lead
  }
}

pairs$LDclumpP <- NA
genes <- unique(pairs[pairs$pvalue < 1e-5, "gene"])
for(gene in genes){
  clump.info <- clump[["P"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpP"] <- lead
  }
}

saveRDS(pairs,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "ratioQTL_right_atrial_appendage_allpairs_clump.RDS",
                          sep="/")))
write.table(pairs,
            file=paste0(get.path("results", local), "/source_data/",
                        "ratioQTL_right_atrial_appendage_allpairs_clump.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")


# add LD clumps res eQTL ----
pairs <- readRDS(file=paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "res_eQTL_right_atrial_appendage_allpairs.RDS",
                                   sep="/")))
clump <- readRDS(file = paste0(get.path("results", local),
                               "imputed/cis/res_eQTL_clump_relaxed.RDS"))

pairs$LDclumpFDR <- NA
genes <- unique(pairs[pairs$FDR < 0.05, "gene"])
for(gene in genes){
  clump.info <- clump[["FDR"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpFDR"] <- lead
  }
}

pairs$LDclumpP <- NA
genes <- unique(pairs[pairs$pvalue < 1e-5, "gene"])
for(gene in genes){
  clump.info <- clump[["P"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpP"] <- lead
  }
}

saveRDS(pairs,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "res_eQTL_right_atrial_appendage_allpairs_clump.RDS",
                          sep="/")))
write.table(pairs,
            file=paste0(get.path("results", local), "/source_data/",
                        "res_eQTL_right_atrial_appendage_allpairs_clump.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")


# add LD clumps res pQTL ----
pairs <- readRDS(file=paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "res_pQTL_right_atrial_appendage_allpairs.RDS",
                                   sep="/")))
clump <- readRDS(file = paste0(get.path("results", local),
                               "imputed/cis/res_pQTL_clump_relaxed.RDS"))

pairs$LDclumpFDR <- NA
genes <- unique(pairs[pairs$FDR < 0.05, "gene"])
for(gene in genes){
  clump.info <- clump[["FDR"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpFDR"] <- lead
  }
}

pairs$LDclumpP <- NA
genes <- unique(pairs[pairs$pvalue < 1e-5, "gene"])
for(gene in genes){
  clump.info <- clump[["P"]][[gene]]
  for(i in 1:dim(clump.info)[1]){
    lead <- clump.info[i, "SNP"]
    snps <- c(lead, clump.info[i, "SP2"][[1]])
    pairs[pairs$gene == gene & pairs$snpid %in% snps, "LDclumpP"] <- lead
  }
}

saveRDS(pairs,
        file=paste0(get.path("results", local),
                    paste("imputed", "cis", "final",
                          "res_pQTL_right_atrial_appendage_allpairs_clump.RDS",
                          sep="/")))
write.table(pairs,
            file=paste0(get.path("results", local), "/source_data/",
                        "res_pQTL_right_atrial_appendage_allpairs_clump.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")

q(save="no")



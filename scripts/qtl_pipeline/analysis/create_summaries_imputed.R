setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")
# setwd("/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/scripts")

library(ggplot2)
library(ggpubr)

source("helper/helper.R")
source("analysis/peer_tables/peer_result.R")
local=F

QTL.res <- readRDS(file="/home/icb/ines.assum/projects/symAtrial_QTL/results/current/peer/QTL_table_imp.RDS")

# Load QTL results ----
# colnames(MatrixEQTLtable)
# [1] "snps"      "gene"      "statistic" "pvalue"    "FDR"       "beta"

if(F){
  eqtl.table <- readRDS(paste0(get.path("results", local),
                               paste("imputed", "cis", "linear", "eqtl", "normalized", "factors_no_cov", "eqtl_nk12.RDS", sep="/")))$cis$eqtls
  eqtl.table$snps <- as.character(eqtl.table$snps)
  eqtl.table$gene <- as.character(eqtl.table$gene)
  eqtl.table$eQTL <- as.numeric(as.logical(eqtl.table$FDR<0.05))
  
  pqtl.table <- readRDS(paste0(get.path("results", local),
                               paste("imputed", "cis", "linear", "pqtl", "normalized", "factors_no_cov", "eqtl_nk10.RDS", sep="/")))$cis$eqtls
  pqtl.table$snps <- as.character(pqtl.table$snps)
  pqtl.table$gene <- as.character(pqtl.table$gene)
  pqtl.table$pQTL <- as.numeric(as.logical(pqtl.table$FDR<0.05))
  
  eqtl.res.table <- readRDS(paste0(get.path("results", local),
                                   paste("imputed", "cis", "linear", "eqtl_res", "normalized", "factors_no_cov", "eqtl_nk08.RDS", sep="/")))$cis$eqtls
  eqtl.res.table$snps <- as.character(eqtl.res.table$snps)
  eqtl.res.table$gene <- as.character(eqtl.res.table$gene)
  eqtl.res.table$reseQTL <- as.numeric(as.logical(eqtl.res.table$FDR<0.05))
  
  pqtl.res.table <- readRDS(paste0(get.path("results", local),
                                   paste("imputed", "cis", "linear", "pqtl_res", "normalized", "factors_no_cov", "eqtl_nk12.RDS", sep="/")))$cis$eqtls
  pqtl.res.table$snps <- as.character(pqtl.res.table$snps)
  pqtl.res.table$gene <- as.character(pqtl.res.table$gene)
  pqtl.res.table$respQTL <- as.numeric(as.logical(pqtl.res.table$FDR<0.05))
  
  ratio.qtl.table <- readRDS(paste0(get.path("results", local),
                                    paste("imputed", "cis", "linear", "ratios", "normalized", "factors_fibro", "eqtl_nk09.RDS", sep="/")))$cis$eqtls
  ratio.qtl.table$snps <- as.character(ratio.qtl.table$snps)
  ratio.qtl.table$gene <- as.character(ratio.qtl.table$gene)
  ratio.qtl.table$ratioQTL <- as.numeric(as.logical(ratio.qtl.table$FDR<0.05))
  
  library(plyr)
  eqtl.table <- rename(eqtl.table, c("pvalue"="pval.eqtl", "FDR"="FDR.eqtl", "beta"="beta.eqtl"))
  pqtl.table <- rename(pqtl.table, c("pvalue"="pval.pqtl", "FDR"="FDR.pqtl", "beta"="beta.pqtl"))
  eqtl.res.table <- rename(eqtl.res.table, c("pvalue"="pval.res_eqtl", "FDR"="FDR.res_eqtl", "beta"="beta.res_eqtl"))
  pqtl.res.table <- rename(pqtl.res.table, c("pvalue"="pval.res_pqtl", "FDR"="FDR.res_pqtl", "beta"="beta.res_pqtl"))
  ratio.qtl.table <- rename(ratio.qtl.table, c("pvalue"="pval.ratios", "FDR"="FDR.ratios", "beta"="beta.ratios"))
  
  QTL.snp <- merge(eqtl.table[, -3],
                   pqtl.table[, -3],
                   by=c("snps", "gene"),
                   all=T, sort=F)
  QTL.snp <- merge(QTL.snp,
                   eqtl.res.table[, -3],
                   by=c("snps", "gene"),
                   all=T, sort=F)
  QTL.snp <- merge(QTL.snp,
                   pqtl.res.table[, -3],
                   by=c("snps", "gene"),
                   all=T, sort=F)
  QTL.snp$rs_id <- NA
  QTL.snp$rs_id[grep("rs", QTL.snp$snps)] <- gsub("\\:.*", "", QTL.snp$snps[grep("rs", QTL.snp$snps)])
  saveRDS(QTL.snp, file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_temp.RDS"))
  
  QTL.snp <- merge(QTL.snp,
                   ratio.qtl.table[, -3],
                   by=c("snps", "gene"),
                   all=T, sort=F)
  QTL.snp$rs_id <- NA
  QTL.snp$rs_id[grep("rs", QTL.snp$snps)] <- gsub("\\:.*", "", QTL.snp$snps[grep("rs", QTL.snp$snps)])
  QTL.snp$snps <- as.character(QTL.snp$snps)
  QTL.snp$gene <- as.character(QTL.snp$gene)
  saveRDS(QTL.snp, file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp.RDS"))
}

# Load merged results: ----

# > colnames(QTL.snp)
# [1] "snps"          "gene"          "pval.eqtl"     "FDR.eqtl"     
# [5] "beta.eqtl"     "eQTL"          "pval.pqtl"     "FDR.pqtl"     
# [9] "beta.pqtl"     "pQTL"          "pval.res_eqtl" "FDR.res_eqtl" 
# [13] "beta.res_eqtl" "reseQTL"       "pval.res_pqtl" "FDR.res_pqtl" 
# [17] "beta.res_pqtl" "respQTL"       "rs_id"         "pval.ratios"  
# [21] "FDR.ratios"    "beta.ratios"   "ratioQTL"

QTL.snp <- readRDS(paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp.RDS"))

locations_p <- readRDS(paste0(get.path("locations"), "protein_locations_ensembl_tss.RDS"))
locations_t <- readRDS(paste0(get.path("locations"), "transcript_locations_ensembl_tss.RDS"))
snp.loc <- read.table(paste0(get.path("genotype", local),
                          "AFHRI-B_imputed.bim"),
                   header=F, sep = "\t", stringsAsFactors = F)
colnames(snp.loc)
snp.loc$V3 <-NULL
colnames(snp.loc) <- c("chr", "snps", "snp.pos", "Allele1", "Allele2")
snp.loc$variant_id <- paste(snp.loc$chr, snp.loc$snp.pos,
                            snp.loc$Allele1, snp.loc$Allele2, sep=":")
locations <- merge(locations_t,
                   locations_p,
                   all=T, sort=F)

QTL.snp2 <- merge(QTL.snp, snp.loc,
                  by="snps",
                  all.x=T, sort=F)
QTL.snp3 <- merge(QTL.snp2, locations,
                  by.x="gene", by.y="name",
                  all.x=T, sort=F)

identical(as.character(QTL.snp3$chr), as.character(QTL.snp3$chromosome))
# colnames(QTL.snp3)
# [1] "gene"           "snps"           "pval.eqtl"      "FDR.eqtl"      
# [5] "beta.eqtl"      "eQTL"           "pval.pqtl"      "FDR.pqtl"      
# [9] "beta.pqtl"      "pQTL"           "pval.res_eqtl"  "FDR.res_eqtl"  
# [13] "beta.res_eqtl"  "reseQTL"        "pval.res_pqtl"  "FDR.res_pqtl"  
# [17] "beta.res_pqtl"  "respQTL"        "rs_id"          "pval.ratios"   
# [21] "FDR.ratios"     "beta.ratios"    "ratioQTL"       "chr"           
# [25] "snp.pos"        "Allele1"        "Allele2"        "variant_id"    
# [29] "chromosome"     "start_location" "end_location"   "strand"        
# [33] "tss"       

library(biomaRt)
genesQTL <- data.frame(symbol=unique(QTL.snp3$gene),
                    stringsAsFactors = F)

#GRCh37.p13
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# filter for desired gene ids
anno <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                           'start_position','end_position','hgnc_symbol',
                           'external_gene_name'),
              filters = 'external_gene_name', values=genesQTL$symbol,
              mart = ensembl)
anno <- anno[anno$chromosome_name %in% as.character(c(1:26)), ]
match <- anno$external_gene_name==anno$hgnc_symbol
table(genesQTL$symbol %in% anno$external_gene_name[match])
miss <- genesQTL$symbol[!(genesQTL$symbol %in% anno$external_gene_name[match])]
anno2 <- anno[anno$external_gene_name %in% miss, ]
anno <- rbind(anno[match, ],
              anno2)
colnames(anno)
QTL.snp4 <- merge(QTL.snp3, anno[, c("ensembl_gene_id", "external_gene_name")],
                  by.x="gene", by.y="external_gene_name",
                  all.x=T, sort=F)
# [1] "gene"            "snps"            "pval.eqtl"       "FDR.eqtl"       
# [5] "beta.eqtl"       "eQTL"            "pval.pqtl"       "FDR.pqtl"       
# [9] "beta.pqtl"       "pQTL"            "pval.res_eqtl"   "FDR.res_eqtl"   
# [13] "beta.res_eqtl"   "reseQTL"         "pval.res_pqtl"   "FDR.res_pqtl"   
# [17] "beta.res_pqtl"   "respQTL"         "rs_id"           "pval.ratios"    
# [21] "FDR.ratios"      "beta.ratios"     "ratioQTL"        "chr"            
# [25] "snp.pos"         "Allele1"         "Allele2"         "variant_id"     
# [29] "chromosome"      "start_location"  "end_location"    "strand"         
# [33] "tss"             "ensembl_gene_id"
QTL.snp5 <- QTL.snp4[, c("snps", "rs_id", "variant_id", "gene", "ensembl_gene_id",
                         "chr", "snp.pos", "Allele1", "Allele2",
                         "start_location", "end_location", "strand", "tss",
                         "pval.eqtl", "FDR.eqtl", "beta.eqtl",
                         "pval.pqtl", "FDR.pqtl", "beta.pqtl",
                         "pval.res_eqtl", "FDR.res_eqtl", "beta.res_eqtl",
                         "pval.res_pqtl", "FDR.res_pqtl", "beta.res_pqtl",
                         "pval.ratios", "FDR.ratios", "beta.ratios",
                         "eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL")]
colnames(QTL.snp5) <- c("snps", "rs_id", "variant_id", "gene", "ensembl_id",
                        "chr", "snp.pos", "Allele1", "Allele2",
                        "gene.start", "gene.end", "strand", "TSS",
                        "pval.eqtl", "FDR.eqtl", "beta.eqtl",
                        "pval.pqtl", "FDR.pqtl", "beta.pqtl",
                        "pval.res_eqtl", "FDR.res_eqtl", "beta.res_eqtl",
                        "pval.res_pqtl", "FDR.res_pqtl", "beta.res_pqtl",
                        "pval.ratios", "FDR.ratios", "beta.ratios",
                        "eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL")
QTL.snp <- QTL.snp5
saveRDS(QTL.snp, file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_anno.RDS"))


QTL.snp$Group1 <- as.numeric(as.logical(QTL.snp$FDR.eqtl<0.05 &
                                          QTL.snp$FDR.pqtl>0.05 &
                                          QTL.snp$FDR.res_eqtl<0.05 &
                                          QTL.snp$FDR.res_pqtl>0.05))
QTL.snp$Group2 <- as.numeric(as.logical(QTL.snp$FDR.eqtl<0.05 &
                                          QTL.snp$FDR.pqtl<0.05 &
                                          QTL.snp$FDR.res_eqtl>0.05 &
                                          QTL.snp$FDR.res_pqtl>0.05))
QTL.snp$Group2b <- as.numeric(as.logical(QTL.snp$FDR.eqtl<0.05 &
                                           QTL.snp$FDR.pqtl<0.05))
QTL.snp$Group3 <- as.numeric(as.logical(QTL.snp$FDR.eqtl<0.05 &
                                          QTL.snp$FDR.pqtl>0.05 &
                                          QTL.snp$FDR.res_eqtl>0.05 &
                                          QTL.snp$FDR.res_pqtl>0.05))
QTL.snp$Group4 <- as.numeric(as.logical(QTL.snp$FDR.eqtl>0.05 &
                                          QTL.snp$FDR.pqtl<0.05 &
                                          QTL.snp$FDR.res_eqtl>0.05 &
                                          QTL.snp$FDR.res_pqtl>0.05))
QTL.snp$Group5 <- as.numeric(as.logical(QTL.snp$FDR.eqtl>0.05 &
                                          QTL.snp$FDR.pqtl<0.05 &
                                          QTL.snp$FDR.res_eqtl>0.05 &
                                          QTL.snp$FDR.res_pqtl<0.05))
QTL.snp$Group6 <- as.numeric(as.logical(QTL.snp$FDR.eqtl>0.05 &
                                          QTL.snp$FDR.pqtl>0.05 &
                                          QTL.snp$FDR.res_eqtl>0.05 &
                                          QTL.snp$FDR.res_pqtl<0.05))
QTL.snp$Group7 <- as.numeric(as.logical(QTL.snp$FDR.res_eqtl<0.05 &
                                          QTL.snp$FDR.res_pqtl<0.05))
QTL.snp$Group8 <- as.numeric(as.logical(QTL.snp$FDR.eqtl>0.05 &
                                          QTL.snp$FDR.pqtl>0.05 &
                                          QTL.snp$FDR.res_eqtl<0.05 &
                                          QTL.snp$FDR.res_pqtl>0.05))
QTL.snp$Group9 <- as.numeric(as.logical(QTL.snp$FDR.eqtl<0.05 &
                                          QTL.snp$FDR.pqtl<0.05 &
                                          QTL.snp$FDR.res_eqtl<0.05 &
                                          QTL.snp$FDR.res_pqtl<0.05))
saveRDS(QTL.snp, file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_group.RDS"))



res.snp.group <- colSums(QTL.snp[, c("eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL",
                               "Group1", "Group2", "Group2b", "Group3", "Group4",
                               "Group5", "Group6", "Group7", "Group8", "Group9")],
                   na.rm = T)

hits.eQTL <- unique(QTL.snp[which(QTL.snp$eQTL==1), "gene"])
hits.pQTL <- unique(QTL.snp[which(QTL.snp$pQTL==1), "gene"])
hits.reseQTL <- unique(QTL.snp[which(QTL.snp$reseQTL==1), "gene"])
hits.respQTL <- unique(QTL.snp[which(QTL.snp$respQTL==1), "gene"])
hits.ratioQTL <- unique(QTL.snp[which(QTL.snp$ratioQTL==1), "gene"])
hits.Group1 <- unique(QTL.snp[which(QTL.snp$Group1==1), "gene"])
hits.Group2 <- unique(QTL.snp[which(QTL.snp$Group2==1), "gene"])
hits.Group2b <- unique(QTL.snp[which(QTL.snp$Group2b==1), "gene"])
hits.Group3 <- unique(QTL.snp[which(QTL.snp$Group3==1), "gene"])
hits.Group4 <- unique(QTL.snp[which(QTL.snp$Group4==1), "gene"])
hits.Group5 <- unique(QTL.snp[which(QTL.snp$Group5==1), "gene"])
hits.Group6 <- unique(QTL.snp[which(QTL.snp$Group6==1), "gene"])
hits.Group7 <- unique(QTL.snp[which(QTL.snp$Group7==1), "gene"])
hits.Group8 <- unique(QTL.snp[which(QTL.snp$Group8==1), "gene"])
hits.Group9 <- unique(QTL.snp[which(QTL.snp$Group9==1), "gene"])
res.gene.group <- c(length(hits.eQTL),
                   length(hits.pQTL),
                   length(hits.reseQTL),
                   length(hits.respQTL),
                   length(hits.ratioQTL),
                   length(hits.Group1),
                   length(hits.Group2),
                   length(hits.Group2b),
                   length(hits.Group3),
                   length(hits.Group4),
                   length(hits.Group5),
                   length(hits.Group6),
                   length(hits.Group7),
                   length(hits.Group8),
                   length(hits.Group9))
names(res.gene.group) <- c("eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL",
                          "Group1", "Group2", "Group2b", "Group3", "Group4",
                          "Group5", "Group6", "Group7", "Group8", "Group9")

res.gene <- list(eQTL = hits.eQTL,
                 pQTL = hits.pQTL,
                 reseQTL = hits.reseQTL,
                 respQTL = hits.respQTL,
                 ratioQTL = hits.ratioQTL,
                 Group1 = hits.Group1,
                 Group2 = hits.Group2,
                 Group2b = hits.Group2b,
                 Group3 = hits.Group3,
                 Group4 = hits.Group4,
                 Group5 = hits.Group5,
                 Group6 = hits.Group6,
                 Group7 = hits.Group7,
                 Group8 = hits.Group8,
                 Group9 = hits.Group9)
saveRDS(res.gene,
        file = paste0(get.path("results", local), "imputed/cis/QTL_res_genes_per_group.RDS"))

res.snp.group
res.gene.group

colSums(!is.na(QTL.snp[, c("pval.eqtl", "pval.pqtl",
                           "pval.res_eqtl", "pval.res_pqtl",
                           "pval.ratios")]))

saveRDS(data.frame(rbind(res.gene.group,res.snp.group),
                   row.names = c("genes", "snps")),
        file = paste0(get.path("results", local), "imputed/cis/QTL_res_group_summary.RDS"))








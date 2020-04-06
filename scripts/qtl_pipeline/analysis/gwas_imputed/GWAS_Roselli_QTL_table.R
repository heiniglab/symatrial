# Table for all GWAS SNPs
# GWAS SNPs: snps with pvalue < 5e-8 in the Roselli data

setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts/analysis/gwas_imputed")
source("../../helper/helper.R")
local=F

get.duplicates <- function(df, cols){
  df2 <- data.frame(table(df[, cols]))
  df2 <- df2[df2$Freq>1, ]
  colnames(df2) <- c(cols, "Freq")
  df <- merge(df,
              df2[cols],
              all.y=T)
}


# create files
if(F){
  if(!file.exists(paste0(get.path("gwas2", local), "Roselli2018_AF_HRC_GWAS_ALLv11.RDS"))){
    GWAS_file <- "/storage/groups/epigenereg01/workspace/public_data/roselli_2018_AF_gwas/AF_HRC_GWAS_ALLv11/AF_HRC_GWAS_ALLv11.txt"
    AF_GWAS <- read.table(GWAS_file,
                          h=T,
                          stringsAsFactors=F)
    saveRDS(AF_GWAS,
            file=paste0(get.path("gwas2", local),
                        "Roselli2018_AF_HRC_GWAS_ALLv11.RDS"))
  }else{
    AF_GWAS <- readRDS(file=paste0(get.path("gwas2", local),
                                   "Roselli2018_AF_HRC_GWAS_ALLv11.RDS"))
  }
  
  map_file <- "/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/AFHRI-B_imputed.bim"
  map = read.table(map_file,
                   stringsAsFactors=F)
  colnames(map) <- c("chr", "snps", "cM", "pos", "A1", "A2")
  map$cM <- NULL
  
  common_snps <- merge(map,
                       AF_GWAS)
  
  QTL <- readRDS(file=paste0(get.path("results", local),
                             "imputed/cis/QTL_res_all_snp_group.RDS"))
  colnames(QTL)
  
  QTLsnps <- unique(QTL$snps)
  common_snps <- common_snps[common_snps$snps %in% QTLsnps, ]
  
  common_snps$allele.match <- NA
  common_snps[which(common_snps$Allele1==common_snps$A1 & common_snps$Allele2==common_snps$A1), "allele.match"] <- -1
  common_snps[which(common_snps$Allele1==common_snps$A2 & common_snps$Allele2==common_snps$A1), "allele.match"] <- 1
  
  saveRDS(common_snps,
          file=paste0(get.path("gwas2", local),
                      "common_snps.RDS"))
  
  table(is.na(common_snps$allele.match))
  table(common_snps$allele.match)
  
  test <- table(common_snps$MarkerName)
  table(test)
  
  common_snps2 <- common_snps[!is.na(common_snps$allele.match), ]
  common_snps2 <- common_snps2[!duplicated(common_snps2), ]
  table(table(common_snps2$MarkerName))
  
  saveRDS(common_snps2,
          file=paste0(get.path("gwas2", local),
                      "common_snps_aligned.RDS"))
  
  AF_snps <- unique(common_snps2$MarkerName)
  QTL_snps <- unique(common_snps2$snps)
  
  QTL.p <- QTL[QTL$snps %in% QTL_snps, ]
  saveRDS(QTL.p, file=paste0(get.path("gwas2", local),
                             "QTL_AF_GWAS_overlap.RDS"))
}
  
common_snps2 <- readRDS(paste0(get.path("gwas2", local),
                               "common_snps_aligned.RDS"))
table(common_snps2$allele.match)

# filter for significance
common_snps <- common_snps2[common_snps2$P.value < 5e-8, ]
colnames(common_snps)
QTL.snp <- readRDS(paste0(get.path("results", local),
                          "imputed/cis/QTL_res_all_snp_group_anno.RDS"))
QTL.snp2 <- QTL.snp[QTL.snp$snps %in% common_snps$snps, ]
colnames(QTL.snp2)

GWAS <- merge(common_snps, QTL.snp2,
              by.x=c("snps", "chr", "A1", "A2"),
              by.y=c("snps", "chr", "Allele1", "Allele2"),
              all.y=T)
colnames(GWAS)
# [1] "snps"          "chr"           "pos"           "A1"            "A2"            "MarkerName"    "Allele1"       "Allele2"      
# [9] "Effect"        "StdErr"        "P.value"       "allele.match"  "rsid"          "gene"          "snp.pos"       "gene.start"   
# [17] "gene.end"      "strand"        "TSS"           "pval.eqtl"     "FDR.eqtl"      "beta.eqtl"     "pval.pqtl"     "FDR.pqtl"     
# [25] "beta.pqtl"     "pval.pqtl_res" "FDR.pqtl_res"  "beta.pqtl_res" "pval.eqtl_res" "FDR.eqtl_res"  "beta.eqtl_res" "eQTL"         
# [33] "pQTL"          "presQTL"       "eresQTL"       "Group1"        "Group2"        "Group2b"       "Group3"        "Group4"       
# [41] "Group5"        "Group6"        "Group7"        "Group8"        "MAF"           "HWE"           "nTSS"          "dist"         
# [49] "in.gene"       "in.transcript" "UTR3"          "UTR5"          "exon"          "splice"        "RBPBS"         "TFBS"         
# [57] "miRBS"         "state"         "VEP"           "hi.ex"         "hi.ex2"  

saveRDS(GWAS,
        file=paste0(get.path("gwas2", local), "GWAS_QTL_overlap_GWAS_significant.RDS"))

write.table(GWAS,
            row.names = F, col.names = T, quote=F, sep = "\t",
            file=paste0(get.path("gwas2", local), "GWAS_QTL_overlap_GWAS_significant.tsv"))


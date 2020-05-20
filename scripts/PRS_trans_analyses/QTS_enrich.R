setwd("~/work/symAtrial_QTL/scripts/polygenic_risk_scores")
source("../helper/helper.R")
local=T
get.duplicates <- function(df, cols){
  df2 <- data.frame(table(df[, cols]))
  df2 <- df2[df2$Freq>1, ]
  colnames(df2) <- c(cols, "Freq")
  df <- merge(df,
              df2[cols],
              all.y=T)
}

library(pROC)   # make ROC curves
library(ROCR)   # compute AUC
library(pscl)   # compute pseudo R^2
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(caret)
library(plyr)
library(data.table)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

HMGU.blue <- "#003E6E"
mygray <- "#C6DDEA"
col.paired <- brewer.pal(n = 11, "Paired")
col.set <- col.paired[c(2,8,9)]


# transcriptomics eQTS ---------------------------------------------------------
res.GPS.trans.file <- paste0(get.path("dataset", local),
                             "risk_score/correlations/",
                             "assoc_transcriptomics_GPS_covs_RIN.RDS")
if(!file.exists(res.GPS.trans.file)){
  texpr <- readRDS(paste0(get.path("dataset", local),
                          "AFHRI_B_transcriptomics_QC_symbol.RDS"))
  tpheno <- readRDS(paste0(get.path("dataset", local),
                           "AFHRI_B_transcriptomics_phenos_symbol.RDS"))
  tpheno$externID2 <- sub(".RAA", "", tpheno$externID)
  covs <- c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "RIN")
  tpheno2 <- merge(tpheno[, c("externID", "externID2", covs)],
                   scores[, c("externID2", "AF.GPS", "percentile.AF")],
                   all.x = T)
  tpheno2 <- tpheno2[!duplicated(tpheno2), ]
  tpheno3 <- tpheno2[complete.cases(tpheno2), ]
  tpheno3$expr <- NA
  tpheno3$res <- NA
  df <- tpheno3
  
  res.GPS.trans <- data.frame(matrix(nrow = dim(texpr)[1], ncol = 13),
                              stringsAsFactors = F,
                              row.names = rownames(texpr))
  colnames(res.GPS.trans) <- c("id",
                               "Estimate", "StdError", "tvalue", "pvalueT", "padjT",
                               "ResDf", "RSS", "Df", "SumOfSq", "F", "pvalueF", "padjF")
  
  for (i in 1:(dim(texpr)[1])){
    #i=1
    id <- rownames(texpr)[i]
    df$expr <- as.numeric(texpr[id, df$externID])
    df$res <- residuals.lm(lm(expr ~ age + sex + BMI + sysBP + CRP + NTproBNP + RIN,
                              data = df))
    lmres <- lm(res ~ percentile.AF,
                data = df)
    res.GPS.trans[id, ] <- c(id,
                             summary(lmres)$coefficients["percentile.AF", ],
                             NA,
                             NA, NA, NA, NA, NA, NA,
                             NA)
    if(i %% 1000 == 0) print(paste0("Gene ", i, " / ", dim(texpr)[1], " done"))
  }
  res.GPS.trans$Estimate <- as.numeric(res.GPS.trans$Estimate)
  res.GPS.trans$StdError <- as.numeric(res.GPS.trans$StdError)
  res.GPS.trans$tvalue <- as.numeric(res.GPS.trans$tvalue)
  res.GPS.trans$pvalueT <- as.numeric(res.GPS.trans$pvalueT)
  res.GPS.trans$padjT <- p.adjust(res.GPS.trans$pvalueT, method = "BH")
  res.GPS.trans$padjF <- p.adjust(res.GPS.trans$pvalueF, method = "BH")
  saveRDS(res.GPS.trans,
          file = res.GPS.trans.file)
} else {
  res.GPS.trans <- readRDS(res.GPS.trans.file)
}

expr <- readRDS(file=paste0(get.path("dataset", local),
                            "AFHRI_B_transcriptomics_QC_symbol.RDS"))
#expr[1:10, 1:5]
genes <- rownames(expr)
median.expr <- data.frame(AFHRIB=apply(expr, 1, median),
                          Description=genes,
                          stringsAsFactors = F,
                          row.names = genes)

# select for expressed genes (GTEx atrial appendage)
tpm <- read.csv(file=paste0(get.path("gtex", local),
                            "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"),
                sep = "\t", stringsAsFactors = F, skip = 2)
tpm.aa <- tpm[, c("gene_id", "Description", "Heart...Atrial.Appendage")]
library(data.table)
tpm.aa <- as.data.table(tpm.aa)
tpm.aa <- tpm.aa[, meanExpr := mean(Heart...Atrial.Appendage), by = Description]
tpm.aa <- as.data.frame(merge(tpm.aa,
                              median.expr,
                              all = T))
cutoff <- 1
tpm.aa$highly.expressed <- as.numeric(log(tpm.aa$meanExpr+1) > cutoff)
tpm.aa2 <- tpm.aa[!duplicated(tpm.aa$Description), c("Description", "meanExpr", "AFHRIB", "highly.expressed")]
tpm.aa2 <- tpm.aa2[!is.na(tpm.aa2$AFHRIB), ]

res.GPS.trans$Description <- rownames(res.GPS.trans)
res.GPS.trans2 <- merge(res.GPS.trans,
                        tpm.aa2[, c("Description", "highly.expressed", "meanExpr")],
                        all.x = T)
res.GPS.trans2$padjT.highly.expressed <- NA
res.GPS.trans2$padjT.highly.expressed[which(res.GPS.trans2$highly.expressed==1)] <- p.adjust(res.GPS.trans2$pvalueT[which(res.GPS.trans2$highly.expressed==1)], method = "BH")
rownames(res.GPS.trans2) <- res.GPS.trans2$id
top.table.trans <- rbind(res.GPS.trans2[order(res.GPS.trans2$pvalueT + as.numeric(res.GPS.trans2$Estimate < 0))[1:10], ],
                         NA,
                         res.GPS.trans2[order(res.GPS.trans2$pvalueT + as.numeric(res.GPS.trans2$Estimate > 0))[1:10], ])
signif(top.table.trans[, c("Estimate", "tvalue", "pvalueT", "padjT", "padjT.highly.expressed")],
       digits = 3)

top.table.trans <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), ]
top.table.trans <- top.table.trans[order(top.table.trans$pvalueT)[1:10], ]
signif(top.table.trans[order(top.table.trans$tvalue, decreasing=T),
                       c("Estimate", "tvalue", "pvalueT")],
       digits = 3)


# transcriptomics GSEA on eQTS -------------------------------------------------
fgsea.trans.file <- paste0(get.path("dataset", local),
                           "risk_score/correlations/",
                           "fgsea_assoc_transcriptomics_GPS_covs_RIN.RDS")
if(!file.exists(fgsea.trans.file)){
  fgsea.trans <- NULL
  library(fgsea)
  
  ## GO biological process
  gmt <- gmtPathways("../../../symAtrial_multiOMICs/data/current/tables/pathways/c5.bp.v6.1.symbols.gmt.txt")
  
  rank <- res.GPS.trans2[, "tvalue"]
  names(rank) <- res.GPS.trans2[, "id"]
  fgsea.trans[["GObp_all"]] <- fgsea(gmt,
                                     rank,
                                     nperm=100000,
                                     minSize = 15,
                                     maxSize=500)
  fgsea.trans[["GObp_all"]] <- fgsea.trans[["GObp_all"]][order(fgsea.trans[["GObp_all"]]$pval), ]
  
  rank <- res.GPS.trans2[, "pvalueT"]
  names(rank) <- res.GPS.trans2[, "id"]
  fgsea.trans[["GObp_all_p"]] <- fgsea(gmt,
                                       rank,
                                       nperm=100000,
                                       minSize = 15,
                                       maxSize=500)
  fgsea.trans[["GObp_all_p"]] <- fgsea.trans[["GObp_all_p"]][order(fgsea.trans[["GObp_all_p"]]$pval), ]
  
  rank <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "tvalue"]
  names(rank) <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "id"]
  fgsea.trans[["GObp_high"]] <- fgsea(gmt,
                                      rank,
                                      nperm=100000,
                                      minSize = 15,
                                      maxSize=500)
  fgsea.trans[["GObp_high"]] <- fgsea.trans[["GObp_high"]][order(fgsea.trans[["GObp_high"]]$pval), ]
  
  rank <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "pvalueT"]
  names(rank) <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "id"]
  fgsea.trans[["GObp_high_p"]] <- fgsea(gmt,
                                        rank,
                                        nperm=100000,
                                        minSize = 15,
                                        maxSize=500)
  fgsea.trans[["GObp_high_p"]] <- fgsea.trans[["GObp_high_p"]][order(fgsea.trans[["GObp_high_p"]]$pval), ]
  
  ## KEGG
  gmt <- gmtPathways("../../../symAtrial_multiOMICs/data/current/tables/pathways/c2.cp.kegg.v6.1.symbols.gmt.txt")
  
  rank <- res.GPS.trans2[, "tvalue"]
  names(rank) <- res.GPS.trans2[, "id"]
  fgsea.trans[["KEGG_all"]] <- fgsea(gmt,
                                     rank,
                                     nperm=100000,
                                     minSize = 15,
                                     maxSize=500)
  fgsea.trans[["KEGG_all"]] <- fgsea.trans[["KEGG_all"]][order(fgsea.trans[["KEGG_all"]]$pval), ]
  
  rank <- res.GPS.trans2[, "pvalueT"]
  names(rank) <- res.GPS.trans2[, "id"]
  fgsea.trans[["KEGG_all_p"]] <- fgsea(gmt,
                                       rank,
                                       nperm=100000,
                                       minSize = 15,
                                       maxSize=500)
  fgsea.trans[["KEGG_all_p"]] <- fgsea.trans[["KEGG_all_p"]][order(fgsea.trans[["KEGG_all_p"]]$pval), ]
  
  rank <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "tvalue"]
  names(rank) <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "id"]
  fgsea.trans[["KEGG_high"]] <- fgsea(gmt,
                                      rank,
                                      nperm=100000,
                                      minSize = 15,
                                      maxSize=500)
  fgsea.trans[["KEGG_high"]] <- fgsea.trans[["KEGG_high"]][order(fgsea.trans[["KEGG_high"]]$pval), ]
  
  rank <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "pvalueT"]
  names(rank) <- res.GPS.trans2[which(res.GPS.trans2$highly.expressed==1), "id"]
  fgsea.trans[["KEGG_high_p"]] <- fgsea(gmt,
                                        rank,
                                        nperm=100000,
                                        minSize = 15,
                                        maxSize=500)
  fgsea.trans[["KEGG_high_p"]] <- fgsea.trans[["KEGG_high_p"]][order(fgsea.trans[["KEGG_high_p"]]$pval), ]
  
  saveRDS(fgsea.trans, file = fgsea.trans.file)
} else {
  fgsea.trans <- readRDS(fgsea.trans.file)
}




# Proteomics pQTS --------------------------------------------------------------
res.GPS.prot.file <- paste0(get.path("dataset", local),
                            "risk_score/correlations/",
                            "assoc_proteomics_GPS_covs_conc.RDS")
if(!file.exists(res.GPS.prot.file)){
  pexpr <- readRDS(paste0(get.path("dataset", local),
                          "AFHRI_B_proteomics_QC_symbol.RDS"))
  ppheno <- readRDS(paste0(get.path("dataset", local),
                           "AFHRI_B_proteomics_phenos_symbol.RDS"))
  ppheno$externID2 <- sub(".RAA", "", ppheno$externID)
  covs <- c("age", "sex", "BMI", "sysBP", "CRP", "NTproBNP", "Protein.c..ug.ul.")
  ppheno2 <- merge(ppheno[, c("externID", "externID2", covs)],
                   scores[, c("externID2", "AF.GPS", "percentile.AF")],
                   all.x = T)
  ppheno2 <- ppheno2[!duplicated(ppheno2), ]
  ppheno3 <- ppheno2[complete.cases(ppheno2), ]
  ppheno3$expr <- NA
  df <- ppheno3
  
  res.GPS.prot <- data.frame(matrix(nrow = dim(pexpr)[1], ncol = 13),
                             stringsAsFactors = F,
                             row.names = rownames(pexpr))
  colnames(res.GPS.prot) <- c("id",
                              "Estimate", "StdError", "tvalue", "pvalueT", "padjT",
                              "ResDf", "RSS", "Df", "SumOfSq", "F", "pvalueF", "padjF")
  for (i in 1:(dim(pexpr)[1])){
    #i=1
    id <- rownames(pexpr)[i]
    df$expr <- as.numeric(pexpr[id, df$externID])
    lmres <- lm(expr ~ AF.GPS + age + sex + BMI + sysBP + CRP + NTproBNP + Protein.c..ug.ul.,
                data = df)
    lmres.null <- lm(expr ~ age + sex + BMI + sysBP + CRP + NTproBNP + Protein.c..ug.ul.,
                     data = df) 
    anova.res <- anova(lmres.null, lmres)
    res.GPS.prot[id, ] <- c(id,
                            summary(lmres)$coefficients["AF.GPS", ],
                            NA,
                            as.data.frame(anova.res[2, ]),
                            NA)
    if(i %% 100 == 0) print(paste0("Gene ", i, " / ", dim(pexpr)[1], " done"))
  }
  res.GPS.prot$padjT <- p.adjust(res.GPS.prot$pvalueT, method = "BH")
  res.GPS.prot$padjF <- p.adjust(res.GPS.prot$pvalueF, method = "BH")
  saveRDS(res.GPS.prot,
          file = res.GPS.prot.file)
} else {
  res.GPS.prot <- readRDS(res.GPS.prot.file)
}
top.table.prot <- rbind(res.GPS.prot[order(res.GPS.prot$pvalueF +
                                             as.numeric(res.GPS.prot$Estimate < 0))[1:10], ],
                        res.GPS.prot[order(res.GPS.prot$pvalueF +
                                             as.numeric(res.GPS.prot$Estimate > 0))[1:10], ])

res.GPS.prot$Description <- rownames(res.GPS.prot)
res.GPS.prot2 <- merge(res.GPS.prot,
                       tpm.aa2[, c("Description", "highly.expressed", "meanExpr")],
                       all.x = T)
res.GPS.prot2$padjT.highly.expressed <- NA
res.GPS.prot2$padjT.highly.expressed[which(res.GPS.prot2$highly.expressed==1)] <- p.adjust(res.GPS.prot2$pvalueT[which(res.GPS.prot2$highly.expressed==1)], method = "BH")
rownames(res.GPS.prot2) <- res.GPS.prot2$id

top.table.prot <- rbind(res.GPS.prot2[order(res.GPS.prot2$pvalueT + as.numeric(res.GPS.prot2$Estimate < 0))[1:10], ],
                        NA,
                        res.GPS.prot2[order(res.GPS.prot2$pvalueT + as.numeric(res.GPS.prot2$Estimate > 0))[1:10], ])
signif(top.table.prot[, c("Estimate", "tvalue", "pvalueT", "padjT", "padjT.highly.expressed")],
       digits = 3)



# Proteomics GSEA pQTS ---------------------------------------------------------
fgsea.prot.file <- paste0(get.path("dataset", local),
                          "risk_score/correlations/",
                          "fgsea_assoc_proteomics_GPS_covs_conc.RDS")
if(!file.exists(fgsea.prot.file)){
  fgsea.prot <- NULL
  library(fgsea)
  
  ## GO biological process
  gmt <- gmtPathways("../../../symAtrial_multiOMICs/data/current/tables/pathways/c5.bp.v6.1.symbols.gmt.txt")
  
  rank <- res.GPS.prot2[, "tvalue"]
  names(rank) <- res.GPS.prot2[, "id"]
  fgsea.prot[["GObp_all"]] <- fgsea(gmt,
                                    rank,
                                    nperm=100000,
                                    minSize = 5,
                                    maxSize=500)
  fgsea.prot[["GObp_all"]] <- fgsea.prot[["GObp_all"]][order(fgsea.prot[["GObp_all"]]$pval), ]
  
  rank <- res.GPS.prot2[, "pvalueT"]
  names(rank) <- res.GPS.prot2[, "id"]
  fgsea.prot[["GObp_all_p"]] <- fgsea(gmt,
                                      rank,
                                      nperm=100000,
                                      minSize = 5,
                                      maxSize=500)
  fgsea.prot[["GObp_all_p"]] <- fgsea.prot[["GObp_all_p"]][order(fgsea.prot[["GObp_all_p"]]$pval), ]
  
  ## KEGG
  gmt <- gmtPathways("../../../symAtrial_multiOMICs/data/current/tables/pathways/c2.cp.kegg.v6.1.symbols.gmt.txt")
  
  rank <- res.GPS.prot2[, "tvalue"]
  names(rank) <- res.GPS.prot2[, "id"]
  fgsea.prot[["KEGG_all"]] <- fgsea(gmt,
                                    rank,
                                    nperm=100000,
                                    minSize = 5,
                                    maxSize=500)
  fgsea.prot[["KEGG_all"]] <- fgsea.prot[["KEGG_all"]][order(fgsea.prot[["KEGG_all"]]$pval), ]
  
  rank <- res.GPS.prot2[, "pvalueT"]
  names(rank) <- res.GPS.prot2[, "id"]
  fgsea.prot[["KEGG_all_p"]] <- fgsea(gmt,
                                      rank,
                                      nperm=100000,
                                      minSize = 5,
                                      maxSize=500)
  fgsea.prot[["KEGG_all_p"]] <- fgsea.prot[["KEGG_all_p"]][order(fgsea.prot[["KEGG_all_p"]]$pval), ]
  
  saveRDS(fgsea.prot, file = fgsea.prot.file)
} else {
  fgsea.prot <- readRDS(fgsea.prot.file)
}

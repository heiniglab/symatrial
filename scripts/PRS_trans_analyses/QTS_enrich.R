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

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(fgsea)

HMGU.blue <- "#003E6E"
mygray <- "#C6DDEA"
col.paired <- brewer.pal(n = 11, "Paired")
col.set <- col.paired[c(2,8,9)]

gmt <- gmtPathways("../../../symAtrial_multiOMICs/data/current/tables/pathways/c5.bp.v6.1.symbols.gmt.txt")
set.seed(1111)

scores <- readRDS(paste0(get.path("dataset", local),
                         "risk_score/GPS_scores_EUR_AFHRI-B.RDS"))
scores$AF_Status <- factor(scores$AF_Status, levels = c(0,1,2), labels = c("ctrl", "preOP AF", "only postOP AF"))
scores$fibro.score.imp.prot.meta <- NULL
scores <- scores[!duplicated(scores), ]
rownames(scores) <- scores$externID2

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



# Transcriptomics eQTS ---------------------------------------------------------
## eQTS lm percentile.AF age/sex/BMI/sysBP/CRP/NTproBNP/RIN --------------------

eQTS.file <- paste0(get.path("results", local),
                    "imputed/trans/eQTS/",
                    "eQTS_percAF_covs_RIN.RDS")
if(!file.exists(eQTS.file)){
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
  
  res.eQTS <- data.frame(matrix(nrow = dim(texpr)[1], ncol = 8),
                         stringsAsFactors = F,
                         row.names = rownames(texpr))
  colnames(res.eQTS) <- c("id", "N", "df",
                          "Estimate", "StdError", "tvalue", "pvalueT", "padjT")
  
  for (i in 1:(dim(texpr)[1])){
    #i=1
    id <- rownames(texpr)[i]
    df$expr <- as.numeric(texpr[id, df$externID])
    lmres <- lm(expr ~ percentile.AF + age + sex + BMI + sysBP + CRP + NTproBNP + RIN,
                data = df)
    res.eQTS[id, ] <- c(id, dim(lmres$model)[1], lmres$df.residual,
                        summary(lmres)$coefficients["percentile.AF", ],
                        NA)
    if(i %% 1000 == 0) print(paste0("Gene ", i, " / ", dim(texpr)[1], " done"))
  }
  res.eQTS$N <- as.numeric(res.eQTS$N)
  res.eQTS$df <- as.numeric(res.eQTS$df)
  res.eQTS$Estimate <- as.numeric(res.eQTS$Estimate)
  res.eQTS$StdError <- as.numeric(res.eQTS$StdError)
  res.eQTS$tvalue <- as.numeric(res.eQTS$tvalue)
  res.eQTS$pvalueT <- as.numeric(res.eQTS$pvalueT)
  res.eQTS$padjT <- p.adjust(res.eQTS$pvalueT, method = "BH")
  
  res.eQTS$Description <- rownames(res.eQTS)
  res.eQTS2 <- merge(res.eQTS,
                     tpm.aa2[, c("Description", "highly.expressed", "meanExpr")],
                     all.x = T)
  res.eQTS2$padjT.highly.expressed <- NA
  res.eQTS2$padjT.highly.expressed[which(res.eQTS2$highly.expressed==1)] <-
    p.adjust(res.eQTS2$pvalueT[which(res.eQTS2$highly.expressed==1)],
             method = "BH")
  rownames(res.eQTS2) <- res.eQTS2$id
  saveRDS(res.eQTS2,
          file = eQTS.file)
} else {
  res.eQTS2 <- readRDS(eQTS.file)
}

top.table.trans <- rbind(res.eQTS2[order(res.eQTS2$pvalueT +
                                           as.numeric(res.eQTS2$Estimate < 0))[1:10], ],
                         NA,
                         res.eQTS2[order(res.eQTS2$pvalueT + 
                                           as.numeric(res.eQTS2$Estimate > 0))[1:10], ])
signif(top.table.trans[, c("Estimate", "tvalue", "pvalueT", "padjT", "padjT.highly.expressed")],
       digits = 3)

top.table.trans <- res.eQTS2[which(res.eQTS2$highly.expressed==1), ]
top.table.trans <- top.table.trans[order(top.table.trans$pvalueT)[1:10], ]
signif(top.table.trans[order(top.table.trans$tvalue, decreasing=T),
                       c("Estimate", "tvalue", "pvalueT")],
       digits = 3)


# GSEA on eQTS -----------------------------------------------------------------
## GSEA on eQTS lm percentile.AF age/sex/BMI/sysBP/CRP/NTproBNP/RIN ------------

fgsea.eQTS.file <- paste0(get.path("results", local),
                          "imputed/trans/eQTS/",
                          "fgsea_eQTS_percAF_covs_RIN.RDS")
if(!file.exists(fgsea.eQTS.file)){
  
  fgsea.eQTS <- NULL
  eQTS.file <- paste0(get.path("results", local),
                      "imputed/trans/eQTS/",
                      "eQTS_percAF_covs_RIN.RDS")
  
  eQTS <- readRDS(eQTS.file)
  
  rank <- eQTS[which(eQTS$highly.expressed==1), "tvalue"]
  names(rank) <- eQTS[which(eQTS$highly.expressed==1), "id"]
  fgsea.eQTS[["GObp_high"]] <- fgsea(gmt,
                                     rank,
                                     nperm=100000,
                                     minSize = 15,
                                     maxSize=500)
  fgsea.eQTS[["GObp_high"]] <- fgsea.eQTS[["GObp_high"]][order(fgsea.eQTS[["GObp_high"]]$pval), ]
  
  fgsea.eQTS[["significant"]] <- fgsea.eQTS[["GObp_high"]][fgsea.eQTS[["GObp_high"]]$padj<0.05, ]
  fgsea.eQTS[["leadingEdge"]] <- data.frame(table(unlist(fgsea.eQTS[["significant"]]$leadingEdge)))
  
  saveRDS(fgsea.eQTS, file = fgsea.eQTS.file)
} else {
  fgsea.eQTS <- readRDS(fgsea.eQTS.file)
}

# Proteomics pQTS --------------------------------------------------------------
## pQTS lm percentile.AF age/sex/BMI/sysBP/CRP/NTproBNP/prot.conc --------------
pQTS.file <- paste0(get.path("results", local),
                    "imputed/trans/pQTS/",
                    "pQTS_percAF_covs_conc.RDS")
if(!file.exists(pQTS.file)){
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
  
  res.pQTS <- data.frame(matrix(nrow = dim(pexpr)[1], ncol = 8),
                         stringsAsFactors = F,
                         row.names = rownames(pexpr))
  colnames(res.pQTS) <- c("id", "N", "df",
                          "Estimate", "StdError", "tvalue", "pvalueT", "padjT")
  
  for (i in 1:(dim(pexpr)[1])){
    #i=1
    id <- rownames(pexpr)[i]
    df$expr <- as.numeric(pexpr[id, df$externID])
    lmres <- lm(expr ~ percentile.AF + age + sex + BMI + sysBP + CRP + NTproBNP + Protein.c..ug.ul.,
                data = df)
    res.pQTS[id, ] <- c(id, dim(lmres$model)[1], lmres$df.residual,
                        summary(lmres)$coefficients["percentile.AF", ],
                        NA)
    if(i %% 100 == 0) print(paste0("Gene ", i, " / ", dim(pexpr)[1], " done"))
  }
  res.pQTS$N <- as.numeric(res.pQTS$N)
  res.pQTS$df <- as.numeric(res.pQTS$df)
  res.pQTS$Estimate <- as.numeric(res.pQTS$Estimate)
  res.pQTS$StdError <- as.numeric(res.pQTS$StdError)
  res.pQTS$tvalue <- as.numeric(res.pQTS$tvalue)
  res.pQTS$pvalueT <- as.numeric(res.pQTS$pvalueT)
  res.pQTS$padjT <- p.adjust(res.pQTS$pvalueT, method = "BH")
  
  res.pQTS$Description <- rownames(res.pQTS)
  res.pQTS2 <- merge(res.pQTS,
                     tpm.aa2[, c("Description", "highly.expressed", "meanExpr")],
                     all.x = T)
  res.pQTS2$padjT.highly.expressed <- NA
  res.pQTS2$padjT.highly.expressed[which(res.pQTS2$highly.expressed==1)] <-
    p.adjust(res.pQTS2$pvalueT[which(res.pQTS2$highly.expressed==1)],
             method = "BH")
  rownames(res.pQTS2) <- res.pQTS2$id
  saveRDS(res.pQTS2,
          file = pQTS.file)
} else {
  res.pQTS2 <- readRDS(pQTS.file)
}

top.table.prot <- rbind(res.pQTS2[order(res.pQTS2$pvalueT + as.numeric(res.pQTS2$Estimate < 0))[1:10], ],
                        NA,
                        res.pQTS2[order(res.pQTS2$pvalueT + as.numeric(res.pQTS2$Estimate > 0))[1:10], ])
signif(top.table.prot[, c("Estimate", "tvalue", "pvalueT", "padjT", "padjT.highly.expressed")],
       digits = 3)


# GSEA on pQTS -----------------------------------------------------------------
## GSEA on pQTS lm percentile.AF age/sex/BMI/sysBP/CRP/NTproBNP/prot.conc ------

fgsea.pQTS.file <- paste0(get.path("results", local),
                          "imputed/trans/pQTS/",
                          "fgsea_pQTS_percAF_covs_conc.RDS")
if(!file.exists(fgsea.pQTS.file)){
  
  fgsea.pQTS <- NULL
  pQTS.file <- paste0(get.path("results", local),
                      "imputed/trans/pQTS/",
                      "pQTS_percAF_covs_conc.RDS")
  
  pQTS <- readRDS(pQTS.file)
  
  rank <- pQTS[which(pQTS$highly.expressed==1), "tvalue"]
  names(rank) <- pQTS[which(pQTS$highly.expressed==1), "id"]
  fgsea.pQTS[["GObp_high"]] <- fgsea(gmt,
                                     rank,
                                     nperm=100000,
                                     minSize = 5,
                                     maxSize=500)
  fgsea.pQTS[["GObp_high"]] <- fgsea.pQTS[["GObp_high"]][order(fgsea.pQTS[["GObp_high"]]$pval), ]
  
  rank <- pQTS[, "tvalue"]
  names(rank) <- pQTS[, "id"]
  fgsea.pQTS[["GObp_all"]] <- fgsea(gmt,
                                    rank,
                                    nperm=100000,
                                    minSize = 5,
                                    maxSize=500)
  fgsea.pQTS[["GObp_all"]] <- fgsea.pQTS[["GObp_all"]][order(fgsea.pQTS[["GObp_all"]]$pval), ]
  
  fgsea.pQTS[["significant"]] <- fgsea.pQTS[["GObp_all"]][fgsea.pQTS[["GObp_all"]]$padj<0.05, ]
  fgsea.pQTS[["leadingEdge"]] <- data.frame(table(unlist(fgsea.pQTS[["significant"]]$leadingEdge)))
  
  saveRDS(fgsea.pQTS, file = fgsea.pQTS.file)
} else {
  fgsea.pQTS <- readRDS(fgsea.pQTS.file)
}


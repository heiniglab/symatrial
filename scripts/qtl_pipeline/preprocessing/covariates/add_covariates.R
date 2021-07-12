source("preprocessing/data_parsing/cut_tables.R")


get.covariates <- function(expression, covariates){
  covariates <- covariates[rownames(covariates) %in% colnames(expression),]
  
  covariates$sex <- as.numeric(levels(covariates$sex))[covariates$sex]
  covariates$overall_AF <- as.numeric(levels(covariates$overall_AF))[covariates$overall_AF]
  
  cov <- t(covariates)
  cov <- order.covariates(expression, cov)
  cov <- t(cov)
  
  return(cov)
}



prepare.covariates <- function(expression, only.fibro=F, batch=F, local=F, pop=F, only.pop=F){
  # expression <- imputed_m_MUC
  # only.fibro = T
  # batch=T
  data_path <- get.path("dataset", local)
  
  covariates <- readRDS(paste0(data_path, "AFHRI_B_phenotypes.RDS"))
  covariates <- covariates[-which(covariates$externID=="SYM0012406"),]
  covariates$externID <- sub(" RAA", "", covariates$externID)
  covariates <- covariates[-grep(" LAA", covariates$externID), ]
  
  covariates <- na.omit(covariates)
  covariates <- covariates[!duplicated(covariates), ]
  #covariates[duplicated(covariates$externID), ]
  
  rownames(covariates) <- covariates$externID
  if(batch){
    tissue <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_tissue.RDS"))
    covariates <- merge(covariates,
                        tissue[, c("id", "batch")],
                        by.x="externID", by.y="id",
                        all.x=T, sort=F)
    covariates$batch <- as.numeric(as.factor(covariates$batch))-1
  }
  
  rownames(covariates) <- covariates$externID
  
  
  
  # covariates$Raucher_aktuell <- NULL
  # covariates$Protein.c..ug.ul.<- NULL
  # covariates$RAA.or.LAA <- NULL
  # covariates$RIN <- NULL
  
  covariates <- covariates[rownames(covariates) %in% colnames(expression),]
  covariates$sex <- as.numeric(levels(covariates$sex))[covariates$sex]
  covariates$overall_AF <- as.numeric(levels(covariates$overall_AF))[covariates$overall_AF]
  
  if(only.fibro){
    covariates <- covariates[, colnames(covariates) %in% c("fibro.score", "batch", "externID")]
    rownames(covariates) <- covariates$externID
  }
  
  if(pop){
    pca.eur <- read.table(paste0(get.path("genotype", local),
                                 "AFHRI-B_EUR_population_pca.eigenvec"),
                          h = F)[,-1]
    colnames(pca.eur) <- c("sample", paste0("EUR_PC", 1:20))
    rownames(pca.eur) <- pca.eur$sample
    covariates[, paste0("EUR_PC", 1:pop)] <- pca.eur[covariates$externID, paste0("EUR_PC", 1:pop)]
  }
  
  if(only.pop){
    covariates[, c("age", "sex", "BMI", "overall_AF", "fibro.score")] <- NULL
  }
  
  covariates[is.na(covariates)] <- 0
  covariates[, "externID"] <- NULL
  
  cov <- t(as.matrix(covariates))
  
  cov <- order.covariates(expression, cov)
  
  cov <- t(cov)
  
  return(cov)
}

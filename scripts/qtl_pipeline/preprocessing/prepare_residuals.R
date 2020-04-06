# calculate residuals

setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")
# setwd("/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/scripts")

source("helper/helper.R")
local=T
data_path <- get.path("dataset", local)
# data_path <- "/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/data/current/AFHRI_B/"
expression_p <- readRDS(paste0(data_path, "AFHRI_B_proteomics_QC_symbol.RDS"))
colnames(expression_p) <- sub(".RAA","", colnames(expression_p))

expression_t <- readRDS(paste0(data_path, "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expression_t) <- sub(".RAA", "", colnames(expression_t))

# residuals without covariates -------------------------------------------------
match_s <- intersect(colnames(expression_t), colnames(expression_p))
match_g <- intersect(rownames(expression_t), rownames(expression_p))
expression_t2 <- expression_t[match_g, match_s]
expression_p2 <- expression_p[match_g, match_s]

if(!file.exists(paste0(data_path, "residuals_trans.RDS"))){
  residuals_trans <- expression_t2
  for(i in 1:length(rownames(residuals_trans))){
    #i=22
    trans <- t(data.frame(expression_t2[i, ]))
    prot <- t(data.frame(expression_p2[i, ]))
    residuals_trans[i, ] <- resid(lm(trans ~ prot, na.action=na.exclude))
  }
  saveRDS(residuals_trans, file=paste0(data_path, "residuals_trans.RDS"))
}else{
  residuals_trans <- readRDS(paste0(data_path, "residuals_trans.RDS"))
}
  
if(!file.exists(paste0(data_path, "residuals_prot.RDS"))){
  residuals_prot <- expression_p2
  for(i in 1:length(rownames(residuals_prot))){
    #i=22
    trans <- t(data.frame(expression_t2[i, ]))
    prot <- t(data.frame(expression_p2[i, ]))
    residuals_prot[i, ] <- resid(lm(prot ~ trans, na.action=na.exclude))
  }
  saveRDS(residuals_prot, file=paste0(data_path, "residuals_prot.RDS"))
}else{
  residuals_prot <- readRDS(paste0(data_path, "residuals_prot.RDS"))
}


# residuals with covariates: ---------------------------------------------------
# c("age", "sex", "BMI", "overall_AF", "fibro_score", "RIN", "Protein.c..ug.ul.")
pheno_t <- readRDS(paste0(data_path, "AFHRI_B_transcriptomics_phenos_symbol.RDS"))
rownames(pheno_t) <- sub(".RAA", "", pheno_t$externID)

pheno_p <- readRDS(paste0(data_path, "AFHRI_B_proteomics_phenos_symbol.RDS"))
rownames(pheno_p) <- sub(".RAA","", pheno_p$externID)

identical(rownames(pheno_t), colnames(expression_t))
identical(rownames(pheno_p), colnames(expression_p))

match_s <- intersect(colnames(expression_t), colnames(expression_p))
match_g <- intersect(rownames(expression_t), rownames(expression_p))
expression_t2 <- expression_t[match_g, match_s]
expression_p2 <- expression_p[match_g, match_s]

phenos <- pheno_t[match_s, c("age", "sex", "BMI", "overall_AF", "fibro.score", "RIN")]
phenos$prot.conc <- pheno_p[match_s, "Protein.c..ug.ul."]

if(!file.exists(paste0(data_path, "residuals_trans_cov.RDS"))){
  residuals_trans <- expression_t2
  for(i in 1:length(rownames(residuals_trans))){
    #i=22
    trans <- t(data.frame(expression_t2[i, ]))
    prot <- t(data.frame(expression_p2[i, ]))
    residuals_trans[i, ] <- resid(lm(data = phenos,
                                     trans ~ prot + age + sex + BMI + overall_AF + fibro.score + prot.conc + RIN,
                                     na.action=na.exclude))
  }
  residuals_trans <- residuals_trans[, colSums(is.na(residuals_trans))<10]
  saveRDS(residuals_trans, file=paste0(data_path, "residuals_trans_cov.RDS"))
}else{
  residuals_trans <- readRDS(paste0(data_path, "residuals_trans_cov.RDS"))
}

if(!file.exists(paste0(data_path, "residuals_prot_cov.RDS"))){
  residuals_prot <- expression_p2
  for(i in 1:length(rownames(residuals_prot))){
    #i=22
    trans <- t(data.frame(expression_t2[i, ]))
    prot <- t(data.frame(expression_p2[i, ]))
    residuals_prot[i, ] <- resid(lm(data = phenos,
                                    prot ~ trans + age + sex + BMI + overall_AF + fibro.score + prot.conc + RIN,
                                    na.action=na.exclude))
  }
  residuals_prot <- residuals_prot[, colSums(is.na(residuals_prot))<10]
  saveRDS(residuals_prot, file=paste0(data_path, "residuals_prot_cov.RDS"))
}else{
  residuals_prot <- readRDS(paste0(data_path, "residuals_prot_cov.RDS"))
}


# adjust locations
path2snplocs <- get.path("snplocs", local=F)

path2genelocs_t <- paste0(get.path("locations", local = F), "transcript_locations_ensembl.txt")
path2genelocs_res <- paste0(get.path("locations", local = F), "residuals_locations_ensembl.txt")

#locs_p <- read.table(path2genelocs_p, h=T, stringsAsFactors = F)
locs_t <- read.table(path2genelocs_t, h=T, stringsAsFactors = F)

identical(rownames(residuals_trans), rownames(residuals_prot))
locs_part <- locs_t[which(locs_t$name %in% rownames(residuals_trans)), ]
identical(intersect(rownames(residuals_trans), locs_part$name), as.character(locs_part$name))

write.table(locs_part, file=path2genelocs_res, quote = F,
            sep = "\t", row.names = F, col.names = T)


# expression corrected for peer factors: -------------------------------------------------
# transcript: 12, protein: 10
peer <- paste0(get.path("peer", local), "normalized/")
peer_t <- read.table(paste0(peer,
                            "eqtl/factors_no_cov/peer_factors_nk12.txt"),
                     h = T)
colnames(peer_t) <- colnames(read.table(paste0(peer,
                                               "eqtl/imputed_expression_t.txt"),
                                        h=T, nrows = 1))

peer_p <- read.table(paste0(peer,
                            "pqtl/factors_no_cov/peer_factors_nk10.txt"),
                     h = T)
colnames(peer_p) <- colnames(read.table(paste0(peer,
                                               "pqtl/imputed_expression_p.txt"),
                                        h=T, nrows = 1))

match_s <- intersect(colnames(peer_t), colnames(peer_p))
match_g <- intersect(rownames(expression_t), rownames(expression_p))
expression_t2 <- expression_t[match_g, match_s]
expression_p2 <- expression_p[match_g, match_s]
peer_t2 <- t(peer_t[, match_s])
peer_p2 <- t(peer_p[, match_s])

if(!file.exists(paste0(data_path, "residuals_trans_peer.RDS"))){
  residuals_trans <- expression_t2
  for(i in 1:length(rownames(residuals_trans))){
    #i=22
    trans <- t(data.frame(expression_t2[i, ]))
    residuals_trans[i, ] <- resid(lm(trans ~ peer_t2,
                                     na.action=na.exclude))
  }
  residuals_trans <- residuals_trans[, colSums(is.na(residuals_trans))<10]
  saveRDS(residuals_trans, file=paste0(data_path, "residuals_trans_peer.RDS"))
}else{
  residuals_trans <- readRDS(paste0(data_path, "residuals_trans_peer.RDS"))
}

if(!file.exists(paste0(data_path, "residuals_prot_peer.RDS"))){
  residuals_prot <- expression_p2
  for(i in 1:length(rownames(residuals_prot))){
    #i=22
    prot <- t(data.frame(expression_p2[i, ]))
    residuals_prot[i, ] <- resid(lm(prot ~ peer_p2,
                                    na.action=na.exclude))
  }
  residuals_prot <- residuals_prot[, colSums(is.na(residuals_prot))<10]
  saveRDS(residuals_prot, file=paste0(data_path, "residuals_prot_peer.RDS"))
}else{
  residuals_prot <- readRDS(paste0(data_path, "residuals_prot_peer.RDS"))
}




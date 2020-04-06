# prepare GenABEL gwaa object

setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")
# setwd("/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/scripts")

source("helper/helper.R")
source("analysis/peer_tables/peer_result.R")

local=F
data_dir <- get.path("data", local)
geno_path <- get.path("genotype")
data_path <- get.path("dataset")

# create GenABEL object for boxplots
library(GenABEL)
# convert the ped/map files to a GenABEL usable format

# map <- read.table(paste0(geno_path,"AFHRI_B.map"), h=F)
# map$V3 <- NULL
# write.table(map, paste0(geno_path,"AFHRI_B_2.map"),
#             row.names = F, col.names = F, quote = F)

convert.snp.ped(pedfile=paste0(geno_path,"AFHRI_B.ped"),
                mapfile=paste0(geno_path,"AFHRI_B_2.map"),
                outfile=paste0(geno_path,"AFHRI_B.raw"),
                mapHasHeaderLine=F)

# load expression data
fam <- read.table(paste0(geno_path,"AFHRI-B.fam"), h=F)

covariates <- readRDS(paste0(data_path, "AFHRI_B_phenotypes.RDS"))
covariates <- covariates[-which(covariates$externID=="SYM0012406"),]
covariates$externID <- sub(" RAA", "", covariates$externID)
covariates <- covariates[-grep(" LAA", covariates$externID), ]
covariates <- na.omit(covariates)
covariates <- covariates[!duplicated(covariates), ]
rownames(covariates) <- covariates$externID
meta_tissue <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_tissue.RDS"))
colnames(meta_tissue)[-1] <- paste0(colnames(meta_tissue)[-1], "_tissue")
meta_MUC <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_tissue_MUC.RDS"))
colnames(meta_MUC)[-1] <- paste0(colnames(meta_MUC)[-1], "_MUC")
meta_UKE <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_tissue_UKE.RDS"))
colnames(meta_UKE)[-1] <- paste0(colnames(meta_UKE)[-1], "_UKE")
meta_serum <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_serum.RDS"))
colnames(meta_serum)[-1] <- paste0(colnames(meta_serum)[-1], "_serum")

metabolomics <- merge(meta_tissue,
                      meta_MUC,
                      all=T, sort=F)
metabolomics <- merge(metabolomics,
                      meta_UKE,
                      all=T, sort=F)
metabolomics <- merge(metabolomics,
                      meta_serum,
                      all=T, sort=F)
rownames(metabolomics) <- metabolomics$id

covariates <- merge(covariates,
                    data.frame(externID=metabolomics[, c("id")], stringsAsFactors = F),
                    all.x=T, sort=F)
rownames(covariates) <- covariates$externID
#covariates$batch <- as.factor(covariates$batch)
covariates$sex <- NULL

expression_p <- readRDS(paste0(data_path, "AFHRI_B_proteomics_QC_symbol.RDS"))
colnames(expression_p) <- sub(".RAA","", colnames(expression_p))

expression_t <- readRDS(paste0(data_path, "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expression_t) <- sub(".RAA", "", colnames(expression_t))

expression_m <- t(metabolomics[,-(1)])

residuals_t <- readRDS(paste0(data_path, "residuals_trans.RDS"))
colnames(residuals_t) <- sub(".RAA","", colnames(residuals_t))

residuals_p <- readRDS(paste0(data_path, "residuals_prot.RDS"))
colnames(residuals_p) <- sub(".RAA", "", colnames(residuals_p))

rownames(expression_t) <- paste0(rownames(expression_t), "_trans")
rownames(expression_p) <- paste0(rownames(expression_p), "_prot")
rownames(residuals_t) <- paste0(rownames(residuals_t), "_t.res")
rownames(residuals_p) <- paste0(rownames(residuals_p), "_p.res")


# merge all expression data
pheno.temp <- merge(data.frame(id=fam$V1,
                               sex=ifelse(fam$V5==2,0,1),
                               stringsAsFactors = F),
                    covariates,
                    by.x = "id", by.y = "externID", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(expression_t), t(expression_t), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(expression_p), t(expression_p), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(expression_m), t(expression_m), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(residuals_p), t(residuals_p), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(residuals_t), t(residuals_t), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(data.frame(id=fam$V1,
                               stringsAsFactors = F),
                    pheno.temp,
                    by = "id", all.x = T, sort=F)

identical(as.character(fam$V1), as.character(pheno.temp$id))
#test <- data.frame(colnames(pheno.temp))
# saveRDS(pheno.temp, file=paste0(data_path, "expression_all.RDS"))
# write.table(pheno.temp,
#             file = paste0(geno_path, "AFHRI_B.pheno"),
#             row.names = F, col.names = T, quote = F,
#             sep="\t")

# create snp.data-class file including phenotypic information
AFHRI_B <- load.gwaa.data(phenofile=paste0(geno_path, "AFHRI_B.pheno"),
                          genofile=paste0(geno_path, "AFHRI_B.raw"))

save(AFHRI_B, file=paste0(geno_path, "AFHRI_B_gwaa_data.RData"))


#### Now for the imputed Genotype file ----

setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")

source("helper/helper.R")
source("analysis/peer_tables/peer_result.R")

local=F
data_dir <- get.path("data", local)

library(GenABEL)

# unfiltered, imputed genotypes: AFHRI-B_imputed.bed/.bim/.fam
# mapfile AFHRI_B_imputed_2.map for filtering
# AFHRI_B_imputed_recode_alleles.map for recoding alleles to A/B (for import)
# create ped & map filte for GenABEL

# system(paste0("cd /home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype \n",
#               "/home/icb/ines.assum/ext_tools/plink/plink --bfile AFHRI-B_imputed --extract AFHRI_B_imputed_2.map --update-alleles AFHRI_B_imputed_recode_alleles.map  --recode --out AFHRI_B_imputed"))

# create GenABEL .raw file
convert.snp.ped(pedfile=paste0(geno_path,"AFHRI_B_imputed.ped"),
                mapfile=paste0(geno_path,"AFHRI_B_imputed_2.map"),
                outfile=paste0(geno_path,"AFHRI_B_imputed.raw"),
                mapHasHeaderLine=T)

# phenofile can be used from before
phenofile=paste0(get.path("genotype"), "AFHRI_B.pheno")

# pheno <- read.table(paste0(get.path("genotype"), "AFHRI_B.pheno"),
#                     h=T)

AFHRI_B_imp <- load.gwaa.data(phenofile=paste0(get.path("genotype"),
                                               "AFHRI_B.pheno"),
                              genofile=paste0(get.path("genotype"),
                                          "AFHRI_B_imputed.raw"))

save(AFHRI_B_imp,
     file=paste0(get.path("genotype"),
                 "AFHRI_B_imputed_gwaa_data.RData"))

saveRDS(AFHRI_B_imp,
        file=paste0(get.path("genotype"),
                 "AFHRI_B_imputed_gwaa_data.RDS"))

# update pheno info ----------------------------------------
AFHRI_B_imp <- readRDS(paste0(get.path("genotype"),
                              "AFHRI_B_imputed_gwaa_data.RDS"))
pheno.cols <- colnames(phdata(AFHRI_B_imp))

# load expression data
fam <- read.table(paste0(geno_path,"AFHRI-B.fam"), h=F)

covariates <- readRDS(paste0(data_path, "AFHRI_B_phenotypes.RDS"))
covariates <- covariates[-which(covariates$externID=="SYM0012406"),]
covariates$externID <- sub(" RAA", "", covariates$externID)
covariates <- covariates[-grep(" LAA", covariates$externID), ]
covariates <- na.omit(covariates)
covariates <- covariates[!duplicated(covariates), ]
rownames(covariates) <- covariates$externID

meta_tissue <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_tissue.RDS"))
colnames(meta_tissue)[-1] <- paste0(colnames(meta_tissue)[-1], "_tissue")
meta_MUC <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_tissue_MUC.RDS"))
colnames(meta_MUC)[-1] <- paste0(colnames(meta_MUC)[-1], "_MUC")
meta_UKE <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_tissue_UKE.RDS"))
colnames(meta_UKE)[-1] <- paste0(colnames(meta_UKE)[-1], "_UKE")
meta_serum <- readRDS(paste0(get.path("dataset"), "AFHRI_B_metabolomics_serum.RDS"))
colnames(meta_serum)[-1] <- paste0(colnames(meta_serum)[-1], "_serum")

metabolomics <- merge(meta_tissue,
                      meta_MUC,
                      all=T, sort=F)
metabolomics <- merge(metabolomics,
                      meta_UKE,
                      all=T, sort=F)
metabolomics <- merge(metabolomics,
                      meta_serum,
                      all=T, sort=F)
rownames(metabolomics) <- metabolomics$id

covariates <- merge(covariates,
                    data.frame(externID=metabolomics[, c("id")], stringsAsFactors = F),
                    all.x=T, sort=F)
rownames(covariates) <- covariates$externID
#covariates$batch <- as.factor(covariates$batch)
covariates$sex <- NULL

expression_t <- readRDS(paste0(data_path, "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expression_t) <- sub(".RAA", "", colnames(expression_t))
rownames(expression_t) <- paste0(rownames(expression_t), "_trans")

expression_p <- readRDS(paste0(data_path, "AFHRI_B_proteomics_QC_symbol.RDS"))
colnames(expression_p) <- sub(".RAA","", colnames(expression_p))
rownames(expression_p) <- paste0(rownames(expression_p), "_prot")

residuals_t <- readRDS(paste0(data_path, "residuals_trans.RDS"))
colnames(residuals_t) <- sub(".RAA","", colnames(residuals_t))
rownames(residuals_t) <- paste0(rownames(residuals_t), "_t.res")

residuals_p <- readRDS(paste0(data_path, "residuals_prot.RDS"))
colnames(residuals_p) <- sub(".RAA", "", colnames(residuals_p))
rownames(residuals_p) <- paste0(rownames(residuals_p), "_p.res")

ratios <- as.data.frame(readRDS(paste0(data_path, "prot_trans_ratios.RDS")))
colnames(ratios) <- sub(".RAA", "", colnames(ratios))
rownames(ratios) <- paste0(rownames(ratios), "_ratios")

expression_m <- t(metabolomics[,-(1)])

# merge all expression data
pheno.temp <- merge(data.frame(id=fam$V1,
                               sex=ifelse(fam$V5==2,0,1),
                               stringsAsFactors = F),
                    covariates,
                    by.x = "id", by.y = "externID", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(expression_t), t(expression_t), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(expression_p), t(expression_p), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(residuals_t), t(residuals_t), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(residuals_p), t(residuals_p), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(ratios), t(ratios), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

colnames(pheno.temp) <- gsub("\\.", "\\-", colnames(pheno.temp))
colnames(pheno.temp) <- gsub("fibro-score", "fibro.score", colnames(pheno.temp))
colnames(pheno.temp) <- gsub("X1060P11-3_trans", "1060P11.3_trans", colnames(pheno.temp))
colnames(pheno.temp) <- gsub("_t-res", "_t.res", colnames(pheno.temp))
colnames(pheno.temp) <- gsub("_p-res", "_p.res", colnames(pheno.temp))

test <- data.frame(colnames(pheno.temp))

pheno.temp <- merge(pheno.temp,
                    data.frame(id=colnames(expression_m), t(expression_m), stringsAsFactors = F),
                    by = "id", all.x = T, sort=F)

pheno.temp <- merge(data.frame(id=fam$V1,
                               stringsAsFactors = F),
                    pheno.temp,
                    by = "id", all.x = T, sort=F)

identical(as.character(fam$V1), as.character(pheno.temp$id))
saveRDS(pheno.temp, file=paste0(data_path, "expression_all.RDS"))
write.table(pheno.temp,
            file = paste0(geno_path, "AFHRI_B.pheno"),
            row.names = F, col.names = T, quote = F,
            sep="\t")


# creation of the gwaa-data object resorts id's!
identical(idnames(AFHRI_B_imp), pheno.temp$id)
pheno.temp <- merge(data.frame(id=idnames(AFHRI_B_imp),
                               stringsAsFactors = F),
                    pheno.temp,
                    by = "id", all.x = T, sort=F)

identical(idnames(AFHRI_B_imp), pheno.temp$id)

# now we can update the gwaa-data object
phdata(AFHRI_B_imp) <- pheno.temp

save(AFHRI_B_imp,
     file=paste0(get.path("genotype"),
                 "AFHRI_B_imputed_gwaa_data.RData"))

saveRDS(AFHRI_B_imp,
        file=paste0(get.path("genotype"),
                    "AFHRI_B_imputed_gwaa_data.RDS"))






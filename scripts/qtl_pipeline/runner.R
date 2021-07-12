#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd("~/work/symAtrial_QTL/scripts")

source("helper/helper.R")
local=F
## uncomment this if you want to run from command-line
#
# if (length(args)==0) {
#   stop(get.message.usage(), call.=FALSE)
# } else if (length(args)==1) {
#   data_dir <- args[1]
# }

## uncomment this if you want to run from RStudio

data_dir <- get.path("data", local)

#################################################################################################

## data generation
#
# this part runs the preprocessing steps for data generation and creates annotation files
# set the if to TRUE if you want to run the corresponding steps

if(F){
  system("plink --bfile AFHRI-B --recodeA --out raw_genotype")
}

if(F){
  source("retrieve_data/gene_annotations/protein_locations_biomart.R")
  source("retrieve_data/gene_annotations/transcript_locations_biomart.R")

  get.protein.locations()
  get.transcript.locations()
}

if(F){
  source("retrieve_data/snp_annotations/snp_locations.R")
  get.snp.locations()
}

if(F){
  source("preprocessing/data_parsing/plink2snp.R")
  get.processed.genotype()
}


#################################################################################################

## the following part runs the whole analysis inclusive plotting
## make sure that the paths (in the helper class) are set right


# load data

data_path <- get.path("dataset", local)

expression_p <- readRDS(paste0(data_path, "AFHRI_B_proteomics_QC_symbol.RDS"))
colnames(expression_p) <- sub(".RAA","", colnames(expression_p))

expression_t <- readRDS(paste0(data_path, "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expression_t) <- sub(".RAA", "", colnames(expression_t))

residuals_t <- readRDS(paste0(data_path, "residuals_trans.RDS"))
colnames(residuals_t) <- sub(".RAA","", colnames(residuals_t))

residuals_p <- readRDS(paste0(data_path, "residuals_prot.RDS"))
colnames(residuals_p) <- sub(".RAA", "", colnames(residuals_p))

ratios <- readRDS(paste0(data_path, "prot_trans_ratios.RDS"))
colnames(ratios) <- sub(".RAA","", colnames(ratios))

locations_p <- readRDS(paste0(get.path("locations"), "protein_locations_ensembl_tss.RDS"))
locations_t <- readRDS(paste0(get.path("locations", local), "transcript_locations_ensembl_tss.RDS"))

# cut locations
source("preprocessing/data_parsing/cut_tables.R")

new_expression_p <- cut.locations(expression_p, locations_p)
new_expression_t <- cut.locations(expression_t, locations_t)
new_expression_res_t <- cut.locations(residuals_t, locations_t)
new_expression_res_p <- cut.locations(residuals_p, locations_p)
new_expression_ratios <- cut.locations(ratios, locations_t)

# drop samples
library(data.table)

processed_genotype_file <- get.path("base genotype imputed", local)
genotype <- fread(processed_genotype_file, sep="\t", header=T)
rownames(genotype) <- genotype[[1]]
genotype[,1] <- NULL

results_t <- common.samples(new_expression_t, genotype)
new_gen_t <- results_t[[1]]
new_exp_t <- results_t[[2]]

results_p <- common.samples(new_expression_p, genotype)
new_gen_p <- results_p[[1]]
new_exp_p <- results_p[[2]]

results_res_t <- common.samples(new_expression_res_t, genotype)
new_gen_res_t <- results_res_t[[1]]
new_exp_res_t <- results_res_t[[2]]

results_res_p <- common.samples(new_expression_res_p, genotype)
new_gen_res_p <- results_res_p[[1]]
new_exp_res_p <- results_res_p[[2]]

results_ratios <- common.samples(new_expression_ratios, genotype)
new_gen_ratios <- results_ratios[[1]]
new_exp_ratios <- results_ratios[[2]]

geno_path <- get.path("genotype", local)

geno_t <- paste0(geno_path, "genotype_imputed_common_samples_t_raa.txt")
geno_p <- paste0(geno_path, "genotype_imputed_common_samples_p_raa.txt")
geno_res_t <- paste0(geno_path, "genotype_imputed_common_samples_res_t_raa.txt")
geno_res_p <- paste0(geno_path, "genotype_imputed_common_samples_res_p_raa.txt")
geno_ratios <- paste0(geno_path, "genotype_imputed_common_samples_ratios.txt")

write.table(new_gen_t, geno_t, sep="\t", quote=F, col.names=NA)
write.table(new_gen_p, geno_p, sep="\t", quote=F, col.names=NA)
write.table(new_gen_res_t, geno_res_t, sep="\t", quote=F, col.names=NA)
write.table(new_gen_res_p, geno_res_p, sep="\t", quote=F, col.names=NA)
write.table(new_gen_ratios, geno_ratios, sep="\t", quote=F, col.names=NA)


# data normalization
if(T) {
  source("preprocessing/correction/normalization/norm.R")
  norm_exp_t <- correct.linewise(new_exp_t)
  norm_exp_p <- correct.linewise(new_exp_p)
  norm_exp_res_t <- correct.linewise(new_exp_res_t)
  norm_exp_res_p <- correct.linewise(new_exp_res_p)
  norm_exp_ratios <- correct.linewise(new_exp_ratios)
}

# impute
source("preprocessing/imputation/run_imputation.R")

imputed_t_norm <- impute.data(norm_exp_t)
imputed_p_norm <- impute.data(norm_exp_p)
imputed_res_t_norm <- impute.data(norm_exp_res_t)
imputed_res_p_norm <- impute.data(norm_exp_res_p)
imputed_ratios_norm <- impute.data(norm_exp_ratios)

imputed_t <- impute.data(new_exp_t)
imputed_p <- impute.data(new_exp_p)
imputed_res_t <- impute.data(new_exp_res_t)
imputed_res_p <- impute.data(new_exp_res_p)
imputed_ratios <- impute.data(new_exp_ratios)

peer_path <- get.path("peer", local)
edir_norm <- paste0(peer_path, "normalized/eqtl")
pdir_norm <- paste0(peer_path, "normalized/pqtl")
edir_res_norm <- paste0(peer_path, "normalized/eqtl_res")
pdir_res_norm <- paste0(peer_path, "normalized/pqtl_res")
dir_ratios_norm <- paste0(peer_path, "normalized/ratios")

edir <- paste0(peer_path, "natural/eqtl")
pdir <- paste0(peer_path, "natural/pqtl")
edir_res <- paste0(peer_path, "natural/eqtl_res")
pdir_res <- paste0(peer_path, "natural/pqtl_res")
dir_ratios <- paste0(peer_path, "natural/ratios")

dir.create(edir_norm, showWarnings=F, recursive=T)
dir.create(pdir_norm, showWarnings=F, recursive=T)
dir.create(edir_res_norm, showWarnings=F, recursive=T)
dir.create(pdir_res_norm, showWarnings=F, recursive=T)
dir.create(dir_ratios_norm, showWarnings=F, recursive=T)

dir.create(edir, showWarnings=F, recursive=T)
dir.create(pdir, showWarnings=F, recursive=T)
dir.create(edir_res, showWarnings=F, recursive=T)
dir.create(pdir_res, showWarnings=F, recursive=T)
dir.create(dir_ratios, showWarnings=F, recursive=T)

write.table(imputed_t, paste0(edir, "/imputed_expression_t.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_p, paste0(pdir, "/imputed_expression_p.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_res_t, paste0(edir_res, "/imputed_expression_t.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_res_p, paste0(pdir_res, "/imputed_expression_p.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_ratios, paste0(dir_ratios, "/imputed_expression_ratios.txt"), sep="\t", quote=F, col.names=NA)

write.table(imputed_t_norm, paste0(edir_norm, "/imputed_expression_t.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_p_norm, paste0(pdir_norm, "/imputed_expression_p.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_res_t_norm, paste0(edir_res_norm, "/imputed_expression_t.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_res_p_norm, paste0(pdir_res_norm, "/imputed_expression_p.txt"), sep="\t", quote=F, col.names=NA)
write.table(imputed_ratios_norm, paste0(dir_ratios_norm, "/imputed_expression_ratios.txt"), sep="\t", quote=F, col.names=NA)

# prepare covariates
source("preprocessing/covariates/add_covariates.R")
covariates_t <- prepare.covariates(imputed_t)
covariates_p <- prepare.covariates(imputed_p)
covariates_res_t <- prepare.covariates(imputed_res_t)
covariates_res_p <- prepare.covariates(imputed_res_p)
covariates_ratios <- prepare.covariates(imputed_ratios, local=local)

fibro_t <- prepare.covariates(imputed_t, only.fibro=T)
fibro_p <- prepare.covariates(imputed_p, only.fibro=T)
fibro_res_t <- prepare.covariates(imputed_res_t, only.fibro=T)
fibro_res_p <- prepare.covariates(imputed_res_p, only.fibro=T)
fibro_ratios <- prepare.covariates(imputed_ratios, only.fibro=T, local=local)

# peer
source("preprocessing/correction/peer/run_peer.R")
run.peer(imputed_t, c(0:30), paste0(edir, "/factors_fibro"), fibro_t)
run.peer(imputed_p, c(0:30), paste0(pdir, "/factors_fibro"), fibro_p)
run.peer(imputed_res_t, c(0:30), paste0(edir_res, "/factors_fibro"), fibro_res_t)
run.peer(imputed_res_p, c(0:30), paste0(pdir_res, "/factors_fibro"), fibro_res_p)
run.peer(imputed_ratios, c(0:30), paste0(dir_ratios, "/factors_fibro"), fibro_ratios)

run.peer(imputed_t_norm, c(0:30), paste0(edir_norm, "/factors_fibro"), fibro_t)
run.peer(imputed_p_norm, c(0:30), paste0(pdir_norm, "/factors_fibro"), fibro_p)
run.peer(imputed_res_t_norm, c(0:30), paste0(edir_res_norm, "/factors_fibro"), fibro_res_t)
run.peer(imputed_res_p_norm, c(0:30), paste0(pdir_res_norm, "/factors_fibro"), fibro_res_p)
run.peer(imputed_ratios_norm, c(0:30), paste0(dir_ratios_norm, "/factors_fibro"), fibro_ratios)

run.peer(imputed_p, c(0:30), paste0(pdir, "/factors_no_cov"))
run.peer(imputed_t, c(0:30), paste0(edir, "/factors_no_cov"))
run.peer(imputed_res_t, c(0:30), paste0(edir_res, "/factors_no_cov"))
run.peer(imputed_res_p, c(0:30), paste0(pdir_res, "/factors_no_cov"))
run.peer(imputed_ratios, c(0:30), paste0(dir_ratios, "/factors_no_cov"))

run.peer(imputed_p_norm, c(0:30), paste0(pdir_norm, "/factors_no_cov"))
run.peer(imputed_t_norm, c(0:30), paste0(edir_norm, "/factors_no_cov"))
run.peer(imputed_res_t_norm, c(0:30), paste0(edir_res_norm, "/factors_no_cov"))
run.peer(imputed_res_p_norm, c(0:30), paste0(pdir_res_norm, "/factors_no_cov"))
run.peer(imputed_ratios_norm, c(0:30), paste0(dir_ratios_norm, "/factors_no_cov"))

run.peer(imputed_p, c(0:30), paste0(pdir, "/factors_cov"), covariates_p)
run.peer(imputed_t, c(0:30), paste0(edir, "/factors_cov"), covariates_t)
run.peer(imputed_res_t, c(0:30), paste0(edir_res, "/factors_cov"), covariates_res_t)
run.peer(imputed_res_p, c(0:30), paste0(pdir_res, "/factors_cov"), covariates_res_p)
run.peer(imputed_ratios, c(0:30), paste0(dir_ratios, "/factors_cov"), covariates_ratios)

run.peer(imputed_t_norm, c(0:30), paste0(edir_norm, "/factors_cov"), covariates_t)
run.peer(imputed_p_norm, c(0:30), paste0(pdir_norm, "/factors_cov"), covariates_p)
run.peer(imputed_res_t_norm, c(0:30), paste0(edir_res_norm, "/factors_cov"), covariates_res_t)
run.peer(imputed_res_p_norm, c(0:30), paste0(pdir_res_norm, "/factors_cov"), covariates_res_p)
run.peer(imputed_ratios_norm, c(0:30), paste0(dir_ratios_norm, "/factors_cov"), covariates_ratios)

run.peer(imputed_p_norm, c(0:30), paste0(pdir_norm, "/factors_all"))
run.peer(imputed_t_norm, c(0:30), paste0(edir_norm, "/factors_all"))
run.peer(imputed_res_t_norm, c(0:30), paste0(edir_res_norm, "/factors_all"))
run.peer(imputed_res_p_norm, c(0:30), paste0(pdir_res_norm, "/factors_all"))
run.peer(imputed_ratios_norm, c(0:30), paste0(dir_ratios_norm, "/factors_all"))

run.peer(imputed_t_norm, c(0:30), paste0(edir_norm, "/factors_pop"), covariates_t)
run.peer(imputed_p_norm, c(0:30), paste0(pdir_norm, "/factors_pop"), covariates_p)
run.peer(imputed_res_t_norm, c(0:30), paste0(edir_res_norm, "/factors_pop"), covariates_res_t)
run.peer(imputed_res_p_norm, c(0:30), paste0(pdir_res_norm, "/factors_pop"), covariates_res_p)
run.peer(imputed_ratios_norm, c(0:30), paste0(dir_ratios_norm, "/factors_pop"), covariates_ratios)


# eQTL
source("eqtl/eqtl.R")

# path2snplocs <- get.path("snplocs")
path2snplocs <- get.path("snplocs imputed", local)

path2genelocs_t <- paste0(get.path("locations", local), "transcript_locations_ensembl.txt")
path2genelocs_p <- paste0(get.path("locations", local), "protein_locations_ensembl.txt")
path2genelocs_res <- paste0(get.path("locations", local), "residuals_locations_ensembl.txt")
path2genelocs_ratios <- paste0(get.path("locations", local), "residuals_locations_ensembl.txt")


cov_subdirs <- c("factors_cov", "factors_fibro", "factors_no_cov",
                 "factors_all", "factors_pop")
norm_subdirs <- c("natural", "normalized")
model_strings <- c("linear")
types <- c("eqtl", "pqtl", "eqtl_res", "pqtl_res", "ratios")

for(cov_sd in cov_subdirs){
  for(norm_sd in norm_subdirs){
    for(model in model_strings){
      for(type in types){
        
        if(type == "eqtl"){
          geno <- geno_t
          path2genelocs <- path2genelocs_t
          i <- "t"
        }
        
        if(type == "pqtl"){
          geno <- geno_p
          path2genelocs <- path2genelocs_p
          i <- "p"
        }
    
        if(type == "eqtl_res"){
          geno <- geno_res_t
          path2genelocs <- path2genelocs_res
          i <- "t"
        }
        
        if(type == "pqtl_res"){
          geno <- geno_res_p
          path2genelocs <- path2genelocs_res
          i <- "p"
        }
        
        if(type == "ratios"){
          geno <- geno_ratios
          path2genelocs <- path2genelocs_ratios
          i <- "ratios"
        }
        if(i=="ratios"){
          file <- paste0("expression_", i , ".txt")
        } else {
          file <- paste0("imputed_expression_", i , ".txt")
        }
        
        eqtls(paste(peer_path, norm_sd, type, file, sep="/"),
              geno,
              path2snplocs,
              path2genelocs,
              model, 
              paste(get.path("results", local), "imputed", "cis", model, type, norm_sd, cov_sd, sep = "/"),
              paste(peer_path, norm_sd, type, cov_sd, sep="/"),
              trans = trans
        )
        
      }
    }
  }
}


# peer table analysis
source("analysis/peer_tables/peer_result.R")

for(cov_sd in cov_subdirs){
  for(norm_sd in norm_subdirs){
    for(model in model_strings){
      for(type in types){
        
        dir <- paste(get.path("results"), "current", "cis", model, type, norm_sd, cov_sd, sep="/")
        outdir <- paste(get.path("results"), "current", "peer", model, type, norm_sd, cov_sd, sep="/")
        
        peer.write.table(dir, outdir)
      }
    }
  }
}

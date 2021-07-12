#!/usr/bin/env Rscript

##
#  runs eqtl analyses for all combinations of model, normalization, type, etc.
#
#
#
#

setwd("~/work/symAtrial_QTL/scripts")

# load libraries
library(batchtools)

# source scripts
source("batchtools_helper.R")

do.stuff <- function(k){
  #k=1
  cov_subdirs <- c("factors_cov", "factors_fibro", "factors_no_cov",
                   "factors_all", "factors_pop")
  norm_subdirs <- c("natural", "normalized")
  model_strings <- c("linear")
  types <- c("eqtl", "pqtl", "eqtl_res", "pqtl_res", "ratios")
  nks <- list(0:5, 6:10, 11:15, 16:20, 21:25, 26:30)
  
  batch.table <- expand.grid(cov_subdirs,
                             norm_subdirs,
                             model_strings,
                             types,
                             nks,
                             stringsAsFactors = F)
  
  nkp <- NULL
  cov_sd <- batch.table[k, 1]
  norm_sd <- batch.table[k, 2]
  model <- batch.table[k, 3]
  type <- batch.table[k, 4]
  nkp <- batch.table[k, 5][[1]]
  print(batch.table[k, ])
  trans <- F
  
  setwd("~work/symAtrial_QTL/scripts")
  
  source("preprocessing/correction/peer/run_peer.R")
  source("helper/helper.R")
  source("eqtl/eqtl.R")
  
  local <- F
  data_path <- get.path("dataset", local)
  peer_path <- get.path("peer", local)
  
  path2snplocs <- get.path("snplocs imputed", local)
  geno_path <- get.path("genotype", local)
  
  geno_p <- paste0(geno_path, "genotype_imputed_common_samples_p_raa.txt")
  geno_t <- paste0(geno_path, "genotype_imputed_common_samples_t_raa.txt")
  geno_res_t <- paste0(geno_path, "genotype_imputed_common_samples_res_t_raa.txt")
  geno_res_p <- paste0(geno_path, "genotype_imputed_common_samples_res_p_raa.txt")
  geno_ratios <- paste0(geno_path, "genotype_imputed_common_samples_ratios.txt")
  
  path2genelocs_p <- paste0(get.path("locations"), "protein_locations_ensembl.txt")
  path2genelocs_t <- paste0(get.path("locations"), "transcript_locations_ensembl.txt")
  path2genelocs_res <- paste0(get.path("locations"), "residuals_locations_ensembl.txt")
  path2genelocs_ratios <- paste0(get.path("locations", local), "residuals_locations_ensembl.txt")

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
    file <- paste0("imputed_expression_", i , ".txt")
  } else {
    file <- paste0("imputed_expression_", i , ".txt")
  }
  
  if(!is.null(nkp)){
    files <- paste0(paste(peer_path, norm_sd, type, cov_sd, sep="/"),
                    paste0("/peer_factors_nk", sprintf("%02d", nkp),".txt"))
    eqtls(expression_file=paste(peer_path, norm_sd, type, file, sep="/"),
          genotype = geno,
          snplocs = path2snplocs,
          genelocs = path2genelocs,
          model = model, 
          outdir = paste(get.path("results"), "imputed", "cis", model, type, norm_sd, cov_sd, sep = "/"),
          covariate_dir = paste(peer_path, norm_sd, type, cov_sd, sep="/"),
          trans=trans,
          files = files
    )
  } else {
    eqtls(expression_file=paste(peer_path, norm_sd, type, file, sep="/"),
          genotype = geno,
          snplocs = path2snplocs,
          genelocs = path2genelocs,
          model = model, 
          outdir = paste(get.path("results"), "imputed", "cis", model, type, norm_sd, cov_sd, sep = "/"),
          covariate_dir = paste(peer_path, norm_sd, type, cov_sd, sep="/"),
          trans=trans
    )
  }
}

cov_subdirs <- c("factors_cov", "factors_fibro", "factors_no_cov",
                 "factors_all", "factors_pop")
norm_subdirs <- c("natural", "normalized")
model_strings <- c("linear")
types <- c("eqtl", "pqtl", "eqtl_res", "pqtl_res", "ratios")
nks <- list(0:5, 6:10, 11:15, 16:20, 21:25, 26:30)

batch.table <- expand.grid(cov_subdirs,
                           norm_subdirs,
                           model_strings,
                           types,
                           nks,
                           stringsAsFactors = F)

names <- paste0(batch.table$Var4,
                sapply(batch.table$Var5,
                       FUN = function(x){x[1]}), "_",
                sapply(batch.table$Var5,
                       FUN = function(x){x[length(x)]}))

run.batchtools(do.stuff, 1:dim(batch.table)[1], more.args=list(),
               names, "tmp_dir_QTL",
               clean.up=F,
               resources=list(partition="icb_cpu",
                              memory="40G",
                              ncpus = 1,
                              measure.memory = TRUE,
                              walltime="12:00:00"))

q(save = "no")

#!/usr/bin/env Rscript

##
#  runs eqtl analyses for all combinations of model, normalization, type, etc.
#
#
#
#

setwd("~/work/symAtrial_QTL/scripts")

# load libraries
library(BatchJobs)

# source scripts
source("batchtools_helper.R")

do.stuff <- function(k){
  #k=1
  cov_subdirs <- c("factors_fibro", "factors_cov", "factors_no_cov")
  norm_subdirs <- c("natural", "normalized")
  model_strings <- c("linear")
  types <- c("eqtl", "pqtl", "eqtl_res", "pqtl_res", "ratios")
  
  batch.table <- expand.grid(cov_subdirs,
                             norm_subdirs,
                             model_strings,
                             types,
                             stringsAsFactors = F)
  
  #k=1
  cov_sd <- batch.table[k, 1]
  norm_sd <- batch.table[k, 2]
  model <- batch.table[k, 3]
  type <- batch.table[k, 4]
  batch.table[k, ]
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
    #cov_file <- "peer_factors_nk12.txt"
  }
  
  if(type == "pqtl"){
    geno <- geno_p
    path2genelocs <- path2genelocs_p
    i <- "p"
    #cov_file <- "peer_factors_nk08.txt"
  }
  
  if(type == "pqtl_trans"){
    geno <- geno_p
    path2genelocs <- path2genelocs_p
    i <- "p"
    trans <- T
    # cov_file <- "peer_factors_nk08.txt"
  }
  
  if(type == "eqtl_res"){
    geno <- geno_res_t
    path2genelocs <- path2genelocs_res
    i <- "t"
    #cov_file <- "peer_factors_nk?.txt"
  }
  
  if(type == "pqtl_res"){
    geno <- geno_res_p
    path2genelocs <- path2genelocs_res
    i <- "p"
    #cov_file <- "peer_factors_nk?.txt"
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

cov_subdirs <- c("factors_fibro", "factors_cov", "factors_no_cov")
norm_subdirs <- c("natural", "normalized")
model_strings <- c("linear")
types <- c("eqtl", "pqtl", "eqtl_res", "pqtl_res", "ratios")

batch.table <- expand.grid(cov_subdirs,
                           norm_subdirs,
                           model_strings,
                           types,
                           stringsAsFactors = F)

run.batchjobs(do.stuff, 1:dim(batch.table)[1], more.args=list(),
              "QTL", "tmp_dir_QTL", clean.up=F,
              resources = list(partition="my_queue",
                               memory="80G",
                               ncpus = 1,
                               measure.memory = TRUE,
                               walltime="48:00:00"))

q(save = "no")

#!/usr/bin/env Rscript

##
#   runs the peer analysis solely 
#
#
#

setwd("~/work/symAtrial_QTL/scripts")

# load libraries
library(BatchJobs)

# source scripts
source("batchtools_helper.R")

do.stuff <- function(k){

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
  #nkp <- batch.table[k, 5]
  batch.table[k, ]
  
  setwd("~/work/symAtrial_QTL/scripts")
  
  source("preprocessing/correction/peer/run_peer.R")
  source("helper/helper.R")
  
  local <- F
  data_path <- get.path("dataset", local)
  peer_path <- get.path("peer", local)
  
  dir <- paste0(peer_path, norm_sd, "/", type)
  if(type %in% c("eqtl", "eqtl_res")){
    imputed <- read.table(paste0(dir, "/imputed_expression_t.txt"),
                          sep="\t", header=T, row.names = 1)
  }else if(type %in% c("pqtl", "pqtl_res")){
    imputed <- read.table(paste0(dir, "/imputed_expression_p.txt"),
                          sep="\t", header=T, row.names = 1)
  }else if(type %in% c("ratios")){
    imputed <- read.table(paste0(dir, "/imputed_expression_ratios.txt"),
                          sep="\t", header=T, row.names = 1)
  }
  
  source("preprocessing/covariates/add_covariates.R")
  source("preprocessing/correction/peer/run_peer.R")
  
  if(cov_sd=="factors_cov"){
    if(type=="mqtl"){
      covariates <- prepare.covariates(imputed, batch = T, local = local)
    }else{
      covariates <- prepare.covariates(imputed, local = local)
    }
    run.peer(imputed, c(0:30), paste0(dir, "/", cov_sd), covariates)
  } else if (cov_sd=="factors_fibro"){
    fibro <- prepare.covariates(imputed, only.fibro=T, local = local)
    run.peer(imputed, c(0:30), paste0(dir, "/", cov_sd), fibro)
  } else if (cov_sd=="factors_no_cov"){
    run.peer(imputed, c(1:30), paste0(dir, "/", cov_sd))
  }
  
  print(paste0(k, " done"))
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
              "peer", "tmp_dir_peer")


q(save="no")

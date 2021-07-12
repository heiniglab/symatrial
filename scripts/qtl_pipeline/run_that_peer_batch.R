#!/usr/bin/env Rscript

##
#   runs the peer analysis solely 
#
#
#

setwd("~/work/symAtrial_QTL/scripts")

# load libraries
library(batchtools)

# source scripts
source("batchtools_helper.R")

do.stuff <- function(k){

  cov_subdirs <- c("factors_cov", "factors_fibro", "factors_no_cov",
                   "factors_all", "factors_pop")
  norm_subdirs <- c("natural", "normalized")
  model_strings <- c("linear")
  types <- c("eqtl", "pqtl", "eqtl_res", "pqtl_res", "ratios")
  nks <- c(0:30)
  
  batch.table <- expand.grid(cov_subdirs,
                             norm_subdirs,
                             model_strings,
                             types,
                             nks,
                             stringsAsFactors = F)
  
  #k=1
  cov_sd <- batch.table[k, 1]
  norm_sd <- batch.table[k, 2]
  model <- batch.table[k, 3]
  type <- batch.table[k, 4]
  nkp <- batch.table[k, 5]
  print(batch.table[k, ])
  
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
    covariates <- prepare.covariates(imputed, local = local)
    run.peer(imputed, nkp, paste0(dir, "/", cov_sd), covariates)
  } else if(cov_sd=="factors_all"){
    covariates <- prepare.covariates(imputed, local = local, pop = 3)
    run.peer(imputed, nkp, paste0(dir, "/", cov_sd), covariates)
  } else if(cov_sd=="factors_pop"){
    covariates <- prepare.covariates(imputed, local = local, pop = 3, only.pop=T)
    run.peer(imputed, nkp, paste0(dir, "/", cov_sd), covariates)
  } else if (cov_sd=="factors_fibro"){
    fibro <- prepare.covariates(imputed, only.fibro=T, local = local)
    run.peer(imputed, c(0:30), paste0(dir, "/", cov_sd), fibro)
  } else if (cov_sd=="factors_no_cov" & nkp>0){
    run.peer(imputed, nkp, paste0(dir, "/", cov_sd))
  }
  
  print(paste0(k, " done"))
}

cov_subdirs <- c("factors_cov", "factors_fibro", "factors_no_cov",
                 "factors_all", "factors_pop")
norm_subdirs <- c("natural", "normalized")
model_strings <- c("linear")
types <- c("eqtl", "pqtl", "eqtl_res", "pqtl_res", "ratios")
nks <- c(0:30)

batch.table <- expand.grid(cov_subdirs,
                           norm_subdirs,
                           model_strings,
                           types,
                           nks,
                           stringsAsFactors = F)

names <- paste0("PEER", gsub("qtl", "", batch.table$Var4, batch.table$Var5))

run.batchtools(do.stuff, 1:dim(batch.table)[1], more.args=list(),
               names, "tmp_dir_PEER",
               clean.up=F,
               resources=list(partition="icb_cpu",
                              memory="8G",
                              ncpus = 1,
                              measure.memory = TRUE,
                              walltime="48:00:00"))

q(save="no")

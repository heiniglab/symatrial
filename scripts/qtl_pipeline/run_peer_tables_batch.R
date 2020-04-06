setwd("~/work/symAtrial_QTL/scripts")

# load libraries
library(batchtools)

# source scripts
source("batchtools_helper.R")

cov_subdirs <- c("factors_cov", "factors_fibro", "factors_no_cov")
norm_subdirs <- c("normalized", "natural")
model_strings <- c("linear")
types <- c("eqtl", "eqtl_res", "pqtl", "pqtl_res", "ratios")

batch.table <- expand.grid(cov_subdirs,
                           norm_subdirs,
                           model_strings,
                           types,
                           stringsAsFactors = F)

do.stuff <- function(k){
# for (k in 1:dim(batch.table)[1]){
  #k=1
  setwd("~/work/symAtrial_QTL/scripts")
  
  source("helper/helper.R")
  
  # peer table analysis
  source("analysis/peer_tables/peer_result.R")
  
  cov_subdirs <- c("factors_cov", "factors_fibro", "factors_no_cov")
  norm_subdirs <- c("normalized", "natural")
  model_strings <- c("linear")
  types <- c("eqtl", "eqtl_res", "pqtl", "pqtl_res", "ratios")
  
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
  print(batch.table[k, ])

  dir <- paste(get.path("results"), "imputed", "cis", model, type, norm_sd, cov_sd, sep="/")
  outdir <- paste(get.path("results"), "peer", "imputed", model, type, norm_sd, cov_sd, sep="/")
  peer.write.table(dir, outdir)
}

run.batchtools(do.stuff, 5:dim(batch.table)[1], more.args=list(),
               "PEER_tables", "tmp_dir_PEER_tables",
               clean.up=F,
               resources=list(partition="icb_cpu",
                                memory="8G",
                                ncpus = 1,
                                measure.memory = TRUE,
                                walltime="01:00:00"))

dir = "~/work/symAtrial_QTL/results/current/peer/imputed/linear"
files <- list.files(path=dir, pattern="*.RDS$", full.names=T, recursive=T)
files <- files[!(files %in% c("~/work/symAtrial_QTL/results/current/peer/imputed/QTL_table.RDS",
                              "~/work/symAtrial_QTL/results/current/peer/imputed/QTL_table_imp.RDS",
                              "~/work/symAtrial_QTL/results/current/peer/imputed/QTL_table_imp_new.RDS"))]
data <- lapply(files, readRDS)
names(data) <- sub("/home/icb/ines.assum/projects/symAtrial_QTL/results/current/peer/imputed/",
                   "",
                   sub("/peer_table.RDS", "", files))
data <- lapply(seq_along(data), function(list.index){
  run <- unlist(strsplit(names(data)[list.index], "/"))
  data.frame(model = run[1], type = run[2], norm_sd = run[3], cov_sd = run[4],
             Nk=as.integer(sub("eqtl_nk", "", sub(".RDS", "", rownames(data[[list.index]])))),
             data[[list.index]])
})
QTL.res <-  do.call(rbind, data)
saveRDS(QTL.res, file="/home/icb/ines.assum/projects/symAtrial_QTL/results/current/peer/QTL_table_imp_new.RDS")

q(save = "no")

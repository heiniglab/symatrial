## runs peer correction for the given expression file
#
# @expression: expression file
# @nk_vector: the nk numbers given to peer
# @directory: output directory

run.peer <- function(expression, nk_vector, directory, covariates = NULL){
  source("preprocessing/correction/peer/peer.R")
  transposed <- t(expression)
  for(nk in nk_vector){
    if(nk==0){
      dir.create(directory, showWarnings = F, recursive = T)
      write.table(t(covariates), paste0(directory, "/peer_factors_nk", sprintf("%02d", nk),".txt"), sep="\t", quote=F, col.names=NA, row.names = T)
    }else{
      # exp_peer <- get.peer.factors(transposed, get.residuals = T, Nk=nk)
      peer_factors <- get.peer.factors(transposed, covariates=covariates,  get.residuals = F, Nk=nk)
      
      # write.table(t(exp_peer), paste0(directory, "/peer_nk", sprintf("%02d", nk),".txt"), sep="\t", quote=F, col.names=NA)
      dir.create(directory, showWarnings = F, recursive = T)
      write.table(t(peer_factors), paste0(directory, "/peer_factors_nk", sprintf("%02d", nk),".txt"), sep="\t", quote=F, col.names=NA, row.names = T)
    }
  }
}
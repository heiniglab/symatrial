# SNP file from plink output

get.processed.genotype <- function() {
  source("helper/helper.R")
  library(data.table)
  
  plink_out <-
    fread(paste0(get.path("data"), "plink/genotype.raw"),
          sep = " ",
          header = T)
  
  rn <- plink_out$IID
  
  
  plink_out <- plink_out[, -1:-6, with = F]
  
  rownames(plink_out) <- rn
  
  snp_file <- transpose(plink_out)
  rownames(snp_file) <- colnames(plink_out)
  colnames(snp_file) <- rownames(plink_out)
  
  rownames(snp_file) <- gsub('.{2}$', '', rownames(snp_file))
  
  write.table(
    snp_file,
    paste0(get.path("data"), "plink/processed_genotype.txt"),
    sep = "\t",
    quote = F,
    col.names = NA
  )
}

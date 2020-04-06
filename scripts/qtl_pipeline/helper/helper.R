
get.project.path <- function() {
  return("/home/icb/ines.assum/projects/symAtrial_QTL")
}

get.path <- function(dest, local = F) {
    
  basepath <- ifelse(local,
                     "/Volumes/Server/symAtrial_QTL",
                     "/home/icb/ines.assum/projects/symAtrial_QTL") 
  basepath2 <- file.path(basepath, "data/current/")
  basepath3 <- file.path(basepath, "results/current/")
  
  if (dest == "dataset") {
    return(file.path(basepath, "data/current/AFHRI_B/"))
  }
  if (dest == "script") {
    return(file.path("/Users/ines/work/symAtrial_QTL/scripts"))
  }
  if (dest == "data") {
    return(basepath2)
  }
  if (dest == "data10") {
    return(basepath2)
  }
  if (dest == "genotype") {
    return(file.path(basepath2, "genotype/"))
  }
  if (dest == "metabolome") {
    return(file.path(basepath2, "metabolome/"))
  }
  if (dest == "gwas") {
    return(file.path(basepath2, "annotation/gwas/"))
  }
  if (dest == "gwas2") {
    return(file.path(basepath2, "annotation/gwas_imputed/"))
  }
  if (dest == "gwas3") {
    return(file.path(basepath2, "annotation/gwas3/"))
  }
  if (dest == "risk") {
    return(file.path(basepath2, "polygenic_risk_scores/"))
  }
  if (dest == "gtex") {
    return(file.path(basepath2, "gtex/v7/"))
  }
  if (dest == "plasma") {
    return(file.path(basepath2, "plasma_pQTL/"))
  }
  if (dest == "peer") {
    return(file.path(basepath2, "peer/"))
  }
  if (dest == "eQTL") {
    return(file.path(basepath2, "eQTL/"))
  }
  if (dest == "pQTL") {
    return(file.path(basepath2, "pQTL/"))
  }
  if (dest == "plots") {
    return(file.path(basepath3, "plots/"))
  }
  if (dest == "locations_ensembl"){
    return(file.path(basepath2, "annotation/biomart_hg19_strand.txt"))
  }
  if (dest == "locations"){
    return(file.path(basepath2, "annotation/"))
  }
  if (dest == "results"){
    return(basepath3)
  }
  if (dest == "base genotype"){
    return(file.path(basepath2, "genotype/processed_genotype.txt"))
  }
  if (dest == "base genotype imputed"){
    return(file.path(basepath2, "genotype/AFHRI_B_imputed_preprocessed_genotypes.txt"))
  }
  if (dest == "snplocs"){
    return(file.path(basepath2, "annotation/snplocs.txt"))
  }
  if (dest == "snplocs imputed"){
    return(file.path(basepath2, "annotation/snplocs_imputed.txt"))
  }
  return(NULL)
}




loaddata <- function() {
  path <- get.path("dataset")
  genotype <- readRDS(paste0(path, "AFHRI_B_genotypes_phenos.RDS"))
  protein <- readRDS(paste0(path, "AFHRI_B_proteomics_phenos_symbol.RDS"))
  protein_expression <- readRDS(paste0(path, "AFHRI_B_proteomics_QC_symbol.RDS"))
  transcript <- readRDS(paste0(path, "AFHRI_B_transcriptomics_phenos_symbol.RDS"))
  transcript_expression <- readRDS(paste0(path, "AFHRI_B_transcriptomics_QC_symbol.RDS"))
}

get.filename <- function(path) {
  filepath <- strsplit(path, "/")[[1]]
  return(tail(filepath, n=1))
}

rds2txt <- function(file) {
  locs <- readRDS(file)
  write.table(locs, paste0(path, ".txt"), sep="\t", row.names = T, quote=F, col.names = NA)
}

sortbynames <- function(names, target_df) {
  return(target_df[order(match(target_df$name, names)),])
}

get.message.usage <- function(){
  return(
    "Usage: Rscript runner.R [data directory] ... "
  )
}

get.metadata <- function(file){
  filename <- get.filename(file)
  name <- strsplit(filename, split = ".", fixed = T)[[1]][1]
  
  data <- tail(strsplit(file, split = "/", fixed = T)[[1]], n = 3)
  type <- data[2]
  model <- data[1]
  
  return(list(name, type, model))
}



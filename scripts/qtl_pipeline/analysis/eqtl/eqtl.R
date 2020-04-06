
eqtl.single <- function(expression_file, genotype, snplocs, genelocs, model, outdir, covariate_file = "no peer", of = F){
  library(MatrixEQTL)
  
  # Genotype file name
  SNP_file_name = genotype
  
  ## Load genotype data
  snps = SlicedData$new()
  
  snps$fileDelimiter = "\t"
  # the TAB character
  snps$fileOmitCharacters = "NA"
  # denote missing values;
  snps$fileSkipRows = 1
  # one row of column labels
  snps$fileSkipColumns = 1
  # one column of row labels
  snps$fileSliceSize = 2000
  # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name)
  
  
  ## Load gene expression data
  gene = SlicedData$new()
  
  gene$fileDelimiter = "\t"
  # the TAB character
  gene$fileOmitCharacters = "NA"
  # denote missing values;
  gene$fileSkipRows = 1
  # one row of column labels
  gene$fileSkipColumns = 1
  # one column of row labels
  gene$fileSliceSize = 2000
  # read file in slices of 2,000 rows
  gene$LoadFile(expression_file)
  
  if (is.na(genelocs)){
    mqtl(gene, snps, snplocs, genelocs, model, outdir, covariate_file, of)
  }else{
    eqtl(gene, snps, snplocs, genelocs, model, outdir, covariate_file, of)
  }
}

eqtls <- function(expression_file, genotype, snplocs, genelocs, model, outdir, covariate_dir, trans=F){
  # expression_file=paste(peer_path, norm_sd, type, file, sep="/")
  # genotype = geno
  # snplocs = path2snplocs
  # genelocs = path2genelocs
  # model = model
  # outdir = paste(get.path("results"), "imputed", "cis", model, type, norm_sd, cov_sd, sep = "/")
  # covariate_dir = paste(peer_path, norm_sd, type, cov_sd, sep="/")
  # trans=trans
  files <- list.files(path=covariate_dir, pattern=".*peer_factors_.*.txt", full.names=T, recursive=FALSE)

  library(MatrixEQTL)
  
  # Genotype file name
  SNP_file_name = genotype
  
  ## Load genotype data
  snps = SlicedData$new()
  
  snps$fileDelimiter = "\t"
  # the TAB character
  snps$fileOmitCharacters = "NA"
  # denote missing values;
  snps$fileSkipRows = 1
  # one row of column labels
  snps$fileSkipColumns = 1
  # one column of row labels
  snps$fileSliceSize = 2000
  # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name)
  
  
  ## Load gene expression data
  gene = SlicedData$new()
  
  gene$fileDelimiter = "\t"
  # the TAB character
  gene$fileOmitCharacters = "NA"
  # denote missing values;
  gene$fileSkipRows = 1
  # one row of column labels
  gene$fileSkipColumns = 1
  # one column of row labels
  gene$fileSliceSize = 2000
  # read file in slices of 2,000 rows
  gene$LoadFile(expression_file)
  
  if (is.na(genelocs)){
    mqtl(gene, snps, model, outdir, "no peer")
  }else{
    eqtl(gene, snps, snplocs, genelocs, model, outdir, "no peer", trans)
  }
  
  #eqtl(gene, snps, snplocs, genelocs, model, outdir, paste0(covariate_dir, "peer_factors_nk01.txt"))
  
  if(is.na(genelocs)){
    lapply(files, function(x){
      try(mqtl(gene, snps, model, outdir, x))
    })
    
  }else{
    lapply(files, function(x){
      tryCatch(eqtl(gene, snps, snplocs, genelocs, model, outdir, x, trans),error=function(e) NULL)
    })
    
  }
  # if(tables=F){
  #   source("analysis/peer_tables/me_analyze.R")
  #   source("helper/helper.R")
  #   dir <- paste(get.path("results"), "current", "cis", model, type, norm_sd, cov_sd, sep="/")
  #   outdir <- paste(get.path("results"), "current", "peer", model, type, norm_sd, cov_sd, sep="/")
  #   
  #   peer.write.table(dir, outdir)
  # }
}

eqtl <- function(gene, genotype, snplocs, genelocs, model, outdir, covariate_file = character(), trans = F, of = F) {
  tryCatch({
  library(MatrixEQTL)
  source("helper/helper.R")

  
  
  if(covariate_file == "no peer") {
    out <- "nk00"
    covariate_file <- character()
  } else {
    out <- get.filename(covariate_file)
    out <- strsplit(out, ".", fixed = T)
    out <- strsplit(out[[1]][1], "_", fixed = T)[[1]][3]
  }
  
  outfile <- paste0(outdir , "/eqtl_", out, ".RDS")
  
  if(of) {
    outfile <- outdir
  } else {
    dir.create(outdir, showWarnings=F, recursive=T)
  }
  
  # Linear model to use
  '%nin%' <- Negate('%in%')
  if(model == "linear"){
    useModel = modelLINEAR  
  }
  if(model == "anova"){
    useModel = modelANOVA
  }
  if(model %nin% c("linear", "anova")){
    stop(paste0("'",model, "'", " is no valid model."))
  }
  
  # Covariates file name
  covariates_file_name = covariate_file
  
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = 1
  
  # pvOutputThreshold_tra = 1e-2; # -> rds 02
  pvOutputThreshold_tra = 0 # -> rds01
  if (trans) {pvOutputThreshold_tra = 1e-2}
  
  # Error covariance matrix
  errorCovariance = numeric()
  
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6
  
  
  # You can also have a look at:
  # library(MatrixEQTL)
  # ?MatrixEQTL_cis_code
  
  
  ## Settings
  

  
  snps_location_file_name = snplocs
  
  
  gene_location_file_name = genelocs
  
  
  # Output file name
  output_file_name_cis = paste0(outdir , "/eqtl_", out, "_cis.txt") # tempfile() tempfile() #
  
  output_file_name_tra = paste0(outdir , "/eqtl_", out, "_trans.txt") # tempfile()
  
  
  ## Genotype data

  snps <- genotype
  
  
  
  ## Load covariates
  cvrt = SlicedData$new()
  
  cvrt$fileDelimiter = "\t"
  # the TAB character
  cvrt$fileOmitCharacters = "NA"
  # denote missing values;
  cvrt$fileSkipRows = 1
  # one row of column labels
  cvrt$fileSkipColumns = 1
  # one column of row labels
  if (length(covariates_file_name) > 0) {
    cvrt$LoadFile(covariates_file_name)
    
  }
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name,
                       header = TRUE,
                       stringsAsFactors = FALSE)
  
  genepos = read.table(gene_location_file_name,
                       header = TRUE,
                       stringsAsFactors = FALSE)
  
  if(trans){
    me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = output_file_name_tra,
      pvOutputThreshold = pvOutputThreshold_tra,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThreshold_cis,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )
    saveRDS(me, outfile)
  }else{
    me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = output_file_name_tra,
      pvOutputThreshold = pvOutputThreshold_tra,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThreshold_cis,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )
    saveRDS(me, outfile)
  }

  
  
  # unlink(output_file_name_tra)
  # 
  unlink(output_file_name_cis)
  
  },
  warning=function(e){
    message(e)
  })
}

mqtl <- function(gene, genotype, model, outdir, covariate_file = character(), trans = F, of = F) {
  try({
    library(MatrixEQTL)
    source("helper/helper.R")
    
    # if(covariate_file == "no peer") {
    #   out <- "nk00"
    #   covariate_file <- character()
    # } else {
    #   out <- get.filename(covariate_file)
    #   out <- strsplit(out, ".", fixed = T)
    #   out <- strsplit(out[[1]][1], "_", fixed = T)[[1]][3]
    # }
    # 
    # outfile <- paste0(outdir , "/eqtl_", out, ".RDS")
    # 
    # if(of) {
    #   outfile <- outdir
    # } else {
    #   dir.create(outdir, showWarnings=F, recursive=T)
    # }
    
    #covariate_file <- "no peer"
    
    if(covariate_file == "no peer") {
      out <- "nk00"
      covariate_file <- character()
    } else {
      out <- get.filename(covariate_file)
      out <- strsplit(out, ".", fixed = T)
      out <- strsplit(out[[1]][1], "_", fixed = T)[[1]][3]
    }
    
    outfile <- paste0(outdir , "/eqtl_", out, ".RDS")
    
    if(of) {
      outfile <- outdir
    } else {
      dir.create(outdir, showWarnings=F, recursive=T)
    }
    
    # Linear model to use
    '%nin%' <- Negate('%in%')
    if(model == "linear"){
      useModel = modelLINEAR  
    }
    if(model == "anova"){
      useModel = modelANOVA
    }
    if(model %nin% c("linear", "anova")){
      stop(paste0("'",model, "'", " is no valid model."))
    }
    
    # Covariates file name
    covariates_file_name = covariate_file
    
    # pvalue output threshold
    pvOutputThreshold = 0.05
    
    # Error covariance matrix
    errorCovariance = numeric()
    
    
    # Distance for local gene-SNP pairs
    cisDist = 0
    
    
    # You can also have a look at:
    # library(MatrixEQTL)
    # ?MatrixEQTL_cis_code
    
    
    ## Settings
    
    
    
    # snps_location_file_name = snplocs
    # 
    # 
    # gene_location_file_name = genelocs
    
    
    # Output file name
    output_file_name = paste0(outdir , "/eqtl_", out, ".txt") # tempfile()
    
    
    ## Genotype data
    
    snps <- genotype
    
    
    
    
    
    ## Load covariates
    cvrt = SlicedData$new()
    
    cvrt$fileDelimiter = "\t"
    # the TAB character
    cvrt$fileOmitCharacters = "NA"
    # denote missing values;
    cvrt$fileSkipRows = 1
    # one row of column labels
    cvrt$fileSkipColumns = 1
    # one column of row labels
    if (length(covariates_file_name) > 0) {
      cvrt$LoadFile(covariates_file_name)
      
    }
    
    ## Run the analysis
    try({
    me = Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      cvrt = cvrt, 
      output_file_name = output_file_name, 
      pvOutputThreshold = pvOutputThreshold, 
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = TRUE,
      noFDRsaveMemory = FALSE
    )
    })
    
    saveRDS(me, outfile)
    #unlink(output_file_name)
  # },
  # warning=function(e){
  #   message(e)
  })
}

# # cis-eQTLs
# show(me$cis$eqtls)

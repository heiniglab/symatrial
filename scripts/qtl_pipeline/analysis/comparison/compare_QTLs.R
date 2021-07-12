## here we do a thorrow comparison of our eqtl results and the gtex data set

library(ggplot2)
library(ggpubr)
library(data.table)

setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")
source("helper/helper.R")
local=F

merged.all = F
threshold <- 1e-5
outfile <- NULL

get.duplicates <- function(df, cols){
  df2 <- data.frame(table(df[, cols]))
  df2 <- df2[df2$Freq>1, ]
  colnames(df2) <- c(cols, "Freq")
  df <- merge(df,
              df2[cols],
              all.y=T)
}

do.stuff <- function(){
  
  ## start the analysis
  outdir = prefix
  dir.create(outdir, recursive=T)
  ## dir.create(paste(prefix, "_compare_to_gtex/", secondary.name, "/", sep=""))
  
  ## load the reference eQTL data set
  primary = read.table(file=primary.file, h=T, stringsAsFactors=F)
  
  ## export all significant QTLs matched with secondary dataset
  if(merged.all){
    if(primary.name=="GTEx_v7"){
      anno <- readRDS(paste0(get.path("locations", local),
                             "GTEx_ensembl_ids_AFHRIB_QTLgenes.RDS"))
      primary$ensembl_gene_id <- gsub("\\..*", "", primary[, primary.gene])
      primary <- merge(primary, anno[, c("ensembl_gene_id", "symbol")],
                       all.x=T, sort=F)
      primary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
        primary[, c(primary.variant, "symbol", primary.pvalue, primary.effect)]
    } else {
      primary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
        primary[, c(primary.variant, primary.gene, primary.pvalue, primary.effect)]
    }
    genes = unique(primary$s.gene)
    if(secondary.name=="GTEx_v7"){
      anno <- readRDS(paste0(get.path("locations", local),
                             "GTEx_ensembl_ids_AFHRIB_QTLgenes.RDS"))
      genes <- anno[anno$symbol %in% genes, "gene_id.v7"]
      genes <- genes[!is.na(genes)]
    }
    gfile = paste(prefix, primary.name, "_egenes.txt", sep="")
    cat(genes, sep="\n", file=gfile)
    secondary = fread(paste("fgrep -f ", gfile, " -w ", secondary.file, sep=""), header=F)
    setnames(secondary, colnames(read.csv(secondary.file, sep="\t", nrows=5)))
    secondary = as.data.frame(secondary)
    if(secondary.name=="GTEx_v7"){
      secondary <- merge(secondary, anno,
                         by.x = "gene_id", by.y = "gene_id.v7",
                         all.x=T, sort=F)
      secondary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
        secondary[, c(secondary.variant, "symbol", secondary.pvalue, secondary.effect)]
    } else {
      secondary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
        secondary[, c(secondary.variant, secondary.gene, secondary.pvalue, secondary.effect)]
    }
    
    merged = merge(primary, secondary,
                   by=c("s.variant", "s.gene"),
                   suffixes=c("_primary", "_secondary"))
    
    ## make sure that the allele encoding is the same
    same = toupper(merged[,"Allele1_primary"]) == toupper(merged[,"Allele1_secondary"])
    if(primary.name=="GTEx_v7" | secondary.name=="GTEx_v7"){
      same = toupper(merged[,"Allele1_primary"]) == toupper(merged[,"Allele2_secondary"])
    }
    
    ## switch the sign of the effect size for the ones that are not matching
    merged[!same,"s.beta_secondary"] = -merged[!same,"s.beta_secondary"]
    
    ## also note in the table whether the effect was switched
    merged = data.frame(merged, switched=!same)
    
    ## first compare the number of genes with eQTL
    secondary.significant = merged[,"s.pvalue_secondary"] < threshold
    
    ## also check on the gene level if there is a significant eQTL
    secondary.genes = unique(secondary[secondary[,"s.pvalue"] < threshold,"s.gene"])
    secondary.gene.level = merged[,"s.gene"] %in% secondary.genes
    
    ## check for concordance of the effect direction
    concordant = sign(merged[,"s.beta_primary"]) == sign(merged[,"s.beta_secondary"])
    
    pal = c("lightgrey", "darkgrey", "red", "black")
    col = pal[as.numeric(factor(interaction(concordant, secondary.significant),
                                levels=c("FALSE.FALSE", "TRUE.FALSE",
                                         "FALSE.TRUE", "TRUE.TRUE")))]
    merged$col <- col
    
    if(is.null(outfile)){
      saveRDS(merged,
              file=paste0(outdir,
                          primary.name, "_to_", secondary.name,
                          "_merged_all_matched.RDS"))
    }else{
      saveRDS(merged,
              file=paste0(outfile,
                          "_merged_all_matched.RDS"))
    }
    
  }
  
  ## extract the corresponding subset from secondary
  ## since this is a big file we filter for the SNPs that are the most
  ## significant per gene in our data set
  best = unlist(tapply(1:nrow(primary), primary[, primary.gene], function(idx) idx[which.min(primary[, c(primary.pvalue)][idx])]))
  best.primary = primary[best,]
  if(primary.name=="plasma_pQTL"){
    best.primary[, primary.pvalue] <- exp(best.primary[, primary.pvalue])
    if(secondary.name=="GTEx_v7"){
      best.primary[, primary.variant] <- gsub("(:[^:]+):.*", "\\1", best.primary[, primary.variant])
    }
  }
  if(primary.name=="GTEx_v7"){
    anno <- readRDS(paste0(get.path("locations", local),
                           "GTEx_ensembl_ids_AFHRIB_QTLgenes.RDS"))
    best.primary$ensembl_gene_id <- gsub("\\..*", "", best.primary[, primary.gene])
    best.primary <- merge(best.primary, anno[, c("ensembl_gene_id", "symbol")],
                          all.x=T, sort=F)
    if(secondary.name=="plasma_pQTL"){
      best.primary[, primary.variant] <- gsub("(_[^_]+)_.*", "\\1", best.primary[, primary.variant])
      best.primary[, primary.variant] <- gsub("_", ":", best.primary[, primary.variant])
    }
    best.primary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
      best.primary[, c(primary.variant, "symbol", primary.pvalue, primary.effect)]
  } else {
    best.primary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
      best.primary[, c(primary.variant, primary.gene, primary.pvalue, primary.effect)]
  }
  genes = unique(best.primary$s.gene)
  
  if(secondary.name=="GTEx_v7"){
    anno <- readRDS(paste0(get.path("locations", local),
                           "GTEx_ensembl_ids_AFHRIB_QTLgenes.RDS"))
    genes <- anno[anno$symbol %in% genes, "gene_id.v7"]
  }
  gfile = paste(prefix, primary.name, "_egenes.txt", sep="")
  cat(genes, sep="\n", file=gfile)
  
  secondary = fread(paste("fgrep -f ", gfile, " -w ", secondary.file, sep=""), header=F)
  setnames(secondary, colnames(read.csv(secondary.file, sep="\t", nrows=5)))
  secondary = as.data.frame(secondary)
  if(secondary.name=="plasma_pQTL"){
    secondary[, secondary.pvalue] <- exp(secondary[, secondary.pvalue])
    if(primary.name=="GTEx_v7"){
      secondary[, secondary.variant] <- gsub("(:[^:]+):.*", "\\1", secondary[, secondary.variant])
    }
  }
  if(secondary.name=="GTEx_v7"){
    if(primary.name=="plasma_pQTL"){
      secondary[, secondary.variant] <- gsub("(_[^_]+)_.*", "\\1", secondary[, secondary.variant])
      secondary[, secondary.variant] <- gsub("_", ":", secondary[, secondary.variant])
    }
    secondary <- merge(secondary, anno,
                       by.x = "gene_id", by.y = "gene_id.v7",
                       all.x=T, sort=F)
    secondary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
      secondary[, c(secondary.variant, "symbol", secondary.pvalue, secondary.effect)]
  } else {
    secondary[, c("s.variant", "s.gene", "s.pvalue", "s.beta")] <-
      secondary[, c(secondary.variant, secondary.gene, secondary.pvalue, secondary.effect)]
  }
  
  ## when a map file is given, the allele encoding is not contained in the eQTL
  ## results and we have to extract it from the map file
  addAlleles <- function(data, map, prefix, name, variant) {
    ## extract the snps
    sfile = paste(prefix, name, "_snps.txt", sep="")
    snps = unique(data[, variant])
    cat(snps, file=sfile, sep="\n")
    
    ## allele encoding
    map = fread(paste("fgrep -f ", sfile, " -w ", map))
    n = colnames(read.table(map, header=T, nrows=5))
    n = gsub("variantid", "ID", gsub("_b37", "", tolower(n)))
    n = gsub("alt", "Allele2", gsub("ref", "Allele1", n))
    setnames(map, n)
    map = as.data.frame(map)
    
    file.remove(sfile)
    
    cols = c("Allele1", "Allele2")
    if ("rs_id_dbsnp142_chg37p13" %in% n) {
      cols = c(cols, "rs_id_dbsnp142_chg37p13")
    }
    secondary = cbind(data, map[match(data$SNP, map$ID),cols])
    if ("rs_id_dbsnp142_chg37p13" %in% n) {
      colnames(data) = gsub("rs_id_dbsnp142_chg37p13", "SNP", gsub("SNP", "gtex_id", colnames(data)))
    }
    return(data)
  }
  
  if (!is.null(secondary.map)) {
    secondary = addAlleles(secondary, secondary.map, prefix, secondary.name, secondary.variant)
  }
  
  if (!is.null(primary.map)) {
    best.primary = addAlleles(best.primary, primary.map, prefix, primary.name, primary.variant)
  }
  
  file.remove(gfile)
  
  
  ## match eQTLs
  merged = merge(best.primary, secondary,
                 by=c("s.variant", "s.gene"),
                 suffixes=c("_primary", "_secondary"))
  
  ## make sure that the allele encoding is the same
  same = toupper(merged[,"Allele1_primary"]) == toupper(merged[,"Allele1_secondary"])
  if(primary.name=="GTEx_v7" | secondary.name=="GTEx_v7"){
    same = toupper(merged[,"Allele1_primary"]) == toupper(merged[,"Allele2_secondary"])
  }
  
  ## switch the sign of the effect size for the ones that are not matching
  merged[!same,"s.beta_secondary"] = -merged[!same,"s.beta_secondary"]
  
  ## also note in the table whether the effect was switched
  merged = data.frame(merged, switched=!same)
  
  ## first compare the number of genes with eQTL
  secondary.significant = merged[,"s.pvalue_secondary"] < threshold
  
  ## also check on the gene level if there is a significant eQTL
  secondary.genes = unique(secondary[secondary[,"s.pvalue"] < threshold,"s.gene"])
  secondary.gene.level = merged[,"s.gene"] %in% secondary.genes
  
  ## check for concordance of the effect direction
  concordant = sign(merged[,"s.beta_primary"]) == sign(merged[,"s.beta_secondary"])
  
  pal = c("lightgrey", "darkgrey", "red", "black")
  col = pal[as.numeric(factor(interaction(concordant, secondary.significant), levels=c("FALSE.FALSE", "TRUE.FALSE", "FALSE.TRUE", "TRUE.TRUE")))]
  merged$col <- col
  
  if(is.null(outfile)){
    saveRDS(merged,
            file=paste0(outdir,
                        primary.name, "_to_", secondary.name,
                        "_merged.RDS"))
  }else{
    saveRDS(merged,
            file=paste0(outfile,
                        "_merged.RDS"))
  }
  
  # merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
  # outdir <- prefix
  
  if(is.null(outfile)){
    pdf(file=paste0(outdir, primary.name, "_to_", secondary.name, "_scatterplot_effect_sizes.pdf"))
    plot(merged[order(merged$s.pvalue_secondary, decreasing = T),
                c("s.beta_primary", "s.beta_secondary")],
         col=col,
         xlab=paste("Effect size", primary.name),
         ylab=paste("Effect size", secondary.name))
    abline(a=0, b=1)
    legend("topleft", pch=1, col=pal, c("discordant | non significant", "concordant | nonsignificant", "discordant | significant", "concordant | significant"))
    dev.off()
    merged$col <- factor(merged$col,
                         levels = c("lightgrey", "darkgrey", "red", "black"),
                         ordered = T)
    gg <- ggplot(merged[order(merged$s.pvalue_secondary, decreasing = T), ],
                 aes(x=s.beta_primary, y=s.beta_secondary, col=col)) +
      geom_abline(intercept = 0, slope = 1, col="grey") +
      geom_point() +
      theme_bw() +
      theme(aspect.ratio = 1) +
      xlim(c(min(c(merged$s.beta_primary, merged$s.beta_secondary)),
             max(c(merged$s.beta_primary, merged$s.beta_secondary)))) +
      ylim(c(min(c(merged$s.beta_primary, merged$s.beta_secondary)),
             max(c(merged$s.beta_primary, merged$s.beta_secondary)))) +
      labs(x=paste("Effect size", primary.name),
           y=paste("Effect size", secondary.name)) +
      scale_color_manual(values = pal,
                         labels = c("discordant | non significant",
                                    "concordant | nonsignificant",
                                    "discordant | significant",
                                    "concordant | significant"),
                         "", drop = F)
    pdf(file=paste0(prefix, primary.name, "_to_", secondary.name, "_scatterplot_effect_sizes_ggplot.pdf"))
    print(gg)
    dev.off()
  }
  
  
  #return(merged)
}

merged.all <- T

## First test: our eQTL vs GTEx -------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "eQTL"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  do.stuff()
}

if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "eQTL"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".1e-5")
  do.stuff()
  outfile <- NULL
}

## GTEx to eQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "eQTL"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "eQTL"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".1e-5_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

### alternative covariate sets ----
#### our eQTL all vs GTEx -------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "eQTL_cov_all"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  do.stuff()
}

if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "eQTL_cov_all"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".1e-5")
  do.stuff()
  outfile <- NULL
}

#### GTEx to eQTL all -----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "eQTL_cov_all"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "eQTL_cov_all"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".1e-5_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

#### our eQTL pop vs GTEx -------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "eQTL_cov_pop"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  do.stuff()
}

if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "eQTL_cov_pop"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".1e-5")
  do.stuff()
  outfile <- NULL
}

#### GTEx to eQTL pop -----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "eQTL_cov_pop"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/alternative/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "eQTL_cov_pop"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".1e-5_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

## our pQTL vs GTEx -------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "pQTL"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  do.stuff()
}

if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "pQTL"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "gtex.variant_id"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".1e-5")
  do.stuff()
  outfile <- NULL
}

## GTEx to pQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/GTEx/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "gtex.variant_id"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".1e-5_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

# # None of our pQTLs hits are measured in plasma -> check!
# ## our pQTL vs plasma_pQTL -------
# if(F){
#   prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL/"))
#   
#   primary.name = "pQTL"
#   secondary.name = "plasma_pQTL"
#   
#   primary.file = paste0(get.path("results", local), "imputed/cis/final/",
#                         primary.name, "_right_atrial_appendage_allpairs.significant.txt")
#   primary.map = NULL
#   primary.gene = "gene"
#   primary.variant = "variant_id"
#   primary.pvalue = "pvalue"
#   primary.effect = "beta"
#   
#   secondary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt")
#   secondary.map = NULL
#   secondary.gene = "gene"
#   secondary.variant = "variant_id"
#   secondary.pvalue = "log.P."
#   secondary.effect = "Effect"
#   
#   do.stuff()
# }



## plasma_pQTL to pQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL/"))
  
  primary.name = "plasma_pQTL"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "variant_id"
  primary.pvalue = "log.P."
  primary.effect = "Effect"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}


## plasma pQTL to GTEx ------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL_GTEx/"))
  
  primary.name = "plasma_pQTL"
  secondary.name = "GTEx_v7"
  
  primary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "variant_id"
  primary.pvalue = "log.P."
  primary.effect = "Effect"
  
  secondary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
  secondary.map = NULL
  secondary.gene = "gene_id"
  secondary.variant = "variant_id"
  secondary.pvalue = "pval_nominal"
  secondary.effect = "slope"
  
  do.stuff()
}

## GTEx to plasma pQTL------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL_GTEx/"))
  
  primary.name = "GTEx_v7"
  secondary.name = "plasma_pQTL"
  
  primary.file = paste0(get.path("gtex", local), "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")
  primary.map = NULL
  primary.gene = "gene_id"
  primary.variant = "variant_id"
  primary.pvalue = "pval_nominal"
  primary.effect = "slope"
  
  secondary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "variant_id"
  secondary.pvalue = "log.P."
  secondary.effect = "Effect"
  
  do.stuff()
}


## eQTL to pQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "eQTL"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "eQTL"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}


## pQTL to eQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "pQTL"
  secondary.name = "eQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "pQTL"
  secondary.name = "eQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

## pQTL to res_pQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "pQTL"
  secondary.name = "res_pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "pQTL"
  secondary.name = "res_pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}


## res_pQTL to pQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "res_pQTL"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "res_pQTL"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}


## ratios to res_pQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "ratios"
  secondary.name = "res_pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "ratios"
  secondary.name = "res_pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

## res_pQTL to ratios-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "res_pQTL"
  secondary.name = "ratios"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "res_pQTL"
  secondary.name = "ratios"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}


## pQTL to ratios -----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "pQTL"
  secondary.name = "ratios"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "pQTL"
  secondary.name = "ratios"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

## ratios to pQTL-----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "ratios"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/AFHRI-B/"))
  
  primary.name = "ratios"
  secondary.name = "pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.FDR.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "snpid"
  primary.pvalue = "FDR"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "snpid"
  secondary.pvalue = "FDR"
  secondary.effect = "beta"
  
  outfile <- paste0(prefix, primary.name, ".FDR_to_", secondary.name, ".FDR")
  threshold <- 0.05
  do.stuff()
  outfile <- NULL
  threshold <- 1e-5
}

## plasma_pQTL to ratios -----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL/"))
  
  primary.name = "plasma_pQTL"
  secondary.name = "ratios"
  
  primary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "variant_id"
  primary.pvalue = "log.P."
  primary.effect = "Effect"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}


## plasma_pQTL to eQTL ------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL/"))
  
  primary.name = "plasma_pQTL"
  secondary.name = "eQTL"
  
  primary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "variant_id"
  primary.pvalue = "log.P."
  primary.effect = "Effect"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

## eQTL to plasma_pQTL -------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL/"))
  
  primary.name = "eQTL"
  secondary.name = "plasma_pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "variant_id"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "variant_id"
  secondary.pvalue = "log.P."
  secondary.effect = "Effect"
  
  do.stuff()
}


## plasma_pQTL to res pQTL -----------
if (T) {
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL/"))
  
  primary.name = "plasma_pQTL"
  secondary.name = "res_pQTL"
  
  primary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "variant_id"
  primary.pvalue = "log.P."
  primary.effect = "Effect"
  
  secondary.file = paste0(get.path("results", local), "imputed/cis/final/",
                          secondary.name, "_right_atrial_appendage_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "variant_id"
  secondary.pvalue = "pvalue"
  secondary.effect = "beta"
  
  do.stuff()
  #merged <- readRDS(file=paste0(prefix, primary.name, "_to_", secondary.name, "_merged.RDS"))
}

## res eQTL vs plasma_pQTL -------
if(T){
  prefix = paste0(get.path("results", local), paste("imputed/cis/final/comparisons/plasma_pQTL/"))
  
  primary.name = "res_eQTL"
  secondary.name = "plasma_pQTL"
  
  primary.file = paste0(get.path("results", local), "imputed/cis/final/",
                        primary.name, "_right_atrial_appendage_allpairs.significant.txt")
  primary.map = NULL
  primary.gene = "gene"
  primary.variant = "variant_id"
  primary.pvalue = "pvalue"
  primary.effect = "beta"
  
  secondary.file = paste0(get.path("plasma", local), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt")
  secondary.map = NULL
  secondary.gene = "gene"
  secondary.variant = "variant_id"
  secondary.pvalue = "log.P."
  secondary.effect = "Effect"
  
  do.stuff()
}

q(save="no")

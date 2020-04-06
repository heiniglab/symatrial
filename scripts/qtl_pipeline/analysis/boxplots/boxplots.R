##
#   makes boxplots for given SNP-gene pair
#
#
#
#
#

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

get.genos <- function(snp, AFHRI_B_imp=NULL, local=F){
  if(is.null(AFHRI_B_imp)){
    load(paste0(get.path("genotype", local), "AFHRI_B_imputed_gwaa_data.RData")) 
  }
  map.file <- paste0(get.path("genotype", local), "AFHRI_B_imputed_recode_alleles.map")
  snp.file <- tempfile()
  cat(snp, file=snp.file, sep="\n")
  write.table(snp, file=snp.file, row.names = F, col.names = F,
              sep = "\t", quote = F)
  map <- read.table(pipe(paste("cat ", map.file, " | fgrep -w -f ", snp.file)),
                    stringsAsFactors=F)
  rownames(map) <- map$V1 
  
  df <- data.frame(id=idnames(AFHRI_B_imp),
                   sex=phdata(AFHRI_B_imp)[, "sex"],
                   sub("/", "", as.character(AFHRI_B_imp[, snp])),
                   stringsAsFactors = F)
  for(id in snp){
    df[, gsub("\\:", "\\.", id)] <- gsub("A", map[id, "V2"], df[, gsub("\\:", "\\.", id)])
    df[, gsub("\\:", "\\.", id)] <- gsub("B", map[id, "V3"], df[, gsub("\\:", "\\.", id)])
  }
  return(df)
}

get.boxplot.data <- function(snp, gene, local=F, meta=F, imputed=T){
  source("~/work/symAtrial_QTL/scripts/helper/helper.R")
  geno_path <- get.path("genotype", local)
  
  library(GenABEL)
  library(reshape2)
  
  if(imputed){
    load(paste0(geno_path, "AFHRI_B_imputed_gwaa_data.RData"))
    gwaa <- AFHRI_B_imp
    map.file <- paste0(geno_path, "AFHRI_B_imputed_recode_alleles.map")
    snp.file <- tempfile()
    cat(snp, file=snp.file, sep="\n")
    write.table(snp, file=snp.file, row.names = F, col.names = F,
                sep = "\t", quote = F)
    map <- read.table(pipe(paste("cat ", map.file, " | fgrep -w -f ", snp.file)),
                      stringsAsFactors=F)
    rownames(map) <- map$V1
  }else{
    load(paste0(geno_path, "AFHRI_B_gwaa_data.RData"))
    gwaa <- AFHRI_B
  }

  genes <- colnames(phdata(gwaa))[grep(paste(gene, collapse = "|"), colnames(phdata(gwaa)))]
  df <- data.frame(id=idnames(gwaa),
                   sex=phdata(gwaa)[, "sex"],
                   sub("/", "", as.character(gwaa[, snp])),
                   phdata(gwaa)[, genes],
                   stringsAsFactors = F)
  if(imputed){
    for(id in snp){
      df[, gsub("\\:", "\\.", id)] <- gsub("A", map[id, "V2"], df[, gsub("\\:", "\\.", id)])
      df[, gsub("\\:", "\\.", id)] <- gsub("B", map[id, "V3"], df[, gsub("\\:", "\\.", id)])
    }
  }
  return(df)
}


eqtl.pqtl.boxplot <- function(data=NULL, SNP, Gene, df){

  library(ggplot2)
  library(RColorBrewer)
  library(ggpubr)
  
  col.paired <- brewer.pal(n = 11, "Paired")
  col.set <- col.paired[c(2,8,9)]
  
  df1 <- df[, c("id", SNP, paste0(Gene, c("_trans", "_prot", "_p.res", "_t.res")))]
  colnames(df1) <- c("externID", "gt", "trans", "prot", "p.res", "t.res")
  df1$gt <- as.factor(df1$gt)
  df2 <- df1[!is.na(df1$gt), ]
  
  g1 <- ggplot(df2, aes(x=gt, y=trans, col=gt, fill=gt)) +
    geom_boxplot(alpha=0.2) +
    geom_point(position = position_jitter(width = 0.2)) +
    scale_color_manual("genotype", values = col.set) +
    scale_fill_manual("genotype", values = col.set) +
    scale_x_discrete(labels = paste0(levels(df2$gt), "\n(N=", table(df2$gt), ")")) +
    labs(list(title=paste0(Gene, " ~ ", SNP
                           # , "\np=", format(anova(lm(data=df, trans ~ SNP))[1,5], digits = 3)
                           ),
              x=SNP, y=paste0(Gene, " transcript"))) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  g1
  
  g2 <- ggplot(df2, aes(x=gt, y=prot, col=gt, fill=gt)) +
    geom_boxplot(alpha=0.2) +
    geom_point(position = position_jitter(width = 0.2)) +
    scale_color_manual("genotype", values = col.set) +
    scale_fill_manual("genotype", values = col.set) +
    scale_x_discrete(labels = paste0(levels(df2$gt), "\n(N=", table(df2$gt), ")")) +
    labs(list(title=paste0(Gene, " ~ ", SNP
                           # , "\np=", format(anova(lm(data=df, prot ~ SNP))[1,5], digits = 3)
                          ),
              x=SNP, y=paste0(Gene, " protein"))) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
              axis.title=element_text(size=16,face="bold"))
  g2
  
  g3 <- ggplot(df2, aes(x=trans, y=prot, col=gt)) +
    geom_point() +
    scale_color_manual("genotype", values = col.set) +
    #scale_color_manual(breaks=levels(df1$gt), values=gg_color_hue(nlevels(df1$gt))) +
    labs(list(title=paste0("Coexpression ", Gene), x=paste0(Gene, " transcript"),
              y=paste0(Gene, " protein"),
              col="genotype")) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  g3
  
  # grid.arrange(g1, g2, g3, ncol=3, nrow=1, top = paste0(gene," ", snp, "\npval = ", pvalue, ", FDR = ", fdr), newpage = T)
  g <- ggarrange(g1, g2, g3,
                 ncol=3, nrow=1,
                 heights = c(1,1,1))#, top = paste0("pval eQTL = ", formatC(pvalue, format = "e", digits = 2), 
                                                        # ", FDR eQTL = ", formatC(fdr, format = "e", digits = 2), 
                                                        # "\npval pQTL = ", formatC(pvalue2, format = "e", digits = 2),
                                                        # ", FDR pQTL = ", formatC(fdr2, format = "e", digits = 2)))
  return(g)
}


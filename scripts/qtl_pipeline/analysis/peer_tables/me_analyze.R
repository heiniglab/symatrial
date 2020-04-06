##  Returns the number of genes and gene-snp-pairs for different FDR /pval cutoffs



analyze <- function(me, trans=F){
  if(trans){
    fdr005 <- me$all$eqtls[me$all$eqtls$FDR < 0.05,]
    fdr02 <- me$all$eqtls[me$all$eqtls$FDR < 0.2,]
    p5x1em5 <- me$all$eqtls[me$all$eqtls$pvalue < 1e-5,]
  }else{
    fdr005 <- me$cis$eqtls[me$cis$eqtls$FDR < 0.05,]
    fdr02 <- me$cis$eqtls[me$cis$eqtls$FDR < 0.2,]
    p5x1em5 <- me$cis$eqtls[me$cis$eqtls$pvalue < 1e-5,]
  }
  fdr005_nr <- nrow(fdr005)
  fdr005_genes <- length(unique(fdr005$gene))
  fdr02_nr <- nrow(fdr02)
  fdr02_genes <- length(unique(fdr02$gene))
  p5x1em5_nr <- nrow(p5x1em5)
  p5x1em5_genes <- length(unique(p5x1em5$gene)) 
  return(c(fdr005_nr, fdr005_genes, fdr02_nr, fdr02_genes, p5x1em5_nr, p5x1em5_genes))
}
  

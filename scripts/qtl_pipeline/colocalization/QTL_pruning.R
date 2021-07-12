# QTL clumping using Plink -----------------------------------------------------

# https://zzz.bwh.harvard.edu/plink/clump.shtml

setwd("~/work/symAtrial_QTL/scripts/revision")
source("../helper/helper.R")
local=F
get.duplicates <- function(df, cols){
  df2 <- data.frame(table(df[, cols]))
  df2 <- df2[df2$Freq>1, ]
  colnames(df2) <- c(cols, "Freq")
  df <- merge(df,
              df2[cols],
              all.y=T)
}

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(qvalue)
library(coloc)
library(readxl)
library(reshape2)

HMGU.blue <- "#003E6E"
mygray <- "#C6DDEA"
col.paired <- brewer.pal(n = 11, "Paired")
col.set <- col.paired[c(2,9,8)]
col.set2 <- brewer.pal(n=11, "BrBG")[7:10]

## Data ------------------------------------------------------------------------

# shared <- readRDS(paste0(get.path("results", local),
#                             "imputed/cis/",
#                             "QTL_res_shared_QTL_genes_group_anno.RDS"))
# eqtl <- readRDS(paste0(get.path("results", local),
#                             "imputed/cis/final/",
#                             "eQTL_right_atrial_appendage_allpairs.RDS"))
# pqtl <- readRDS(paste0(get.path("results", local),
#                             "imputed/cis/final/",
#                             "pQTL_right_atrial_appendage_allpairs.RDS"))
# e.genes <- unique(eqtl[eqtl$FDR < 0.05, "gene"])
# p.genes <- unique(pqtl[pqtl$FDR < 0.05, "gene"])
# s.genes <- unique(shared$gene)

## eQTL/pQTL integrated clumping -----------------------------------------------

f.clump <- paste0(get.path("results", local),
                            "imputed/cis/eQTL_pQTL_overlap_clump.RDS")
if(!file.exists(f.clump)){
  
  eqtls <- readRDS(file=paste0(get.path("results", local),
                               "/imputed/cis/final/",
                               "eQTL_right_atrial_appendage_allpairs.RDS"))
  pqtls <- readRDS(file=paste0(get.path("results", local),
                               "/imputed/cis/final/",
                               "pQTL_right_atrial_appendage_allpairs.RDS"))
  
  QTLs <- merge(eqtls, pqtls,
                by = c("gene", "snpid",
                       "chr", "variant_pos", "Allele1", "Allele2",
                       "variant_id", "rs_id", "gtex.variant_id",
                       "gtex.match_alleles", "gene_id"),
                suffixes = c(".eqtl", ".pqtl"))
  colnames(QTLs)
  
  s.genes <- unique(QTLs[QTLs$FDR.eqtl<0.05 | QTLs$FDR.pqtl<0.05, "gene"])
  clump <- list()
  tmp <- tempdir()
  clump[["FDR"]] <- list()
  clump[["P"]] <- list()
  dir.create(tmp)
  print(paste0("Temporary files: ", tmp))
  
  for(gene in s.genes){
    # gene <- "MYOZ1"
    
    # FDR ----
    clump[["FDR"]][[gene]] <- list()
    stats <- QTLs[QTLs$gene == gene, ]
    stats1 <- stats[, c("snpid", "FDR.eqtl", "pvalue.eqtl", "FDR.eqtl")]
    colnames(stats1) <- c("SNP", "P", "pvalue.eqtl", "FDR.eqtl")
    stats2 <- stats[, c("snpid", "FDR.pqtl", "pvalue.pqtl", "FDR.pqtl")]
    colnames(stats2) <- c("SNP", "P", "pvalue.pqtl", "FDR.pqtl")
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_eqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(stats2,
                file = paste0(tmp, "/", gene, "_pqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_eqtl.assoc", ",", gene, "_pqtl.assoc",
                  " --clump-p1 0.05", " --clump-p2 0.8",
                  #" --clump-best --clump-annotate pval.eqtl,FDR.eqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_epqtl"))
    
    clump[["FDR"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                         "_epqtl.clumped"),
                                          stringsAsFactors = F, h = T)

    clump[["FDR"]][[gene]]$SP2 <- sapply(
      clump[["FDR"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      })
    
    print(paste0("Gene ", gene, " done ( ",
                 which(s.genes==gene), " / ",
                 length(s.genes), " )"))
  }
  
  dir.create(tmp)
  s.genes <- unique(QTLs[QTLs$pvalue.eqtl<1e-5 | QTLs$pvalue.pqtl<1e-5, "gene"])  
  for(gene in s.genes){
  
    # P value ----
    clump[["P"]][[gene]] <- list()
    stats <- QTLs[QTLs$gene == gene, ]
    stats1 <- stats[, c("snpid", "pvalue.eqtl", "pvalue.eqtl", "FDR.eqtl")]
    colnames(stats1) <- c("SNP", "P", "pvalue.eqtl", "FDR.eqtl")
    stats2 <- stats[, c("snpid", "pvalue.pqtl", "pvalue.pqtl", "FDR.pqtl")]
    colnames(stats2) <- c("SNP", "P", "pvalue.pqtl", "FDR.pqtl")
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_eqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(stats2,
                file = paste0(tmp, "/", gene, "_pqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_eqtl.assoc", ",", gene, "_pqtl.assoc",
                  " --clump-p1 1e-5", " --clump-p2 0.05",
                  #" --clump-best --clump-annotate pval.eqtl,FDR.eqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_epqtl"))
    
    clump[["P"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                         "_epqtl.clumped"),
                                          stringsAsFactors = F, h = T)

    clump[["P"]][[gene]]$SP2 <- sapply(
      clump[["P"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      })

    
    print(paste0("Gene ", gene, " done ( ",
                 which(s.genes==gene), " / ",
                 length(s.genes), " )"))
  }
  saveRDS(clump, file = f.clump)
  unlink(tmp, recursive = T)
}else{
  ep.clump <- readRDS(file = f.clump)
}

## Clumping per QTL type -------------------------------------------------------
### eQTL -----------------------------------------------------------------------

# eqtl <- readRDS(paste0(get.path("results", local),
#                        "imputed/cis/final/",
#                        "eQTL_right_atrial_appendage_allpairs.RDS"))
# e.genes <- unique(eqtl[eqtl$FDR < 0.05, "gene"])

f.clump <- paste0(get.path("results", local),
                  "imputed/cis/eQTL_clump_relaxed.RDS")
if(!file.exists(f.clump)){
  eqtl <- readRDS(paste0(get.path("results", local),
                         "imputed/cis/final/",
                         "eQTL_right_atrial_appendage_allpairs.RDS"))
  clump <- list()
  tmp <- tempdir()
  dir.create(tmp)
  print(paste0("Temporary files: ", tmp))
  
  # FDR cutoff
  e.genes <- unique(eqtl[eqtl$FDR < 0.05, "gene"])
  for(gene in e.genes){
    # gene <- "MYOZ1"
    clump[["FDR"]][[gene]] <- list()
    stats <- eqtl[eqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    
    ## FDR < 0.05
    stats1 <- stats[, c("snpid", "FDR", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_eqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_eqtl.assoc",
                  " --clump-p1 0.05", " --clump-p2 0.8",
                  #" --clump-r2 0.25",
                  #" --clump-kb 500",
                  #" --clump-best --clump-annotate pval.eqtl,FDR.eqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_eqtl"))
    clump[["FDR"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_eqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["FDR"]][[gene]]$SP2 <- sapply(
      clump[["FDR"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    
    print(paste0("Gene ", gene, " done ( ",
                 which(e.genes==gene), " / ",
                 length(e.genes), " )  (FDR)"))
  }
    
  # P value cutoff
  e.genes <- unique(eqtl[eqtl$pvalue < 1e-5, "gene"])
  for(gene in e.genes){
    # gene <- "MYOZ1"
    stats <- eqtl[eqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    
    ## P < 1e-5
    clump[["P"]][[gene]] <- list()
    stats1 <- stats[, c("snpid", "pvalue", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_eqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_eqtl.assoc",
                  " --clump-p1 1e-5", " --clump-p2 0.05",
                  #" --clump-r2 0.25",
                  #" --clump-kb 500",
                  #" --clump-best --clump-annotate pval.eqtl,FDR.eqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_eqtl"))
    clump[["P"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_eqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["P"]][[gene]]$SP2 <- sapply(
      clump[["P"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )

    print(paste0("Gene ", gene, " done ( ",
                 which(e.genes==gene), " / ",
                 length(e.genes), " )  (P value)"))
  }
  
  saveRDS(clump, file = f.clump)
  unlink(tmp, recursive = T)
}else{
  e.clump <- readRDS(file = f.clump)
}

### pQTL -----------------------------------------------------------------------

# pqtl <- readRDS(paste0(get.path("results", local),
#                        "imputed/cis/final/",
#                        "pQTL_right_atrial_appendage_allpairs.RDS"))
# p.genes <- unique(pqtl[pqtl$FDR < 0.05, "gene"])

f.clump <- paste0(get.path("results", local),
                  "imputed/cis/pQTL_clump_relaxed.RDS")
if(!file.exists(f.clump)){
  pqtl <- readRDS(paste0(get.path("results", local),
                         "imputed/cis/final/",
                         "pQTL_right_atrial_appendage_allpairs.RDS"))
  clump <- list()
  tmp <- tempdir()
  dir.create(tmp)
  print(paste0("Temporary files: ", tmp))
  
  # FDR cutoff
  p.genes <- unique(pqtl[pqtl$FDR < 0.05, "gene"])
  for(gene in p.genes){
    # gene <- "MYOZ1"
    clump[["FDR"]][[gene]] <- list()
    stats <- pqtl[pqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "FDR", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_pqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_pqtl.assoc",
                  " --clump-p1 0.05", " --clump-p2 0.8",
                  #" --clump-best --clump-annotate pval.pqtl,FDR.pqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_pqtl"))
    clump[["FDR"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_pqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["FDR"]][[gene]]$SP2 <- sapply(
      clump[["FDR"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(p.genes==gene), " / ",
                 length(p.genes), " )  (FDR)"))
  }
  
  # P value cutoff
  p.genes <- unique(pqtl[pqtl$pvalue < 1e-5, "gene"])
  for(gene in p.genes){
    # gene <- "MYOZ1"
    clump[["P"]][[gene]] <- list()
    stats <- pqtl[pqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "pvalue", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_pqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_pqtl.assoc",
                  " --clump-p1 1e-5", " --clump-p2 0.05",
                  #" --clump-best --clump-annotate pval.pqtl,FDR.pqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_pqtl"))
    clump[["P"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_pqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["P"]][[gene]]$SP2 <- sapply(
      clump[["P"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(p.genes==gene), " / ",
                 length(p.genes), " )  (P value)"))
  }
  
  saveRDS(clump, file = f.clump)
  unlink(tmp, recursive = T)
}else{
  p.clump <- readRDS(file = f.clump)
}

### res eQTL -------------------------------------------------------------------

# res.eqtl <- readRDS(paste0(get.path("results", local),
#                        "imputed/cis/final/",
#                        "res_eQTL_right_atrial_appendage_allpairs.RDS"))
# res.genes <- unique(res.eqtl[res.eqtl$FDR < 0.05, "gene"])

f.clump <- paste0(get.path("results", local),
                  "imputed/cis/res_eQTL_clump_relaxed.RDS")
if(!file.exists(f.clump)){
  res.eqtl <- readRDS(paste0(get.path("results", local),
                             "imputed/cis/final/",
                             "res_eQTL_right_atrial_appendage_allpairs.RDS"))
  clump <- list()
  tmp <- tempdir()
  dir.create(tmp)
  print(paste0("Temporary files: ", tmp))
  
  # FDR cutoff
  res.genes <- unique(res.eqtl[res.eqtl$FDR < 0.05, "gene"])
  for(gene in res.genes){
    # gene <- "MYOZ1"
    clump[["FDR"]][[gene]] <- list()
    stats <- res.eqtl[res.eqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "FDR", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_res.eqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_res.eqtl.assoc",
                  " --clump-p1 0.05", " --clump-p2 0.8",
                  #" --clump-best --clump-annotate pval.res.eqtl,FDR.res.eqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_res.eqtl"))
    clump[["FDR"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_res.eqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["FDR"]][[gene]]$SP2 <- sapply(
      clump[["FDR"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(res.genes==gene), " / ",
                 length(res.genes), " )  (FDR)"))
  }
  
  # P value cutoff
  res.genes <- unique(res.eqtl[res.eqtl$pvalue < 1e-5, "gene"])
  for(gene in res.genes){
    # gene <- "MYOZ1"
    clump[["P"]][[gene]] <- list()
    stats <- res.eqtl[res.eqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "pvalue", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_res.eqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_res.eqtl.assoc",
                  " --clump-p1 1e-5", " --clump-p2 0.05",
                  #" --clump-best --clump-annotate pval.res.eqtl,FDR.res.eqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_res.eqtl"))
    clump[["P"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_res.eqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["P"]][[gene]]$SP2 <- sapply(
      clump[["P"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(res.genes==gene), " / ",
                 length(res.genes), " )  (P value)"))
  }
  
  saveRDS(clump, file = f.clump)
  unlink(tmp, recursive = T)
}else{
  rese.clump <- readRDS(file = f.clump)
}

### res pQTL -------------------------------------------------------------------

# res.pqtl <- readRDS(paste0(get.path("results", local),
#                        "imputed/cis/final/",
#                        "res_pQTL_right_atrial_appendage_allpairs.RDS"))
# res.genes <- unique(res.pqtl[res.pqtl$FDR < 0.05, "gene"])

f.clump <- paste0(get.path("results", local),
                  "imputed/cis/res_pQTL_clump_relaxed.RDS")
if(!file.exists(f.clump)){
  res.pqtl <- readRDS(paste0(get.path("results", local),
                             "imputed/cis/final/",
                             "res_pQTL_right_atrial_appendage_allpairs.RDS"))
  clump <- list()
  tmp <- tempdir()
  dir.create(tmp)
  print(paste0("Temporary files: ", tmp))
  
  # FDR cutoff
  res.genes <- unique(res.pqtl[res.pqtl$FDR < 0.05, "gene"])
  for(gene in res.genes){
    # gene <- "MYOZ1"
    clump[["FDR"]][[gene]] <- list()
    stats <- res.pqtl[res.pqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "FDR", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_res.pqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_res.pqtl.assoc",
                  " --clump-p1 0.05", " --clump-p2 0.8",
                  #" --clump-best --clump-annotate pval.res.pqtl,FDR.res.pqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_res.pqtl"))
    clump[["FDR"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_res.pqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["FDR"]][[gene]]$SP2 <- sapply(
      clump[["FDR"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(res.genes==gene), " / ",
                 length(res.genes), " )  (FDR)"))
  }
  
  # P value cutoff
  res.genes <- unique(res.pqtl[res.pqtl$pvalue < 1e-5, "gene"])
  for(gene in res.genes){
    # gene <- "MYOZ1"
    clump[["P"]][[gene]] <- list()
    stats <- res.pqtl[res.pqtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "pvalue", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_res.pqtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_res.pqtl.assoc",
                  " --clump-p1 1e-5", " --clump-p2 0.05",
                  #" --clump-best --clump-annotate pval.res.pqtl,FDR.res.pqtl",
                  #" --clump-verbose",
                  " --out ", gene, "_res.pqtl"))
    clump[["P"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_res.pqtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["P"]][[gene]]$SP2 <- sapply(
      clump[["P"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(res.genes==gene), " / ",
                 length(res.genes), " )  (P value)"))
  }
  
  saveRDS(clump, file = f.clump)
  unlink(tmp, recursive = T)
}else{
  resp.clump <- readRDS(file = f.clump)
}

### ratioQTL -----------------------------------------------------------------------

# ratio.qtl <- readRDS(paste0(get.path("results", local),
#                        "imputed/cis/final/",
#                        "ratios_right_atrial_appendage_allpairs.RDS"))
# ratio.genes <- unique(ratio.qtl[ratio.qtl$FDR < 0.05, "gene"])

f.clump <- paste0(get.path("results", local),
                  "imputed/cis/ratioQTL_clump_relaxed.RDS")
if(!file.exists(f.clump)){
  ratio.qtl <- readRDS(paste0(get.path("results", local),
                             "imputed/cis/final/",
                             "ratios_right_atrial_appendage_allpairs.RDS"))
  clump <- list()
  tmp <- tempdir()
  dir.create(tmp)
  print(paste0("Temporary files: ", tmp))
  
  # FDR cutoff
  ratio.genes <- unique(ratio.qtl[ratio.qtl$FDR < 0.05, "gene"])
  for(gene in ratio.genes){
    # gene <- "MYOZ1"
    clump[["FDR"]][[gene]] <- list()
    stats <- ratio.qtl[ratio.qtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "FDR", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_ratio.qtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_ratio.qtl.assoc",
                  " --clump-p1 0.05", " --clump-p2 0.8",
                  #" --clump-best --clump-annotate pval.ratio.qtl,FDR.ratio.qtl",
                  #" --clump-verbose",
                  " --out ", gene, "_ratio.qtl"))
    clump[["FDR"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_ratio.qtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["FDR"]][[gene]]$SP2 <- sapply(
      clump[["FDR"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(ratio.genes==gene), " / ",
                 length(ratio.genes), " )  (FDR)"))
  }
  
  # P value cutoff
  ratio.genes <- unique(ratio.qtl[ratio.qtl$pvalue < 1e-5, "gene"])
  for(gene in ratio.genes){
    # gene <- "MYOZ1"
    clump[["P"]][[gene]] <- list()
    stats <- ratio.qtl[ratio.qtl$gene == gene, ]
    chr <- unique(stats$chr)
    bfile <- paste0(get.path("genotype", local),
                    "AFHRI_B_all_trial")
    stats1 <- stats[, c("snpid", "pvalue", "pvalue", "FDR")]
    colnames(stats1) <- c("SNP", "P", "pvalue", "FDR")
    write.table(stats1,
                file = paste0(tmp, "/", gene, "_ratio.qtl.assoc"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    system(paste0("cd ", tmp, "\n",
                  "plink --bfile ",  bfile,
                  " --clump ", gene, "_ratio.qtl.assoc",
                  " --clump-p1 1e-5", " --clump-p2 0.05",
                  #" --clump-best --clump-annotate pval.ratio.qtl,FDR.ratio.qtl",
                  #" --clump-verbose",
                  " --out ", gene, "_ratio.qtl"))
    clump[["P"]][[gene]] <- read.table(file = paste0(tmp, "/", gene,
                                                        "_ratio.qtl.clumped"),
                                          stringsAsFactors = F, h = T)
    clump[["P"]][[gene]]$SP2 <- sapply(
      clump[["P"]][[gene]]$SP2,
      function(x){
        x <- strsplit(gsub("\\(1\\)", "", x), ",")
      }
    )
    print(paste0("Gene ", gene, " done ( ",
                 which(ratio.genes==gene), " / ",
                 length(ratio.genes), " )  (P value)"))
  }
  
  saveRDS(clump, file = f.clump)
  unlink(tmp, recursive = T)
}else{
  ratio.clump <- readRDS(file = f.clump)
}

### Summary -----------------------------------------------------------------------

library(plyr)
e.clump <- readRDS(file = paste0(get.path("results", local),
                                 "imputed/cis/eQTL_clump_relaxed.RDS"))
p.clump <- readRDS(file = paste0(get.path("results", local),
                                 "imputed/cis/pQTL_clump_relaxed.RDS"))
rese.clump <- readRDS(file = paste0(get.path("results", local),
                                    "imputed/cis/res_eQTL_clump_relaxed.RDS"))
resp.clump <- readRDS(file = paste0(get.path("results", local),
                                    "imputed/cis/res_pQTL_clump_relaxed.RDS"))
ratio.clump <- readRDS(file = paste0(get.path("results", local),
                                     "imputed/cis/ratioQTL_clump_relaxed.RDS"))
  
e.clump2 <- ldply(e.clump[["FDR"]], data.frame)
p.clump2 <- ldply(p.clump[["FDR"]], data.frame)
rese.clump2 <- ldply(rese.clump[["FDR"]], data.frame)
resp.clump2 <- ldply(resp.clump[["FDR"]], data.frame)
ratio.clump2 <- ldply(ratio.clump[["FDR"]], data.frame)
clump <- rbind(data.frame(omic = "eQTL",
                          gene = e.clump2$.id,
                          e.clump2,
                          stringsAsFactors = F),
               data.frame(omic = "pQTL",
                          gene = p.clump2$.id,
                          p.clump2,
                          stringsAsFactors = F),
               data.frame(omic = "res_eQTL",
                          gene = rese.clump2$.id,
                          rese.clump2,
                          stringsAsFactors = F),
               data.frame(omic = "res_pQTL",
                          gene = resp.clump2$.id,
                          resp.clump2,
                          stringsAsFactors = F),
               data.frame(omic = "ratioQTL",
                          gene = ratio.clump2$.id,
                          ratio.clump2,
                          stringsAsFactors = F))
sentinels <- data.frame(table(clump$gene, clump$omic))

ggplot(sentinels[sentinels$Freq>0 & sentinels$Var2 %in% c("eQTL", "pQTL"), ],
       aes(x=Freq, col=Var2, fill=Var2)) +
  geom_histogram(binwidth = 1) +
  scale_color_manual("omic", values = col.paired[c(2,10)]) +
  scale_fill_manual("omic", values = col.paired[c(1,9)]) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~Var2,
             scales = "free_y") +
  labs(title = "Number of independent cis SNP-clumps per gene",
       x = "Number of independent clumps")


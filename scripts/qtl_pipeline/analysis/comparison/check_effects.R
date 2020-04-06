# Figures for Poster / Paper

path <- "/home/icb/ines.assum/projects/symAtrial_QTL/scripts/"
#path <- "/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/scripts/"
setwd(path)
source("helper/helper.R")
source("analysis/boxplots/boxplots.R")
local=F

library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

HMGU.blue <- "#003E6E"
mygray <- "#C6DDEA"
col.paired <- brewer.pal(n = 11, "Paired")
col.set <- col.paired[c(2,8,9)]
col.set2 <- col.paired[c(9,5,7,3,8,2)]
col.set3 <- brewer.pal(n = 11, "PRGn")[c(4,6,8)]
col.set3 <- c(col.set[3], "#F7F7F7", col.set[1])
col.set4 <- c(col.set[1], col.set[3])

merged <- readRDS(file="/home/icb/ines.assum/projects/symAtrial_QTL/results/current/imputed/cis_old/final/comparisons/GTEx/eQTL_to_GTEx_v7_merged.RDS")
crit <- merged[!(abs(merged$s.beta_primary)<2),]

snps <- crit$snpid
genes <- crit$gene
df <- get.boxplot.data(snps, genes)

df2 <- df[, c(1,3:22)]
apply(df2[, -1], 2, table)

snps
genes

geno.file <- paste0(get.path("genotype", local), "genotype_imputed_common_samples_t_raa.txt")
snp.file <- tempfile()
write.table(snps, file=snp.file,
            row.names = F, col.names = F, sep="\n", quote=F)
snp.data.final = fread(paste("fgrep -f ", snp.file, " -w ", geno.file, sep=""), header=F)
setnames(snp.data.final, colnames(read.csv(geno.file, sep="\t", nrows=5)))
snp.data.final = as.data.frame(snp.data.final)
snp.data.final$X <- paste0(snp.data.final$X, "_recoded")
rownames(snp.data.final) <- snp.data.final$X
snp.data.final$X <- NULL

geno.file.impute <- paste0(get.path("genotype", local), "AFHRI_B_imputed.txt")
snp.data.impute = fread(paste("fgrep -f ", snp.file, " -w ", geno.file.impute, sep=""), header=F)
setnames(snp.data.impute, colnames(read.csv(geno.file.impute, sep="\t", nrows=5)))
snp.data.impute = as.data.frame(snp.data.impute)
snp.data.impute$imputed_marker_id <- paste0(snp.data.impute$imputed_marker_id, "_dosage")
rownames(snp.data.impute) <- snp.data.impute$imputed_marker_id
snp.data.impute[,1:6] <- NULL

unlink(snp.file)

df.all <- merge(df, data.frame(id=colnames(snp.data.final),
                               t(snp.data.final),
                               stringsAsFactors = F),
                by=c("id"), all=T)
df.all <- merge(df.all, data.frame(id=colnames(snp.data.impute),
                                   t(snp.data.impute),
                                   stringsAsFactors = F),
                by=c("id"), all=T)

make.plot <- function(df.all, snps, genes, crit, i){
  #i=17
  df3 <- df.all[, c("id",
                    gsub(":", ".", snps[i]),
                    paste0(gsub(":", ".", snps[i]), "_dosage"),
                    paste0(gsub(":", ".", snps[i]), "_recoded"),
                    paste0(gsub("-", ".", genes[i]), "_trans"))]
  colnames(df3) <- c("id", "genotype", "dosage", "recoded", "gene")
  
  g1 <- ggplot(df3, aes(x=genotype, y=gene, col=genotype)) +
    geom_point(position = position_jitter(width = 0.3), size=3) +
    geom_boxplot(alpha=0.2) +
    theme_bw() +
    labs(x=snps[i],
         y=genes[i],
         title=paste0("beta: ", crit$beta[i], ", pvalue: ", crit$pvalue[i]))
  
  g2 <- ggplot(df3, aes(x=recoded, y=gene, col=genotype)) +
    geom_point(size=3) +
    theme_bw() +
    labs(x=snps[i],
         y=genes[i],
         title=paste0("beta: ", crit$beta[i], ", pvalue: ", crit$pvalue[i]))
  
  gg <- ggarrange(g1, g2)
  return(gg)
}

make.plot(df.all, snps, genes, crit, 1)
make.plot(df.all, snps, genes, crit, 2)
make.plot(df.all, snps, genes, crit, 3)
make.plot(df.all, snps, genes, crit, 4)
make.plot(df.all, snps, genes, crit, 5)
make.plot(df.all, snps, genes, crit, 6)
make.plot(df.all, snps, genes, crit, 7)
make.plot(df.all, snps, genes, crit, 8)
make.plot(df.all, snps, genes, crit, 9)
make.plot(df.all, snps, genes, crit, 10)
make.plot(df.all, snps, genes, crit, 11)
make.plot(df.all, snps, genes, crit, 12)
make.plot(df.all, snps, genes, crit, 13)
make.plot(df.all, snps, genes, crit, 14)
make.plot(df.all, snps, genes, crit, 15)
make.plot(df.all, snps, genes, crit, 16)
make.plot(df.all, snps, genes, crit, 17)
make.plot(df.all, snps, genes, crit, 18)
make.plot(df.all, snps, genes, crit, 19)
make.plot(df.all, snps, genes, crit, 20)

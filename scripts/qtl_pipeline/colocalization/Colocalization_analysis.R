# Colocalization analysis ------------------------------------------------------

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

## Independent pQTLs with GTEx eQTLs -------------------------------------------

gtex.pqtl <- readRDS(file=paste0(get.path("results", local),
                                 "imputed/cis/final/comparisons/GTEx/",
                                 "pQTL.FDR_to_GTEx_v7.1e-5_merged_all_matched.RDS"))

eqtl.pqtl <- readRDS(file=paste0(get.path("results"),
                                 "imputed/cis/QTL_res_cis_regions.RDS"))
eqtl.pqtl <- eqtl.pqtl[which(eqtl.pqtl$FDR.eqtl < 0.05 |
                               eqtl.pqtl$FDR.pqtl <0.05), ]
eqtl.pqtl$type <- "eQTL"
eqtl.pqtl$type[which(eqtl.pqtl$FDR.pqtl < 0.05)] <- "pQTL"
eqtl.pqtl$type[which(eqtl.pqtl$FDR.eqtl < 0.05 &
                       eqtl.pqtl$FDR.pqtl < 0.05)] <- "eQTL/pQTL"
eqtl.pqtl$type[eqtl.pqtl$Group2 == 1] <- "shared eQTL/pQTL"
eqtl.pqtl$type[eqtl.pqtl$Group5 == 1] <- "independent pQTL"
eqtl.pqtl$type[eqtl.pqtl$Group1 == 1] <- "independent eQTL"

eqtl.pqtl$type <- factor(eqtl.pqtl$type,
                         levels = c("shared eQTL/pQTL",
                                    "independent eQTL",
                                    "independent pQTL",
                                    "eQTL/pQTL", "eQTL", "pQTL"),
                         ordered = T)

indp <- eqtl.pqtl[eqtl.pqtl$type == "independent pQTL", ]

indp <- merge(indp, gtex.pqtl,
              by.x = c("snps", "gene"), by.y = c("snpid", "gene"),
              all.x = T,
              suffixes = c(".afhri", ".comp"))
# dim(indp)
# colnames(indp)
# table(indp$pval_nominal<1e-5, useNA = "always")
indp$same <- sign(indp$beta.pqtl) == sign(indp$s.beta_secondary)
# plot(indp$beta.pqtl, indp$s.beta_secondary, col=indp$col)
# plot(indp$beta.pqtl, indp$s.beta_secondary, col=indp$same)

indp$legend <- indp$gene
indp$legend[indp$pval_nominal>1e-5] <- NA

library(ggrepel)

indp$gtex.1e3 <- "SNP-gene pair not evaluated in GTEx"
indp$gtex.1e3[which(indp$pval_nominal>1e-3)] <- "GTEx eQTL > 1e-3"
indp$gtex.1e3[which(indp$pval_nominal<1e-4)] <- "GTEx eQTL < 1e-4 (AAMDC)"
indp$gtex.1e3[which(indp$pval_nominal<1e-5)] <- "GTEx eQTL < 1e-5 (YBX3)"
indp$gtex.1e3[which(indp$pval_nominal<1e-12)] <- "GTEx eQTL < 1e-12 (CRYZ)"
indp$gtex.1e3[which(indp$pval_nominal<1e-18)] <- "GTEx eQTL < 1e-18 (ANXA5)"
indp$gtex.1e3[which(indp$pval_nominal<1e-19)] <- "GTEx eQTL < 1e-19 (NIPSNAP1)"

indp$gtex.1e3 <- factor(indp$gtex.1e3,
                        levels = c("GTEx eQTL < 1e-19 (NIPSNAP1)",
                                   "GTEx eQTL < 1e-18 (ANXA5)",
                                   "GTEx eQTL < 1e-12 (CRYZ)",
                                   "GTEx eQTL < 1e-5 (YBX3)",
                                   "GTEx eQTL < 1e-4 (AAMDC)",
                                   "GTEx eQTL > 1e-3",
                                   "SNP-gene pair not evaluated in GTEx"),
                        ordered = T)
gtex.top <- indp[order(indp$pval_nominal), ]
gtex.top <- gtex.top[!duplicated(gtex.top$gene), ]
gtex.top <- gtex.top[which(gtex.top$pval_nominal<1e-3), ]

col.paired <- brewer.pal(n = 11, "Paired")
col.set <- col.paired[c(2,9,8)]
col.set2 <- brewer.pal(n=11, "BrBG")[c(10:7, 4)]

ggarrange(
  ggplot(indp[order(indp$pval_nominal, decreasing = T), ],
         aes(x=beta.pqtl, y=s.beta_secondary,
             col=gtex.1e3)) +
    geom_abline(intercept = 0, slope = 1, col="darkgrey") +
    geom_hline(yintercept = 0, col="darkgrey") +
    geom_vline(xintercept = 0, col="darkgrey") +
    geom_point(size=3) +
    geom_label_repel(data = gtex.top,
                     aes(x=beta.pqtl, y=s.beta_secondary,
                         label=gene, fill=gtex.1e3),
                     nudge_y = -0.2, col="black") +
    scale_color_manual(values = c(col.set2, col.set)) +
    scale_fill_manual(values = c(col.set2, col.set)) +
    guides(fill=FALSE) +
    theme_bw() +
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          #legend.position = "bottom",
          legend.direction="vertical") +
    xlim(c(min(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T),
           max(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T))) +
    ylim(c(min(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T),
           max(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T))) +
    labs(title="All independent pQTLs:\nAFHRI pQTLs compared to GTEx eQTLs",
         x="Effect size AFHRI pQTL",
         y="Effect size GTEx eQTL"),
  ggplot(indp[order(indp$gtex.1e3, decreasing = T), ],
                  aes(x=beta.pqtl, y=beta.eqtl,
                      col=gtex.1e3)) +
    geom_abline(intercept = 0, slope = 1, col="darkgrey") +
    geom_hline(yintercept = 0, col="darkgrey") +
    geom_vline(xintercept = 0, col="darkgrey") +
    geom_point(size=3) +
    # geom_label_repel(data = gtex.top,
    #                  aes(x=beta.pqtl, y=beta.eqtl,
    #                      label=gene, fill=gtex.1e3),
    #                  nudge_y = -0.2, col="black") +
    scale_color_manual(values = c(col.set2, col.set)) +
    scale_fill_manual(values = c(col.set2, col.set)) +
    guides(fill=FALSE) +
    theme_bw() +
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          #legend.position = "bottom",
          legend.direction="vertical") +
    xlim(c(min(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T),
           max(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T))) +
    ylim(c(min(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T),
           max(c(indp$s.beta_secondary, indp$beta.eqtl, indp$beta.pqtl),
               na.rm = T))) +
    labs(title="All independent pQTLs:\nAFHRI pQTLs compared to AFHRI eQTLs",
         x="Effect size AFHRI pQTL",
         y="Effect size AFHRI eQTL"),
  nrow = 2, ncol = 1,
  common.legend = T, legend = "bottom"
)

## Run coloc -------------------------------------------------------------------

if(!file.exists(paste0(get.path("results", local),
                       "/imputed/cis/",
                       "eQTL_pQTL_coloc_results.RDS"))){
  eqtls <- readRDS(file=paste0(get.path("results", local),
                               "/imputed/cis/final/",
                               "eQTL_right_atrial_appendage_allpairs.RDS"))
  pqtls <- readRDS(file=paste0(get.path("results", local),
                               "/imputed/cis/final/",
                               "pQTL_right_atrial_appendage_allpairs.RDS"))
  qtls.group <- readRDS(file=paste0(get.path("results", local),
                                    "/imputed/cis/",
                                    "QTL_res_all_snp_group_anno.RDS"))
  colnames(qtls.group)
   
  
  QTLs <- merge(eqtls, pqtls,
                by = c("gene", "snpid",
                       "chr", "variant_pos", "Allele1", "Allele2",
                       "variant_id", "rs_id", "gtex.variant_id",
                       "gtex.match_alleles", "gene_id"),
                suffixes = c(".eqtl", ".pqtl"))
  QTLres <- QTLs[!is.na(QTLs$pvalue.eqtl) &
                   !is.na(QTLs$pvalue.pqtl), ]
  QTLres2 <- merge(QTLres, qtls.group[c("snps", "gene",
                                        "eQTL", "pQTL",
                                        "Group2", "Group1", "Group5")],
                  by.x = c("snpid", "gene"), by.y = c("snps", "gene"),
                  all.x = T)
  res.genes <- unique(QTLres$gene)
  
  expr.t <- readRDS(paste0(get.path("dataset", local),
                           "AFHRI_B_transcriptomics_QC_symbol.RDS"))
  expr.p <- readRDS(paste0(get.path("dataset", local),
                           "AFHRI_B_proteomics_QC_symbol.RDS"))
  expr.t.scale  <- read.table(paste0(get.path("peer", local),
                                     "normalized/eqtl/",
                                     "imputed_expression_t.txt"),
                              row.names = 1)
  expr.p <- readRDS(paste0(get.path("dataset", local),
                           "AFHRI_B_proteomics_QC_symbol.RDS"))
  expr.p.scale  <- read.table(paste0(get.path("peer", local),
                                     "normalized/pqtl/",
                                     "imputed_expression_p.txt"),
                              row.names = 1)
  coloc.abf.res <- list()
  coloc.summary <- data.frame(symbol = res.genes,
                              matrix(NA,
                                     nrow = length(res.genes),
                                     ncol = 1+5+5),
                              stringsAsFactors = F,
                              row.names = res.genes)
  colnames(coloc.summary) <- c("symbol",
                               "N.snps",
                               "PP.noQTL.abf", "PP.eQTL.abf", "PP.pQTL.abf",
                               "PP.indQTL.abf", "PP.sharedQTL.abf",
                               "eQTL", "pQTL",
                               "sharedQTL", "indeQTL", "indpQTL")
  
  # H0: neither trait has a genetic association in the region
  # H1: only trait 1 has a genetic association in the region
  # H2: only trait 2 has a genetic association in the region
  # H3: both traits are associated, but with different causal variants
  # H4: both traits are associated and share a single causal variant

  for(gene in res.genes){
    QTL <- QTLres2[QTLres2$gene == gene, ]
    eqtl <- list(N = 75,
                 beta = QTL$beta.eqtl,
                 varbeta = (QTL$beta.eqtl/QTL$statistic.eqtl)^2,
                 type = "quant",
                 sdY = sd(expr.t.scale[gene, ]),
                 snp = QTL$snpid)
    pqtl <- list(N = 75,
                 beta = QTL$beta.pqtl,
                 varbeta = (QTL$beta.pqtl/QTL$statistic.pqtl)^2,
                 type = "quant",
                 sdY = sd(expr.p.scale[gene, ]),
                 snp = QTL$snpid)
    print(paste0("Gene ", gene,
                 " (", which(res.genes==gene),
                 "/", length(res.genes), "):"))
    coloc.abf.res[[gene]] <- coloc.abf(eqtl, pqtl)
    coloc.summary[gene, "symbol"] <- gene
    coloc.summary[gene, -1] <- c(coloc.abf.res[[gene]]$summary,
                                 colSums(QTL[, c("eQTL", "pQTL",
                                                 "Group2", "Group1", "Group5")]))
  }
  saveRDS(list(coloc.abf.res = coloc.abf.res,
               coloc.summary = coloc.summary),
        file = paste0(get.path("results", local),
                      "/imputed/cis/",
                      "eQTL_pQTL_coloc_results.RDS"))
} else {
  coloc.eqtl.pqtl <- readRDS(paste0(get.path("results", local),
                                    "/imputed/cis/",
                                    "eQTL_pQTL_coloc_results.RDS"))
}
coloc.summary <- coloc.eqtl.pqtl[["coloc.summary"]]


ggplot() +
  geom_histogram(data=coloc.summary[coloc.summary$eQTL+coloc.summary$pQTL>0 &
                                      coloc.summary$sharedQTL+
                                      coloc.summary$indeQTL+
                                      coloc.summary$indpQTL==0, ],
                 aes(x=PP.sharedQTL.abf,
                     fill="d", col="d"),
                 alpha = 0.5, bins = 20, size=0.3) +
  geom_histogram(data=coloc.summary[coloc.summary$sharedQTL>0, ],
                 aes(x=PP.sharedQTL.abf,
                     fill="a", col="a"),
                 alpha = 0.5, bins = 20, size=0.3) +
  geom_histogram(data=coloc.summary[coloc.summary$indeQTL>0, ],
                 aes(x=PP.sharedQTL.abf,
                     fill="b", col="b"),
                 alpha = 0.5, bins = 20, size=0.3) +
  geom_histogram(data=coloc.summary[coloc.summary$indpQTL>0, ],
                 aes(x=PP.sharedQTL.abf,
                     fill="c", col="c"),
                 alpha = 0.5, bins = 20, size=0.3) +
  theme_bw() +
  scale_color_manual("functional QTL category",
                    values = c(col.set, "darkgrey"),
                    labels = c("shared eQTL/pQTL",
                               "independent eQTL",
                               "independent pQTL",
                               "all other")) +
  scale_fill_manual("functional QTL category",
                    values = c(col.set, "darkgrey"),
                    labels = c("shared eQTL/pQTL",
                               "independent eQTL",
                               "independent pQTL",
                               "all other")) +
  labs(title = "Posterior probability for same causal variant")

coloc.summary[coloc.summary$sharedQTL>0 & coloc.summary$indpQTL>0, ]
coloc.summary[coloc.summary$indpQTL>0 & coloc.summary$PP.sharedQTL.abf>0.5, ]
coloc.summary[coloc.summary$sharedQTL>0, ]

print("CDH13 is the only gene with a shared QTL, that does not colocalize at a PP > 0.5")


library(plyr)
ep.clump <- readRDS(paste0(get.path("results", local),
                           "imputed/cis/eQTL_pQTL_overlap_clump.RDS"))
ep.clump2 <- ldply(ep.clump[["FDR"]], data.frame)

test <- ep.clump2[1:3, ]
test1 <- gsub("\\(.*", "", unlist(ep.clump2[1, "SP2"]))
test2 <- gsub("\\(.*", "", unlist(ep.clump2[2, "SP2"]))
test3 <- gsub("\\(.*", "", unlist(ep.clump2[3, "SP2"]))


## coloc for eQTL/pQTL clumps --------------------------------------------------

# setwd("~/work/symAtrial_QTL/scripts/revision")
# source("../helper/helper.R")
# local=F

if(!file.exists(paste0(get.path("results", local),
                       "/imputed/cis/",
                       "eQTL_pQTL_clumps_coloc_results.RDS"))){
  library(plyr)
  ep.clump <- readRDS(paste0(get.path("results", local),
                             "imputed/cis/eQTL_pQTL_overlap_clump.RDS"))
  ep.clump2 <- ldply(ep.clump[["FDR"]], data.frame)
  dim(ep.clump2)
  rownames(ep.clump2) <- paste0(ep.clump2$.id, "_", ep.clump2$SNP)
  
  if(!file.exists(paste0(get.path("results", local),
                         "/imputed/cis/",
                         "eQTL_pQTL_coloc_input.RDS"))){
    eqtls <- readRDS(file=paste0(get.path("results", local),
                                 "/imputed/cis/final/",
                                 "eQTL_right_atrial_appendage_allpairs.RDS"))
    pqtls <- readRDS(file=paste0(get.path("results", local),
                                 "/imputed/cis/final/",
                                 "pQTL_right_atrial_appendage_allpairs.RDS"))
    qtls.group <- readRDS(file=paste0(get.path("results", local),
                                      "/imputed/cis/",
                                      "QTL_res_all_snp_group_anno.RDS"))
    colnames(qtls.group)
    
    QTLs <- merge(eqtls, pqtls,
                  by = c("gene", "snpid",
                         "chr", "variant_pos", "Allele1", "Allele2",
                         "variant_id", "rs_id", "gtex.variant_id",
                         "gtex.match_alleles", "gene_id"),
                  suffixes = c(".eqtl", ".pqtl"))
  
    QTLs <- merge(QTLs, qtls.group[c("snps", "gene",
                                     "eQTL", "pQTL",
                                     "reseQTL", "respQTL", "ratioQTL",
                                     "Group2", "Group1", "Group5")],
                  by.x = c("snpid", "gene"), by.y = c("snps", "gene"),
                  all.x = T)
    saveRDS(QTLs,
            file = paste0(get.path("results", local),
                          "/imputed/cis/",
                          "eQTL_pQTL_coloc_input.RDS"))
  }else{
    QTLs <- readRDS(paste0(get.path("results", local),
                           "/imputed/cis/",
                           "eQTL_pQTL_coloc_input.RDS"))
  }
  
  expr.t.scale  <- read.table(paste0(get.path("peer", local),
                                     "normalized/eqtl/",
                                     "imputed_expression_t.txt"),
                              row.names = 1)
  expr.p.scale  <- read.table(paste0(get.path("peer", local),
                                     "normalized/pqtl/",
                                     "imputed_expression_p.txt"),
                              row.names = 1)
  coloc.abf.res <- list()
  coloc.summary <- data.frame(leadSNP = ep.clump2$SNP,
                              gene = ep.clump2$.id,
                              matrix(NA,
                                     nrow = dim(ep.clump2)[1],
                                     ncol = 1+5+5+3+5+3),
                              stringsAsFactors = F,
                              row.names = paste0(ep.clump2$.id, "_", ep.clump2$SNP))
  colnames(coloc.summary) <- c("leadSNP", "gene",
                               "N.clump",
                               "PP.noQTL.abf", "PP.eQTL.abf", "PP.pQTL.abf",
                               "PP.indQTL.abf", "PP.sharedQTL.abf",
                               "lead.eQTL", "lead.pQTL",
                               "lead.reseQTL", "lead.respQTL", "lead.ratioQTL",
                               "lead.sharedQTL", "lead.indeQTL", "lead.indpQTL",
                               "eQTLs", "pQTLs",
                               "reseQTLs", "respQTLs", "ratioQTLs",
                               "sharedQTLs", "indeQTLs", "indpQTLs")
  
  # H0: neither trait has a genetic association in the region
  # H1: only trait 1 has a genetic association in the region
  # H2: only trait 2 has a genetic association in the region
  # H3: both traits are associated, but with different causal variants
  # H4: both traits are associated and share a single causal variant

  for(clump in paste0(ep.clump2$.id, "_", ep.clump2$SNP)){
    gene <- coloc.summary[clump, "gene"]
    lead <- coloc.summary[clump, "leadSNP"]
    snps <- c(ep.clump2[clump, "SNP"],
              gsub("\\(.*", "", unlist(ep.clump2[clump, "SP2"])))
    QTL <- QTLs[QTLs$gene == gene & QTLs$snpid %in% snps, ]
    eqtl <- list(N = 75,
                 beta = QTL$beta.eqtl,
                 varbeta = (QTL$beta.eqtl/QTL$statistic.eqtl)^2,
                 type = "quant",
                 sdY = sd(expr.t.scale[gene, ]),
                 snp = QTL$snpid)
    pqtl <- list(N = 75,
                 beta = QTL$beta.pqtl,
                 varbeta = (QTL$beta.pqtl/QTL$statistic.pqtl)^2,
                 type = "quant",
                 sdY = sd(expr.p.scale[gene, ]),
                 snp = QTL$snpid)
    
    print(paste0("Clump ", clump,
                 " (", which(rownames(ep.clump2)==clump),
                 "/", length(rownames(ep.clump2)), "):"))
    coloc.abf.res[[clump]] <- coloc.abf(eqtl, pqtl)
    
    coloc.summary[clump, c("N.clump",
                           "PP.noQTL.abf", "PP.eQTL.abf", "PP.pQTL.abf",
                           "PP.indQTL.abf", "PP.sharedQTL.abf")] <-
      coloc.abf.res[[clump]]$summary
    coloc.summary[clump, c("lead.eQTL", "lead.pQTL",
                           "lead.reseQTL", "lead.respQTL", "lead.ratioQTL",
                           "lead.sharedQTL", "lead.indeQTL", "lead.indpQTL")] <-
      QTL[QTL$snpid==lead, c("eQTL", "pQTL",
                             "reseQTL", "respQTL", "ratioQTL",
                             "Group2", "Group1", "Group5")]
    coloc.summary[clump, c("eQTLs", "pQTLs",
                           "reseQTLs", "respQTLs", "ratioQTLs",
                           "sharedQTLs", "indeQTLs", "indpQTLs")] <-
      colSums(QTL[QTL$snpid %in% snps, c("eQTL", "pQTL",
                                         "reseQTL", "respQTL", "ratioQTL",
                                         "Group2", "Group1", "Group5")])
    
  }
  saveRDS(list(coloc.abf.res = coloc.abf.res,
               coloc.summary = coloc.summary),
        file = paste0(get.path("results", local),
                      "/imputed/cis/",
                      "eQTL_pQTL_clump_coloc_results.RDS"))
} else {
  coloc.eqtl.pqtl.clump <- readRDS(paste0(get.path("results", local),
                                          "/imputed/cis/",
                                          "eQTL_pQTL_clump_coloc_results.RDS"))
}
coloc.ep <- coloc.eqtl.pqtl.clump$coloc.summary

dim(coloc.ep)
length(unique(coloc.ep$gene))
table(coloc.ep$PP.sharedQTL.abf>0.5)
length(unique(coloc.ep$gene[coloc.ep$PP.sharedQTL.abf>0.5]))

ep.ov <- ep.qtl.ov[coloc.ep$lead.eQTL==1 & coloc.ep$lead.pQTL==1, ]
dim(ep.ov)
unique(ep.ov$gene)


ep.ov.shared <- coloc.ep[coloc.ep$lead.sharedQTL==1, ]
ep.ov.inde <- coloc.ep[coloc.ep$lead.indeQTL==1, ]
ep.ov.indp <- coloc.ep[coloc.ep$lead.indpQTL==1, ]

ep.ov.ep <- coloc.ep[coloc.ep$lead.eQTL==1 & coloc.ep$lead.pQTL==1, ]
ep.ov.e <- coloc.ep[coloc.ep$lead.eQTL==1 & coloc.ep$lead.pQTL==0, ]
ep.ov.p <- coloc.ep[coloc.ep$lead.eQTL==0 & coloc.ep$lead.pQTL==1, ]

func.qtls <- data.frame(eQTL_pQTL = c(dim(ep.ov.ep)[1],
                                      NA,
                                      sum(ep.ov.ep$lead.sharedQTL==1),
                                      sum(ep.ov.ep$lead.indeQTL==1),
                                      sum(ep.ov.ep$lead.indpQTL==1),
                                      sum(ep.ov.ep$lead.eQTL==1 &
                                            ep.ov.ep$lead.pQTL==1 &
                                            ep.ov.ep$lead.reseQTL==1 &
                                            ep.ov.ep$lead.respQTL==1),
                                      NA,
                                      sum(ep.ov.ep$PP.sharedQTL.abf>0.5),
                                      sum(ep.ov.ep$PP.eQTL.abf>0.5),
                                      sum(ep.ov.ep$PP.pQTL.abf>0.5),
                                      sum(ep.ov.ep$PP.indQTL.abf>0.5)),
                        eQTL_nopQTL = c(dim(ep.ov.e)[1],
                                      NA,
                                      sum(ep.ov.e$lead.sharedQTL==1),
                                      sum(ep.ov.e$lead.indeQTL==1),
                                      sum(ep.ov.e$lead.indpQTL==1),
                                      sum(ep.ov.e$lead.eQTL==1 &
                                            ep.ov.e$lead.pQTL==1 &
                                            ep.ov.e$lead.reseQTL==1 &
                                            ep.ov.e$lead.respQTL==1),
                                      NA,
                                      sum(ep.ov.e$PP.sharedQTL.abf>0.5),
                                      sum(ep.ov.e$PP.eQTL.abf>0.5),
                                      sum(ep.ov.e$PP.pQTL.abf>0.5),
                                      sum(ep.ov.e$PP.indQTL.abf>0.5)),
                        noeQTL_pQTL = c(dim(ep.ov.p)[1],
                                      NA,
                                      sum(ep.ov.p$lead.sharedQTL==1),
                                      sum(ep.ov.p$lead.indeQTL==1),
                                      sum(ep.ov.p$lead.indpQTL==1),
                                      sum(ep.ov.p$lead.eQTL==1 &
                                            ep.ov.p$lead.pQTL==1 &
                                            ep.ov.p$lead.reseQTL==1 &
                                            ep.ov.p$lead.respQTL==1),
                                      NA,
                                      sum(ep.ov.p$PP.sharedQTL.abf>0.5),
                                      sum(ep.ov.p$PP.eQTL.abf>0.5),
                                      sum(ep.ov.p$PP.pQTL.abf>0.5),
                                      sum(ep.ov.p$PP.indQTL.abf>0.5)),
                        row.names = c("Loci",
                                      "Residual_approach",
                                      "Shared_eQTL_pQTL",
                                      "Independent_eQTL",
                                      "Independent_pQTL",
                                      "Independent_eQTL_pQTL",
                                      "Coloc_posterior_probability",
                                      "Shared_eQTL_pQTL.abf",
                                      "Independent_eQTL.abf",
                                      "Independent_pQTL.abf",
                                      "Independent_eQTL_pQTL.abf"))

func.qtls

# Out of 1 243 genes with transcriptomics and proteomics measurement, 190 genes
# had at least one significant eQTL or pQTL (FDR<0.05). Performing LD clumping
# (p1 FDR 0.05, p2 FDR 0.8) we identified 314 clumps for 180 genes,
# of which 30 clumps (for 21 genes) formally colocalized
# (posterior probability for a shared causal variant >0.5).
# Compared to 642 SNP-gene-pairs (for 21 genes), with an overlapping
# significant (FDR<0.05) eQTL and pQTL, for 19 clumps including 12 genes,
# the lead SNP was a significant eQTL and pQTL (FDR<0.05) and of those,
# all but one (lead SNP rs7202898 for gene CDH13) formally colocalize (PP>0.5).

# new table: eQTL / pQTL overlap

coloc.ep <- readRDS(paste0(get.path("results", local),
                           "/imputed/cis/",
                           "eQTL_pQTL_clump_coloc_results.RDS"))$coloc.summary

ep.ov.ep <- coloc.ep[coloc.ep$lead.eQTL==1 & coloc.ep$lead.pQTL==1, ]
ep.ov.e <- coloc.ep[coloc.ep$lead.eQTL==1 & coloc.ep$lead.pQTL==0, ]
ep.ov.p <- coloc.ep[coloc.ep$lead.eQTL==0 & coloc.ep$lead.pQTL==1, ]

func.qtls <- data.frame(eQTL_pQTL = c(dim(ep.ov.ep)[1],
                                      NA,
                                      sum(ep.ov.ep$lead.sharedQTL==1),
                                      sum(ep.ov.ep$lead.indeQTL==1),
                                      sum(ep.ov.ep$lead.indpQTL==1),
                                      sum(ep.ov.ep$lead.eQTL==1 &
                                            ep.ov.ep$lead.pQTL==1 &
                                            ep.ov.ep$lead.reseQTL==1 &
                                            ep.ov.ep$lead.respQTL==1),
                                      NA,
                                      sum(ep.ov.ep$PP.sharedQTL.abf>0.5),
                                      sum(ep.ov.ep$PP.eQTL.abf>0.5),
                                      sum(ep.ov.ep$PP.pQTL.abf>0.5),
                                      sum(ep.ov.ep$PP.indQTL.abf>0.5)),
                        eQTL_nopQTL = c(dim(ep.ov.e)[1],
                                      NA,
                                      sum(ep.ov.e$lead.sharedQTL==1),
                                      sum(ep.ov.e$lead.indeQTL==1),
                                      sum(ep.ov.e$lead.indpQTL==1),
                                      sum(ep.ov.e$lead.eQTL==1 &
                                            ep.ov.e$lead.pQTL==1 &
                                            ep.ov.e$lead.reseQTL==1 &
                                            ep.ov.e$lead.respQTL==1),
                                      NA,
                                      sum(ep.ov.e$PP.sharedQTL.abf>0.5),
                                      sum(ep.ov.e$PP.eQTL.abf>0.5),
                                      sum(ep.ov.e$PP.pQTL.abf>0.5),
                                      sum(ep.ov.e$PP.indQTL.abf>0.5)),
                        noeQTL_pQTL = c(dim(ep.ov.p)[1],
                                      NA,
                                      sum(ep.ov.p$lead.sharedQTL==1),
                                      sum(ep.ov.p$lead.indeQTL==1),
                                      sum(ep.ov.p$lead.indpQTL==1),
                                      sum(ep.ov.p$lead.eQTL==1 &
                                            ep.ov.p$lead.pQTL==1 &
                                            ep.ov.p$lead.reseQTL==1 &
                                            ep.ov.p$lead.respQTL==1),
                                      NA,
                                      sum(ep.ov.p$PP.sharedQTL.abf>0.5),
                                      sum(ep.ov.p$PP.eQTL.abf>0.5),
                                      sum(ep.ov.p$PP.pQTL.abf>0.5),
                                      sum(ep.ov.p$PP.indQTL.abf>0.5)),
                        row.names = c("Loci",
                                      "Residual_approach",
                                      "Shared_eQTL_pQTL",
                                      "Independent_eQTL",
                                      "Independent_pQTL",
                                      "Independent_eQTL_pQTL",
                                      "Coloc_posterior_probability",
                                      "Shared_eQTL_pQTL.abf",
                                      "Independent_eQTL.abf",
                                      "Independent_pQTL.abf",
                                      "Independent_eQTL_pQTL.abf"))

func.qtls

# xtable(func.qtls,
#        digits = 0, display = rep("fg", 4))

library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv",
                               package = "UpSetR"),
                   header = T, sep = ";")

upset(movies, sets = c("Action", "Adventure", "Comedy", "Drama", "Mystery",
                       "Thriller", "Romance", "War", "Western"),
      mb.ratio = c(0.55, 0.45), order.by = "freq")

library(cowplot)
library(gridExtra)
library(UpSetR)

coloc.ep <- readRDS(paste0(get.path("results", local),
                           "/imputed/cis/",
                           "eQTL_pQTL_clump_coloc_results.RDS"))$coloc.summary

coloc.ep$all.eQTL.or.pQTL <- 1

coloc.ep$lead.eQTL.pQTL <- (coloc.ep$lead.eQTL==1 & coloc.ep$lead.pQTL==1)*1
coloc.ep$lead.eQTL.no_pQTL <- (coloc.ep$lead.eQTL==1 & coloc.ep$lead.pQTL==0)*1
coloc.ep$lead.no_eQTL.pQTL <- (coloc.ep$lead.eQTL==0 & coloc.ep$lead.pQTL==1)*1

coloc.ep$coloc.sharedQTL <- (coloc.ep$PP.sharedQTL.abf>0.5)*1
coloc.ep$coloc.onlyeQTL <- (coloc.ep$PP.eQTL.abf>0.5)*1
coloc.ep$coloc.onlypQTL <- (coloc.ep$PP.pQTL.abf>0.5)*1

upset1 <-upset(coloc.ep,
               sets = rev(c("lead.eQTL.pQTL",
                            "lead.sharedQTL",
                            "coloc.sharedQTL")),
               mb.ratio = c(0.55, 0.45),
               order.by = "freq",
               #group.by = "sets",
               main.bar.color = brewer.pal(8, "Dark2")[3],
               keep.order = TRUE)
upset1

upset2 <-upset(coloc.ep,
               sets = rev(c("lead.eQTL.no_pQTL", "lead.indeQTL", "coloc.onlyeQTL",
                            # "lead.no_eQTL.pQTL", "lead.indpQTL", "coloc.onlypQTL",
                            "all.eQTL.or.pQTL")),
               mb.ratio = c(0.55, 0.45),
               order.by = "freq",
               #group.by = "sets",
               main.bar.color = "#1F78B4",
               keep.order = TRUE)
upset2

upset3 <-upset(coloc.ep,
               sets = rev(c("lead.no_eQTL.pQTL", "lead.indpQTL", "coloc.onlypQTL",
                            "all.eQTL.or.pQTL")),
               mb.ratio = c(0.55, 0.45),
               order.by = "freq",
               #group.by = "sets",
               main.bar.color = "#6A3D9A",
               keep.order = TRUE)
upset3
as.ggplot(upset1$Main_bar)
ggarrange(NULL, as.ggplot(upset1$Main_bar) +
            theme(plot.margin = unit(c(0,0,0,2), "cm")), #trbl
          as.ggplot(upset1$Sizes) +
            theme(plot.margin = unit(c(0,0,0,0), "cm")),
          as.ggplot(upset1$Matrix) +
            theme(plot.margin = unit(c(0,0,0,0), "cm")),
          nrow = 2, ncol = 2,
          align = "hv")

uu_1 <- cowplot::plot_grid(NULL, upset1$Main_bar, upset1$Sizes, upset1$Matrix,
                           nrow=2,
                           #axis = "lr",
                           align='hv',
                           rel_heights = c(3,2),
                           rel_widths = c(2,3))
uu_2 <- cowplot::plot_grid(NULL, upset2$Main_bar, upset2$Sizes, upset2$Matrix,
                           nrow=2, align='hv',
                           rel_heights = c(3,2),
                           rel_widths = c(2,3))
uu_3 <- cowplot::plot_grid(NULL, upset3$Main_bar, upset3$Sizes, upset3$Matrix,
                           nrow=2, align='hv',
                           rel_heights = c(3,2),
                           rel_widths = c(2,3))

grid.arrange(uu_1, uu_2, uu_3,
             nrow = 3)





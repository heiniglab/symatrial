## start here!
setwd("~/work/symAtrial_QTL/scripts/analysis/gwas_imputed")
source("../../helper/helper.R")
local=F

get.duplicates <- function(df, cols){
  df2 <- data.frame(table(df[, cols]))
  df2 <- df2[df2$Freq>1, ]
  colnames(df2) <- c(cols, "Freq")
  df <- merge(df,
              df2[cols],
              all.y=T)
}

gwas.traits <- read.csv(paste0(get.path("gwas3", local),
                               "gwas-efo-trait-mappings_r2019-11-26.tsv"),
                        sep="\t", h=T, quote="", fill=FALSE,
                        stringsAsFactors = F)

cvm = "EFO_0004298"
measures = gwas.traits[grep(cvm, gwas.traits$Parent.URI), c("Disease.trait", "EFO.URI")]

heart.disease = "EFO_0000319" #EFO_0003777 Cardiovascular disease	
diseases = gwas.traits[grep(heart.disease, gwas.traits$Parent.URI), c("Disease.trait", "EFO.URI")]

ids <- rbind(measures, diseases)
heart.traits <- gwas.traits[gwas.traits$EFO.URI %in% ids$EFO.URI,]
write.table(heart.traits,
            file=paste0(get.path("gwas3", local),
                        "GWAS_catalog_trait_mapping_heart.tsv"),
            quote = F, row.names = F, col.names = T, sep="\t")

gwas = read.csv(paste0(get.path("gwas3", local), "gwas-catalog-associations_ontology-annotated_r2019-11-26.tsv"),
                sep="\t", h=T, quote="", fill=FALSE,
                stringsAsFactors = F)
colnames(gwas)
dim(gwas)

heart.gwas = gwas[grep(paste(unique(ids$EFO.URI), collapse = "|"), gwas$MAPPED_TRAIT_URI), ]
dim(heart.gwas)

arrythmia <- c("atrial fibrillation",
               "cardiac arrhythmia",
               "sudden cardiac arrest",
               "supraventricular ectopy",
               "early cardiac repolarization measurement",
               "heart rate",
               "heart rate variability measurement",
               "P wave duration",
               "P wave terminal force measurement",
               "PR interval",
               "PR segment",
               "QRS amplitude",
               "QRS complex",
               "QRS duration",
               "QT interval",
               "R wave amplitude",
               "resting heart rate",
               "RR interval",
               "S wave amplitude",
               "T wave amplitude")

rythm.traits <- heart.traits[heart.traits[, "EFO.term"] %in% arrythmia, ]
AF.traits <- heart.traits[heart.traits[,"EFO.term"] %in% c("atrial fibrillation", "QT interval"), ]
cat(heart.gwas$SNPS, file=paste0(get.path("gwas3", local), "heart.gwas.snps.txt"), sep="\n")

## annotate with LD
library(data.table)

# # gwas snps -> alle LD snps 0.8
if(!file.exists(paste0(get.path("gwas3", local), "snipa_proxies.RDS"))){
  source("snipe.R")
  proxies <- snipa.get.ld.by.snp(heart.gwas[grep("rs.*", heart.gwas$SNPS), "SNPS"],
                                 rsquare=0.8,
                                 population=c('eur'))
  saveRDS(proxies, file=paste0(get.path("gwas3", local), "snipa_proxies.RDS"))
} else {
  proxies <- readRDS(paste0(get.path("gwas3", local), "snipa_proxies.RDS"))
}
colnames(proxies)

snap <- proxies[, c("QRSID", "RSID", "RSALIAS", "DIST", "R2", "DPRIME")]
colnames(snap) <- c("SNP", "Proxy", "Alias", "Distance", "RSquared", "DPrime")

heart.gwas = merge(heart.gwas, snap, by.x="SNPS", by.y="SNP", all.x=T)

## set SNPs that were not found as their own proxies
heart.gwas[is.na(heart.gwas$Proxy),"Proxy"] = heart.gwas[is.na(heart.gwas$Proxy),"SNPS"]
cat(setdiff(unique(heart.gwas$Proxy), c(NA, "")),
    file=paste0(get.path("gwas3", local), "snps-for-eqtls.txt"),
    sep="\n")
write.table(heart.gwas, file=paste0(get.path("gwas3", local),
                                    "heart.gwas.txt"),
            sep="\t", quote=F, row.names=F, col.names = T)

rythm.gwas <- heart.gwas[grep(paste(unique(rythm.traits$EFO.URI), collapse = "|"), heart.gwas$MAPPED_TRAIT_URI), ]
write.table(rythm.gwas, file=paste0(get.path("gwas3", local), "rythm.gwas.txt"),
            sep="\t", quote=F, row.names=F, col.names = T)

AF.gwas <- heart.gwas[grep(paste(unique(AF.traits$EFO.URI), collapse = "|"), heart.gwas$MAPPED_TRAIT_URI), ]
write.table(AF.gwas, file=paste0(get.path("gwas3", local), "af.gwas.txt"),
            sep="\t", quote=F, row.names=F, col.names = T)




source("../../helper/helper.R")
local=F

heart.gwas <- read.csv(paste0(get.path("gwas3", local), "heart.gwas.txt"),
                       sep="\t", stringsAsFactors = F, h = T)
# rythm.gwas <- read.csv(paste0(get.path("gwas3", local), "rythm.gwas.txt"),
#                        sep="\t", stringsAsFactors = F, h = T)

## ratioQTL ----
our.eqtls <- readRDS(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "ratios_right_atrial_appendage_allpairs.RDS",
                                  sep="/")))
our.eqtls <- our.eqtls[our.eqtls$rs_id %in% heart.gwas$Proxy, ]

## find pairs of GWA and eQTL SNPs significant in at least one of the two sets
our.gwa = heart.gwas
## remove proxies for SNPs that were actually tested for eQTLs
tested = intersect(heart.gwas$SNPS, our.eqtls$snpid)
keep = our.gwa$SNPS == our.gwa$Proxy | !our.gwa$SNPS %in% tested
our.gwa = our.gwa[which(keep),]
our.gwa = our.gwa[!is.na(our.gwa$Proxy),]

our.gwa = merge(our.gwa, our.eqtls, by.x="Proxy", by.y="snpid")
our.gwa$GWAS.FDR <- p.adjust(our.gwa$pvalue, method = "BH")
our.gwa.ratios.t <- our.gwa
write.table(our.gwa.ratios.t,
            file=paste0(get.path("gwas3", local), "our.gwa.ratios.temp.txt"),
            row.names = F, col.names = T, quote=F, sep="\t")

# pvalue
our.pairs = our.gwa[which(our.gwa[,"pvalue"] < 1e-5),
                    c("SNPS", "Proxy", "gene", "pvalue")]
## select the best proxy
best <- NULL
best = tapply(1:nrow(our.pairs), paste(our.pairs$SNPS, our.pairs$gene), function(idx) {
  return(idx[which.min(our.pairs$pvalue[idx])])
})
our.pairs = our.pairs[best,]
our.gwa2 <- our.gwa[which(our.gwa[,"pvalue"] < 1e-5), ]
our.gwa2 <- our.gwa2[best,]
our.pairs.ratios <- our.pairs
our.gwa.ratios <- our.gwa2

write.table(our.gwa.ratios, file=paste0(get.path("gwas3", local), "our.gwa.ratios.txt"),
            row.names = F, col.names = T, quote=F, sep="\t")

# FDR  
our.pairs = our.gwa[which(our.gwa[,"FDR"] < 0.05),
                    c("SNPS", "Proxy", "gene", "pvalue")]
## select the best proxy
best <- NULL
best = tapply(1:nrow(our.pairs), paste(our.pairs$SNPS, our.pairs$gene), function(idx) {
  return(idx[which.min(our.pairs$pvalue[idx])])
})
our.pairs = our.pairs[best,]
our.gwa2 <- our.gwa[which(our.gwa[,"FDR"] < 0.05), ]
our.gwa2 <- our.gwa2[best,]

write.table(our.gwa2, file=paste0(get.path("gwas3", local),
                                  "our.gwa.ratios.FDR.txt"),
            row.names = F, col.names = T, quote=F, sep="\t")

# GWAS FDR  
our.pairs = our.gwa[which(our.gwa[,"GWAS.FDR"] < 0.05),
                    c("SNPS", "Proxy", "gene", "pvalue")]
## select the best proxy
best <- NULL
best = tapply(1:nrow(our.pairs), paste(our.pairs$SNPS, our.pairs$gene), function(idx) {
  return(idx[which.min(our.pairs$pvalue[idx])])
})
our.pairs = our.pairs[best,]
our.gwa2 <- our.gwa[which(our.gwa[,"GWAS.FDR"] < 0.05), ]
our.gwa2 <- our.gwa2[best,]

write.table(our.gwa2, file=paste0(get.path("gwas3", local),
                                  "our.gwa.ratios.GWAS.FDR.txt"),
            row.names = F, col.names = T, quote=F, sep="\t")

rm(our.eqtls)


# GWAS enrichment on AFHRI-B data ----
setwd("~/work/symAtrial_QTL/scripts/analysis/gwas_imputed")
source("../../helper/helper.R")
local=F

heart.gwas <- read.csv(paste0(get.path("gwas3", local),
                              "heart.gwas.txt"),
                       sep="\t", h = T, stringsAsFactors = F)

rythm.gwas <- read.csv(paste0(get.path("gwas3", local),
                              "rythm.gwas.txt"),
                       sep="\t", h = T, stringsAsFactors = F)

AF.gwas <- read.csv(paste0(get.path("gwas3", local),
                           "af.gwas.txt"),
                    sep="\t", h = T, stringsAsFactors = F)

heart.traits <- read.csv(file=paste0(get.path("gwas3", local),
                                     "GWAS_catalog_trait_mapping_heart.tsv"),
                         h = T, sep="\t", stringsAsFactors = F)

arrythmia <- c("atrial fibrillation",
               "cardiac arrhythmia",
               "sudden cardiac arrest",
               "supraventricular ectopy",
               "early cardiac repolarization measurement",
               "heart rate",
               "heart rate variability measurement",
               "P wave duration",
               "P wave terminal force measurement",
               "PR interval",
               "PR segment",
               "QRS amplitude",
               "QRS complex",
               "QRS duration",
               "QT interval",
               "R wave amplitude",
               "resting heart rate",
               "RR interval",
               "S wave amplitude",
               "T wave amplitude")

rythm.traits <- heart.traits[heart.traits[, "EFO.term"] %in% arrythmia, ]
AF.traits <- heart.traits[heart.traits[,"EFO.term"] %in% c("atrial fibrillation", "QT interval"), ]

if (F){
  QTL.snp <- readRDS(file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_2.RDS"))
  colnames(QTL.snp)
  our.qtls <- QTL.snp[QTL.snp$rs_id %in% heart.gwas$Proxy, ]
  saveRDS(our.qtls, file=paste0(get.path("gwas2", local), "our-qtls-heart.RDS"))
}

## enrichment for eQTLs ----
our.qtls <- readRDS(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "eQTL_right_atrial_appendage_allpairs.RDS",
                                  sep="/")))
our.gwa.qtl <- read.csv(paste0(get.path("gwas3", local),
                                "our.gwa.eqtl.temp.txt"),
                         h = T, stringsAsFactors=F, sep="\t")

our.gwa.qtl.heart <- our.gwa.qtl
saveRDS(our.gwa.qtl.heart, paste0(get.path("gwas3", local), "our_gwa_eqtl_heart.RDS"))

our.gwa.qtl.rythm <- our.gwa.qtl[grep(paste(unique(rythm.traits$EFO.URI), collapse = "|"),
                                        our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.rythm, paste0(get.path("gwas3", local), "our_gwa_eqtl_arrythmia.RDS"))

our.gwa.qtl.AF <- our.gwa.qtl[grep(paste(unique(AF.traits$EFO.URI), collapse = "|"),
                                     our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.AF, paste0(get.path("gwas3", local), "our_gwa_eqtl_AF.RDS"))

## also reconstruct a contingency table to test for enrichment of gwa hits
tested.snps = unique(our.qtls$rs_id)

# our.gwa.eqtl <- read.csv(paste0(get.path("gwas2", local), "our.gwa.eqtl.txt"),
#             h=T, sep="\t")
sig.qtls = read.csv(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "eQTL_right_atrial_appendage_allpairs.significant.txt",
                                  sep="/")),
                     sep = "\t", stringsAsFactors = F, h = T)

tab.heart = table(gwa.proxy=tested.snps %in% our.gwa.qtl.heart$Proxy, eqtl=tested.snps %in% sig.qtls$rs_id)
tab.heart

tab.rythm = table(gwa.proxy=tested.snps %in% our.gwa.qtl.rythm$Proxy, eqtl=tested.snps %in% sig.qtls$rs_id)
tab.rythm

tab.AF = table(gwa.proxy=tested.snps %in% our.gwa.qtl.AF$Proxy, eqtl=tested.snps %in% sig.qtls$rs_id)
tab.AF

sink(file="new/gwa-snps-vs-sign-eqtls_heart.txt")
cat("Observed:\n")
print(tab.heart)
cat("Expected:\n")
print(chisq.test(tab.heart)$expected)
print(fisher.test(tab.heart))
sink()

sink(file="new/gwa-snps-vs-sign-eqtls_arrythmia.txt")
cat("Observed:\n")
print(tab.rythm)
cat("Expected:\n")
print(chisq.test(tab.rythm)$expected)
print(fisher.test(tab.rythm))
sink()

sink(file="new/gwa-snps-vs-sign-eqtls_AF.txt")
cat("Observed:\n")
print(tab.AF)
cat("Expected:\n")
print(chisq.test(tab.AF)$expected)
print(fisher.test(tab.AF))
sink()

## enrichment for pQTLs ----
our.qtls <- readRDS(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "pQTL_right_atrial_appendage_allpairs.RDS",
                                  sep="/")))
our.gwa.qtl <- read.csv(paste0(get.path("gwas3", local), "our.gwa.pqtl.temp.txt"),
                        h = T, stringsAsFactors=F, sep="\t")


our.gwa.qtl.heart <- our.gwa.qtl
saveRDS(our.gwa.qtl.heart, paste0(get.path("gwas3", local), "our_gwa_pqtl_heart.RDS"))

our.gwa.qtl.rythm <- our.gwa.qtl[grep(paste(unique(rythm.traits$EFO.URI), collapse = "|"),
                                      our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.rythm, paste0(get.path("gwas3", local), "our_gwa_pqtl_arrythmia.RDS"))

our.gwa.qtl.AF <- our.gwa.qtl[grep(paste(unique(AF.traits$EFO.URI), collapse = "|"),
                                   our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.AF, paste0(get.path("gwas3", local), "our_gwa_pqtl_AF.RDS"))

## also reconstruct a contingency table to test for enrichment of gwa hits
tested.snps = unique(our.qtls$rs_id)

sig.qtls = read.csv(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "pQTL_right_atrial_appendage_allpairs.significant.txt",
                                  sep="/")),
                     sep = "\t", stringsAsFactors = F, h = T)

tab.heart = table(gwa.proxy=tested.snps %in% our.gwa.qtl.heart$Proxy, pqtl=tested.snps %in% sig.qtls$rs_id)
tab.heart

tab.rythm = table(gwa.proxy=tested.snps %in% our.gwa.qtl.rythm$Proxy, pqtl=tested.snps %in% sig.qtls$rs_id)
tab.rythm

tab.AF = table(gwa.proxy=tested.snps %in% our.gwa.qtl.AF$Proxy, pqtl=tested.snps %in% sig.qtls$rs_id)
tab.AF

sink(file="new/gwa-snps-vs-sign-pqtls_heart.txt")
cat("Observed:\n")
print(tab.heart)
cat("Expected:\n")
print(chisq.test(tab.heart)$expected)
print(fisher.test(tab.heart))
sink()

sink(file="new/gwa-snps-vs-sign-pqtls_arrythmia.txt")
cat("Observed:\n")
print(tab.rythm)
cat("Expected:\n")
print(chisq.test(tab.rythm)$expected)
print(fisher.test(tab.rythm))
sink()

sink(file="new/gwa-snps-vs-sign-pqtls_AF.txt")
cat("Observed:\n")
print(tab.AF)
cat("Expected:\n")
print(chisq.test(tab.AF)$expected)
print(fisher.test(tab.AF))
sink()


## enrichment for res eQTLs ----
our.qtls <- readRDS(paste0(get.path("results", local),
                           paste("imputed", "cis", "final",
                                 "res_eQTL_right_atrial_appendage_allpairs.RDS",
                                 sep="/")))
our.gwa.qtl <- read.csv(paste0(get.path("gwas3", local),
                               "our.gwa.res_eqtl.temp.txt"),
                        h = T, stringsAsFactors=F, sep="\t")

our.gwa.qtl.heart <- our.gwa.qtl
saveRDS(our.gwa.qtl.heart, paste0(get.path("gwas3", local), "our_gwa_res_eqtl_heart.RDS"))

our.gwa.qtl.rythm <- our.gwa.qtl[grep(paste(unique(rythm.traits$EFO.URI), collapse = "|"),
                                      our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.rythm, paste0(get.path("gwas3", local), "our_gwa_res_eqtl_arrythmia.RDS"))

our.gwa.qtl.AF <- our.gwa.qtl[grep(paste(unique(AF.traits$EFO.URI), collapse = "|"),
                                   our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.AF, paste0(get.path("gwas3", local), "our_gwa_res_eqtl_AF.RDS"))

## also reconstruct a contingency table to test for enrichment of gwa hits
tested.snps = unique(our.qtls$rs_id)

sig.qtls = read.csv(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "res_eQTL_right_atrial_appendage_allpairs.significant.txt",
                                  sep="/")),
                     sep = "\t", stringsAsFactors = F, h = T)

tab.heart = table(gwa.proxy=tested.snps %in% our.gwa.qtl.heart$Proxy, res_qtl=tested.snps %in% sig.qtls$rs_id)
tab.heart

tab.rythm = table(gwa.proxy=tested.snps %in% our.gwa.qtl.rythm$Proxy, res_eqtl=tested.snps %in% sig.qtls$rs_id)
tab.rythm

tab.AF = table(gwa.proxy=tested.snps %in% our.gwa.qtl.AF$Proxy, res_eqtl=tested.snps %in% sig.qtls$rs_id)
tab.AF

sink(file="new/gwa-snps-vs-sign-res_eqtl_heart.txt")
cat("Observed:\n")
print(tab.heart)
cat("Expected:\n")
print(chisq.test(tab.heart)$expected)
print(fisher.test(tab.heart))
sink()

sink(file="new/gwa-snps-vs-sign-res_eqtl_arrythmia.txt")
cat("Observed:\n")
print(tab.rythm)
cat("Expected:\n")
print(chisq.test(tab.rythm)$expected)
print(fisher.test(tab.rythm))
sink()

sink(file="new/gwa-snps-vs-sign-res_eqtl_AF.txt")
cat("Observed:\n")
print(tab.AF)
cat("Expected:\n")
print(chisq.test(tab.AF)$expected)
print(fisher.test(tab.AF))
sink()



## enrichment for res pQTLs ----
our.qtls <- readRDS(paste0(get.path("results", local),
                           paste("imputed", "cis", "final",
                                 "res_pQTL_right_atrial_appendage_allpairs.RDS",
                                 sep="/")))
our.gwa.qtl <- read.csv(paste0(get.path("gwas3", local),
                               "our.gwa.res_pqtl.temp.txt"),
                        h = T, stringsAsFactors=F, sep="\t")

our.gwa.qtl.heart <- our.gwa.qtl
saveRDS(our.gwa.qtl.heart, paste0(get.path("gwas3", local), "our_gwa_res_pqtl_heart.RDS"))

our.gwa.qtl.rythm <- our.gwa.qtl[grep(paste(unique(rythm.traits$EFO.URI), collapse = "|"),
                                      our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.rythm, paste0(get.path("gwas3", local), "our_gwa_res_pqtl_arrythmia.RDS"))

our.gwa.qtl.AF <- our.gwa.qtl[grep(paste(unique(AF.traits$EFO.URI), collapse = "|"),
                                   our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.AF, paste0(get.path("gwas3", local), "our_gwa_res_pqtl_AF.RDS"))

## also reconstruct a contingency table to test for enrichment of gwa hits
tested.snps = unique(our.qtls$rs_id)

sig.qtls = read.csv(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "res_pQTL_right_atrial_appendage_allpairs.significant.txt",
                                  sep="/")),
                     sep = "\t", stringsAsFactors = F, h = T)

tab.heart = table(gwa.proxy=tested.snps %in% our.gwa.qtl.heart$Proxy, res_pqtl=tested.snps %in% sig.qtls$rs_id)
tab.heart

tab.rythm = table(gwa.proxy=tested.snps %in% our.gwa.qtl.rythm$Proxy, res_pqtl=tested.snps %in% sig.qtls$rs_id)
tab.rythm

tab.AF = table(gwa.proxy=tested.snps %in% our.gwa.qtl.AF$Proxy, res_pqtl=tested.snps %in% sig.qtls$rs_id)
tab.AF

sink(file="new/gwa-snps-vs-sign-res_pqtl_heart.txt")
cat("Observed:\n")
print(tab.heart)
cat("Expected:\n")
print(chisq.test(tab.heart)$expected)
print(fisher.test(tab.heart))
sink()

sink(file="new/gwa-snps-vs-sign-res_pqtl_arrythmia.txt")
cat("Observed:\n")
print(tab.rythm)
cat("Expected:\n")
print(chisq.test(tab.rythm)$expected)
print(fisher.test(tab.rythm))
sink()

sink(file="new/gwa-snps-vs-sign-res_pqtl_AF.txt")
cat("Observed:\n")
print(tab.AF)
cat("Expected:\n")
print(chisq.test(tab.AF)$expected)
print(fisher.test(tab.AF))
sink()


## enrichment for ratio QTLs ----
our.qtls <- readRDS(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "ratios_right_atrial_appendage_allpairs.RDS",
                                  sep="/")))
our.gwa.qtl <- read.csv(paste0(get.path("gwas3", local), "our.gwa.ratios.temp.txt"),
                         h = T, stringsAsFactors=F, sep="\t")

our.gwa.qtl.heart <- our.gwa.qtl
saveRDS(our.gwa.eqtl.heart, paste0(get.path("gwas3", local), "our_gwa_ratio_QTL_heart.RDS"))

our.gwa.qtl.rythm <- our.gwa.qtl[grep(paste(unique(rythm.traits$EFO.URI), collapse = "|"),
                                      our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.rythm, paste0(get.path("gwas3", local), "our_gwa_ratio_QTL_arrythmia.RDS"))

our.gwa.qtl.AF <- our.gwa.qtl[grep(paste(unique(AF.traits$EFO.URI), collapse = "|"),
                                   our.gwa.qtl$MAPPED_TRAIT_URI), ]
saveRDS(our.gwa.qtl.AF, paste0(get.path("gwas3", local), "our_gwa_ratio_QTL_AF.RDS"))

## also reconstruct a contingency table to test for enrichment of gwa hits
tested.snps = unique(our.eqtls$rs_id)

# our.gwa.eqtl <- read.csv(paste0(get.path("gwas2", local), "our.gwa.eqtl.txt"),
#             h=T, sep="\t")
sig.eqtls = read.csv(paste0(get.path("results", local),
                            paste("imputed", "cis", "final",
                                  "ratios_right_atrial_appendage_allpairs.significant.txt",
                                  sep="/")),
                     sep = "\t", stringsAsFactors = F, h = T)

tab.heart = table(gwa.proxy=tested.snps %in% our.gwa.qtl.heart$Proxy, ratio_qtl=tested.snps %in% sig.qtls$rs_id)
tab.heart

tab.rythm = table(gwa.proxy=tested.snps %in% our.gwa.qtl.rythm$Proxy, ratio_qtl=tested.snps %in% sig.qtls$rs_id)
tab.rythm

tab.AF = table(gwa.proxy=tested.snps %in% our.gwa.qtl.AF$Proxy, ratio_qtl=tested.snps %in% sig.qtls$rs_id)
tab.AF

sink(file="new/gwa-snps-vs-sign-ratio_qtl_heart.txt")
cat("Observed:\n")
print(tab.heart)
cat("Expected:\n")
print(chisq.test(tab.heart)$expected)
print(fisher.test(tab.heart))
sink()

sink(file="new/gwa-snps-vs-sign-ratio_qtl_arrythmia.txt")
cat("Observed:\n")
print(tab.rythm)
cat("Expected:\n")
print(chisq.test(tab.rythm)$expected)
print(fisher.test(tab.rythm))
sink()

sink(file="new/gwa-snps-vs-sign-ratio_qtl_AF.txt")
cat("Observed:\n")
print(tab.AF)
cat("Expected:\n")
print(chisq.test(tab.AF)$expected)
print(fisher.test(tab.AF))
sink()


# Rheumatoid arthritis GWAS for revision ----
if(T){
  source("../../helper/helper.R")
  local=F
  
  gwas.traits <- read.csv(paste0(get.path("gwas3", local),
                                 "gwas-efo-trait-mappings_r2019-11-26.tsv"),
                          sep="\t", h=T, quote="", fill=FALSE,
                          stringsAsFactors = F)
  
  gwas = read.csv(paste0(get.path("gwas3", local), "gwas-catalog-associations_ontology-annotated_r2019-11-26.tsv"),
                  sep="\t", h=T, quote="", fill=FALSE,
                  stringsAsFactors = F)
  
  # Rheumatoid arthritis
  RA <- "EFO_0000685"
  RA.traits <- gwas.traits[grep(RA, gwas.traits$EFO.URI),
                           c("Disease.trait", "EFO.URI")]
  RA.gwas <- gwas[grep(paste(unique(RA.traits$EFO.URI),
                             collapse = "|"),
                       gwas$MAPPED_TRAIT_URI), ]
  
  source("snipe.R")
  RA.snps <- RA.gwas[grep("rs.*", RA.gwas$SNPS), "SNPS"]
  #test <- data.frame(RA.snps)
  proxies.RA <- snipa.get.ld.by.snp(RA.snps,
                                    rsquare=0.8,
                                    population=c('eur'))
  saveRDS(proxies.RA, file=paste0(get.path("gwas3", local),
                                  "snipa_proxies_RA.RDS"))
  
  snap <- proxies.RA[, c("QRSID", "RSID", "RSALIAS", "DIST", "R2", "DPRIME")]
  colnames(snap) <- c("SNP", "Proxy", "Alias", "Distance", "RSquared", "DPrime")
  
  RA.gwas = merge(RA.gwas, snap,
                  by.x="SNPS", by.y="SNP", all.x=T)
  
  ## set SNPs that were not found as their own proxies
  RA.gwas[is.na(RA.gwas$Proxy),"Proxy"] = RA.gwas[is.na(RA.gwas$Proxy),"SNPS"]
  cat(setdiff(unique(RA.gwas$Proxy), c(NA, "")),
      file=paste0(get.path("gwas3", local), "snps-for-eqtls-RA.txt"),
      sep="\n")
  write.table(RA.gwas, file=paste0(get.path("gwas3", local),
                                   "RA.gwas.txt"),
              sep="\t", quote=F, row.names=F, col.names = T)
  
  ## eQTL ----
  our.eqtls <- readRDS(paste0(get.path("results", local),
                              paste("imputed", "cis", "final",
                                    "eQTL_right_atrial_appendage_allpairs.RDS",
                                    sep="/")))
  our.eqtls <- our.eqtls[our.eqtls$rs_id %in% RA.gwas$Proxy, ]
  
  #library(plyr)
  # colnames(our.eqtls)
  # [1] "gene"               "snpid"              "chr"                "variant_pos"       
  # [5] "Allele1"            "Allele2"            "variant_id"         "rs_id"             
  # [9] "statistic"          "pvalue"             "FDR"                "beta"              
  # [13] "gtex.variant_id"    "gtex.match_alleles" "gene_id"    
  
  ## find pairs of GWA and eQTL SNPs significant in at least one of the two sets
  our.gwa = RA.gwas
  ## remove proxies for SNPs that were actually tested for eQTLs
  tested = intersect(RA.gwas$SNPS, our.eqtls$snpid)
  keep = our.gwa$SNPS == our.gwa$Proxy | !our.gwa$SNPS %in% tested
  our.gwa = our.gwa[which(keep),]
  our.gwa = our.gwa[!is.na(our.gwa$Proxy),]
  
  our.gwa = merge(our.gwa, our.eqtls, by.x="Proxy", by.y="snpid")
  our.gwa.eqtl.t <- our.gwa
  write.table(our.gwa.eqtl.t,
              file=paste0(get.path("gwas3", local),
                          "our.gwa.RA.eqtl.temp.txt"),
              row.names = F, col.names = T, quote=F, sep="\t")
  
  # FDR  
  our.pairs = our.gwa[which(our.gwa[,"FDR"] < 0.05),
                      c("SNPS", "Proxy", "gene", "pvalue")]
  ## select the best proxy
  best <- NULL
  try(
    best = tapply(1:nrow(our.pairs),
                  paste(our.pairs$SNPS, our.pairs$gene),
                  function(idx) {
                    return(idx[which.min(our.pairs$pvalue[idx])])
                  })
  )
  our.pairs = our.pairs[best,]
  our.gwa2 <- our.gwa[which(our.gwa[,"FDR"] < 0.05), ]
  our.gwa2 <- our.gwa2[best,]
  
  write.table(our.gwa2, file=paste0(get.path("gwas3", local),
                                    "our.gwa.RA.eqtl.FDR.txt"),
              row.names = F, col.names = T, quote=F, sep="\t")
  
  rm(our.eqtls)
  
  
  
  
  ## pQTL ----
  our.pqtls <- readRDS(paste0(get.path("results", local),
                              paste("imputed", "cis", "final",
                                    "pQTL_right_atrial_appendage_allpairs.RDS",
                                    sep="/")))
  our.pqtls <- our.pqtls[our.pqtls$rs_id %in% RA.gwas$Proxy, ]
  
  ## find pairs of GWA and eQTL SNPs significant in at least one of the two sets
  our.gwa = RA.gwas
  ## remove proxies for SNPs that were actually tested for eQTLs
  tested = intersect(RA.gwas$SNPS, our.pqtls$snpid)
  keep = our.gwa$SNPS == our.gwa$Proxy | !our.gwa$SNPS %in% tested
  our.gwa = our.gwa[which(keep),]
  our.gwa = our.gwa[!is.na(our.gwa$Proxy),]
  
  our.gwa = merge(our.gwa, our.pqtls, by.x="Proxy", by.y="snpid")
  our.gwa.pqtl.t <- our.gwa
  write.table(our.gwa.pqtl.t, file=paste0(get.path("gwas3", local),
                                          "our.gwa.RA.pqtl.temp.txt"),
              row.names = F, col.names = T, quote=F, sep="\t")
  
  # FDR  
  our.pairs = our.gwa[which(our.gwa[,"FDR"] < 0.05),
                      c("SNPS", "Proxy", "gene", "pvalue")]
  ## select the best proxy
  best <- NULL
  try(
    best = tapply(1:nrow(our.pairs),
                  paste(our.pairs$SNPS, our.pairs$gene),
                  function(idx) {
                    return(idx[which.min(our.pairs$pvalue[idx])])
                  })
  )
  
  our.pairs = our.pairs[best,]
  our.gwa2 <- our.gwa[which(our.gwa[,"FDR"] < 0.05), ]
  our.gwa2 <- our.gwa2[best,]
  
  write.table(our.gwa2, file=paste0(get.path("gwas3", local),
                                    "our.gwa.RA.pqtl.FDR.txt"),
              row.names = F, col.names = T, quote=F, sep="\t")
  
  rm(our.pqtls)
  
  
  ## enrichment for eQTLs ----
  our.qtls <- readRDS(paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "eQTL_right_atrial_appendage_allpairs.RDS",
                                   sep="/")))
  our.gwa.qtl <- read.csv(paste0(get.path("gwas3", local),
                                 "our.gwa.RA.eqtl.temp.txt"),
                          h = T, stringsAsFactors=F, sep="\t")
  
  our.gwa.RA.qtl <- our.gwa.qtl
  saveRDS(our.gwa.RA.qtl,
          paste0(get.path("gwas3", local),
                 "our_gwa_RA_eqtl.RDS"))
  
  ## also reconstrRAt a contingency table to test for enrichment of gwa hits
  tested.snps = unique(our.qtls$rs_id)
  
  # our.gwa.eqtl <- read.csv(paste0(get.path("gwas2", local), "our.gwa.eqtl.txt"),
  #             h=T, sep="\t")
  sig.qtls = read.csv(paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "eQTL_right_atrial_appendage_allpairs.significant.txt",
                                   sep="/")),
                      sep = "\t", stringsAsFactors = F, h = T)
  
  tab.RA = table(gwa.proxy=tested.snps %in% our.gwa.RA.qtl$Proxy,
                 eqtl=tested.snps %in% sig.qtls$rs_id)
  tab.RA
  
  sink(file="new/gwa-RA-snps-vs-sign-eqtls.txt")
  cat("Observed:\n")
  print(tab.RA)
  cat("Expected:\n")
  print(chisq.test(tab.RA)$expected)
  print(fisher.test(tab.RA))
  sink()
  
  cat("Observed:\n")
  print(tab.RA)
  cat("Expected:\n")
  print(chisq.test(tab.RA)$expected)
  print(fisher.test(tab.RA))
  
  ## enrichment for pQTLs ----
  our.qtls <- readRDS(paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "pQTL_right_atrial_appendage_allpairs.RDS",
                                   sep="/")))
  our.gwa.qtl <- read.csv(paste0(get.path("gwas3", local), "our.gwa.RA.pqtl.temp.txt"),
                          h = T, stringsAsFactors=F, sep="\t")
  
  our.gwa.RA.qtl <- our.gwa.qtl
  saveRDS(our.gwa.RA.qtl,
          paste0(get.path("gwas3", local),
                 "our_gwa_RA_pqtl.RDS"))
  
  ## also reconstrRAt a contingency table to test for enrichment of gwa hits
  tested.snps = unique(our.qtls$rs_id)
  
  sig.qtls = read.csv(paste0(get.path("results", local),
                             paste("imputed", "cis", "final",
                                   "pQTL_right_atrial_appendage_allpairs.significant.txt",
                                   sep="/")),
                      sep = "\t", stringsAsFactors = F, h = T)
  
  tab.RA = table(gwa.proxy=tested.snps %in% our.gwa.RA.qtl$Proxy,
                 pqtl=tested.snps %in% sig.qtls$rs_id)
  tab.RA
  
  sink(file="new/gwa-RA-snps-vs-sign-pqtls.txt")
  cat("Observed:\n")
  print(tab.RA)
  cat("Expected:\n")
  print(chisq.test(tab.RA)$expected)
  print(fisher.test(tab.RA))
  sink()
  
  cat("Observed:\n")
  print(tab.RA)
  cat("Expected:\n")
  print(chisq.test(tab.RA)$expected)
  print(fisher.test(tab.RA))
}


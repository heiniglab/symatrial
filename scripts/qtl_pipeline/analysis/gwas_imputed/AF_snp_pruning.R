snps.t <- read.table(paste0(get.path("results", local),
                            "imputed/trans/AF_GWAS_catalogue_snps.txt"),
                     stringsAsFactors = F)
snps.p <- read.table(paste0(get.path("results", local),
                            "imputed/trans/AF_GWAS_catalogue_snps_prot.txt"),
                     stringsAsFactors = F)
af.snps <- unique(c(snps.t$V1, snps.p$V1))



snplocs <- read.csv(get.path("snplocs imputed", local),
                    sep = "\t", h=T, stringsAsFactors = F)
colnames(snplocs) <- c("snps", "chr", "pos")

snplocs.af <- snplocs[snplocs$snps %in% af.snps, ]
snplocs.af$rsID <- gsub("\\:.*", "", snplocs.af$snps)
write.table(snplocs.af$rsID,
            file = paste0(get.path("gwas2", local), "AF_snps_unpruned.txt"),
            row.names = F, col.names = F, sep = "\n", quote = F)

gwas <- readRDS(paste0(get.path("gwas2", local), "Roselli2018_AF_HRC_GWAS_ALLv11.RDS"))
gwas.af <- gwas[gwas$MarkerName %in% snplocs.af$rsID, ]
gwas.af <- gwas.af[order(gwas.af$P.value), ]

LD <- read.table(paste0(get.path("gwas2", local),
                        "proxySearch.results/proxySearch.results.csv"),
                 sep = "\t", h = T, stringsAsFactors = F)

gwas.af.LD <- merge(gwas.af, LD,
                    by.x="MarkerName", by.y="QRSID",
                    all.x=T, sort=F)
#saveRDS(gwas.af.LD, file="gwas.af.LD.RDS")

gwas.af.LD <- gwas.af.LD[order(gwas.af.LD$P.value, gwas.af.LD$R2), ]
gwas.af.prune <- gwas.af.LD
i=0
gwas.af.prune <- gwas.af.prune[gwas.af.prune$MAF>=0.1, ]
Nmax <- dim(gwas.af.prune)[1]
while(i < dim(gwas.af.prune)[1] & i < Nmax){
  i <- i+1
  sentinel <- gwas.af.prune$MarkerName[i]
  print(sentinel)
  proxies <- LD[LD$QRSID==sentinel, "RSID"]
  proxies <- proxies[!proxies==sentinel]
  print(proxies)
  gwas.af.prune <- gwas.af.prune[!(gwas.af.prune$MarkerName %in% proxies), ]
}
gwas.af.prune <- gwas.af.prune[!duplicated(gwas.af.prune$MarkerName), ]
test <- snplocs.af[snplocs.af$rsID %in% gwas.af.prune$MarkerName, ]

write.table(snplocs.af[snplocs.af$rsID %in% gwas.af.prune$MarkerName, "snps"],
            file = paste0(get.path("gwas2", local), "AF_snps_pruned.txt"),
            row.names = F, col.names = F, sep = "\n", quote = F)
write.table(snplocs.af[snplocs.af$rsID %in% gwas.af.prune$MarkerName, "rsID"],
            file = paste0("AF_snps_pruned.txt"),
            row.names = F, col.names = F, sep = "\n", quote = F)

af.proxies <- list(rs9481842.G = "rs9481842-G",
                   rs11658168.A = c("rs9675122-C",
                                    "rs8073937-G"),
                   rs11588763.T = c("rs34292822-C",
                                    "rs13376333-T",
                                    "rs6666258-C"))
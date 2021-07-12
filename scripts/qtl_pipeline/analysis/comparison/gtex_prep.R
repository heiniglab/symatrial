##
#  prepares data, to make the GTEx data comparable to our analyses
#
#

# first, we need to match the ids
##
#  reduces GTEx data to only contain SNPs, we analyzed in our study
#  (match SNPs)
#
#
#

path <- "/home/icb/ines.assum/projects/symAtrial_QTL/scripts"
setwd(path)
source("helper/helper.R")
local=F

get.duplicates <- function(df, cols){
  df2 <- data.frame(table(df[, cols]))
  df2 <- df2[df2$Freq>1, ]
  colnames(df2) <- c(cols, "Freq")
  df <- merge(df,
              df2[cols],
              all.y=T)
}

snps <- read.table(paste0(get.path("genotype", local),
                          "AFHRI-B_imputed.bim"),
                   header=F, sep = "\t", stringsAsFactors = F)
snps$V3 <-NULL
colnames(snps) <- c("chr", "snpid", "variant_pos", "alt", "ref")
snps$variant_id <- paste(snps$chr, snps$variant_pos,
                         snps$ref, snps$alt, "b37", sep="_")
snps$variant_id2 <- paste(snps$chr, snps$variant_pos,
                          snps$alt, snps$ref, "b37", sep="_")
snps$rsid <- NA
snps$rsid[grep("rs", snps$snpid)] <- gsub("\\:.*", "", snps$snpid[grep("rs", snps$snpid)])

gtex.snps <- read.table(paste0(get.path("gtex", local),
                               "GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt"),
                        header=T, sep = "\t", stringsAsFactors = F)

sum(gtex.snps$variant_id %in% snps$variant_id)
# 4240961

sum(gtex.snps$variant_id %in% snps$variant_id2)
# 1520824

sum(gtex.snps$rs_id_dbSNP147_GRCh37p13 %in% snps$rsid)
# 5514965


snps$same_id <- as.numeric(snps$variant_id %in% gtex.snps$variant_id)
snps$rev_id <- as.numeric(snps$variant_id2 %in% gtex.snps$variant_id)
snps$match_rsid <- as.numeric(snps$rsid %in% gtex.snps$rs_id_dbSNP147_GRCh37p13)
snps$match_rsid[is.na(snps$rsid)] <- NA
snps$match_alleles <- NA
snps$match_alleles[snps$rev_id==1] <- -1
snps$match_alleles[snps$same_id==1] <- 1

table(is.na(snps$match_alleles), is.na(snps$match_rsid))
#       FALSE    TRUE
# FALSE 5751117   10288
# TRUE   816896  104287

dim(snps)[1]
print("104 287 out of 6 682 588 could not be identified by either position or rs-id")
print("816 896 could not be identified by rs-id, but not by their position")
print("5 751 117 can be identified by rs-id and position, but we have to see if the two match")

print("Do we have duplicated SNPs?")
table(snps$same_id, snps$rev_id)
test <- snps[snps$same_id==1 & snps$rev_id==1, ]
test2 <- get.duplicates(test, c("chr", "variant_pos"))

table(snps$match_alleles, snps$match_rsid)
#          0       1
# -1   73356 1446042
# 1   166737 4064982

# snps 6 682 588
snps2 <- merge(snps, gtex.snps,
               all.x = T)
saveRDS(snps2,
        file=paste0(get.path("gtex", local),
                    "tmp/match_ids1.RDS"))
saveRDS(gtex.snps,
        file=paste0(get.path("gtex", local),
                    "tmp/gtex_snps.RDS"))

snps3 <- merge(snps2, gtex.snps,
               by.x=c("chr", "variant_pos", "alt", "ref", "variant_id2"),
               by.y=c("chr", "variant_pos", "ref", "alt", "variant_id"),
               prefix.y=gtex.2, all.x = T)
snps3$num_alt_per_site.x <- NULL
snps3$num_alt_per_site.y <- NULL

saveRDS(snps3,
        file=paste0(get.path("gtex", local),
                    "tmp/match_ids2.RDS"))
# snps3 <- readRDS(paste0(get.path("gtex", local),
#                         "tmp/match_ids2.RDS"))

snps4 <- merge(snps3[which((snps3$match_rsid==1) & is.na(snps3$match_alleles)), ], gtex.snps,   #[!is.na(snps3$rsid), ]
               by.x=c("rsid", "chr", "variant_pos"),
               by.y=c("rs_id_dbSNP147_GRCh37p13", "chr", "variant_pos"),
               all.x = T)
saveRDS(snps4,
        file=paste0(get.path("gtex", local),
                    "tmp/match_ids3.RDS"))
snps4$match_alleles_rsid <- NA
snps4$match_alleles_rsid[which(snps4$ref.x==snps4$ref.y &
                                 snps4$alt.x==snps4$alt.y)] <- -1
snps4$match_alleles_rsid[which(snps4$alt.x==snps4$ref.y &
                                 snps4$ref.x==snps4$alt.y)] <- 1
table(is.na(snps4$match_alleles_rsid))

snps3$gtexid <- snps3$variant_id
saveRDS(snps3[complete.cases(snps3[, c("snpid", "gtexid", "match_alleles")]), c("snpid", "gtexid", "match_alleles")],
        file=paste0(get.path("gtex", local),
                    "tmp/match_gtex_AFHRI_B_SNPs.RDS"))
write.table(data.frame(variant_id=snps3[complete.cases(snps3[, c("snpid", "gtexid", "match_alleles")]), "gtexid"],
                       stringsAsFactors = F),
            file=paste0(get.path("gtex", local),
                        "gtex_snps_in_afhrib.txt"),
            row.names = F, col.names = T, quote = F)

# snps.write <-readRDS(paste0(get.path("gtex", local),
#                     "tmp/match_gtex_AFHRI_B_SNPs.RDS"))
# write.table(data.frame(variant_id=snps.write$gtexid, stringsAsFactors = F),
#             file=paste0(get.path("gtex", local),
#                         "gtex_snps_in_afhrib.txt"),
#             row.names = F, col.names = T, quote = F)


texpr <- readRDS(paste0(get.path("dataset", local),
                        "AFHRI_B_transcriptomics_QC_symbol.RDS"))
pexpr <- readRDS(paste0(get.path("dataset", local),
                        "AFHRI_B_proteomics_QC_symbol.RDS"))
genes <- unique(rownames(texpr), rownames(pexpr))


gtex.genes1a <- read.csv(paste0(get.path("gtex", local),
                                "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"),
                         header = T, stringsAsFactors = F, skip = 2, sep = "\t")[, c("gene_id", "Description")]
gtex.genes1a$gene_id_sub <- gsub("\\..*", "", gtex.genes1a$gene_id)
gtex.genes1b <- read.csv(paste0(get.path("gtex", local),
                                "Heart_Atrial_Appendage.v7.normalized_expression.bed"),
                         header = T, stringsAsFactors = F, sep = "\t")[, c("gene_id", "gene_id")]
colnames(gtex.genes1b) <- c("gene_id", "gene_id_sub")
gtex.genes1b$gene_id_sub <- gsub("\\..*", "", gtex.genes1b$gene_id)
gtex.genes1 <- merge(gtex.genes1a, gtex.genes1b,
                     all = T)

gtex.genes2a <- read.csv(paste0("~/work/symAtrial_QTL/data/current//gtex/v8/",
                                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"),
                         header = T, stringsAsFactors = F, skip = 2, sep = "\t")[, c("Name", "Description")]
gtex.genes2a$gene_id <- gtex.genes2a$Name
gtex.genes2a$gene_id_sub <- gsub("\\..*", "", gtex.genes2a$Name)
gtex.genes2a$Name <- NULL
gtex.genes2b <- read.csv(paste0("~/work/symAtrial_QTL/data/current//gtex/v8/",
                                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"),
                         header = T, stringsAsFactors = F, sep = "\t", skip = 2)[, c("Name", "Description")]
gtex.genes2b$gene_id <- gtex.genes2b$Name
gtex.genes2b$gene_id_sub <- gsub("\\..*", "", gtex.genes2b$Name)
gtex.genes2b$Name <- NULL
gtex.genes2 <- merge(gtex.genes2a, gtex.genes2b,
                     all = T)

gtex.genes <- merge(gtex.genes1[, c("gene_id_sub", "gene_id")],
                    gtex.genes2[, c("gene_id_sub", "gene_id")],
                    by = c("gene_id_sub"),
                    suffixes = c(".v7", ".v8"),
                    all = T)
gtex.genes <- merge(gtex.genes,
                    gtex.genes2[!is.na(gtex.genes2$Description), c("gene_id_sub", "Description")],
                    by = c("gene_id_sub"),
                    all = T)
gtex.genes <- merge(gtex.genes,
                    gtex.genes1[!is.na(gtex.genes1$Description), c("gene_id_sub", "Description")],
                    by = c("gene_id_sub", "Description"),
                    all = T)
gtex.genes <- gtex.genes[!duplicated(gtex.genes), ]

# gtex.anno <- merge(data.frame(symbol=genes,
#                                stringsAsFactors = F),
#                    gtex.genes,
#                     by.x = "symbol", by.y = "Description",
#                     all = T)

library(biomaRt)
genesQTL <- readRDS("genes_QTL_temp.RDS")

#GRCh37.p13
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# filter for desired gene ids
anno <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                           'start_position','end_position','hgnc_symbol',
                           'external_gene_name'),
              filters = 'external_gene_name', values=genesQTL$symbol,
              mart = ensembl)
anno <- anno[anno$chromosome_name %in% as.character(c(1:26)), ]
match <- anno$external_gene_name==anno$hgnc_symbol
table(genesQTL$symbol %in% anno$external_gene_name[match])
miss <- genesQTL$symbol[!(genesQTL$symbol %in% anno$external_gene_name[match])]
anno2 <- anno[anno$external_gene_name %in% miss, ]
anno <- rbind(anno[match, ],
              anno2)
anno$hgnc_symbol <- NULL
saveRDS(anno,
        file=paste0(get.path("locations", local),
                    "ensembl_ids_AFHRIB_QTLgenes.RDS"))
# anno <- readRDS(paste0(get.path("locations", local),
#                        "ensembl_ids_AFHRIB_QTLgenes.RDS"))
anno2 <- merge(anno[complete.cases(anno[, c("ensembl_gene_id", "external_gene_name")]),
                    c("ensembl_gene_id", "external_gene_name")],
               gtex.genes[complete.cases(gtex.genes[, c("gene_id_sub", "Description")]),
                          c("gene_id_sub", "Description")],
               by.x=c("ensembl_gene_id", "external_gene_name"),
               by.y=c("gene_id_sub", "Description"),
               all = T)
anno2 <- merge(anno2,
               gtex.genes[, c("gene_id_sub", "gene_id.v7", "gene_id.v8")],
               by.x=c("ensembl_gene_id"),
               by.y=c("gene_id_sub"),
               all = T)
anno2 <- anno2[!duplicated(anno2), ]
colnames(anno2) <- c("ensembl_gene_id", "symbol", "gene_id.v7", "gene_id.v8")
test <- genesQTL$symbol[!(genesQTL$symbol %in% anno2$symbol)]

anno3 <- anno2[anno2$symbol %in% genesQTL$symbol, ]

saveRDS(anno3,
        file=paste0(get.path("locations", local),
                    "GTEx_ensembl_ids_AFHRIB_QTLgenes.RDS"))
write.table(data.frame(gene_id=unique(anno2$gene_id.v7[anno2$symbol %in% genesQTL$symbol]),
                       stringsAsFactors = F),
            file=paste0(get.path("gtex", local),
                        "gtex_genes_in_afhrib.txt"),
            row.names = F, col.names = T, quote = F)

# filter gtex file for AFHRI-B SNPs and then genes
gtex.file <- paste0(get.path("gtex", local),
                    "Heart_Atrial_Appendage.allpairs.txt")
afhribgenes <- paste0(get.path("gtex", local),
                      "gtex_genes_in_afhrib.txt")
afhribSNPs <- paste0(get.path("gtex", local),
                     "gtex_snps_in_afhrib.txt")
gtex.cut1 <- paste0(get.path("gtex", local),
                    "Heart_Atrial_Appendage.allpairs.cut.genes.txt")
tmp <- tempfile()
gtex.cut2 <- paste0(get.path("gtex", local),
                    "Heart_Atrial_Appendage.allpairs.cut.genes.snps.txt")
gtex.cut3 <- paste0(get.path("gtex", local),
                    "Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt")

system(paste("awk 'FNR==NR { a[$1]; next } $1 in a'", afhribgenes, gtex.file,
             ">", gtex.cut1, sep=" "))
system(paste("awk 'FNR==NR { a[$1]; next } $2 in a'", afhribSNPs, gtex.cut1,
             ">", tmp, sep=" "))
system(paste("awk 'NR==1 {print $0 \"\tAllele1\tAllele2\"} NR!=1 {split($2,a,\"_\"); print $0 \"\t\" a[3] \"\t\" a[4]}'",
             tmp, ">", gtex.cut2, sep=" "))
unlink(tmp)

#awk 'NR==1 {print $0 "\t Allele1 \t Allele2"} NR!=1 {split($2,a,"_"); print $0 "\t" a[3] "\t" a[4]}' Heart_Atrial_Appendage.allpairs.cut.genes.txt > Heart_Atrial_Appendage.allpairs.cut.genes.txt
#awk 'NR==1 {print $0 "\t Allele1 \t Allele2"} NR!=1 {split($2,a,"_"); print $0 "\t" a[3] "\t" a[4]}' Heart_Atrial_Appendage.allpairs.cut.genes.snps.significant.txt > testfile.txt
# gtex.cut = read.table(pipe(paste("cat ", gtex.file, " | fgrep -w -f ", afhribSNPs)),
#                       stringsAsFactors=F)
# awk '(NR==1) || ($9 < 0.05) ' myfile.qassoc > myfile_subset.qassoc
system(paste("awk 'NR==1 {print $0} NR!=1 { if ($7 <= 1e-5) { print } }'", gtex.cut2,
             ">", gtex.cut3, "&", sep=" "))

gtex.cut <- read.table(gtex.cut2,
                       sep="\t", h=T,
                       stringsAsFactors = F)
#colnames(gtex.cut) <- colnames(read.table(gtex.file, nrows=1, sep="\t", h=T))
saveRDS(gtex.cut,
        file=paste0(get.path("gtex", local),
                    "tmp/GTEx_Heart_Atrial_Appendage.allpairs.cut.RDS"))

colnames(gtex.cut)
# "gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"




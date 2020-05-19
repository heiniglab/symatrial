# polygenic risk scores for heart disease ----

setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts/polygenic_risk_scores/")
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

# prepare 1000 genomes SNP data ----
score_file <- paste0(get.path("risk", local), "AF/AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt")
AF = read.table(score_file,
                #skip = 15,
                stringsAsFactors=F, h=T)

AF_snp_file <- paste0(get.path("risk", local),
                      "AF/AtrialFibrillation_PRS_LDpred_rho0.003_v3_SNP_positions.txt")
write.table(AF[, c("chr", "position_hg19", "position_hg19", "A1", "A2", "variant")],
            file = AF_snp_file,
            row.names = F, col.names = F,
            sep = "\t", quote = F)

score_file <- paste0(get.path("risk", local), "CAD/CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt")
CAD = read.table(score_file,
                 stringsAsFactors=F, h=T)

CAD_snp_file <- paste0(get.path("risk", local),
                       "CAD/CoronaryArteryDisease_PRS_LDpred_rho0.001_v3_SNP_positions.txt")
write.table(CAD[, c("chr", "position_hg19", "position_hg19", "A1", "A2", "variant")],
            file = CAD_snp_file,
            row.names = F, col.names = F,
            sep = "\t", quote = F)

AF$disease <- "AF"
CAD$disease <- "CAD"
AF_CAD <- rbind(AF, CAD)
AF_CAD$disease <- factor(AF_CAD$disease, levels = c("CAD", "AF"), ordered = T)

saveRDS(AF_CAD[, c("variant", "effect_weight", "effect_allele", "disease")],
        file=paste0(get.path("risk", local),
                    "AF_CAD_weights.RDS"))

snp_file <- paste0(get.path("risk", local),
                   "AF_CAD_SNP_positions.txt")
write.table(AF_CAD[!duplicated(AF_CAD[, c("chr", "position_hg19", "A1", "A2")]),
                   c("chr", "position_hg19", "position_hg19", "A1", "A2")],
            file = snp_file,
            row.names = F, col.names = F,
            sep = "\t", quote = F)

# ## List of variants and weights comprising a polygenic score for Atrial Fibrillation
# ## Optimal score derived via LDpred with rho=0.003
# ##
# ## From: Amit V. Khera, Mark Chaffin, Krishna G. Aragam, Mary E. Haas, Carolina Roselli, Seung Hoan Choi, Pradeep Natarajan, Eric S. Lander, Steven A. Lubitz, Patrick T. Ellinor, and Sekar Kathiresan. Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. (2018) Nature Genetics, in press.
# ##
# ## Created: 16 July 2018
# ##
# ## Column definitions:
# ## variant=Variant ID, as chromosome:position_hg19:allele1:allele2
# ## effect_allele=Allele to which effect weight is relevant to
# ## effect_weight=Weight of the effect estimate for effect_allele
# ## chr=Chromosome
# ## position_hg19=Position based on hg19
# ## A1=Allele1
# ## A2=Allele2
# variant effect_allele   effect_weight   chr     position_hg19   A1      A2
# 1:30923:G:T     G       4.7853e-06      1       30923   G       T
# 1:52238:T:G     G       2.7000e-06      1       52238   G       T
# 1:54490:A:G     A       1.5657e-06      1       54490   A       G
# 1:58814:A:G     G       2.4831e-05      1       58814   G       A
# 1:59040:C:T     T       1.9251e-05      1       59040   T       C
# 1:63671:A:G     G       3.5490e-05      1       63671   G       A
# 1:77874:A:G     G       5.5080e-06      1       77874   G       A

for (i in 1:22){
  chr <- paste0("chr", i)
  vcf <- paste0("ALL.", chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")
  system(paste0("cd /storage/groups/epigenereg01/workspace/public_data/1000genomes/Phase3/ \n",
                "plink --vcf ", vcf, ".vcf.gz --double-id --make-bed --out plink/", vcf))
}

info <- "/storage/groups/epigenereg01/workspace/public_data/1000genomes/Phase3/integrated_call_samples_v2.20130502.ALL.ped"
sample_info <- read.csv(info,
                        stringsAsFactors=F, h=T, sep="\t")
samples <- sample_info[sample_info$Population %in% c("CEU", "TSI", "FIN", "GBR", "IBS"), ]
samples <- samples[samples$Relationship %in% c("unrel", "father", "mother"), ]
write.table(samples[, c("Individual.ID", "Individual.ID")],
            file="/storage/groups/epigenereg01/workspace/public_data/1000genomes/Phase3/plink/EUR_founders.txt",
            row.names = F, col.names = F, quote = F)


for (i in 1:22){
  chr <- paste0("chr", i)
  bfile <- paste0("ALL.", chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")
  out <- paste0(get.path("risk", local),
                "1000genomes/riskSNPs_AF_CAD_EUR_", chr)
  system(paste0("cd /storage/groups/epigenereg01/workspace/public_data/1000genomes/Phase3/plink/ \n",
                "plink --bfile ", bfile, " --extract range ", snp_file, " --keep EUR_founders.txt --make-bed --out ", out))
}

write.table(paste0(get.path("risk", local), "1000genomes/riskSNPs_AF_CAD_EUR_", paste0("chr", 1:22)),
            file=paste0(get.path("risk", local), "1000genomes/merge_list.txt"),
            sep="\n", col.names = F, row.names = F, quote = F)

out <- "riskSNPs_AF_CAD_EUR_temp"
system(paste0("cd ", get.path("risk", local), "1000genomes/ \n",
              "plink --merge-list merge_list.txt --make-bed --out ", out))

out <- paste0(get.path("risk", local),
              "1000genomes/riskSNPs_AF_CAD_EUR")
for (i in 1:22){
  chr <- paste0("chr", i)
  system(paste0("cd /storage/groups/epigenereg01/workspace/public_data/1000genomes/Phase3/plink/ \n",
                "plink --bfile ", out, "_", chr, " --exclude ", out, "_temp-merge.missnp --keep EUR_founders.txt --make-bed --out ", out, "_", chr))
}
system(paste0("cd ", get.path("risk", local), "1000genomes/ \n",
              "plink --merge-list merge_list.txt --make-bed --out ", out))

map_EUR = read.table(paste0(get.path("risk", local), "1000genomes/riskSNPs_AF_CAD_EUR.bim"),
                     stringsAsFactors=F)
colnames(map_EUR) <- c("chr", "snps", "cM", "pos", "Allele1", "Allele2")
map_EUR$cM <- NULL

## for 1000genomes, we don't find those duplicated SNPs, probably artefact from imputation of the AFHRI-B cohort
# doubles_map_EUR <- get.duplicates(map_EUR, c("chr", "pos"))
# doubles_map_EUR$id <- paste(doubles_map_EUR$chr, doubles_map_EUR$pos, doubles_map_EUR$Allele2, doubles_map_EUR$Allele1, sep=":")
# doubles_map_EUR$id2 <- paste(doubles_map_EUR$chr, doubles_map_EUR$pos, doubles_map_EUR$Allele1, doubles_map_EUR$Allele2, sep=":")
# flipped_alleles <- intersect(doubles_map_EUR$id, doubles_map_EUR$id2)

# prepare AFHRI-B data
map_afhri = read.table(paste0(get.path("genotype", local), "AFHRI-B_imputed.bim"),
                       stringsAsFactors=F)
colnames(map_afhri) <- c("chr", "snps", "cM", "pos", "Allele1", "Allele2")
map_afhri$cM <- NULL

doubles_map_afhri <- get.duplicates(map_afhri, c("chr", "pos"))
doubles_map_afhri$id <- paste(doubles_map_afhri$chr, doubles_map_afhri$pos, doubles_map_afhri$Allele2, doubles_map_afhri$Allele1, sep=":")
doubles_map_afhri$id2 <- paste(doubles_map_afhri$chr, doubles_map_afhri$pos, doubles_map_afhri$Allele1, doubles_map_afhri$Allele2, sep=":")
flipped_alleles <- intersect(doubles_map_afhri$id, doubles_map_afhri$id2)

write.table(doubles_map_afhri[doubles_map_afhri$id %in% flipped_alleles |
                                doubles_map_afhri$id2 %in% flipped_alleles, "snps"],
            file=paste0(get.path("genotype", local), "AFHRI-B_imputed_doubles.txt"),
            sep = "\n", row.names = F, col.names = F, quote = F)
afhri <- paste0(get.path("genotype", local), "AFHRI-B_imputed")
system(paste0("cd ", get.path("genotype", local), " \n",
              "plink --bfile ", afhri, " --extract AFHRI-B_imputed_doubles.txt --make-bed --recode --missing --out ", afhri, "_doubles"))

lmiss_afhri <- read.table(paste0(get.path("genotype", local),
                                 "AFHRI-B_imputed_doubles.lmiss"),
                          h=T)

doubles_map_afhri2 <- merge(doubles_map_afhri,
                            lmiss_afhri,
                            by.x=c("chr", "snps"),
                            by.y=c("CHR", "SNP"),
                            all.y=T, sort=F)

doubles_map_afhri2 <- doubles_map_afhri2[order(doubles_map_afhri2$F_MISS, decreasing = T), ]
doubles_map_afhri3 <- doubles_map_afhri2[!duplicated(doubles_map_afhri2[, c("chr", "pos")]), ]

write.table(doubles_map_afhri3$snps,
            file=paste0(get.path("genotype", local), "AFHRI-B_imputed_remove_doubles.txt"),
            sep = "\n", row.names = F, col.names = F, quote = F)
afhri <- paste0(get.path("genotype", local), "AFHRI-B_imputed")
system(paste0("cd ", get.path("genotype", local), " \n",
              "plink --bfile ", afhri, " --exclude AFHRI-B_imputed_remove_doubles.txt --make-bed --out ", afhri, "2"))

map_afhri = read.table(paste0(get.path("genotype", local), "AFHRI-B_imputed2.bim"),
                       stringsAsFactors=F)
colnames(map_afhri) <- c("chr", "snps", "cM", "pos", "Allele1", "Allele2")
map_afhri$cM <- NULL

map <- merge(map_EUR,
             map_afhri,
             suffixes = c(".EUR", ".afhri"),
             by = c("chr", "pos"),
             all=T)

map$rsid <- NA
map$rsid[grep("rs", map$snps.EUR)] <- map$snps.EUR[grep("rs", map$snps.EUR)]
map$rsid2 <- NA
map$rsid2[grep("rs", map$snps.afhri)] <- gsub("*\\:rs", "rs", map$snps.afhri[grep("rs", map$snps.afhri)])
map$rsid2[grep("rs", map$rsid2)] <- gsub("\\:.*", "", map$rsid2[grep("rs", map$rsid2)])
map[is.na(map$rsid), "rsid"] <- map[is.na(map$rsid), "rsid2"]
table(is.na(map$rsid))
map$rsid2 <- NULL

map$allele.match <- NA
map[which(map$Allele1.EUR==map$Allele1.afhri & map$Allele2.EUR==map$Allele2.afhri), "allele.match"] <- 1
map[which(map$Allele1.EUR==map$Allele2.afhri & map$Allele2.EUR==map$Allele1.afhri), "allele.match"] <- -1
map[is.na(map$snps.EUR) | is.na(map$snps.afhri), "allele.match"] <- 0

saveRDS(map,
        file=paste0(get.path("risk", local),
                    "1000genomes/temp_merge_map_EUR_AFHRI_B.RDS"))
rm(map_EUR, map_afhri)
map <- readRDS(paste0(get.path("risk", local),
                      "1000genomes/temp_merge_map_EUR_AFHRI_B.RDS"))

# map complete: c("chr", "snps", "cM", "pos", "Allele1", "Allele2")
map1 <- map[which(map$allele.match==1), ]
mapm1 <- map[which(map$allele.match==-1), ]
map0EUR <- map[is.na(map$snps.EUR), ]
map0afhri <- map[is.na(map$snps.afhri), ]
mapmiss <- map[is.na(map$allele.match), ]

map.complete <- rbind(data.frame(chr=map1$chr, rsid=map1$rsid, cM=0, pos=map1$pos,
                                 A1=map1$Allele1.EUR, A2=map1$Allele2.EUR,
                                 match=map1$allele.match, id.EUR=map1$snps.EUR, id.afhri=map1$snps.afhri,
                                 stringsAsFactors=F),
                      data.frame(chr=mapm1$chr, rsid=mapm1$rsid, cM=0, pos=mapm1$pos,
                                 A1=mapm1$Allele1.EUR, A2=mapm1$Allele2.EUR,
                                 match=mapm1$allele.match, id.EUR=mapm1$snps.EUR, id.afhri=mapm1$snps.afhri,
                                 stringsAsFactors=F),
                      data.frame(chr=map0EUR$chr, rsid=map0EUR$rsid, cM=0, pos=map0EUR$pos,
                                 A1=map0EUR$Allele1.afhri, A2=map0EUR$Allele2.afhri,
                                 match=map0EUR$allele.match, id.EUR=map0EUR$snps.EUR, id.afhri=map0EUR$snps.afhri,
                                 stringsAsFactors=F),
                      data.frame(chr=map0afhri$chr, rsid=map0afhri$rsid, cM=0, pos=map0afhri$pos,
                                 A1=map0afhri$Allele1.EUR, A2=map0afhri$Allele2.EUR,
                                 match=map0afhri$allele.match, id.EUR=map0afhri$snps.EUR, id.afhri=map0afhri$snps.afhri,
                                 stringsAsFactors=F),
                      data.frame(chr=mapmiss$chr, rsid=mapmiss$rsid, cM=0, pos=mapmiss$pos,
                                 A1=mapmiss$Allele1.EUR, A2=mapmiss$Allele2.EUR,
                                 match=mapmiss$allele.match, id.EUR=mapmiss$snps.EUR, id.afhri=NA,
                                 stringsAsFactors=F),
                      data.frame(chr=mapmiss$chr, rsid=mapmiss$rsid, cM=0, pos=mapmiss$pos,
                                 A1=mapmiss$Allele1.afhri, A2=mapmiss$Allele2.afhri,
                                 match=mapmiss$allele.match, id.EUR=NA, id.afhri=mapmiss$snps.afhri,
                                 stringsAsFactors=F))

map.complete$id <- paste(map.complete$chr, map.complete$pos,
                         map.complete$A2, map.complete$A1,
                         map.complete$rsid, sep = ":")
map.complete$id2 <- paste(map.complete$chr, map.complete$pos,
                          map.complete$A2, map.complete$A1,
                          sep = ":")

saveRDS(map.complete,
        file=paste0(get.path("risk", local),
                    "1000genomes/temp_merge_map_EUR_AFHRI_B_complete.RDS"))
rm(map, map1, mapm1, map0afhri, map0EUR, mapmiss)

map.complete <- readRDS(paste0(get.path("risk", local),
                               "1000genomes/temp_merge_map_EUR_AFHRI_B_complete.RDS"))

# get common names for SNPs (chr:pos:A1:A2) ----------
## prepare 1000genomes data
rename.EUR <- map.complete[, c("id.EUR", "id2")]
rename.EUR <- rename.EUR[complete.cases(rename.EUR), ]
dim(rename.EUR)
rename.EUR <- rename.EUR[!duplicated(rename.EUR), ]
dim(rename.EUR)
test <- get.duplicates(rename.EUR, "id.EUR")
test2 <- get.duplicates(rename.EUR, "id2")
# good: no duplicate IDs!

write.table(rename.EUR,
            file=paste0(get.path("risk", local), "1000genomes/update_snp_names_EUR.txt"),
            sep="\t", col.names = F, row.names = F, quote = F)
EUR <- "riskSNPs_AF_CAD_EUR"
system(paste0("cd ", get.path("risk", local), "1000genomes/ \n",
              "plink --bfile ", EUR, " --update-name update_snp_names_EUR.txt --make-bed --out ", EUR, "_merge"))
rm(rename.EUR, test, test2)

## prepare AFHRI-B data
rename.afhri <- map.complete[complete.cases(map.complete[, c("id.afhri", "id2")]), ]
rename.afhri <- rename.afhri[complete.cases(rename.afhri[, c("id.afhri", "id2")]), ]
dim(rename.afhri)
rename.afhri <- rename.afhri[!duplicated(rename.afhri[, c("id.afhri", "id2")]), ]
dim(rename.afhri)
test <- get.duplicates(rename.afhri, "id.afhri")
test2 <- rename.afhri[(rename.afhri$id2 %in% test$id2 &
                         rename.afhri$id.afhri %in% test$id.afhri &
                         is.na(rename.afhri$match)), ]
rename.afhri <- rename.afhri[!(rename.afhri$id2 %in% test$id2 &
                                 rename.afhri$id.afhri %in% test$id.afhri &
                                 is.na(rename.afhri$match)), ]

test3 <- get.duplicates(rename.afhri, "id.afhri")
test4 <- get.duplicates(rename.afhri, "id2")
test5 <- get.duplicates(rename.afhri, "id")

rename.afhri$id3 <- rename.afhri$id2
rename.afhri$id3[rename.afhri$id2 %in% test4$id2] <- rename.afhri$id[rename.afhri$id2 %in% test4$id2]

write.table(rename.afhri[, c("id.afhri", "id3")],
            file=paste0(get.path("genotype", local), "update_snp_names_AFHRI-B.txt"),
            sep="\t", col.names = F, row.names = F, quote = F)

afhri <- "AFHRI-B_imputed2"
system(paste0("cd ", get.path("genotype", local), " \n",
              "plink --bfile ", afhri, " --update-name update_snp_names_AFHRI-B.txt --make-bed --out ", afhri, "_rename"))
rm(test, test2, test3, test4)

########################
# flip reference alleles in AFHRI-B data to match EUR
map_EUR = read.table(paste0(get.path("risk", local), "1000genomes/riskSNPs_AF_CAD_EUR_merge.bim"),
                     stringsAsFactors=F)
colnames(map_EUR) <- c("chr", "snps", "cM", "pos", "Allele1", "Allele2")
map_EUR$cM <- NULL

map_afhri = read.table(paste0(get.path("genotype", local), "AFHRI-B_imputed2_rename.bim"),
                       stringsAsFactors=F)
colnames(map_afhri) <- c("chr", "snps", "cM", "pos", "Allele1", "Allele2")
map_afhri$cM <- NULL

test <- get.duplicates(map_afhri, "snps")
map <- merge(map_afhri,
                map_EUR,
                by=c("chr", "snps", "pos"),
                all.x=T, sort=F,
                suffixes=c(".afhri", ".EUR"))

map$allele.match <- NA
map[which(map$Allele1.EUR==map$Allele1.afhri & map$Allele2.EUR==map$Allele2.afhri), "allele.match"] <- 1
map[which(map$Allele1.EUR==map$Allele2.afhri & map$Allele2.EUR==map$Allele1.afhri), "allele.match"] <- -1
map$ref <- map$Allele1.afhri
map$ref[which(map$allele.match=="-1")] <- map$Allele1.EUR[which(map$allele.match=="-1")]

table(map$allele.match)
table(is.na(map$allele.match))

write.table(map[, c("snps", "ref")],
            file = paste0(get.path("genotype", local), "AFHRI-B_imputed2_rename_reference.txt"),
            row.names = F,col.names = F, quote=F)

# reference allele == A1 allele!!!
# mylist.txt=
# rs10001 A
# rs10002 T
# rs10003 T
# plink --bfile mydata --reference-allele mylist.txt --assoc

system(paste0("cd ", get.path("genotype", local), " \n",
              "plink --bfile AFHRI-B_imputed2_rename --reference-allele AFHRI-B_imputed2_rename_reference.txt --make-bed --out AFHRI-B_imputed_merge"))


afhri <- paste0(get.path("genotype", local), "AFHRI-B_imputed_merge")
EUR <- paste0(get.path("risk", local), "1000genomes/riskSNPs_AF_CAD_EUR_merge")
out <- paste0(get.path("genotype", local),
              "riskSNPs_AF_CAD_EUR_AFHRI_B")
system(paste0("cd ", get.path("genotype", local), " \n",
              "plink --bfile ", afhri, " --bmerge ", EUR, " --make-bed --out ", out))


# Calculate risk scores on merged data-------------

setwd("~work/symAtrial_QTL/scripts/polygenic_risk_scores/")
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

map_file <- paste0(get.path("genotype", local),
                   "riskSNPs_AF_CAD_EUR_AFHRI_B.bim")
map = read.table(map_file, stringsAsFactors=F)
colnames(map) <- c("chr", "snps", "cM", "pos", "Allele1", "Allele2")
map$cM <- NULL


# AF-risk score ----
score_file <- paste0(get.path("risk", local), "AF/AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt")
AF = read.table(score_file,
                stringsAsFactors=F, h=T)
colnames(AF)
dim(AF)

merged <- merge(AF, map,
                by.x = c("chr", "position_hg19"),
                by.y = c("chr", "pos"),
                all.x = T)
dim(merged)

merged$allele.match <- NA
merged[which(merged$Allele1==merged$A1 & merged$Allele2==merged$A2), "allele.match"] <- 1
merged[which(merged$Allele1==merged$A2 & merged$Allele2==merged$A1), "allele.match"] <- -1

merged <- merged[!is.na(merged$allele.match), ]
colnames(merged)
dim(merged)

test <- get.duplicates(merged, c("variant"))

ref_risk_file <- paste0(get.path("risk", local),
                        "AF/AtrialFibrillation_PRS_LDpred_rho0.003_v3_reference_EUR_AFHRI-B.txt")
write.table(merged[, c("snps", "effect_allele", "effect_weight")],
            file = ref_risk_file,
            row.names = F, col.names = F,
            sep = "\t", quote = F)

bfile <- paste0(get.path("genotype", local), "riskSNPs_AF_CAD_EUR_AFHRI_B")
outfile <- paste0(get.path("genotype", local), "AF_EUR_AFHRI-B_risk_score")
system(paste0("plink --bfile ", bfile,
              " --score ", ref_risk_file, " sum ",
              " --out ", outfile))

saveRDS(merged,
        paste0(get.path("risk", local),
               "AF/AtrialFibrillation_PRS_LDpred_rho0.003_v3_reference_EUR_AFHRI-B.RDS"))

outfile <- paste0(get.path("genotype", local), "AF_EUR_AFHRI-B_risk_score")
AF_score <- read.table(paste0(outfile, ".profile"),
                       h=T, stringsAsFactors=F)
phenos <- readRDS(paste0(get.path("dataset", local),
                         "AFHRI_B_phenotypes.RDS"))
phenos$IID <- gsub(" LAA", "", gsub(" RAA", "", phenos$externID))
AF_score2 <- merge(AF_score,
                   phenos,
                   by = "IID",
                   all.x = T)

plot(AF_score2$SCORESUM ~AF_score2$overall_AF)
wilcox.test(AF_score2$SCORESUM ~AF_score2$overall_AF)

library(ggplot2)
ggplot(AF_score2, aes(x=overall_AF, y=SCORESUM, col=overall_AF,  fill=overall_AF)) +
  geom_point(position = position_jitter(width = 0.3), size=1) +
  geom_boxplot(alpha=0.2) +
  labs(list(title="Genome-wide polygenic score for Atrial Fibrillation",
            y="AF GPS")) +
  theme_bw()

ggplot(AF_score2, aes(x=overall_AF, y=rank(SCORESUM)/dim(AF_score2)[1], col=overall_AF,  fill=overall_AF)) +
  geom_point(position = position_jitter(width = 0.3), size=1) +
  geom_boxplot(alpha=0.2) +
  labs(list(title="Genome-wide polygenic score for Atrial Fibrillation",
            y="percentiles of AF GPS")) +
  theme_bw()


# CAD ----
score_file <- paste0(get.path("risk", local), "CAD/CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt")
CAD = read.table(score_file,
                 stringsAsFactors=F, h=T)
colnames(CAD)
dim(CAD)

merged_CAD <- merge(CAD, map,
                by.x = c("chr", "position_hg19"),
                by.y = c("chr", "pos"),
                all.x = T)
dim(merged_CAD)

merged_CAD$allele.match <- NA
merged_CAD[which(merged_CAD$Allele1==merged_CAD$A1 &
                   merged_CAD$Allele2==merged_CAD$A2), "allele.match"] <- 1
merged_CAD[which(merged_CAD$Allele1==merged_CAD$A2 &
                   merged_CAD$Allele2==merged_CAD$A1), "allele.match"] <- -1

merged_CAD <- merged_CAD[!is.na(merged_CAD$allele.match), ]
colnames(merged_CAD)
dim(merged_CAD)

test <- get.duplicates(merged, c("variant"))

ref_risk_file <- paste0(get.path("risk", local),
                        "CAD/CoronaryArteryDisease_PRS_LDpred_rho0.001_v3_reference_EUR_AFHRI_B.txt")

write.table(merged_CAD[, c("snps", "effect_allele", "effect_weight")],
            file = ref_risk_file,
            row.names = F, col.names = F,
            sep = "\t", quote = F)

bfile <- paste0(get.path("genotype", local), "riskSNPs_AF_CAD_EUR_AFHRI_B")
outfile <- paste0(get.path("genotype", local), "CAD_EUR_AFHRI-B_risk_score")
system(paste0("plink --bfile ", bfile,
              # " --recodeA --recode-allele ", ref_risk_file,
              # " --extract ", ref_risk_file,
              " --score ", ref_risk_file, " sum ",
              " --out ", outfile))

table(table(merged_CAD$snps))

saveRDS(merged_CAD,
        paste0(get.path("risk", local),
               "CAD/CoronaryArteryDisease_PRS_LDpred_rho0.001_v3_reference_EUR_AFHRI-B.RDS"))

outfile <- paste0(get.path("genotype", local), "CAD_EUR_AFHRI-B_risk_score")
CAD_score <- read.table(paste0(outfile, ".profile"),
                        h=T, stringsAsFactors=F)

phenos <- readRDS(paste0(get.path("dataset", local),
                         "AFHRI_B_phenotypes.RDS"))
phenos$IID <- gsub(" LAA", "", gsub(" RAA", "", phenos$externID))
CAD_score2 <- merge(CAD_score,
                   phenos,
                   by = "IID",
                   all.x = T)
CAD_score2$cohort <- "AFHRI-B"
CAD_score2$cohort[is.na(CAD_score2$externID)] <- "1000 genomes"

library(ggplot2)
ggplot(CAD_score2, aes(x=cohort, y=SCORESUM, col=cohort, fill=cohort)) +
  geom_point(position = position_jitter(width = 0.3), size=1) +
  geom_boxplot(alpha=0.2) +
  labs(list(title="Genome-wide polygenic score for CAD",
            y="CAD GPS")) +
  theme_bw()

ggplot(CAD_score2, aes(x=cohort, y=rank(SCORESUM)/dim(CAD_score2)[1], col=cohort, fill=cohort)) +
  geom_point(position = position_jitter(width = 0.3), size=1) +
  geom_boxplot(alpha=0.2) +
  labs(list(title="Genome-wide polygenic score for CAD",
            y="percentiles of CAD GPS")) +
  theme_bw()

ggplot(CAD_score2, aes(x=cohort, y=rank(SCORESUM)/dim(CAD_score2)[1], col=cohort, fill=cohort)) +
  geom_point(position = position_jitter(width = 0.3), size=1) +
  geom_violin(alpha=0.2) +
  labs(list(title="Genome-wide polygenic score for CAD",
            y="percentiles of CAD GPS")) +
  theme_bw()
setwd("~/work/symAtrial_QTL/scripts")

library(VariantAnnotation)
library(GenomicRanges)

#http://feb2014.archive.ensembl.org/biomart/martview/
source("helper/helper.R")

local=F
#data_dir <- get.path("data")

# bcftools view -i 'VT ~"SNP"' /Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/data_local/20180427/annotation/lappalainen/snps.vcf.bgz > filtered.vcf
# bcftools view -i 'EUR_AF >0.01' filtered.vcf > filtered2.vcf
# bcftools view -i 'EUR_AF <0.99' filtered2.vcf > filtered3.vcf

map.seq.levels <- function(gr){
  newStyle <- mapSeqlevels(seqlevels(gr), "UCSC")
  newStyle <- newStyle[!is.na(newStyle)]
  gr <- renameSeqlevels(gr, newStyle)
  return(gr)
}


# exons ----

# # exon_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/martquery_exon.txt"), header = T, sep = ",", stringsAsFactors = F)
# # exon <- GRanges(
# #   seqnames = Rle(exon_data$Chromosome.Name),
# #   ranges = IRanges(exon_data$Exon.Chr.Start..bp.,
# #                    end = exon_data$Exon.Chr.End..bp.),
# #   strand = Rle(exon_data$Strand),
# #   rank = exon_data$Exon.Rank.in.Transcript,
# #   symbol = exon_data$Associated.Gene.Name,
# #   feature = "exon"
# # )
# # exon <- map.seq.levels(exon)
# # saveRDS(exon, paste0(get.path("locations", local), "functional_annotation/exon_GRanges.RDS"))
# exon <- readRDS(paste0(get.path("locations", local), "functional_annotation/exon_GRanges.RDS"))
# exon

# exon_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/martquery_exon2.txt"), header = T, sep = ",", stringsAsFactors = F)
# exon_data <- exon_data[!duplicated(exon_data[,c("Chromosome.Name", "Exon.Chr.Start..bp.",
#                                                 "Exon.Chr.End..bp.", "Ensembl.Exon.ID",
#                                                 "Strand", "Ensembl.Transcript.ID",
#                                                 "Exon.Rank.in.Transcript", "Associated.Gene.Name")]),]
# exon_data <- exon_data[exon_data$Chromosome.Name %in% c(1:30, "X", "Y", "M", "MT"), ]
# exon <- GRanges(
#   seqnames = Rle(exon_data$Chromosome.Name),
#   ranges = IRanges(exon_data$Exon.Chr.Start..bp.,
#                    end = exon_data$Exon.Chr.End..bp.,
#                    names = exon_data$Ensembl.Exon.ID),
#   strand = Rle(exon_data$Strand),
#   transcript = exon_data$Ensembl.Transcript.ID,
#   rank = exon_data$Exon.Rank.in.Transcript,
#   symbol = exon_data$Associated.Gene.Name,
#   feature = "exon",
#   gene_start = exon_data$Gene.Start..bp.,
#   gene_end = exon_data$Gene.End..bp.
# )
# exon <- map.seq.levels(exon)
# saveRDS(exon, paste0(get.path("locations", local), "functional_annotation/exon_ids_GRanges.RDS"))
exon <- readRDS(paste0(get.path("locations", local), "functional_annotation/exon_ids_GRanges.RDS"))
exon



# splice sites ----

# spl.start <- GRanges(
#   seqnames = Rle(exon_data$Chromosome.Name),
#   ranges = IRanges(exon_data$Exon.Chr.Start..bp.-1,
#                    end = exon_data$Exon.Chr.Start..bp.),
#   strand = Rle(exon_data$Strand),
#   transcript = exon_data$Ensembl.Transcript.ID,
#   rank = exon_data$Exon.Rank.in.Transcript,
#   symbol = exon_data$Associated.Gene.Name,
#   feature = "splice_site",
#   position = "start"
# )
# spl.end <- GRanges(
#   seqnames = Rle(exon_data$Chromosome.Name),
#   ranges = IRanges(exon_data$Exon.Chr.End..bp.,
#                    end = exon_data$Exon.Chr.End..bp.+1),
#   strand = Rle(exon_data$Strand),
#   transcript = exon_data$Ensembl.Transcript.ID,
#   rank = exon_data$Exon.Rank.in.Transcript,
#   symbol = exon_data$Associated.Gene.Name,
#   feature = "splice_site",
#   position = "end"
# )
# splice <- c(spl.start, spl.end)
# splice <- map.seq.levels(splice)
# saveRDS(splice, paste0(get.path("locations", local), "functional_annotation/splicesites_GRanges.RDS"))
splice <- readRDS(paste0(get.path("locations", local), "functional_annotation/splicesites_GRanges.RDS"))
splice



# UTR ----

# UTR_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/martquery_UTR.txt"), header = T, sep = "\t", stringsAsFactors = F)
# UTR_data <- UTR_data[, c("Chromosome.Name", "Strand", "Associated.Gene.Name",
#                          "Ensembl.Transcript.ID",
#                          "X5..UTR.Start", "X5..UTR.End", "X3..UTR.Start", "X3..UTR.End")]
# UTR_data <- UTR_data[UTR_data$Chromosome.Name %in% c(1:30, "X", "Y", "M", "MT"), ]
# UTR_data_5p <- UTR_data[which(!is.na(UTR_data$X5..UTR.Start)),
#                         c("Chromosome.Name", "Strand", "Ensembl.Transcript.ID", "Associated.Gene.Name",
#                           "X5..UTR.Start", "X5..UTR.End")]
# UTR_data_3p <- UTR_data[which(!is.na(UTR_data$X3..UTR.Start)),
#                         c("Chromosome.Name", "Strand", "Ensembl.Transcript.ID", "Associated.Gene.Name",
#                           "X3..UTR.Start", "X3..UTR.End")]
# UTR_data_5p <- UTR_data_5p[!duplicated(UTR_data_5p), ]
# UTR_data_3p <- UTR_data_3p[!duplicated(UTR_data_3p), ]
# UTR5p <- GRanges(
#   seqnames = Rle(UTR_data_5p$Chromosome.Name),
#   ranges = IRanges(UTR_data_5p$X5..UTR.Start,
#                    end = UTR_data_5p$X5..UTR.End),
#   strand = Rle(UTR_data_5p$Strand),
#   transcript = UTR_data_5p$Ensembl.Transcript.ID,
#   symbol = UTR_data_5p$Associated.Gene.Name,
#   feature = "utr5"
# )
# UTR3p <- GRanges(
#   seqnames = Rle(UTR_data_3p$Chromosome.Name),
#   ranges = IRanges(UTR_data_3p$X3..UTR.Start,
#                    end = UTR_data_3p$X3..UTR.End),
#   strand = Rle(UTR_data_3p$Strand),
#   transcript = UTR_data_3p$Ensembl.Transcript.ID,
#   symbol = UTR_data_3p$Associated.Gene.Name,
#   feature = "utr3"
# )
# UTR <- c(UTR5p, UTR3p)
# UTR <- map.seq.levels(UTR)
# saveRDS(UTR, paste0(get.path("locations", local), "functional_annotation/UTR_GRanges.RDS"))
UTR <- readRDS(paste0(get.path("locations", local), "functional_annotation/UTR_GRanges.RDS"))
UTR


# transcript locations ----

# locs_data <- read.csv(paste0(get.path("locations", local),
#                              "functional_annotation/martquery_transcript_locs_type.txt"),
#                       header = T, sep = "\t", stringsAsFactors = F)
# locs.data <- locs_data[, c("Chromosome.Name", "Transcript.Start..bp.", "Transcript.End..bp.",
#                            "Strand", "Ensembl.Transcript.ID", "Associated.Gene.Name",
#                            "Gene.Start..bp.", "Gene.End..bp.", "Transcript.Biotype")]
# locs.data <- locs.data[locs.data$Chromosome.Name %in% c(1:30, "X", "Y", "MT"), ]
# locs.data <- locs.data[!duplicated(locs.data), ]
# t_locs <- GRanges(
#     seqnames = Rle(locs.data$Chromosome.Name),
#     ranges = IRanges(locs.data$Transcript.Start..bp.,
#                      end = locs.data$Transcript.End..bp.),
#     strand = Rle(locs.data$Strand),
#     transcript = locs.data$Ensembl.Transcript.ID,
#     symbol = locs.data$Associated.Gene.Name,
#     gene.start = locs.data$Gene.Start..bp.,
#     gene.end = locs.data$Gene.End..bp.,
#     feature = locs.data$Transcript.Biotype,
#     type = "transcript"
#   )
# # locs_data <- read.csv(paste0(get.path("locations", local),
# #                              "functional_annotation/martquery_gene_locs_type.txt"),
# #                       header = T, sep = "\t", stringsAsFactors = F)
# locs.data <- locs_data[, c("Chromosome.Name", "Gene.Biotype",
#                            "Strand", "Associated.Gene.Name",
#                            "Gene.Start..bp.", "Gene.End..bp.")]
# locs.data <- locs.data[!duplicated(locs.data), ]
# #locs.data <- locs.data[order(locs.data$Associated.Gene.Name), ]
# locs.data <- locs.data[locs.data$Chromosome.Name %in% c(1:30, "X", "Y", "MT"), ]
# locs.data <- locs.data[!duplicated(locs.data), ]
# g_locs <- GRanges(
#   seqnames = Rle(locs.data$Chromosome.Name),
#   ranges = IRanges(locs.data$Gene.Start..bp.,
#                    end = locs.data$Gene.End..bp.),
#   strand = Rle(locs.data$Strand),
#   transcript = NA,
#   symbol = locs.data$Associated.Gene.Name,
#   gene.start = locs.data$Gene.Start..bp.,
#   gene.end = locs.data$Gene.End..bp.,
#   feature = locs.data$Gene.Biotype,
#   type = "gene"
# )
# g_locs <- g_locs[!duplicated(g_locs)]
# locs <- c(t_locs, g_locs)
# locs <- map.seq.levels(locs)
# saveRDS(locs, paste0(get.path("locations", local), "functional_annotation/locs_GRanges.RDS"))
locs <- readRDS(paste0(get.path("locations", local), "functional_annotation/locs_GRanges.RDS"))
locs


# RBPBS -----

# RBP1 <- readRDS(paste0(get.path("locations", local), "functional_annotation/RBP/significant.eClip.calls.K562.RDS"))
# RBP2 <- readRDS(paste0(get.path("locations", local), "functional_annotation/RBP/significant.eClip.calls.HepG2.RDS"))
# RBPBS <- c(RBP1, RBP2)
# saveRDS(RBPBS, paste0(get.path("locations", local), "functional_annotation/RBP/RBPBS_GRanges.RDS"))
# RBPBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/RBP/RBPBS_GRanges.RDS"))
# load(paste0(get.path("locations", local), "functional_annotation/RBP/interactome.rda"))
# symbols <- gsub("_.*", "", RBPBS$target)
# RBPBS_heart <- RBPBS[grep(paste(interactome$gene_symbol, collapse = "|"), symbols, ignore.case=TRUE)]
# saveRDS(RBPBS_heart, paste0(get.path("locations", local), "functional_annotation/RBP/RBPBS_heart_GRanges.RDS"))
RBPBS_heart <- readRDS(paste0(get.path("locations", local), "functional_annotation/RBP/RBPBS_heart_GRanges.RDS"))

# TFBS ----

# TF_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/TFBS/remap2018_nr_macs2_hg19_v1_2.bed"),
#                     header = F, sep = "\t", stringsAsFactors = F)
# TFBS <- GRanges(
#   seqnames = Rle(TF_data$V1),
#   ranges = IRanges(TF_data$V2+1,
#                    end = TF_data$V3),
#   name = TF_data$V4,
#   score = TF_data$V5,
#   thickStart = TF_data$V7,
#   thickEnd = TF_data$V8,
#   itemRgb = TF_data$V9
# )
# saveRDS(TFBS, paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_GRanges.RDS"))
TFBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_GRanges.RDS"))

# TF_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/TFBS/remap2018_all_macs2_hg19_v1_2.bed"),
#                     header = F, sep = "\t", stringsAsFactors = F)
# saveRDS(TF_data, paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_single.RDS"))
# 
# TFBS <- GRanges(
#   seqnames = Rle(TF_data$V1),
#   ranges = IRanges(TF_data$V2+1,
#                    end = TF_data$V3),
#   name = TF_data$V4,
#   score = TF_data$V5,
#   thickStart = TF_data$V7,
#   thickEnd = TF_data$V8,
#   itemRgb = TF_data$V9
# )
# saveRDS(TFBS, paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_single_GRanges.RDS"))
TFBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_single_GRanges.RDS"))

# # library(plyr)
# # TF_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/TFBS/remap2018_all_macs2_hg19_v1_2.bed"),
# #                     header = F, sep = "\t", stringsAsFactors = F)
# # write.table(TF_data$V4,
# #             file=paste0(get.path("locations", local),
# #                         "functional_annotation/TFBS/temp.txt"),
# #             row.names = F, col.names = F, quote = F)
# # temp <- read.csv(file=paste0(get.path("locations", local),
# #                              "functional_annotation/TFBS/temp.txt"),
# #                  h=F, sep = ".")
# # colnames(temp) <- c("dataset", "name", "cell")
# # print("strsplit performed")
# # write.table(cbind(TF_data[, 1:3], temp, TF_data[, 5:9]),
# #             file=paste0(get.path("locations", local),
# #                         "functional_annotation/TFBS/remap2018_all_split_macs2_hg19_v1_2.bed"),
# #                   row.names = F, col.names = T, quote = F, sep = "\t")
# # rm(temp)
# TF_data <- read.csv(file=paste0(get.path("locations", local),
#                                 "functional_annotation/TFBS/remap2018_all_split_macs2_hg19_v1_2.bed"),
#                     h = T, sep = "\t", stringsAsFactors=F)
# TFBS <- GRanges(
#   seqnames = Rle(TF_data$V1),
#   ranges = IRanges(TF_data$V2+1,
#                    end = TF_data$V3),
#   dataset = TF_data$dataset,
#   name = TF_data$name,
#   cell = TF_data$cell,
#   score = TF_data$V5,
#   thickStart = TF_data$V7,
#   thickEnd = TF_data$V8,
#   itemRgb = TF_data$V9
# )
# saveRDS(TFBS, paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_single_split_GRanges.RDS"))
TFBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_single_split_GRanges.RDS"))


# # miRNA BS
# miR_data <- read.csv(file=paste0(get.path("locations", local),
#                                 "functional_annotation/targetscan/Predicted_Target_Locations.default_predictions.hg19.bed"),
#                      stringsAsFactors = F, h = F, sep = "\t")
# miRBS <- GRanges(
#   seqnames = Rle(miR_data$V1),
#   ranges = IRanges(miR_data$V2+1,
#                    end = miR_data$V3),
#   strand = Rle(miR_data$V6),
#   name = miR_data$V4,
#   score = miR_data$V5,
#   thickStart = miR_data$V7,
#   thickEnd = miR_data$V8,
#   itemRgb = miR_data$V9
# )
# saveRDS(miRBS, paste0(get.path("locations", local), "functional_annotation/miRBS_GRanges.RDS"))
miRBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/miRBS_GRanges.RDS"))


# Roadmap 15 state ----
# # /home/icb/ines.assum/projects/symAtrial_QTL/data/current/annotation/functional_annotation/roadmap/E104_15/E104_15_coreMarks_dense.bed
# states_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/roadmap/E104_15/E104_15_coreMarks_dense.bed"),
#                     header = F, skip = 1, sep = "\t", stringsAsFactors = F)
# states <- GRanges(
#   seqnames = Rle(states_data$V1),
#   ranges = IRanges(states_data$V2+1,
#                    end = states_data$V3),
#   state = states_data$V4,
#   thickStart = states_data$V7,
#   thickEnd = states_data$V8,
#   itemRgb = states_data$V9
# )
# mean(width(states))
# median(width(states))
# saveRDS(states, paste0(get.path("locations", local), "functional_annotation/15_states_GRanges.RDS"))
states <- readRDS(paste0(get.path("locations", local), "functional_annotation/15_states_GRanges.RDS"))


# VEP
# VEP_data <- read.csv(paste0(get.path("locations", local), "functional_annotation/VEP/martquery_VEP.txt"),
#                      sep = "\t", stringsAsFactors = F, #skip = 1,
#                      header = T)
# VEP_data <- VEP_data[VEP_data$Consequence.to.transcript %in%
#                        c("stop_gained","frameshift_variant", "stop_lost",
#                          "inframe_insertion", "inframe_deletion",
#                          "missense_variant", "NMD_transcript_variant"),
#                      ]
# 
# QTL.snp <- readRDS(file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_group.RDS"))
# # test <- head(QTL.snp$snps)
# # test2 <- test[grep("rs", test)]
# # test3 <- gsub(":.*", "", test2)
# # test3
# 
# QTL.snp2 <- QTL.snp[grep("rs", QTL.snp$snps), ]
# QTL.snp2$rs.id <- gsub(":.*", "", QTL.snp2$snps)
# 
# pairs <- QTL.snp2[, c("rs.id", "gene")]
# VEP_data2 <- merge(VEP_data, pairs,
#                    by.x=c("Variation.Name", "Associated.Gene.Name"),
#                    by.y=c("rs.id", "gene"))
# VEP_data2 <- VEP_data2[!duplicated(VEP_data2), ]
# write.table(VEP_data2,
#             file=paste0(get.path("locations", local), "functional_annotation/VEP/VEP_data_pairs_imputed.txt"),
#             row.names = F, col.names = T, quote = F, sep = "\t")
# VEP_data2 <- read.table(file=paste0(get.path("locations", local), "functional_annotation/VEP/VEP_data_pairs_imputed.txt"),
#                         h = T, stringsAsFactors = F, sep = "\t")
# 
# VEP_data3 <- VEP_data2[!duplicated(VEP_data2[, c("Variation.Name", "Consequence.to.transcript", "Associated.Gene.Name")]), ]
# table(table(VEP_data3$Variation.Name))
# name.triples <- unique(VEP_data3$Variation.Name)[table(VEP_data3$Variation.Name)==3]
# doubles <- VEP_data3[which(duplicated(VEP_data3[, c("Variation.Name", "Associated.Gene.Name")])),
#                      c("Variation.Name", "Associated.Gene.Name")]
# doubles <- doubles[!(doubles$Variation.Name %in% name.triples), ]
# triples <- VEP_data3[VEP_data3$Variation.Name==name.triples, ]
# triples <- triples[order(triples$Consequence.to.transcript), ]
# rownames(triples) <- NULL
# temp <- VEP_data3[VEP_data3$Variation.Name %in% doubles$Variation.Name, ]
# temp <- temp[order(temp$Consequence.to.transcript), ]
# temp2 <- data.frame(Variation.Name=temp[!duplicated(temp$Variation.Name), "Variation.Name"],
#                     Associated.Gene.Name=temp[!duplicated(temp$Variation.Name), "Associated.Gene.Name"],
#                     VEP1=temp[!duplicated(temp$Variation.Name), "Consequence.to.transcript"],
#                     VEP2=temp[duplicated(temp$Variation.Name), "Consequence.to.transcript"],
#                     VEP3=NA,
#                     Consequence.to.transcript=paste(temp[!duplicated(temp$Variation.Name), "Consequence.to.transcript"],
#                                   temp[duplicated(temp$Variation.Name), "Consequence.to.transcript"],
#                                     sep="/"))
# temp3 <- rbind(temp2,
#                data.frame(Variation.Name=triples[!duplicated(triples$Variation.Name), "Variation.Name"],
#                           Associated.Gene.Name=triples[!duplicated(triples$Variation.Name), "Associated.Gene.Name"],
#                           VEP1=triples[1, "Consequence.to.transcript"],
#                           VEP2=triples[2, "Consequence.to.transcript"],
#                           VEP3=triples[3, "Consequence.to.transcript"],
#                           Consequence.to.transcript=paste(triples[1, "Consequence.to.transcript"],
#                                                           triples[2, "Consequence.to.transcript"],
#                                                           triples[3, "Consequence.to.transcript"],
#                                                           sep="/")))
# temp4 <- merge(VEP_data3[!(VEP_data3$Variation.Name %in% temp3$Variation.Name),
#                          c("Variation.Name", "Consequence.to.transcript", "Associated.Gene.Name")],
#                temp3[,c("Variation.Name", "Consequence.to.transcript", "Associated.Gene.Name")],
#                  all=T, sort=F, stringsAsFactors=F)
# colnames(temp4) <- c("SNP", "VEP", "gene")
# temp5 <- merge(temp4,
#                QTL.snp2[, c("snps", "rs.id")],
#                by.x="SNP", by.y = "rs.id",
#                all.x = T, sort=F)
# temp5 <- temp5[!duplicated(temp5), ]
# colnames(temp5) <- c("rs.id", "VEP", "gene", "SNP")
# 
# write.table(temp5[, c("SNP", "VEP", "gene", "rs.id")],
#             file=paste0(get.path("locations", local), "functional_annotation/VEP/VEP_data_imputed_agg.txt"),
#             row.names = F, col.names = T, quote = F, sep = "\t")
# 
# VEP <- GRanges(
#   seqnames = Rle(VEP_data$Chromosome.name),
#   ranges = IRanges(VEP_data$Position.on.Chromosome..bp.,
#                    end = VEP_data$Position.on.Chromosome..bp.),
#   strand = Rle(VEP_data$Transcript.strand),
#   consequence = VEP_data$Consequence.to.transcript,
#   gene = VEP_data$Ensembl.Gene.ID,
#   transcript = VEP_data$Ensembl.Transcript.ID,
#   symbol = VEP_data$Associated.Gene.Name,
#   rs.id = VEP_data$Variation.Name
# )
# saveRDS(VEP, paste0(get.path("locations", local), "functional_annotation/VEP_GRanges.RDS"))
VEP <- readRDS(paste0(get.path("locations", local), "functional_annotation/VEP_GRanges.RDS"))


# snp-gene pair annotations ----

# # get all cis QTL combinations
# QTL.snp <- readRDS(file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_group.RDS"))
# 
# snp.data <- QTL.snp[, c("snps", "chr", "snp.pos", "gene", "gene.start", "gene.end", "strand", "TSS")]
# 
# 
# # Add MAF, HWE pvalue and VEP
# load(paste0(get.path("genotype", local), "AFHRI_B_imputed_gwaa_data.RData"))
# stats <- summary(AFHRI_B_imp)
# VEP_data <- read.table(file=paste0(get.path("locations", local), "functional_annotation/VEP/VEP_data_imputed_agg.txt"),
#                        h = T, sep = "\t", stringsAsFactors = F)
# 
# snp.data <- merge(snp.data,
#                   data.frame(snps=rownames(stats), stats[, c("Q.2", "Pexact")]),
#                   all.x=T, sort=F)
# 
# snp.data <- merge(snp.data,
#                   VEP_data,
#                   by.x=c("snps", "gene"),
#                   by.y=c("SNP", "gene"),
#                   all.x=T, sort=F)
# 
# snps <- GRanges(
#   seqnames = Rle(snp.data$chr),
#   ranges = IRanges(snp.data$snp.pos,
#                    end = snp.data$snp.pos),
#   SNP = snp.data$snps,
#   symbol = snp.data$gene,
#   gene.start = snp.data$gene.start,
#   gene.end = snp.data$gene.end,
#   TSS = snp.data$TSS,
#   MAF = snp.data$Q.2,
#   HWE = snp.data$Pexact,
#   VEP = snp.data$VEP
# )
# snps <- map.seq.levels(snps)
# 
# saveRDS(snps, paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_raw_temp.RDS"))
# 
# # Add Roadmap Right Atrium 15 state coremarks
# snps$state <- NA
# states <- readRDS(paste0(get.path("locations", local), "functional_annotation/15_states_GRanges.RDS"))
# ovs <- findOverlaps(snps, states)
# if(any(duplicated(queryHits(ovs)))){
#   table(snps$state)
#   test <- data.frame(SNP=snps$SNP[queryHits(ovs)],
#                      state=states[subjectHits(ovs)]$state,
#                      stringsAsFactors = F)
#   test <- test[!duplicated(test), ]
#   saveRDS(test, file=paste0(get.path("locations", local), "functional_annotation/snp_state_link_imputed.RDS"))
#   doubles <- test[duplicated(test$SNP),"SNP"]
#   test2 <- test[test$SNP %in% doubles, ]
#   test3 <- data.frame(SNP=test2[!duplicated(test2$SNP), "SNP"],
#                       state1=test2[!duplicated(test2$SNP), "state"],
#                       state2=test2[duplicated(test2$SNP), "state"],
#                       agg_state=paste(test2[!duplicated(test2$SNP), "state"],
#                                       test2[duplicated(test2$SNP), "state"],
#                                       sep="/"))
#   test4 <- merge(test,
#                  test3[,c("SNP", "agg_state")],
#                  all=T, sort=F, stringsAsFactors=F)
#   test4$agg_state <- as.character(test4$agg_state)
#   test4[is.na(test4$agg_state), "agg_state"] <- test4[is.na(test4$agg_state), "state"]
# 
#   snpstate <- test4[, c("SNP", "agg_state")]
#   snpstate <- snpstate[!duplicated(snpstate), ]
# 
#   saveRDS(snpstate, paste0(get.path("locations", local), "functional_annotation/snp_states_imputed.RDS"))
#   anno.data <- as.data.frame(mcols(snps)[, c("SNP", "symbol")])
#   anno.data <- merge(anno.data, snpstate,
#                     by=c("SNP"),
#                     sort=F, all.x=T)
#   if(identical(as.character(snps$SNP), as.character(anno.data$SNP))){
#     snps$state <- anno.data$agg_state
#   }
# 
# }else{
#   snps$state[queryHits(ovs)] <- states[subjectHits(ovs)]$state
# }
# snps
# saveRDS(snps, paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_raw.RDS"))

snps.raw <- readRDS(paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_raw.RDS"))
genes <- unique(snps$symbol)
length(genes)

# # system(paste("cd projects/symAtrial_QTL/scripts/functional_analysis/",
# #                      "R CMD BATCH build_annotation_matrix_batch_imputed.R &", sep="\n"))

snps <- readRDS(paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_all.RDS"))

hist(snps$dist)
genes <- unique(snps$symbol)

anno.data <- as.data.frame(mcols(snps)[, c("SNP", "symbol", "MAF", "HWE", "nTSS", "dist",
                                           "in.gene", "in.transcript", "UTR3", "UTR5",
                                           "exon", "splice", "RBPBS", "TFBS", "miRBS",
                                           "state", "VEP")])

if(F){
  info <- data.frame(mcols(readRDS(paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_raw.RDS"))))
  
  anno.data[duplicated(anno.data[, c("SNP", "symbol")]), ]
  anno.data[c(355, 356), ]
  anno.data <- anno.data[-356, ]
  
  identical(as.character(info$SNP), as.character(anno.data$SNP))
  anno.data$VEP <- NULL
  anno.data <- merge(anno.data,
                     info[, c("SNP", "symbol", "VEP")],
                     all=T, sort=F,
                     by=c("SNP", "symbol"))
  anno.data$state <- anno.data$state.y
  anno.data$VEP <- anno.data$VEP.y
  anno.data[, c("state.x", "state.y", "VEP.x", "VEP.y")] <- NULL
  saveRDS(anno.data,
          file=paste0(get.path("locations", local),
                      "functional_annotation/snp_imputed_final_anno.RDS"))
}

rm(snps)

QTL.snp <- readRDS(file=paste0(get.path("results", local),
                               "imputed/cis/QTL_res_all_snp_group.RDS"))

QTL.snp <- merge(QTL.snp, anno.data,
                  by.x=c("snps", "gene"), by.y=c("SNP", "symbol"),
                  sort=F, all.x=T)

saveRDS(QTL.snp, file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_group_anno.RDS"))


# add highly expressed genes annotation ----------------------------------------------------------

# get median per sample TPM from GTEx Heart Atrial Appendage tissue
tpm <- read.csv(file=paste0(get.path("gtex", local),
                            "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"),
                sep = "\t", stringsAsFactors = F, skip = 2)
tpm.aa <- tpm[, c("gene_id", "Description", "Heart...Atrial.Appendage")]
library(data.table)
tpm.aa <- as.data.table(tpm.aa)
tpm.aa <- tpm.aa[, meanExpr := mean(Heart...Atrial.Appendage), by = Description]
cutoff <- 1
tpm.aa$highly.expressed <- as.numeric(log(tpm.aa$meanExpr+1) > cutoff)
tpm.aa <- tpm.aa[!duplicated(tpm.aa$Description), c("Description", "meanExpr", "highly.expressed")]
colnames(tpm.aa) <- c("gene", "meanExpr", "hi.ex")

QTL.snp <- readRDS(paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_group_anno.RDS"))

# We want to add: highly.expressed 0/1
QTL.snp2 <- data.frame(merge(QTL.snp, tpm.aa[, c("gene", "hi.ex")],
                             all.x = T),
                       stringsAsFactors = F)
QTL.snp2$hi.ex2 <- as.numeric(QTL.snp2$hi.ex==1 | !is.na(QTL.snp2$pval.pqtl))

saveRDS(QTL.snp2, file=paste0(get.path("results", local), "imputed/cis/QTL_res_all_snp_group_anno_hiex.RDS"))


table(is.na(QTL.snp2$hi.ex))
# # FALSE     TRUE 
# 56263065   187272 

table(QTL.snp2$hi.ex)
# 0        1 
# 22713705 33549360 

table(is.na(QTL.snp2$hi.ex2))
# FALSE     TRUE 
# 56269436   180901 
# > table(QTL.snp2$hi.ex2)

# 0        1 
# 22248366 34021070 


pQTL <- QTL.snp[which(QTL.snp$pQTL==1), ]
pQTL <- pQTL[order(pQTL$pval.pqtl), ]
pQTL2 <- pQTL[!duplicated(pQTL$gene), ]

table(pQTL2$exon)
table(pQTL2$in.gene)
table(pQTL2$in.transcript)
table(pQTL2$TFBS)
table(pQTL2$RBPBS)

table(pQTL2$in.transcript, pQTL2$TFBS)


# filter TFBS and add NKX2-5 annotations ----------
setwd("~/work/symAtrial_QTL/scripts")

library(rtracklayer)
library(GenomicRanges)
library(biomaRt)

source("helper/helper.R")

states <- rtracklayer::import(paste0(get.path("locations", local),
                                     "functional_annotation/roadmap/E104_15/",
                                     "E104_15_coreMarks_dense.bed"))
TFBS <- rtracklayer::import(paste0(get.path("locations", local),
                                   "functional_annotation/TFBS/remap2018_nr_macs2_hg19_v1_2.bed"))
NKX <- rtracklayer::import(paste0(get.path("locations", local),
                                  "functional_annotation/NKX2-5/GSE133833/GSE133833_NKX2-5.bed"))

tpm <- read.csv(file=paste0(get.path("gtex", local),
                            "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"),
                sep = "\t", stringsAsFactors = F, skip = 2)
tpm.aa <- tpm[, c("gene_id", "Description", "Heart...Atrial.Appendage")]
tpm.aa$gene_id.cut <- gsub("\\..*", "", tpm.aa$gene_id)
tpm.aa$high <- as.numeric(tpm.aa$Heart...Atrial.Appendage>1)
table(tpm.aa$high)

tfs <- data.frame(table(TFBS$name), stringsAsFactors = F)

# get the basic ensembl annotations based on GRCh37
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes <- tfs$Var1
# filter for desired gene ids
anno <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                           'start_position','end_position','hgnc_symbol'),
              # filters = 'ensembl_gene_id', values=genes, mart = ensembl)
              filters = 'hgnc_symbol', values=genes, mart = ensembl)
anno$patch <- as.numeric(!(anno$chromosome_name %in% c(as.character(1:22), "X", "Y", "MT")))
table(duplicated(anno$hgnc_symbol), anno$patch)

sum(tfs$Var1 %in% tpm.aa$Description)
sum(tfs$Var1 %in% tpm.aa$Description[tpm.aa$high==1])
sum(anno$ensembl_gene_id %in% tpm.aa$gene_id.cut)
sum(anno$ensembl_gene_id %in% tpm.aa$gene_id.cut[tpm.aa$high==1])

NKX$name <- "NKX2-5"
TF <- c(NKX,
        TFBS[TFBS$name %in% tpm.aa$Description[tpm.aa$high==1], c("name")])

active <- c("1_TssA", "2_TssAFlnk", "10_TssBiv",
            "6_EnhG", "7_Enh", "11_BivFlnk", "12_EnhBiv")

TF.open <- subsetByOverlaps(TF, states[states$name %in% active],
                            ignore.strand=T, minoverlap = 25)

length(TF)
length(TF.open)

# saveRDS(TF.open,
#         paste0(get.path("locations", local),
#                "functional_annotation/TFBS/TFBS_remap_and_NKX2-5_openChromatin_GRanges.RDS"))

TF.open <- readRDS(paste0(get.path("locations", local),
                          "functional_annotation/TFBS/TFBS_remap_and_NKX2-5_openChromatin_GRanges.RDS"))
TF.open

snplocs <- read.csv(get.path("snplocs imputed", local),
                    sep = "\t", h=T, stringsAsFactors = F)
head(snplocs)

snps <- GRanges(
  seqnames = Rle(paste0("chr", snplocs$chrom)),
  ranges = IRanges(snplocs$position, width = 1,
                   names = snplocs$imputed_marker_id),
  name = snplocs$imputed_marker_id,
  snps = snplocs$imputed_marker_id)
snps

snps$TFBS2 <- as.numeric(overlapsAny(snps, TF.open,
                                     maxgap = 0,
                                     type="any",
                                     ignore.strand=T))
snps2 <- as.data.frame(snps)
colnames(snps2)
saveRDS(snps2,
        file = paste0(get.path("locations", local),
                     "functional_annotation/TFBS/SNP_anno_TFBS_remap_and_NKX2-5_openChromatin_GRanges.RDS"))



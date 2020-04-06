# plasma pQTL compairison ----

# compairison files: https://app.box.com/s/u3flbp13zjydegrxjb2uepagp1vb6bj2/
# csv overview: https://app.box.com/s/u3flbp13zjydegrxjb2uepagp1vb6bj2/file/296358772531

setwd("~/work/symAtrial_QTL/scripts")

# library(VariantAnnotation)
# library(GenomicRanges)
#library(plyr)
library(dplyr)
library(tidyr)
library(R.utils)
library(GenomicRanges)

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

if(F){
  # get overview table for the measured plasma proteins
  info <- read.csv(paste0(get.path("plasma"), "SOMALOGIC_GWAS_protein_info.csv"),
                   sep=",", h=T, stringsAsFactors = F, na.strings = c("", " "))
  info$symbol <- gsub("\\.[0-9]+", "", info$SOMAMER_ID)
  
  # get measured proteins in the symatrial AFHRI-B cohort
  load(paste0(get.path("plasma"), "protein_annotation.RData"))
  
  # compute overlaps between the two datasets depending on protein/gene annotations
  proteins <- intersect(unlist(strsplit(as.character(info$UniProt), ",")), as.character(anno$Row))
  genes <- intersect(unlist(strsplit(as.character(info$symbol), ",")), as.character(anno$Gene.name))
  
  # get plasma protein candidates for validation
  val.set <- info[unique(grep(paste(proteins, collapse="|"), info$UniProt),
                         grep(paste(genes, collapse="|"), info$symbol)), ]
  
  # split for different proteins/genes captured by the same aptamer
  val.set2 <- val.set%>%
      mutate(UniProts=strsplit(UniProt,","))%>%
      mutate(symbols=strsplit(symbol,"\\."))%>%
      unnest(UniProts,symbols)
  
  # get chromosome annotations for selecting files necessary for cis pQTL validation
  library(biomaRt)
  #GRCh37.p13
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  
  # chromosome annotation based on gene name
  chr <- getBM(attributes=c("chromosome_name", "external_gene_name"),
               values=genes,
               mart=ensembl,
               bmHeader=T)
  colnames(chr) <- c("Chromosome", "symbol")
  chr <- chr[chr$Chromosome %in% c(1:30, "X", "Y", "M", "MT"), ]
  
  # chromosome annotation based on Uniprot IDs
  chr2 <- getBM(attributes=c("chromosome_name", "uniprotswissprot"),
                values=proteins,
                mart=ensembl,
                bmHeader=T)
  colnames(chr2) <- c("Chromosome", "uniprot_id")
  chr2 <- chr2[chr2$Chromosome %in% c(1:30, "X", "Y", "M", "MT"), ]
  
  # merge with our validation candidates
  val.set3 <- merge(val.set2, chr,
                    by.x="symbols", by.y="symbol",
                    all.x=T, sort=F)
  val.set4 <- merge(val.set3, chr2,
                    by.x="UniProts", by.y="uniprot_id",
                    all.x=T, sort=F)
  
  # check, if annotations are complete
  #crit <- val.set4[is.na(val.set4$Chromosome.x) | is.na(val.set4$Chromosome.y), ]
  
  # get final dataset with all information
  val.set5 <- val.set4
  val.set5[is.na(val.set5$Chromosome.x), "Chromosome.x"] <- val.set5[is.na(val.set5$Chromosome.x), "Chromosome.y"]
  val.set5[is.na(val.set5$Chromosome.y), "Chromosome.y"] <- val.set5[is.na(val.set5$Chromosome.y), "Chromosome.x"]
  val.set5 <- val.set5[!duplicated(val.set5), ]
  
  # write it out
  write.table(val.set5, file=paste0(get.path("plasma"), "validation_set.tsv"),
              row.names = F, col.names = T, sep = "\t", quote = F)
  
  # create the full list of required files
  files <- rbind(val.set5[, c("SOMAMER_ID", "UniProt", "symbol", "TargetFullName", "Target", "Chromosome.y")],
                 val.set5[, c("SOMAMER_ID", "UniProt", "symbol", "TargetFullName", "Target", "Chromosome.y")])
  files[1:dim(val.set5)[1], ] <- val.set5[, c("SOMAMER_ID", "UniProt", "symbol", "TargetFullName", "Target", "Chromosome.x")]
  colnames(files) <- c("SOMAMER_ID", "UniProt", "symbol", "TargetFullName", "Target", "Chromosome")
  files <- files[!duplicated(files), ]
  
  # results only available for Chromosomes 1:22
  apt.x <- files[files$Chromosome=="X", "SOMAMER_ID"]
  files[files$SOMAMER_ID %in% apt.x, ]
  
  apt.non.x <- files[files$Chromosome!="X", "SOMAMER_ID"]
  
  files.final <- files[files$SOMAMER_ID %in% apt.non.x, ]
  files.final$filepath <- paste0(files.final$SOMAMER_ID, "/", files.final$SOMAMER_ID, "_chrom_", files.final$Chromosome, "_meta_final_v1.tsv")
  files.final$filename <- paste0(files.final$SOMAMER_ID, "_chrom_", files.final$Chromosome, "_meta_final_v1.tsv")
  files.final <- files.final[order(files.final$SOMAMER_ID, as.numeric(as.character(files.final$Chromosome))),]
  
  # final set of validation candidates
  dim(files.final)
  write.table(files.final, file=paste0(get.path("plasma"), "filelist.tsv"),
              row.names = F, col.names = T, sep = "\t", quote = F)
  write.table(files.final, file=paste0(get.path("plasma"), "filelist.csv"),
              row.names = F, col.names = T, sep = ",", quote = T)
}

plasma <- read.csv(file=paste0(get.path("plasma"),
                                 "filelist.tsv"),
                     h = T, sep = "\t", stringsAsFactors = F)
plasma2 <- plasma[, c("SOMAMER_ID", "UniProt", "symbol", "Chromosome", "filepath", "filename")]%>%
  mutate(UniProts=strsplit(UniProt,","))%>%
  mutate(symbols=strsplit(symbol,"\\."))%>%
  unnest(UniProts,symbols)

# Strategy
# 1. unzip all files
# 2. merge files: only cis-chromosome per SOMAMER-ID
# 2. b) SNPs in cis range only
# 3. create SNP-id that can be matched
# 4. filter lists for gene symbols that are in data set
# 5. merge gene symbol to somamer data

if(F){
  # 1. unzip all files
  zipF <- list.files(path = paste0(get.path("plasma"), "plasma_pQTL/"), pattern = "*.zip", full.names = TRUE)
  
  # unzip all your files
  ldply(.data = zipF, .fun = unzip, exdir = paste0(get.path("plasma"), "plasma_pQTL/"))
  apply(plasma, 1,
        FUN= function(x) gunzip(filename = paste0(get.path("plasma"), "plasma_pQTL/", x[7], ".gz"),
                                destname = paste0(get.path("plasma"), "plasma_pQTL/cis/", x[8]),
                                skip=T, remove = F, overwrite = T))
  zipF.rm <- sub(".zip", "", list.files(path = paste0(get.path("plasma"), "plasma_pQTL/"), pattern = "*.zip", full.names = T))
  unlink(zipF.rm, recursive=T)
}


snps <- read.table(paste0(get.path("genotype", local),
                          "AFHRI-B_imputed.bim"),
                   header=F, sep = "\t", stringsAsFactors = F)
snps$V3 <-NULL
colnames(snps) <- c("chr", "snpid", "variant_pos", "alt", "ref")
snps$variant_id <- paste(snps$chr, snps$variant_pos,
                         snps$alt, snps$ref, sep=":")
snps$rsid <- NA
snps$rsid[grep("rs", snps$snpid)] <- gsub("\\:.*", "", snps$snpid[grep("rs", snps$snpid)])



##### add gene positions
library(biomaRt)

# get the basic ensembl annotations based on GRCh37
ensembl = useEnsembl(biomart="ensembl",
                     dataset="hsapiens_gene_ensembl",
                     GRCh=37)

# define genes of interest with ensembl ids
genes <- unique(plasma2$symbols)

# filter for desired gene ids
anno <- getBM(attributes=c('hgnc_symbol',
                           'chromosome_name',
                           'start_position',
                           'end_position'),
              filters = 'hgnc_symbol', values=genes,
              mart = ensembl)
#chr6:112,429,134-112,576,141
anno <- anno[anno$chromosome_name %in% as.character(1:22), ]

anno <- GRanges(seqnames = Rle(anno$chromosome_name),
                ranges = IRanges(anno$start_position, end = anno$end_position),
                hgnc_symbol = anno$hgnc_symbol,
                chromosome_name = anno$chromosome_name)
#######


file <- read.csv(paste0(get.path("plasma", local),
                        "plasma_pQTL/cis/", plasma$filename[1]),
                 sep="\t", h=T, stringsAsFactors = F)
file$variant_id <- paste(file$chromosome, file$position,
                         toupper(file$Allele2), toupper(file$Allele1), sep=":")
file$variant_id2 <- paste(file$chromosome, file$position,
                          toupper(file$Allele1), toupper(file$Allele2), sep=":")
file$plasma.match_alleles <- NA
file$plasma.match_alleles[file$variant_id2 %in% snps$variant_id] <- -1
file$plasma.match_alleles[file$variant_id %in% snps$variant_id] <- 1
file <- file[!is.na(file$plasma.match_alleles), ]
file.GR <- GRanges(seqnames = Rle(file$chromosome),
                   ranges = IRanges(file$position, end = file$position),
                   VARIANT_ID = file$VARIANT_ID,
                   plasma.match_alleles = file$plasma.match_alleles)
variants <- subsetByOverlaps(file.GR,
                            anno[anno$hgnc_symbol %in% plasma2[plasma2$SOMAMER_ID==plasma$SOMAMER_ID[1], "symbol"]],
                            maxgap = 1.01e6,
                            type = "any",
                            ignore.strand=T)$VARIANT_ID
file <- file[file$VARIANT_ID %in% variants, ]

# write.table(data.frame(SOMAMER_ID=plasma$SOMAMER_ID[1],
#                        file, stringsAsFactors = F),
#             file = paste0(get.path("plasma"), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt"),
#             row.names = F, col.names = T, sep = "\t", quote = F)

for(i in 244:dim(plasma)[1]){     # i in 2:length(cisF)
  file <- read.csv(paste0(get.path("plasma", local),
                          "plasma_pQTL/cis/", plasma$filename[i]),
                   sep="\t", h=T, stringsAsFactors = F)
  file$variant_id <- paste(file$chromosome, file$position,
                           toupper(file$Allele2), toupper(file$Allele1), sep=":")
  file$variant_id2 <- paste(file$chromosome, file$position,
                            toupper(file$Allele1), toupper(file$Allele2), sep=":")
  file$plasma.match_alleles <- NA
  file$plasma.match_alleles[file$variant_id2 %in% snps$variant_id] <- -1
  file$plasma.match_alleles[file$variant_id %in% snps$variant_id] <- 1
  file <- file[!is.na(file$plasma.match_alleles), ]
  file.GR <- GRanges(seqnames = Rle(file$chromosome),
                     ranges = IRanges(file$position, end = file$position),
                     VARIANT_ID = file$VARIANT_ID,
                     plasma.match_alleles = file$plasma.match_alleles)
  variants <- subsetByOverlaps(file.GR,
                               anno[anno$hgnc_symbol %in% plasma2[plasma2$SOMAMER_ID==plasma$SOMAMER_ID[i], "symbols"]],
                               maxgap = 1.01e6,
                               type = "any",
                               ignore.strand=T)$VARIANT_ID
  file <- file[file$VARIANT_ID %in% variants, ]
  if(dim(file)[1]>0){
    write.table(data.frame(SOMAMER_ID=plasma$SOMAMER_ID[i],
                           file, stringsAsFactors = F),
                file = paste0(get.path("plasma"), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt"),
                append = T, row.names = F, col.names = F, sep = "\t", quote = F)
  }
  print(paste0("File ", i, " / ", dim(plasma)[1]))
}

ppQTL <- read.table(paste0(get.path("plasma"), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt"),
                    h = T, sep = "\t")

pexpr <- readRDS(paste0(get.path("dataset", local),
                        "AFHRI_B_proteomics_QC_symbol.RDS"))
genes <- rownames(pexpr)

plasma.genes <- plasma2[plasma2$symbols %in% genes, c("SOMAMER_ID", "symbols")]
plasma.genes <- plasma.genes[!duplicated(plasma.genes), ]
colnames(plasma.genes) <- c("SOMAMER_ID", "gene")

ppQTL <- merge(ppQTL, plasma.genes,
               all=T)
write.table(ppQTL,
            file = paste0(get.path("plasma"), "plasma_pQTL/cis/plasma_pQTL_allpairs.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)

write.table(ppQTL[ppQTL$log.P. < log(1e-5), ],
            file = paste0(get.path("plasma"),
                          "plasma_pQTL/cis/plasma_pQTL_allpairs.significant.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)

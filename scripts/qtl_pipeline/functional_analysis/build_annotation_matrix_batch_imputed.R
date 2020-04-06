setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")
# setwd("/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/scripts")

# load libraries
library(BatchJobs, lib.loc = "/storage/groups/epigenereg01/tools/2017/R/3.4/")

# source scripts
source(paste0("/home/icb/ines.assum/projects/symAtrial_QTL/scripts/BatchJobs.R"))

source("helper/helper.R")

do.stuff <- function(k){

  setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")
  # setwd("/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/scripts")
  
  library(VariantAnnotation)
  library(GenomicRanges)
  
  source("helper/helper.R")
  
  local=F
  
  exon <- readRDS(paste0(get.path("locations", local), "functional_annotation/exon_ids_GRanges.RDS"))
  splice <- readRDS(paste0(get.path("locations", local), "functional_annotation/splicesites_GRanges.RDS"))
  UTR <- readRDS(paste0(get.path("locations", local), "functional_annotation/UTR_GRanges.RDS"))
  locs <- readRDS(paste0(get.path("locations", local), "functional_annotation/locs_GRanges.RDS"))
  RBPBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/RBP/RBPBS_heart_GRanges.RDS"))
  TFBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/TFBS/TFBS_GRanges.RDS"))
  miRBS <- readRDS(paste0(get.path("locations", local), "functional_annotation/miRBS_GRanges.RDS"))
  
  # snp-gene pair annotations ----
  snps <- readRDS(paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_raw.RDS"))
  genes <- unique(snps$symbol)
  
  # k=1
  k1 <- (k-1)*500 + 1
  k2 <- k1+499
  if (k==33){
    k2 <- length(genes)
  }
  
  # extract necessary SNPs
  snps <- snps[snps$symbol %in% genes[k1:k2]]
  
  # SNP within gene?
  snps$in.gene <- as.numeric(NA)
  
  # SNP within any transcript annotated to the gene symbol?
  snps$in.transcript <- as.numeric(NA)
  
  # SNP within any exon annotated to a transcript of this gene symbol?
  snps$exon <- as.numeric(NA)
  
  # SNP within a 10bp region of a splice site of a transcript of this gene symbol?
  snps$splice <- as.numeric(NA)
  
  # SNP located in any transcript annotated to the gene symbol AND overlapping with RBP BS?
  snps$RBPBS <- as.numeric(NA)
  
  # SNP overlapping with a TF BS?
  snps$TFBS <- as.numeric(NA)
  
  # SNP overlapping with a miRNA BS & 3p UTR?
  snps$miRBS <- as.numeric(NA)
  
  # SNP located in the 3'/5' UTR of a transcript of this gene symbol?
  snps$UTR3 <- as.numeric(NA)
  snps$UTR5 <- as.numeric(NA)
  
  # Genomic location in bp of the closest TSS belonging to a transcript of this gene
  snps$nTSS <- as.numeric(NA)
  
  # strand-matched distance to the closest TSS (in bp, - = upstream, + = downstream )
  snps$dist <- as.numeric(NA)
  
  for (i in k1:k2){
    #i=k1+1
    # select all snp-gene pairs with gene i
    snps.g <- snps[snps$symbol==genes[i]]
    
    # subset locations annotated to gene i:
    # transcript start and stop positions in bp
    locs.t.g <- locs[locs$symbol==genes[i] & locs$type=="transcript"]
    # gene start and stop positions in bp
    locs.g.g <- locs[locs$symbol==genes[i] & locs$type=="gene"]
    locs.g.g <- locs.g.g[!duplicated(locs.g.g)]
    
    # TSS locations for transcripts of gene i (i.e. start/end = start of transcript)
    locs.tss <- locs.t.g
    end(locs.tss) <- start(locs.tss)
    
    # get annotations :-)
    # in.gene
    snps.g$in.gene <- as.numeric(overlapsAny(snps.g,
                                             locs.g.g,
                                             maxgap = 0,
                                             type="any",
                                             ignore.strand=T))

    # in.transcript
    snps.g$in.transcript <- as.numeric(overlapsAny(snps.g,
                                                   locs.t.g,
                                                   maxgap = 0,
                                                   type="any",
                                                   ignore.strand=T))

    # exon
    snps.g$exon <- as.numeric(overlapsAny(snps.g,
                                          exon[exon$symbol==genes[i]],
                                          maxgap = 0,
                                          type="any",
                                          ignore.strand=T))

    # splice
    snps.g$splice <- as.numeric(overlapsAny(snps.g,
                                            splice[splice$symbol==genes[i]],
                                            maxgap = 9,
                                            type="any",
                                            ignore.strand=T))

    # RBP BS
    snps.g$RBPBS <- as.numeric(overlapsAny(snps.g,
                                           subsetByOverlaps(RBPBS,
                                                            locs.t.g,
                                                            maxgap = 0,
                                                            type = "any",
                                                            ignore.strand=F),
                                           maxgap = 0,
                                           type="any",
                                           ignore.strand=T))

    # TF BS
    snps.g$TFBS <- as.numeric(overlapsAny(snps.g,
                                          TFBS,
                                          maxgap = 0,
                                          type="any",
                                          ignore.strand=T))

    # UTR3
    snps.g$UTR3 <- as.numeric(overlapsAny(snps.g,
                                          UTR[UTR$symbol==genes[i] & UTR$feature=="utr3"],
                                          maxgap = 0,
                                          type="any",
                                          ignore.strand=T))
    # UTR5
    snps.g$UTR5 <- as.numeric(overlapsAny(snps.g,
                                          UTR[UTR$symbol==genes[i] & UTR$feature=="utr5"],
                                          maxgap = 0,
                                          type="any",
                                          ignore.strand=T))
    
    # miRNA BS
    snps.g$miRBS <- as.numeric(overlapsAny(snps.g,
                                           subsetByOverlaps(miRBS,
                                                            UTR[UTR$symbol==genes[i] & UTR$feature=="utr3"],
                                                            maxgap = 0,
                                                            type = "any",
                                                            ignore.strand=F),
                                           maxgap = 0,
                                           type="any",
                                           ignore.strand=T))
    
    # nearest TSS
    nTSS.transcript <- locs.tss[nearest(snps.g,
                                        locs.tss,
                                        select = "arbitrary",
                                        ignore.strand=T)]
    snps.g$nTSS <- as.numeric(start(nTSS.transcript))
    
    # dist to nearest TSS
    snps.g$dist <- (-1) * (start(snps.g) - start(nTSS.transcript)) *
                            ifelse(strand(nTSS.transcript)=="-", -1, 1)
    
    # return results for gene i
    snps[snps$symbol==genes[i]] <- snps.g
    
    print(paste0(i-k1+1, " / ", k2-k1+1))
  }
  # save results for this part
  saveRDS(snps, paste0(get.path("locations", local), "functional_annotation/temp_imputed/snp_imputed_anno_", k1, "_", k2, ".RDS"))
}

local=F
snps <- readRDS(paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_raw.RDS"))
genes <- unique(snps$symbol)
run.batchjobs(do.stuff, 1:33, more.args=list(), "anno", "tmp_dir_anno_imputed", clean.up=T,
              resources = list(queue="long_fed25", memory="40G"))

snps <- readRDS(paste0(get.path("locations", local), "functional_annotation/temp_imputed/snp_imputed_anno_1_500.RDS"))
for (k in 2:33){
  k1 <- (k-1)*500 + 1
  k2 <- k1+499
  if (k==33){
    k2 <- length(genes)
  }
  temp <- readRDS(paste0(get.path("locations", local), "functional_annotation/temp_imputed/snp_imputed_anno_", k1, "_", k2, ".RDS"))
  snps <- c(snps, temp)
}
saveRDS(snps, paste0(get.path("locations", local), "functional_annotation/snp_imputed_anno_all.RDS"))


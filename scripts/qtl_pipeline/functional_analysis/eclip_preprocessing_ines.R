###################################################### FUNCTIONS #################################################
##################################################################################################################
###
# This function creates a comprehensive gene annotation file
###
load.gene.annot <- function(){

  target.file <- paste(INPUT, GENEANNOT, "Gencode/comprehensive_gene_annotation_hg19.RDS", sep="")
  
  # Build list of annotations by annotation type
  annotations <- vector("list", 5)
  names(annotations) <- c("five_prime", 
                          "exons", 
                          "introns", 
                          "coding.exons", 
                          "three_prime")

  # Create transcript database object from genecode gff file
  TxDb <- makeTxDbFromGFF(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.annotation.gff3", sep =""))
  genecode.annot <- import(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.annotation.gff3", sep =""))

  # Retrieve genes
  genes.data <- genes(TxDb, columns='tx_name')
  # Change naming of col "tx_name" to "transcript id
  names(mcols(genes.data))[1] <- "transcript_id"
  # Append gene_id column
  mcols(genes.data)["gene_id"] <- names(genes.data)
  # Remove version numbers from ids
  mcols(genes.data)["transcript_id"] <- strip.version(mcols(genes.data)[["transcript_id"]], names = T)
  mcols(genes.data)["gene_id"] <-  strip.version(mcols(genes.data)[["gene_id"]])
  # Append gene type
  mcols(genes.data)[, "gene_type"] <- mcols(genecode.annot)[ match(names(genes.data), mcols(genecode.annot)[, "gene_id"]),"gene_type"]
  # Reformat missing gene type annotations (NA to "none")
  genes.data$gene_type[which(is.na(genes.data$gene_type))] <- "none"
  # Append gene and transcript names
  mcols(genes.data)[, "gene_name"] <-  mcols(genecode.annot)[ match(names(genes.data), mcols(genecode.annot)[, "gene_id"]),"gene_name"]
  mcols(genes.data)[, "transcript_name"] <-  mcols(genecode.annot)[ match(names(genes.data), mcols(genecode.annot)[, "gene_id"]),"transcript_name"]
  # Retrieve transcripts
  transcripts.data <- transcripts(TxDb, columns=c("tx_name", "gene_id"))
  # Change naming of col "tx_name" to "transcript id
  names(mcols(transcripts.data))[1] <- "transcript_id"
  # Append transcript type
  mcols(transcripts.data)[, "transcript_type"] <- mcols(genecode.annot)[ match(transcripts.data$transcript_id, mcols(genecode.annot)[, "transcript_id"]),"transcript_type"]
  # Append gene and transcript names
  mcols(transcripts.data)[, "gene_name"] <-  mcols(genecode.annot)[ match(transcripts.data$transcript_id, mcols(genecode.annot)[, "transcript_id"]),"gene_name"]
  mcols(transcripts.data)[, "transcript_name"] <-  mcols(genecode.annot)[ match(transcripts.data$transcript_id, mcols(genecode.annot)[, "transcript_id"]),"transcript_name"]
  # Append gene ids to names attribute
  #names(transcripts.data) <- mcols(transcripts.data)[, "gene_id"]
  # Remove version numbers from ids
  mcols(transcripts.data)["transcript_id"] <- strip.version(mcols(transcripts.data)[["transcript_id"]])
  mcols(transcripts.data)["gene_id"] <- strip.version(mcols(transcripts.data)[["gene_id"]], names = T)
  # Mapping of transcripts to genes
  tx2gene <- unlist(elementMetadata(transcripts.data)$gene_id)
  names(tx2gene) <- elementMetadata(transcripts.data)$transcript_id
  ## Extract relevant annotation types
  # Retrieve 3' UTR
  threeUTRs.data <- threeUTRsByTranscript(TxDb, use.names=T)
  names(threeUTRs.data) <- strip.version(names(threeUTRs.data))
  #names(threeUTRs.data) <- tx2gene[names(threeUTRs.data)]
  annotations$three_prime <- threeUTRs.data
  # Retrieve 5' UTR
  fiveUTRs.data <- fiveUTRsByTranscript(TxDb, use.names=T)
  names(fiveUTRs.data) <- strip.version(names(fiveUTRs.data))
  #names(fiveUTRs.data) <- tx2gene[names(fiveUTRs.data)]
  annotations$five_prime <- fiveUTRs.data
  ## Retrieve exons
  exons.data <- exonsBy(TxDb, by='tx', use.names=T)
  names(exons.data) <- strip.version(names(exons.data))
  #names(exons.data) <- tx2gene[names(exons.data)]
  annotations$exons <- exons.data
  ## Retrieve introns
  introns.data <- intronsByTranscript(TxDb, use.names=T)
  ## Append intron identifiers and ranks (analogous to exon_name)
  num.introns.per.tx <- lengths(introns.data)
  intron.ranks <- unlist(lapply(num.introns.per.tx, function(num.introns){
    if(num.introns == 0){return(0)}
    1:num.introns
  }))
  # Exclude transcripts with no intron
  intron.ranks <- intron.ranks[intron.ranks !=0]
  # Append intron id (Analogous to exon id)
  temp <- unlist(introns.data)
  mcols(temp)["intron_id"] <- 1:length(temp)
  mcols(temp)["intron_name"] <- paste(names(temp), intron.ranks, sep = ":" )
  mcols(temp)["intron_rank"] <- intron.ranks
  # Split intron data by transcripts
  introns.data <- split(temp, names(temp))
  # Strip version number from transcripts ids
  names(introns.data) <- strip.version(names(introns.data))
  #names(introns.data) <- tx2gene[names(introns.data)]
  annotations$introns <- introns.data
  ## Retrieve coding sequences
  cds.data <- cdsBy(TxDb, by='tx', use.names=T)
  names(cds.data) <- strip.version(names(cds.data))
  #names(cds.data) <- tx2gene[names(cds.data)]
  annotations$coding.exons <- cds.data
  # Reformat annotations to GRanges
  annotations <- lapply(annotations, function(annot){
    annot.reformatted <- unlist(annot, use.names=FALSE)
    names(annot.reformatted) <- rep(names(annot), elementNROWS(annot))
    annot.reformatted
  })
  # Strip off version numbers from sub region annotations
  annotations <- lapply(annotations, function(annot){
    mcols(annot)[, 2] <- strip.version(mcols(annot)[, 2], pattern = "\\.[0-9]+:", replacement = ":", names = F)
    annot
  })
  # Retrieve genes
  annotations$genes <- genes.data
  # Retrieve transcripts
  annotations$transcripts <- transcripts.data
  # Remove version numbers from all gene ids in each GRanges object
  annotations <- lapply(annotations, function(annot){
    names(annot) <- strip.version(names(annot))
    annot
  })
  names(annotations$transcripts) <- NULL
  names(annotations$genes) <- NULL
  
  ## Create mapping of transcript ids to gene ids (along with ref seq ids)
  # Mapping of transcripts to genes
  id.map <- unlist(elementMetadata(annotations$transcripts)$gene_id)
  names(id.map) <- elementMetadata(annotations$transcripts)$transcript_id
  id.map <- as.data.frame(id.map)
  # Append column for transcript id
  id.map$ensembl_transcript_id <- rownames(id.map)
  # Append column for protein id
  id.map$ensembl_protein_id <- genecode.annot$protein_id[match(id.map$ensembl_transcript_id, strip.version(genecode.annot$transcript_id))]
  id.map$ensembl_protein_id <- strip.version(id.map$ensembl_protein_id)
  # Load refseq ids
  ensemble.refseq.mapping <- read.table(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.metadata.RefSeq",
                                              sep =""), sep = "\t")
  # Delete version numbers from ids
  ensemble.refseq.mapping <- as.data.frame(apply(ensemble.refseq.mapping, 2, function(ids){gsub("\\..*", "", ids)}))
  colnames(ensemble.refseq.mapping) <- c("ensembl_transcript_id", "refseq_transcript_id", "refseq_protein_id")
  ## Extend target df by refseq annotaitons
  id.map$refseq_transcript_id <- ensemble.refseq.mapping$refseq_transcript_id[match(id.map$ensembl_transcript_id, ensemble.refseq.mapping$ensembl_transcript_id)]
  id.map$refseq_protein_id <-ensemble.refseq.mapping$refseq_protein_id[match(id.map$ensembl_transcript_id, ensemble.refseq.mapping$ensembl_transcript_id)]
  colnames(id.map) <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_protein_id", "refseq_transcript_id", "refseq_protein_id")
  # Retrieve and append hgnc symbols
  ensembl2hgnc <- read.table(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.metadata.HGNC", sep =""), sep = "\t")
  colnames(ensembl2hgnc) <- c("ensemble_transcript_id", "hgnc")
  ensembl2hgnc$ensemble_transcript_id <- strip.version(ensembl2hgnc$ensemble_transcript_id)
  id.map$hgnc_symbol <- ensembl2hgnc$hgnc[match(id.map$ensembl_transcript_id, ensembl2hgnc$ensemble_transcript_id)]
  # Save id mapping to annotations
  annotations$id_map <- id.map
  
  # Append list of refseq transcript identifier for quick access
  refseq.transcripts <- GENE.ANNOT$id_map$ensembl_transcript_id[!is.na(GENE.ANNOT$id_map$refseq_transcript_id)]
  refseq.supported.ensbl.transcripts <- as.character(refseq.transcripts)
  annotations$refseq_transcript_support <- refseq.supported.ensbl.transcripts

  ## Save data sets
  # Define output directory for sub data sets
  out.dir <- paste(INPUT, GENEANNOT, "Gencode", sep="")
  lapply(names(annotations), function(annot){
    saveRDS(annotations[[annot]], paste(out.dir, "/", annot, ".RDS", sep = ""))
    return()
  })
  saveRDS(annotations, target.file)
  return(annotations)
}

###
# Removes a pattern in a character vector and optionally from names attribute
##
strip.version <- function(id, pattern = "\\..*", replacement = "", names = F){
  if(names & (!is.null(names(id)))){
    names(id)  <- gsub(pattern, replacement, names(id))
    return(gsub(pattern, replacement = replacement, id))
  }else{
    return(gsub(pattern, replacement = replacement, id))
  }
}

###
# This function loads the ENCODE Clipper peak calls from the eCLIP-seq bed-file data sets
###
load.eCLIPseq.data <- function(){

  target.file <- paste(INPUT, eCLIP, "/eClip.calls", ".RDS", sep = "")

  eClip.calls <- lapply(CARGS$cell.line, function(cline) {
    
    # Get eClip file paths for both replicates per cell line per rbp
    files <- list(list.files(paste(INPUT, eCLIP, cline, "/non_controls/bed/rep1", sep = ""), full.names = T),
                  list.files(paste(INPUT, eCLIP, cline, "/non_controls/bed/rep2", sep = ""), full.names = T))
    
    ## Parse data
    cnames <- c("chrom", "start", "end", "target", 
                "score", "strand", "l2_eClip_VS_SMinput_enrichment", 
                "l10_eClip_VS_SMinput_fisher", "V9", "V10")
    clusterExport(WORKERS$high, list("cnames"), environment())
    eClip.calls <- lapply(files, function(file.paths){
      parLapply(WORKERS$high, file.paths, function(pData, pCnames = cnames){
        eClip.data <- read.table(gzfile(pData), header = F)
        #! columns "V9" and "V10" do not convey any information and 
        #  exist to maintain the bed format and thus are ignored
        colnames(eClip.data) <- pCnames
        GRanges(seqnames = eClip.data$chrom, 
                strand = eClip.data$strand, 
                ranges = IRanges(start = eClip.data$start, end = eClip.data$end), 
                target = eClip.data$target,
                score = eClip.data$score,
                enrichment = eClip.data$l2_eClip_VS_SMinput_enrichment,
                significance = eClip.data$l10_eClip_VS_SMinput_fisher)
      })
    })
    # Append gene symbols to list name attribute
    eClip.calls <- lapply(seq_along(eClip.calls), function(index){
      names(eClip.calls[[index]]) <- gsub(".bed.gz", "", basename(files[[index]]))
      eClip.calls[[index]]
    })
    names(eClip.calls) <- c("rep1", "rep2")
    eClip.calls
  })
  names(eClip.calls) <- CARGS$cell.line
  saveRDS(eClip.calls, target.file)
  return(eClip.calls)
}

###
# This function takes as input the ENCODE Clipper peak calls from the eCLIP-seq data
# and filters the calls by signficant, enriched and overlapping binding sites between replicates
###
filter.eClip.calls <- function(clip.calls, alpha, ratio){
  target.file <- paste(INPUT, eCLIP, "/merged.significant.eClip.calls", ".RDS", sep = "")
  merged.sig.eClip.calls <- lapply(names(clip.calls), function(cline){
    
    ## Filter clip calls by signal strength, i.e. enrichment and significance of the observed
    # clip signals relative to the size-matched-input control (SMinput)
    sig.eClip.calls <- lapply(clip.calls[[cline]], function(replicate.calls){
      lapply(replicate.calls, function(call){
        call[mcols(call)[ ,"enrichment"] > log2(ratio) & 
               mcols(call)[ ,"significance"] > -log10(alpha), ]
      })
    })
    ## Merge clip calls by overlapping replicates at a minimum overlap threshold
    # Get names of target rbp to iterate over
    targets <- names(sig.eClip.calls$rep1)
    merged.sig.eClip.calls <- lapply(targets, function(rbp){
      # Get peak calls for a specific target from both replicates
      rep1.calls <- sig.eClip.calls$rep1[[rbp]]
      rep2.calls <- sig.eClip.calls$rep2[[rbp]]
      # Overlap peak calls and retain those with high overlap
      reduce(reduceGRanges(rep1.calls, rep2.calls, min = 0.5, intersect = T, both = F))
    })
    names(merged.sig.eClip.calls) <- targets
    # Append gene ids and gene name metadata columns
    merged.sig.eClip.calls <- lapply(names(merged.sig.eClip.calls), function(rbp){
      mcols(merged.sig.eClip.calls[[rbp]])["ensembl_id"] <- GENE.ANNOT$id_map$ensembl_gene_id[which(GENE.ANNOT$id_map$hgnc == rbp, arr.ind = T)[1]]
      mcols(merged.sig.eClip.calls[[rbp]])["hgnc_symbol"] <- rbp
      merged.sig.eClip.calls[[rbp]]
    })
    # Merge, save and return data
    merged.sig.eClip.calls <- do.call(c, merged.sig.eClip.calls)
    ## Reformat name attribute, i.e. remove replicate info and seperate out cell info
    mcols(merged.sig.eClip.calls)["cell_line"] <- cline
    return(merged.sig.eClip.calls)
  })
  names(merged.sig.eClip.calls) <- CARGS$cell.line
  # Save intermediate input data
  saveRDS(merged.sig.eClip.calls, target.file)
  return(merged.sig.eClip.calls)
}

###
# As opposed to the GRanges 'reduce' or 'findOverlaps' function, this function 
# allows relative overlapping and subsequent merging of overlapping ranges
##
reduceGRanges <- function(x, y, min=0.6, intersect = T, both = F) {

  hits <- findOverlaps(x, y)
  xhits <- x[queryHits(hits)]
  yhits <- y[subjectHits(hits)]
  # Fragments from both replicates must meet the "min" criterion
  if(both){
    overlap <- width(pintersect(xhits, yhits))
    frac1 <- overlap/ width(xhits)
    frac2 <- overlap/ width(yhits)
    frac1 <- frac1 >= min
    frac2 <- frac2 >= min
    merge <- frac1 & frac2
  # Only the smaller fragment must fulfill the min criterion
  }else{
    frac <- width(pintersect(xhits, yhits)) / pmin(width(xhits), width(yhits))
    merge <- frac >= min
  }
  # Return merged fragments
  if(!intersect){
    return(reduce(c(xhits[merge], yhits[merge])))
  }else{
  # Return only overlapping region of fragments
    return(pintersect(xhits, yhits)[merge])
  }
}

###################################################### MAIN ######################################################
##################################################################################################################
library(GenomicRanges)
library(GenomicFeatures)

CARGS <- list(cell.line = c("K562", "HepG2"), 
              new       = T,
              workers   = setNames(c(6,12,18),c("low", "avg", "high")))

#! Set paths, must contain suffix "/"
INPUT = ""
OUTPUT = ""
GENEANNOT <- "GeneAnnotations/"
eCLIP <- "eClip/"
dir.create(paste0(INPUT, GENEANNOT, "Gencode"), recursive = T)
dir.create(paste0(INPUT, eCLIP, CARGS$cline[1], "/non_controls/bed/rep1"), recursive = T)
dir.create(paste0(INPUT, eCLIP, CARGS$cline[2], "/non_controls/bed/rep2"), recursive = T)

GENE.ANNOT <- load.gene.annot()
eCLIPseq <- load.eCLIPseq.data()
eCLIPseq <- filter.eClip.calls(eCLIPseq, 0.05, 1)



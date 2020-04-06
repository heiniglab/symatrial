# ------------------------------------------------------------------------------
#' Functional enrichment of regulatory/posttranscriptional elements
#' for cis eQTL/pQTL hits (mQTL hits) based on annotated snp-gene pairs
#' TODO: including LD info
#'
#' @author Ines Assum
#'
# ------------------------------------------------------------------------------

# Let's start

setwd("~/work/symAtrial_QTL/scripts")

library(ggplot2)
library(GenomicRanges)

source("functional_analysis/enrichment_analysis/enrichment_functions.R")
source("helper/helper.R")

local=F

# enrichment cases:
# FDR.pqtl, FDR.res_pqtl -> Group5/NULL -> hit/bg -> all annotations
# FDR.eqtl, FDR.res_eqtl -> Group1/NULL -> hit/bg -> all annotations

input.path <- paste0(get.path("results", local),
                     "imputed/cis/enrichment/input")
dir.create(input.path, showWarnings = F, recursive = T)
result.path <- paste0(get.path("results", local),
                      "imputed/cis/enrichment/results")
dir.create(result.path, showWarnings = F, recursive = T)

types1 <- c("FDR.eqtl", "FDR.res_eqtl")
types2 <- c("FDR.pqtl", "FDR.res_pqtl")
types3 <- c("FDR.eqtl", "FDR.pqtl")
types4 <- c("FDR.ratios")
cond1 <- c("Group1", NA)
cond2 <- c("Group5", NA)
cond3 <- c("Group2")
cond4 <- c(NA)

batch.table <- rbind(expand.grid(types1, cond1,
                                 stringsAsFactors = F),
                     expand.grid(types2, cond2,
                                 stringsAsFactors = F),
                     expand.grid(types3, cond3,
                                 stringsAsFactors = F),
                     expand.grid(types4, cond4,
                                 stringsAsFactors = F))

# create input for single enrichments ----
if (F) {
  # load data
  QTL.snp <- readRDS(file=paste0(get.path("results", local),
                                 "imputed/cis/QTL_res_all_snp_group_anno.RDS"))

  batch.table$file <- paste0(input.path, "/", batch.table$Var1, "_", batch.table$Var2, ".RDS")
  
  for (k in 1:dim(batch.table)[1]){
    #k=11
    print(paste0(k, " / ", dim(batch.table)[1]))
    type <- batch.table[k, 1]
    cond <- batch.table[k, 2]
    file <- paste0(input.path, "/", type, "_", cond, ".RDS")
    
    
    if (is.na(cond)) {cond <- NULL}
    
    hits <- get.hits(QTL.snp, "gene", type, cond)
    bg <- get.bg.match(QTL.snp, hits, "gene", cond, N = 100)
    
    res <- list(hits=hits, bg=bg[[1]], index=bg[[2]])
    saveRDS(res, file)
  }
  # states <- sort(unique(QTL.snp$state))
  # states
}

# create input for top5 enrichments ----
if (F) {
  # load data
  QTL.snp <- readRDS(file=paste0(get.path("results", local),
                                 "imputed/cis/QTL_res_all_snp_group_anno.RDS"))

  batch.table$file <- paste0(input.path, "/", batch.table$Var1, "_", batch.table$Var2, "_5.RDS")
  
  for (k in 1:dim(batch.table)[1]){
    #k=1
    print(paste0(k, " / ", dim(batch.table)[1]))
    type <- batch.table[k, 1]
    cond <- batch.table[k, 2]
    file <- paste0(input.path, "/", type, "_", cond, "_5.RDS")
    
    
    if (is.na(cond)) {cond <- NULL}
    
    hits <- rbind(get.hits(QTL.snp, "gene", type, cond),
                  get.hits(QTL.snp, "gene", type, cond, k=2),
                  get.hits(QTL.snp, "gene", type, cond, k=3),
                  get.hits(QTL.snp, "gene", type, cond, k=4),
                  get.hits(QTL.snp, "gene", type, cond, k=5))
    print(paste0("Hits for ", batch.table[k, ], " done"))
    
    bg <- get.bg.match(QTL.snp, hits, "gene", cond, N = 100)
    print(paste0("Background for ", batch.table[k, ], " done"))
    
    res <- list(hits=hits, bg=bg[[1]], index=bg[[2]])
    saveRDS(res, file)
  }
  # states <- sort(unique(QTL.snp$state))
  # print(states)
}

# Do the enrichments ----
if(T){
  tf.anno <- readRDS(file = paste0(get.path("locations", local),
                        "functional_annotation/TFBS/",
                        "SNP_anno_TFBS_remap_and_NKX2-5_openChromatin_GRanges.RDS"))[, c("snps", "TFBS2")]
  
  annos <- c("in.gene", "in.transcript",
             "UTR3", "UTR5", "exon", "splice",
             "TFBS", "TFBS2", "miRBS", "RBPBS")
  
  states <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk",
              "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv",
              "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
  
  states_simple <- c("TssActive", "Transcript", "Enhancer", "Heterochromatin",
                     "TssBivalent", "Repressive", "Quiescent")
  
  map = c(
    "1_TssA"      = "TssActive",
    "2_TssAFlnk"  = "TssActive",
    "3_TxFlnk"    = "Transcript",
    "4_Tx"        = "Transcript",
    "5_TxWk"      = "Transcript",
    "6_EnhG"      = "Enhancer",
    "7_Enh"       = "Enhancer",
    "8_ZNF/Rpts"  = "Heterochromatin",
    "9_Het"       = "Heterochromatin",
    "10_TssBiv"   = "TssBivalent",
    "11_BivFlnk"  = "TssBivalent",
    "12_EnhBiv"   = "Enhancer",
    "13_ReprPC"   = "Repressive",
    "14_ReprPCWk" = "Repressive",
    "15_Quies"    = "Quiescent"
  )
  
  states2 <- data.frame(state=states,
                        state_simple=map[states],
                        stringsAsFactors = F)
  
  VEPs <- c("missense", "NMD")
  
  annos_all <- c(annos, states, states_simple, VEPs)
  cols <- c("N", "type", "cond", "variable",
            "bg0", "bg1", "hit0", "hit1",
            "fish.pval", "fish.or", "diff")
  N <- 100
  enrichment <-  data.frame(matrix(as.numeric(NA),
                                   nrow = N*dim(batch.table)[1]*length(annos_all),
                                   ncol = length(cols)))
  colnames(enrichment) <- cols
  
  # calculate single enrichments -----
  if (T){
    batch.table$file <- paste0(input.path, "/", batch.table$Var1, "_", batch.table$Var2, ".RDS")
    
    t <- 0
    for (k in 1:dim(batch.table)[1]){
      # k=1
      res <- readRDS(file = batch.table$file[k])
      print(paste0("res = ", batch.table$file[k]))
      print(paste0("Start ", batch.table[k, 1:2]))
      
      print("Add new TFBS annotations:")
      res[[1]]$TFBS2 <- tf.anno[res[[1]]$snps, "TFBS2"]
      res[[2]]$TFBS2 <- tf.anno[res[[2]]$snps, "TFBS2"]
      
      for (m in 1:length(states)) {
        # m=1
        anno <- states[m]
        res[[1]][, anno] <- 0
        res[[1]][, anno] <- as.numeric(as.logical(res[[1]][, "state"]==anno))
        res[[2]][, anno] <- 0
        res[[2]][, anno] <- as.numeric(as.logical(res[[2]][, "state"]==anno))
      }
      
      for (m in 1:length(states_simple)) {
        # m=1
        anno <- states_simple[m]
        targets <- unique(states2[states2$state_simple==anno, "state"])
        res[[1]][, anno] <- 0
        res[[1]][grep(paste(targets, collapse = "|"), res[[1]]$state), anno] <- 1
        res[[2]][, anno] <- 0
        res[[2]][grep(paste(targets, collapse = "|"), res[[2]]$state), anno] <- 1
      }
      
      anno <- "missense"
      targets <- c("stop_gained","frameshift_variant", "stop_lost",
                   "start_lost",
                   "inframe_insertion", "inframe_deletion",
                   "missense_variant", "protein_altering_variant")
      res[[1]][, anno] <- 0
      res[[1]][grep(paste(targets, collapse = "|"), res[[1]]$VEP), anno] <- 1
      res[[2]][, anno] <- 0
      res[[2]][grep(paste(targets, collapse = "|"), res[[2]]$VEP), anno] <- 1
      
      anno <- "NMD"
      targets <- c("NMD_transcript_variant")
      res[[1]][, anno] <- 0
      res[[1]][grep(paste(targets, collapse = "|"), res[[1]]$VEP), anno] <- 1
      res[[2]][, anno] <- 0
      res[[2]][grep(paste(targets, collapse = "|"), res[[2]]$VEP), anno] <- 1
    
      for(l in 1:N){
        print(paste0("Iteration: ", l, " / ", N))
        # l=1
        for (m in 1:length(annos_all)) {
          # m=1
          t <- t+1
          anno <- annos_all[m]
          enrichment[t, c("N", "type", "cond", "variable")] <-
            c(l, batch.table[k, "Var1"], batch.table[k, "Var2"], anno)
    
          er <- enrich(res[[1]], res[[2]], res[[3]], l, anno)
          
          enrichment[t, c("bg0", "bg1", "hit0", "hit1")] <- 
            as.numeric(c(er[["bg0"]], er[["bg1"]], er[["hit0"]], er[["hit1"]]))
          
          enrichment[t, "diff"] <- 
            sign(enrichment[t, "hit1"]-enrichment[t, "bg1"])
          
          enrichment[t, c("fish.pval", "fish.or")] <- 
            c(er[["fish.pval"]], er[["fish.or"]])
        }
      }
      saveRDS(enrichment, file=paste0(result.path, "/enrichment_results_temp.RDS"))
      print(paste0("Saved results: ", k, " / ", dim(batch.table)[1]))
    }
    saveRDS(enrichment, file=paste0(result.path, "/enrichment_results.RDS"))
  }
  
  annos_all <- c(annos, states, states_simple, VEPs)
  cols <- c("N", "type", "cond", "variable",
            "bg0", "bg1", "hit0", "hit1",
            "fish.pval", "fish.or", "diff")
  N <- 100
  enrichment5 <-  data.frame(matrix(as.numeric(NA),
                                   nrow = N*dim(batch.table)[1]*length(annos_all),
                                   ncol = length(cols)))
  colnames(enrichment5) <- cols
  
  # calculate top5 enrichments -----
  if (T){
    batch.table$file <- paste0(input.path, "/", batch.table$Var1, "_", batch.table$Var2, "_5.RDS")
    
    t <- 0
    for (k in 1:dim(batch.table)[1]){#2){#
      # k=1
      res <- readRDS(file = batch.table$file[k])
      print(paste0("res = ", batch.table$file[k]))
      print(paste0("Start ", batch.table[k, 1:2]))
      
      print("Add new TFBS annotations:")
      res[[1]]$TFBS2 <- tf.anno[res[[1]]$snps, "TFBS2"]
      res[[2]]$TFBS2 <- tf.anno[res[[2]]$snps, "TFBS2"]
      
      for (m in 1:length(states)) {
        # m=1
        anno <- states[m]
        res[[1]][, anno] <- 0
        res[[1]][, anno] <- as.numeric(as.logical(res[[1]][, "state"]==anno))
        res[[2]][, anno] <- 0
        res[[2]][, anno] <- as.numeric(as.logical(res[[2]][, "state"]==anno))
      }
      
      for (m in 1:length(states_simple)) {
        # m=1
        anno <- states_simple[m]
        targets <- unique(states2[states2$state_simple==anno, "state"])
        res[[1]][, anno] <- 0
        res[[1]][grep(paste(targets, collapse = "|"), res[[1]]$state), anno] <- 1
        res[[2]][, anno] <- 0
        res[[2]][grep(paste(targets, collapse = "|"), res[[2]]$state), anno] <- 1
      }
      
      anno <- "missense"
      targets <- c("stop_gained","frameshift_variant", "stop_lost",
                   "start_lost",
                   "inframe_insertion", "inframe_deletion",
                   "missense_variant", "protein_altering_variant")
      res[[1]][, anno] <- 0
      res[[1]][grep(paste(targets, collapse = "|"), res[[1]]$VEP), anno] <- 1
      res[[2]][, anno] <- 0
      res[[2]][grep(paste(targets, collapse = "|"), res[[2]]$VEP), anno] <- 1
      
      anno <- "NMD"
      targets <- c("NMD_transcript_variant")
      res[[1]][, anno] <- 0
      res[[1]][grep(paste(targets, collapse = "|"), res[[1]]$VEP), anno] <- 1
      res[[2]][, anno] <- 0
      res[[2]][grep(paste(targets, collapse = "|"), res[[2]]$VEP), anno] <- 1
      
      for(l in 1:N){
        print(paste0("Iteration: ", l, " / ", N))
        # l=1
        for (m in 1:length(annos_all)) {
          # m=1
          t <- t+1
          anno <- annos_all[m]
          enrichment5[t, c("N", "type", "cond", "variable")] <-
            c(l, batch.table[k, "Var1"], batch.table[k, "Var2"], anno)
          
          er <- enrich(res[[1]], res[[2]], res[[3]], l, anno)
          
          enrichment5[t, c("bg0", "bg1", "hit0", "hit1")] <- 
            as.numeric(c(er[["bg0"]], er[["bg1"]], er[["hit0"]], er[["hit1"]]))
          
          enrichment5[t, "diff"] <- 
            sign(enrichment5[t, "hit1"]-enrichment5[t, "bg1"])
          
          enrichment5[t, c("fish.pval", "fish.or")] <- 
            c(er[["fish.pval"]], er[["fish.or"]])
        
        }
      }
      saveRDS(enrichment5, file=paste0(result.path, "/enrichment_results_5_temp.RDS"))
      print(paste0("Saved results: ", k, " / ", dim(batch.table)[1]))
    }
    saveRDS(enrichment5, file=paste0(result.path, "/enrichment_results_5.RDS"))
  }
  
  
  # create plots -----
  input.path <- paste0(get.path("results", local),
                       "imputed/cis/enrichment/input")
  dir.create(input.path, showWarnings = F, recursive = T)
  result.path <- paste0(get.path("results", local),
                        "imputed/cis/enrichment/results")
  dir.create(result.path, showWarnings = F, recursive = T)
  pics <- paste0(result.path, "/plots/")
  dir.create(pics, recursive = T)
  
  # single enrichment plots ----
  if (F) {
    enrichment <- readRDS(file=paste0(result.path, "/enrichment_results.RDS"))
    library(ggplot2)
    library(RColorBrewer)
    library(lemon)
    HMGU.blue <- "#003E6E"
    mygray <- "#C6DDEA"
    col.paired <- brewer.pal(n = 11, "Paired")
    col.set <- col.paired[c(2,8,9)]
    col.set2 <- col.paired[c(9,5,7,3,8,2)]
    col.set3 <- brewer.pal(n = 11, "PRGn")[c(4,6,8)]
    col.set3 <- c(col.set[3], "#F7F7F7", col.set[1])
    col.set4 <- c(col.set[1], col.set[3])
    
    df1 <- enrichment[(enrichment$type=="FDR.eqtl" | enrichment$type=="FDR.res_eqtl"), ]
    gg_eQTL <- ggplot(data=df1, aes(x=variable, fill=as.factor(diff))) +
      geom_bar(stat="count")+
      geom_hline(yintercept = 95, col=col.set[2], size = 0.5) +
      geom_hline(yintercept = 5, col=col.set[2], size = 0.5) +
      scale_y_continuous(breaks = c(5, 95),
                         labels = c("depleted", "enriched")) +
      scale_fill_manual("Enrichment",
                        values=col.set3,
                        labels=c("depleted in QTLs", "equal counts", "enriched in QTLs"))+
      facet_rep_grid(type~cond, scales = "free_x", space = "free_x",
                     repeat.tick.labels = T) +
      labs(list(x="Annotation", y="Enrichment",
                title="Enrichment for functional element annotations between single QTL / non QTL variants - imputed eQTL data")) +
      theme_minimal() + 
      theme(plot.title = element_text(size=12,
                                      face="bold",
                                      colour=HMGU.blue,
                                      hjust=0.5),
            axis.title=element_text(size=18,
                                    face="bold",
                                    colour=HMGU.blue),
            axis.text.x=element_text(size=10,
                                     face="bold",
                                     colour=HMGU.blue,
                                     hjust = 1,
                                     angle = 90,
                                     vjust = 0.5),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(size = 1, linetype = 'solid',
                                        colour = HMGU.blue,
                                 fill = NA))
    
    gg_eQTL
    
    
    pdf(paste0(pics, "enrichment_results_eQTL.pdf"), width = 11, height = 8)
    print(gg_eQTL)
    dev.off()
    
    df2 <- enrichment[(enrichment$type=="FDR.pqtl" | enrichment$type=="FDR.res_pqtl"), ]
    gg_pQTL <- ggplot(data=df2, aes(x=variable, fill=as.factor(diff))) +
      geom_bar(stat="count")+
      geom_hline(yintercept = 95, col=col.set[2], size = 0.5) +
      geom_hline(yintercept = 5, col=col.set[2], size = 0.5) +
      scale_y_continuous(breaks = c(5, 95),
                         labels = c("depleted", "enriched")) +
      scale_fill_manual("Enrichment",
                        values=col.set3,
                        labels=c("depleted in QTLs", "equal counts", "enriched in QTLs"))+
      facet_rep_grid(type~cond, scales = "free_x", space = "free_x",
                     repeat.tick.labels = T) +
      labs(list(x="Annotation", y="Enrichment",
                title="Enrichment for functional element annotations between single QTL / non QTL variants - imputed pQTL data")) +
      theme_minimal() + 
      theme(plot.title = element_text(size=12,
                                      face="bold",
                                      colour=HMGU.blue,
                                      hjust=0.5),
            axis.title=element_text(size=18,
                                    face="bold",
                                    colour=HMGU.blue),
            axis.text.x=element_text(size=10,
                                     face="bold",
                                     colour=HMGU.blue,
                                     hjust = 1,
                                     angle = 90,
                                     vjust = 0.5),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(size = 1, linetype = 'solid',
                                        colour = HMGU.blue,
                                        fill = NA))
    pdf(paste0(pics, "enrichment_results_pQTL.pdf"), width = 11, height = 8)
    print(gg_pQTL)
    dev.off()
  }
  
  # top5 enrichment plots ----
  if (F) {
    enrichment5 <- readRDS(file=paste0(result.path, "/enrichment_results_5_temp.RDS"))
    library(ggplot2)
    
    df3 <- enrichment5[(enrichment5$type=="FDR.eqtl" | enrichment5$type=="FDR.res_eqtl"), ]
    gg_eQTL5 <- ggplot(data=df3, aes(x=variable, fill=as.factor(diff))) +
      geom_bar(stat="count")+
      geom_hline(yintercept = 95, col=col.set[2], size = 0.5) +
      geom_hline(yintercept = 5, col=col.set[2], size = 0.5) +
      scale_y_continuous(breaks = c(5, 95),
                         labels = c("depleted", "enriched")) +
      scale_fill_manual("Enrichment",
                        values=col.set3,
                        labels=c("depleted in QTLs", "equal counts", "enriched in QTLs"))+
      facet_rep_grid(type~cond, scales = "free_x", space = "free_x",
                     repeat.tick.labels = T) +
      labs(list(x="Annotation", y="Enrichment",
                title="Enrichment for functional element annotations between top5 QTL / non QTL variants - imputed eQTL data")) +
      theme_minimal() + 
      theme(plot.title = element_text(size=12,
                                      face="bold",
                                      colour=HMGU.blue,
                                      hjust=0.5),
            axis.title=element_text(size=18,
                                    face="bold",
                                    colour=HMGU.blue),
            axis.text.x=element_text(size=10,
                                     face="bold",
                                     colour=HMGU.blue,
                                     hjust = 1,
                                     angle = 90,
                                     vjust = 0.5),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(size = 1, linetype = 'solid',
                                        colour = HMGU.blue,
                                        fill = NA))
    pdf(paste0(pics, "enrichment_results_top5_eQTL.pdf"), width = 11, height = 8)
    print(gg_eQTL5)
    dev.off()
    
    df4 <- enrichment5[(enrichment5$type=="FDR.pqtl" | enrichment5$type=="FDR.res_pqtl"), ]
    gg_pQTL5 <- ggplot(data=df4, aes(x=variable, fill=as.factor(diff))) +
      geom_bar(stat="count")+
      geom_hline(yintercept = 95, col=col.set[2], size = 0.5) +
      geom_hline(yintercept = 5, col=col.set[2], size = 0.5) +
      scale_y_continuous(breaks = c(5, 95),
                         labels = c("depleted", "enriched")) +
      scale_fill_manual("Enrichment",
                        values=col.set3,
                        labels=c("depleted in QTLs", "equal counts", "enriched in QTLs"))+
      facet_rep_grid(type~cond, scales = "free_x", space = "free_x",
                     repeat.tick.labels = T) +
      labs(list(x="Annotation", y="Enrichment",
                title="Enrichment for functional element annotations between top5 QTL / non QTL variants - imputed pQTL data")) +
      theme_minimal() + 
      theme(plot.title = element_text(size=12,
                                      face="bold",
                                      colour=HMGU.blue,
                                      hjust=0.5),
            axis.title=element_text(size=18,
                                    face="bold",
                                    colour=HMGU.blue),
            axis.text.x=element_text(size=10,
                                     face="bold",
                                     colour=HMGU.blue,
                                     hjust = 1,
                                     angle = 90,
                                     vjust = 0.5),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(size = 1, linetype = 'solid',
                                        colour = HMGU.blue,
                                        fill = NA))
    pdf(paste0(pics, "enrichment_results_top5_pQTL.pdf"), width = 11, height = 8)
    print(gg_pQTL5)
    dev.off()
  }
}

# heatmap plots ----------------------------------------------------------------
if(F){
  result.path <- paste0(get.path("results", local),
                        "imputed/cis/enrichment/results")
  pics <- paste0(result.path, "/plots/")
  
  ## single enrichments --------------------------------------------------------
  enrichment.one <- readRDS(file=paste0(result.path, "/enrichment_results.RDS"))
  enrichment.one <- enrichment.one[enrichment.one$variable!="in.transcript", ]
  levels <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk",
              "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv",
              "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies",
              "in.gene", "UTR3", "UTR5", "exon", "splice",
              "TFBS", "miRBS", "RBPBS", "missense", "NMD")
  enrichment.one$variable <- factor(enrichment.one$variable, levels = levels)
  groups <- c("eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL",
              "sharedQTL", "only_eQTL", "only_pQTL")
  # groups <- c("eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL",
  #             "sharedQTL", "independentQTL", "only_eQTL", "only_pQTL")
  df.one <- matrix(NA,
                   nrow = nlevels(enrichment.one$variable),
                   ncol = length(groups),
                   dimnames = list(levels(enrichment.one$variable), groups))
  
  for (group in groups){
    if(group=="eQTL"){
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.eqtl" & is.na(enrichment.one$cond)), ]
    } else if (group=="pQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.pqtl" & is.na(enrichment.one$cond)), ]
    } else if (group=="res_eQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.res_eqtl" & is.na(enrichment.one$cond)), ]
    } else if (group=="res_pQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.res_pqtl" & is.na(enrichment.one$cond)), ]
    } else if (group=="ratioQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.ratios" & is.na(enrichment.one$cond)), ]
    } else if (group=="sharedQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.pqtl" & enrichment.one$cond=="Group2"), ]
    } else if (group=="independentQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.res_pqtl" & enrichment.one$cond=="xxx"), ]
    } else if (group=="only_eQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.res_eqtl" & enrichment.one$cond=="Group1"), ]
    } else if (group=="only_pQTL") {
      df1 <- enrichment.one[which(enrichment.one$type=="FDR.res_pqtl" & enrichment.one$cond=="Group5"), ]
    }
    temp <- as.data.frame.matrix(table(df1$variable, df1$diff))
    temp$counts <- apply(temp[, c("-1","1")], 1, max)
    temp$counts.sgn <- 1
    temp[temp$`-1`>temp$`1`, "counts.sgn"] <- -temp[temp$`-1`>temp$`1`, "counts.sgn"]
    df.one[, group] <- temp$counts*temp$counts.sgn
  }
  
  ## top5 enrichments ----------------------------------------------------------
  enrichment.five <- readRDS(file=paste0(result.path, "/enrichment_results_5.RDS"))
  enrichment.five <- enrichment.five[enrichment.five$variable!="in.transcript", ]
  levels <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk",
              "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv",
              "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies",
              "in.gene", "UTR3", "UTR5", "exon", "splice",
              "TFBS", "miRBS", "RBPBS", "missense", "NMD")
  enrichment.five$variable <- factor(enrichment.five$variable, levels = levels)
  groups <- c("eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL",
              "sharedQTL", "only_eQTL", "only_pQTL")
  # groups <- c("eQTL", "pQTL", "reseQTL", "respQTL", "ratioQTL",
  #             "sharedQTL", "independentQTL", "only_eQTL", "only_pQTL")
  df.five <- matrix(NA, nrow = nlevels(enrichment.five$variable), ncol = length(groups),
                    dimnames = list(levels(enrichment.five$variable), groups))
  
  for (group in groups){
    if(group=="eQTL"){
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.eqtl" & is.na(enrichment.five$cond)), ]
    } else if (group=="pQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.pqtl" & is.na(enrichment.five$cond)), ]
    } else if (group=="res_eQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.res_eqtl" & is.na(enrichment.five$cond)), ]
    } else if (group=="res_pQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.res_pqtl" & is.na(enrichment.five$cond)), ]
    } else if (group=="ratioQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.ratios" & is.na(enrichment.five$cond)), ]
    } else if (group=="sharedQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.pqtl" & enrichment.five$cond=="Group2"), ]
    } else if (group=="independentQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.res_pqtl" & enrichment.five$cond=="xxx"), ]
    } else if (group=="only_eQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.res_eqtl" & enrichment.five$cond=="Group1"), ]
    } else if (group=="only_pQTL") {
      df1 <- enrichment.five[which(enrichment.five$type=="FDR.res_pqtl" & enrichment.five$cond=="Group5"), ]
    }
    temp <- as.data.frame.matrix(table(df1$variable, df1$diff))
    temp$counts <- apply(temp[, c("-1","1")], 1, max)
    temp$counts.sgn <- 1
    temp[temp$`-1`>temp$`1`, "counts.sgn"] <- -temp[temp$`-1`>temp$`1`, "counts.sgn"]
    df.five[, group] <- temp$counts*temp$counts.sgn
  }
  
  colnames(df.five) <- paste0(colnames(df.five), "_top5")
  df <- cbind(df.one, df.five)
  
  saveRDS(df, paste0(result.path, "/enrichment_heatmap_data_new.RDS"))
  library(pheatmap)
  pheatmap(t(df),
           cluster_rows = F,
           cluster_cols = T,
           #filename = paste0(pics, "heatmap_enrichments.pdf"),
           main="enrichment of functional elements (QTL vs. nonQTL)",
           legend_breaks = c(-100, 0, 100),
           legend_labels = c("- emp.pval < 0.05 \n  depletion",
                             "- emp.pval = 1",
                             "  enrichment \n- emp.pval < 0.05"),
           # legend_breaks = c(-95, 0, 95),
           legend = T,
           gaps_row = length(groups),
           cellwidth = 20,
           cellheight = 20)
  pheatmap(t(df),
           cluster_rows = F,
           cluster_cols = T,
           filename = paste0(pics, "heatmap_enrichments.pdf"),
           main="enrichment of functional elements (QTL vs. nonQTL)",
           legend_breaks = c(-100, 0, 100),
           legend_labels = c("- emp.pval < 0.05 \n  depletion",
                             "- emp.pval = 1",
                             "  enrichment \n- emp.pval < 0.05"),
           # legend_breaks = c(-95, 0, 95),
           legend = T,
           gaps_row = length(groups),
           cellwidth = 20,
           cellheight = 20)
  dev.off()
  pheatmap(t(df),
           cluster_rows = F,
           cluster_cols = T,
           filename = paste0(pics, "heatmap_enrichments.png"),
           main="enrichment of functional elements (QTL vs. nonQTL)",
           legend_breaks = c(-100, 0, 100),
           legend_labels = c("- emp.pval < 0.05 \n  depletion",
                             "- emp.pval = 1",
                             "  enrichment \n- emp.pval < 0.05"),
           # legend_breaks = c(-95, 0, 95),
           legend = T,
           gaps_row = length(groups),
           cellwidth = 20,
           cellheight = 20)
  dev.off()
  
}



q(save="no")



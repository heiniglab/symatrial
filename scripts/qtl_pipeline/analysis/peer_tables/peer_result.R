##
#   loads MatrixEQTL results and makes a summary table
#
#

peer.write.table <- function(base_dir, outdir = NULL, trans = 0) {
  #base_dir <- dir
  
  source("analysis/peer_tables/me_analyze.R")
  source("helper/helper.R")
  
  files <- list.files(path = base_dir, pattern = "*.RDS", full.names = T, recursive = T)

  data <- sapply(files, function(x) {
      d <- readRDS(x)
      return(analyze(d, trans))
    })
    
    data <- data.frame(t(data))
    colnames(data) <-
      c(
        "FDR<0.05",
        "FDR<0.05_Genes",
        "FDR<0.2",
        "FDR<0.2_Genes",
        "pval<1e-5",
        "pval<1e-5_Genes"
      )
    rownames(data) <-
      sapply(rownames(data), function(x) {
        # paste0(get.filename(x), "_", substr(x, 120, 125))
        get.filename(x)
      })

    outfile <- "peer_table.RDS"
    
    if(is.null(outdir)){
      return(data)
    }
    
    dir.create(outdir, showWarnings = F, recursive = T)
    saveRDS(data, paste(outdir, outfile, sep="/"))
  }

plot.peer.table <- function(QTL.res, model, cutoff, types=NULL) {
  # model <- "linear"
  # type <- "pqtl_res"
  # cutoff <- "FDR.0.05_Genes"
  library(ggplot2)
  library(RColorBrewer)
  col.paired <- brewer.pal(n = 11, "Paired")
  col.set <- col.paired[c(2,8,9)]
  
  if(is.null(types)){
    types <- unique(QTL.res$type)
  }
  df <- QTL.res[QTL.res$model==model & QTL.res$type %in% types, ]
  df$val <- df[, cutoff]
  df$type <- factor(df$type, levels = types)
  g <- ggplot(df, aes(x=Nk, y=val, col=cov_sd)) +
    geom_point() +
    geom_line() +
    scale_color_manual("covariates:",
                       values = col.set,
                       labels = c("PEER factors only",
                                  "PEER factors and fibrosis score",
                                  "PEER factors, age, sex, BMI, disease status and fibrosis score")) +
    facet_grid(norm_sd~type,scales = "free_x", space = "free_x") +
    labs(title=cutoff,
              x="number of PEER factors",
              y="count",
              col="covariates") +
    theme_bw() +
    theme(legend.position = "bottom")#+
        #theme(axis.text.x = element_text(angle = 45, hjust = 1))
        return(g)
}


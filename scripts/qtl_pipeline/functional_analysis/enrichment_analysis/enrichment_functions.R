# ------------------------------------------------------------------------------
#' Functional enrichment of regulatory/posttranscriptional elements
#' for cis eQTL/pQTL hits (mQTL hits) based on annotated snp-gene pairs
#'
#' @author Ines Assum
#'
# ------------------------------------------------------------------------------



get.hits <- function(data, hit, value, condition=NULL, cutoff = 0.05, k=1) {
  # data <- QTL.snp
  # hit <- "gene"
  # value <- "FDR.pqtl"
  # cutoff = 0.05
  # k=2
  # condition="Group5"
  hits <- data[which(data[, value]<cutoff), ]
  if(!is.null(condition)){
    hits <- hits[which(hits[, condition]==1), ]
  }
  hits <- hits[order(hits[, hit], hits[, value]), ]
  table <- table(hits[, hit])
  hits <- hits[hits[, hit] %in% names(table[table>=k]), ]
  hits <- hits[which(!duplicated(hits[, hit]))+(k-1), ]
  return(hits)
}

get.bg.match <- function(data, hits, value, condition=NULL, cutoff = 0.05, MAF.cutoff = 0.05, dist.cutoff = 1000, N=100) {
  # data <- QTL.snp
  # hits <- hits.eqtl
  # value <- "FDR.eqtl"
  # # condition="Group5"
  # condition=NULL
  # cutoff = 0.05
  # MAF.cutoff = 0.05
  # dist.cutoff = 1000
  # N=100
  bg <- data[which(data[, value] > cutoff), ]
  if(!is.null(condition)){
    bg <- bg[which(bg[, condition]==0), ]
  }
  rownames(bg) <- NULL
  index <- data.frame(matrix(NA, nrow=dim(hits)[1], ncol=N))
  for (i in 1:length(hits$snps)){
    cand <- which((abs(bg$MAF-hits$MAF[i]) < MAF.cutoff) &
                    (abs(bg$dist-hits$dist[i]) < dist.cutoff))
    if (length(cand)<N){
      index[i, 1:N] <- sample(cand, N, replace=T)
      print(paste0("too little samples, replace=T, ", length(cand), "/", N))
    }else{
      index[i, 1:N] <- sample(cand, N, replace=F)
    }
  }
  return(list(bg, index))
}


enrich <- function(hits, bg, index, col, element){
  # hits <- hits.pqtl
  # bg <- bg.pqtl[[1]]
  # index <- bg.pqtl[[2]]
  # col <- 1
  # element <- "RBPBS"
  
  # construct table:
  #        |     annotation
  #        |      0       1
  # --------------------------
  # QTL  0 |  bg==0   bg==1
  #      1 | hit==0  hit==1
  
  bg0 <- sum(bg[index[, col], element]==0)
  bg1 <- sum(bg[index[, col], element]==1)
  hit0 <- sum(hits[, element]==0)
  hit1 <- sum(hits[, element]==1)
  
  tab <- matrix(c(bg0, bg1, hit0, hit1),
                nrow = 2, ncol = 2,
                byrow = T)
  # tab
  fish <- tryCatch(
    {
      fisher.test(tab)
    },
    error=function(cond) {
      # Choose a return value in case of error
      return(list(bg0=bg0, bg1=bg1, hit0=hit0, hit1=hit1, tab=tab, fish.pval=NA, fish.or=NA))
    })
  # fish
  return(list(bg0=bg0, bg1=bg1, hit0=hit0, hit1=hit1, tab=tab, fish.pval=fish$p.value, fish.or=fish$estimate))
}

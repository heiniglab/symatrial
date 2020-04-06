library(peer)

#' Calculate peer factors for given data and covariates
#' 
#' @param data nxg matrix (n=samples, g=genes/variables)
#' @param covariates nxc matrix (n=samples, c=covariates)
#' @param Nk Number of factors to estimate. Default: N/4
#' 
#' @return Matrix of peer factors calculated
#' 
#' @author Johann Hawe, slightly modified by Fabian Denbsky
#' 
#' @date 20170328
#' 
get.peer.factors <- function(data, 
                             covariates=NULL, 
                             get.residuals=F,
                             Nk=ceiling(nrow(data)*0.25)) {
     
  # create model
  model <- PEER();
  # input has to be a matrix!
  PEER_setPhenoMean(model, as.matrix(data));
  
  # add the mean estimation as default since it is recommended in the tutorial
  # of peer. will return Nk+1 factors
  PEER_setAdd_mean(model, TRUE)
  
  # set number of hidden factors to identify. If unknown, a good measure is N/4 (see howto)
  PEER_setNk(model, Nk);
  
  # should not be neccessary but increase anyways
  PEER_setNmax_iterations(model, 5000);
  
  if(!is.null(covariates)){
    # set matrix of known and important covariates, since we want to acknowledge their effect
    PEER_setCovariates(model, as.matrix(covariates));
  }
  
  # learn 
  PEER_update(model);

  # directly return the residuals if wanted
  if(get.residuals) {
    re <- PEER_getResiduals(model)
    colnames(re) <- colnames(data)
    rownames(re) <- rownames(data)
    return(re)
  }
  
  # get identified factors, contains design.matrix in the first few columns!
  factors <- PEER_getX(model);

  # if(!is.null(covariates)){
  #   # return only the calculated factors, ignoring the original design components  
  #   factors <- factors[,-c(1:ncol(covariates))]
  # }
  
  if(!is.null(covariates)){
    # return not the first constant "factor"
    factors <- factors[,-c(ncol(covariates) + 1)]
  }
  else {
    factors <- as.data.frame(factors[, -c(1)])
  }
  
  colnames(factors) <- paste0("f", seq(1:ncol(factors)))
  return(factors);
}

#'
#' @param data The data for which to get the residuals
#' @param factors Matrix of identified peer factors (with get.peer.factors)
#'
#' @return The residual matrix after using lm with the data and factors
#' 
#' @author Johann Hawe
#' 
#' @date 20170328
#' 
get.peer.residuals <- function(data, factors) {           
  if(!is.matrix(factors)){
          stop("Peer factors need to be in matrix format.")
  }
      
  n <- ncol(factors)
  cat("Applying linear model for ", n, "factor(s)...");
  
  result <- apply(data, 1, function(y) { 
          dat <- cbind.data.frame(y=y,factors)
          #colnames(dat) <- c("y", colnames(factors));
          lm(data=dat, 
             paste("y", "~", sep=" ", paste(colnames(factors), 
                                            collapse=" + ")))
  });
  
  residual.mat <- matrix(data=unlist(lapply(result, resid)), 
                         nrow=nrow(data), 
                         byrow=T)
  colnames(residual.mat) <- colnames(data)                
  rownames(residual.mat) <- rownames(data)
  
  return(residual.mat);
}


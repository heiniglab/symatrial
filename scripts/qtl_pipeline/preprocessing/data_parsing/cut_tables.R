#' cuts down given expression table to only contain rows where locations are known
#' 
#' @param expression expression table
#' @param locations location table
#' 
#' 
#' @author Fabian Denbsky
#' 
#' @date 20171214
#' 
cut.locations <- function(expression, locations) {
  return(expression[rownames(expression) %in% locations$name,])
}



#' gives datasets with common samples
#' 
#' @param expression expression table
#' @param genotype genotype table
#' 
#' 
#' @author Fabian Denbsky
#' 
#' @date 20171008
#' 
common.samples <- function(expression, genotype) {
  which_gen <- which(colnames(genotype) %in% colnames(expression))
  which_exp <- which(colnames(expression) %in% colnames(genotype))
  
  new_gen <- genotype[, which_gen, with=F]
  rownames(new_gen) <- rownames(genotype)
  new_exp <- expression[, which_exp]
  
  new_gen <- order.genotype(new_exp, new_gen)
  rownames(new_gen) <- rownames(genotype)
  
  return(list(new_gen, new_exp))
}



#' orderes samples of genotype table to to match the order of the samples in the expression table
#' 
#' @param expression expression table
#' @param genotype genotype table
#' 
#' 
#' @author Fabian Denbsky
#' 
#' @date 20171214
#' 
order.genotype <- function(expression, genotype){
  genotype <- as.data.frame(genotype)
  return(genotype[, order(match(colnames(genotype), colnames(expression)))])
}


#' orderes samples of covariate table to match the order of the samples in the expression table
#' 
#' @param expression expression table
#' @param covariates covariate table
#' 
#' 
#' @author Fabian Denbsky
#' 
#' @date 20171214
#' 
order.covariates <- function(expression, covariates){
  covariates <- as.data.frame(covariates)
  return(covariates[, order(match(colnames(covariates), colnames(expression)))])
}



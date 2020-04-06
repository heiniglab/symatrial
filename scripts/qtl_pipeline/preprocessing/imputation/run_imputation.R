## Run Imputation
#
#
#
#
#


impute.data <- function(expression){
  library(impute)
  exp <- data.matrix(expression)
  exp_imp <- t(impute.knn(exp))
  return(as.data.frame(exp_imp[[1]]))
}


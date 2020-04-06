# transform the scaled counts to std normal per gene

stdnorm <- function(x) {
  r = rank(x, ties.method = "random")
  qnorm(r / (length(x) + 1))
}

correct.linewise <- function(expression) {
  transformed = apply(expression, 1, stdnorm)
  return(t(transformed))
}
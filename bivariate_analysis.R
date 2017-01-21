

get.correlated <- function(df, target.name){
  
  targetIndex <- which(colnames(df) == target.name)
  indp.indx <- which(colnames(df) != target.name)
  
  pvalue <- unname(sapply(indp.indx, function(i) chisq.test(x = df[,i], y = df[,targetIndex], correct = F)$p.value))

  named.pvalue <- setNames(pvalue, colnames(df[-targetIndex]))
  
  correlated <- sort(named.pvalue[named.pvalue <= .05])
  
  return (correlated)

}

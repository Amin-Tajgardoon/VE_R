
## prepares FLU dataset

data.preprocess <-  function(data){
  
  ## list all unique values in the dataframe
  #as.character(unique(unlist(flu)))
  
  ##  replace Nulls with 'N'
  data[is.na(data)] <- 'N'
  
  ##  replace 'ND' & 'MC' with 'N'
  data[data %in% c('ND', 'MC', 'MISSING CONFIRMED')] <- 'N'
  
  ##  remove columns with no variance
  data <-  data[, sapply(data, function(x) length(unique(x))) > 1]
  
  ## convert columns to factor 
  data <- data.frame(lapply(data, as.factor), stringsAsFactors=FALSE)
  
  return (data)

}
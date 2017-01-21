library(caret)
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/resample.R")


cross_validation <- function(data, target, pos_class="MATCH", balance = FALSE, blacklist, kfold=10, 
                             algo, score="k2", maxp=5, estimator="bayes",
                             typeBalance="ubSMOTE", percOver=200, percUnder=200,
                             kBalance=5, percBalance=50, methodBalance="percPos", wBalance=NULL) {
  
  flds <-  createFolds(y = data[ ,target], k = kfold, list = TRUE, returnTrain = TRUE)
  
  y_test <- vector(length = nrow(data))
  match_probs <- vector(length = nrow(data))
  
  start.time <- Sys.time()

  i <-  1;
  
  for (fld in flds){
    
    print(paste("fold: ",i))
    i <- i + 1
    
    train_set <- data[fld,]
    
    if (balance == TRUE){
      targetIndex <- which(colnames(data) == target)
      train_set <- data.resample(X= train_set[, -targetIndex], Y= train_set[, targetIndex], positive= pos_class, 
                                 type=typeBalance, percOver=percOver, percUnder=percUnder, verbose=TRUE,
                                 k=kBalance, perc=percBalance, method=methodBalance, w=wBalance)
      train_set <- set.target.name(train_set, targetIndex, target)
      
    }
    
    
    test_set <- data[-fld,]
    
    if(algo == "hc"){
      bn <-  hc(x = train_set, blacklist = blacklist, score = score, maxp = maxp, optimized = TRUE)
    }
    
    else if(algo == "tabu"){
      bn <-  tabu(x = train_set, blacklist = blacklist, score = score, maxp = maxp, optimized = TRUE,
                  tabu = 10, max.tabu = 10)
    }
    else {
      stop(paste0("unknown algorithm: ", algo))
    }
    
    
    fit <-  bn.fit(x = bn, data = train_set, method = estimator)
    
    pred <- predict(object = fit, data = test_set, node = target, method = "parents", prob = TRUE, debug = FALSE)
    pred_attr <- attributes(pred)
    
    match_probs[-fld] <- pred_attr$prob[pos_class,]
    #y_test[-fld] <- test_set[ ,target]
    
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0(algo, "_cross-validation time: "))
  print(time.taken)
  
  
  result <- list("y_test" = data[ ,target], "probs" = match_probs)
  
  return (result)
}



set.target.name <-  function (data, targetIndex, targetName){
  colnames(data)[targetIndex] <- target
  return (data)
}

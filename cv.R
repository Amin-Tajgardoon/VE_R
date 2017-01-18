library(caret)


cross_validation <- function(data, target, pos_class="MATCH", blacklist, k=10, score="k2", maxp=5, estimator="bayes" ) {
  
  flds <-  createFolds(y = data[ ,target], k = k, list = TRUE, returnTrain = TRUE)
  
  y_test <- vector(length = nrow(data))
  match_probs <- vector(length = nrow(data))
  
  start.time <- Sys.time()

  i <-  1;
  
  for (fld in flds){
    
    print(paste("fold: ",i))
    i <- i + 1
    
    train_set <- data[fld,]
    test_set <- data[-fld,]
    
    hc = hc(x = train_set, blacklist = blacklist, score = score, maxp = maxp, optimized = TRUE)
    fit <-  bn.fit(x = hc, data = train_set, method = estimator)
    
    pred <- predict(object = fit, data = test_set, node = target, method = "parents", prob = TRUE, debug = FALSE)
    pred_attr <- attributes(pred)
    
    match_probs[-fld] <- pred_attr$prob[pos_class,]
    #y_test[-fld] <- test_set[ ,target]
    
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("croos validation time: ", time.taken, "\n", sep = "")
  
  result <- list("probs" = match_probs, "y_test" = data[ ,target])
  
  return (result)
}

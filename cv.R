## performs cross-validation
## returns evaluation results

library(caret)
library(bnlearn)
library(smbinning)

source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/hoslem.test.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/calibration_curve.R")


cross_validation <- function(df, target, pos_class, balance = FALSE, blacklist = NULL, kfold= 10, seed = NULL,
                             algo, score= NULL, maxp= NULL, estimator= NULL, restart = NULL, optimized = optimized,
                             test = NULL, alpha = NULL, undirected = FALSE, restrict = NULL, maximize = NULL, max.iter = Inf,
                             typeBalance="ubSMOTE", percOver=200, percUnder=200, 
                             kBalance=5, percBalance=50, methodBalance="percPos", wBalance=NULL) {
  
  set.seed(seed)
  
  print(paste("Algorithm","score","test","estimator"))
  print(str(c(algo, score, test, estimator)))
  
  print(paste("seed number:",seed))
  
  flds <-  createFolds(y = df[ ,target], k = kfold, list = TRUE, returnTrain = TRUE)
  
  y_test <- rep("", nrow(df)) #vector(length = nrow(df))
  match_probs <-rep(-1, nrow(df)) # vector(length = nrow(df))
  row.index = c()
  
  local_score = list()
  
  if(algo == "glm"){
    
    pos.indx <- df$ANTIGEN_MATCHED == pos_class
    neg.indx <- df$ANTIGEN_MATCHED != pos_class
    levels(df$ANTIGEN_MATCHED)[match("MATCH",levels(df$ANTIGEN_MATCHED))] = "1"
    levels(df$ANTIGEN_MATCHED)[match("NMATCH",levels(df$ANTIGEN_MATCHED))] = "0"
    df$ANTIGEN_MATCHED[pos.indx] = 1
    df$ANTIGEN_MATCHED[neg.indx] = 0
    
    print(table(df[,target]))
  }
  
  start.time <- Sys.time()
  
  for (i in 1:length(flds)){
    
    possibleError <- NULL
    
    print(paste("fold: ",i))
    
    fld <- flds[[i]]
    train_set <- df[fld,]
    test_set <- df[-fld,]
    
    # Discretize AGE_wdemog
    smbin = smbinning(df = train_set , x = "AGE_wdemog", y = "ANTIGEN_MATCHED", p = 0.05)
    cuts = c(min(df$AGE_wdemog), smbin$cuts, max(df$AGE_wdemog))
    
    train_set$AGE_wdemog = cut(x = train_set$AGE_wdemog, breaks = cuts, include.lowest = TRUE)
    test_set$AGE_wdemog = cut(x = test_set$AGE_wdemog, breaks = cuts, include.lowest = TRUE)
    train_set$ANTIGEN_MATCHED = as.factor(train_set$ANTIGEN_MATCHED)
    test_set$ANTIGEN_MATCHED = as.factor(test_set$ANTIGEN_MATCHED)
    
    
    if (balance == TRUE){
      targetIndex <- which(colnames(df) == target)
      train_set <- df.resample(X= train_set[, -targetIndex], Y= train_set[, targetIndex], positive= pos_class, 
                                 type=typeBalance, percOver=percOver, percUnder=percUnder, verbose=TRUE,
                                 k=kBalance, perc=percBalance, method=methodBalance, w=wBalance)
      train_set <- set.target.name(train_set, targetIndex, target)
      
    }
    
    
   
    
    if(algo == "hc"){
      bn <-  hc(x = train_set, blacklist = blacklist, score = score, maxp = maxp, restart = restart, optimized = optimized)
    }
    
    else if(algo == "tabu"){
      bn <-  tabu(x = train_set, blacklist = blacklist, score = score, maxp = maxp, optimized = optimized,
                   max.iter = max.iter) # tabu = 10, max.tabu = 10,
    }
    
    else if(algo == "gs"){
      bn <- gs(x = train_set, blacklist = blacklist, test = test, alpha,
         B = 3, undirected = undirected, optimized = optimized)
    }
    
    else if(algo == "iamb"){
      bn <- iamb(x = train_set, blacklist = blacklist, test = test, alpha ,
         B = 2, undirected = undirected, optimized = optimized, strict = F, debug = F)
    }
    
    else if(algo == "fast.iamb"){
      bn <- fast.iamb(x = train_set, blacklist = blacklist, test = test, alpha ,
         B = NULL, undirected = undirected, optimized = optimized)
    }
    
    else if(algo == "inter.iamb"){
      bn <- inter.iamb(x = train_set, blacklist = blacklist, test = test, alpha ,
         B = NULL, undirected = undirected, optimized = optimized)
    }
    
    else if(algo == "mmpc"){
      bn <- mmpc(x = train_set, blacklist = blacklist, test = test, alpha ,
         B = NULL, undirected = undirected, optimized = optimized)
    }
    
    else if(algo == "si.hiton.pc"){
      bn <- si.hiton.pc(x = train_set, blacklist = blacklist, test = test, alpha ,
         B = NULL, undirected = undirected, optimized = optimized)
    }
    
    else if(algo == "rsmax2"){
      bn <- rsmax2(x = train_set, blacklist = blacklist, restrict = restrict, maximize = maximize,
             test = test, score = score, alpha , B = NULL,
             maximize.args = list(restart = restart, maxp = maxp),
             optimized = optimized, strict = FALSE) # , tabu = 10, max.tabu = 10
    }
    
    else if(algo == "mmhc"){
      bn <- mmhc(x = train_set, blacklist = blacklist, test = test, score = score, 
                 alpha , B = NULL, restart = restart, optimized = optimized, strict = FALSE)
    }
    
    else if(algo == "glm"){
      l <- sapply(train_set, function(col){
        (is.factor(col) & length(unique(col)) > 1)
      })
      train_set = train_set[ , match(names(l[l==TRUE]), colnames(train_set))]
      test_set = test_set[, colnames(test_set) %in% colnames(train_set)]
      glm <- glm(formula = ANTIGEN_MATCHED ~ .,family = binomial(link='logit'), data = train_set)
      
      glm$xlevels$RACE_wdemog = union(glm$xlevels$RACE_wdemog, levels(test_set$RACE_wdemog))
      glm$xlevels$diag_HEPATOBILIARY = union(glm$xlevels$diag_HEPATOBILIARY, levels(test_set$diag_HEPATOBILIARY))
      
      p <- predict(glm, newdata = test_set, type = "response")
      match_probs[-fld] <- p
      y_test[-fld] <- test_set[ , target]
      next
    }
    
    else if(algo == "NB"){
      bn = naive.bayes(train_set, target)
      pred = predict(bn, test_set, prob = TRUE)
      pred_attr <- attributes(pred)
      
      match_probs[-fld] <- pred_attr$prob[pos_class,]
      y_test[-fld] <- test_set[ , target]
      next
    }
    
    else {
      stop(paste0("unknown algorithm: ", algo))
    }
    
    if(directed(bn) == FALSE){
      print("Not a DAG. Extending it to a DAG ...")
      possibleError <- tryCatch(
        bn <- cextend(bn, strict = TRUE),
        error=function(e) {
          print(e)
          #graphviz.plot(bn, layout = "fdp", main = algo)
        }
      )
                  
    }
    if(inherits(possibleError, "error")){
      print(paste("skip loop", i))
      next
    }
    
    print(bn$nodes$ANTIGEN_MATCHED$parents)
    
    fit <-  bn.fit(x = bn, data = train_set, method = estimator)
    
    pred <- predict(object = fit, data = test_set, node = target, method = "parents", prob = TRUE, debug = FALSE)
    pred_attr <- attributes(pred)
    
    match_probs[-fld] <- pred_attr$prob[pos_class,]
    y_test[-fld] <- test_set[ , target]
    
    rows = 1:length(y_test)
    row.index = c(row.index, rows[-fld])
    #print(sample(1:10^4,1))
    
    ls = getLocalScore(fitted = fit, target = target, DTE = test_set)
    print(local_score)
    local_score$z = c(local_score$z, ls$z)
    local_score$pvalue = c(local_score$pvalue, ls$pvalue)
    
  }
  
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  print(paste0(algo, "_cross-validation time: ",time.taken))
  
  temp <- data.frame("y_test" = y_test, "probs" = match_probs)
  temp = temp[temp$probs != -1 & temp$y_test != "", ]
  temp$y_test <- droplevels(temp$y_test)
  print(paste("result size:", dim(temp)))
  print(table(temp$y_test))
  
  result <- list("y_test" = temp$y_test, "probs" = temp$probs, "cv_time" = time.taken,
                 "alg" = algo, "score" = score, "test" = test, "estimator" = estimator,
                 "row.index" = row.index, "local_score"=local_score)
  
  return (result)
}



set.target.name <-  function (data, targetIndex, targetName){
  colnames(data)[targetIndex] <- target
  return (data)
}

cv <- function(df, df.name, target, pos_class, algorithm, max.iter, score = NULL, test = NULL,
                  estimator, restrict, maximize, blacklist, optimized, results, maxp = NULL,
                  restart = NULL, seed = NULL){
  
  
  if(!(is.null(score))){
    for (algo in algorithm){
      
      for(sc in score){
        
        for (est in estimator){
          print("test")
          result <-  cross_validation(balance = FALSE, df= df, target= target, pos_class= pos_class,
                                      algo = algo, blacklist= blacklist, kfold= 10, score= sc, maxp,
                                      estimator= est, restart, optimized = optimized, max.iter=max.iter, seed = seed )
          
          yt  <-  as.numeric(result$y_test)
          yt[yt == 1]  <-  0
          yt[yt == 2]  <-  1
          pt = result$probs
          ls = result$local_score
          write.table(t(c(algo, sc, est, "",
                          mean(ls$z),mean(ls$pvalue))),
                      "local_score_score-based.csv", 
                      append = TRUE,sep=",", col.names = FALSE)
        
          
        }
      }
    }
  }
  
  else if(!(is.null(test))){
    
    for (algo in algorithm){
      
      for(t in test){
        
        for (est in estimators){
          result <-  cross_validation(balance = FALSE, df= df, target= target, pos_class= pos_class, 
                                      algo = algo, test = t, kfold= 10, estimator = est, optimized = optimized,
                                      seed = seed)
          
          yt  <-  as.numeric(result$y_test)
          yt[yt == 1]  <-  0
          yt[yt == 2]  <-  1
          pt = result$probs
          ls = result$local_score
          write.table(t(c(algo, "", est, t,
                          mean(ls$z),mean(ls$pvalue))),
                      "local_score_const-based.csv", 
                      append = TRUE,sep=",", col.names = FALSE)
          

        }
      }
    }
  }
  
  else if(algorithm == "glm"){
    result = cross_validation(balance = FALSE, df= df, target= target, pos_class= pos_class, 
                              algo = algorithm, kfold= 10, seed = seed)
    
    y_test <- as.numeric(result$y_test)
    #y_test  <-  as.numeric(levels(y_test))[y_test]
    y_test[y_test == 2] = 0
    print(table(y_test))
    measures <- getMeasures(y_test, result$probs)
    
    results[nrow(results)+1, ] <- 
      c(df.name, length(result$y_test), algorithm, "NA", "NA", "NA",
        unlist(measures, use.names = F), result$cv_time)
    
  }
  
  return(results)
}

cv.nb <- function(df, df.name, target, pos_class, results, seed = NULL){
  
  

          result <-  cross_validation(balance = FALSE, df= df, target= target, pos_class= pos_class,
                                      algo = "NB", kfold= 10, seed = seed )
          
          yt  <-  as.numeric(result$y_test)
          yt[yt == 2]  <-  0
          pt = result$probs
          measures <- getMeasures(yt, pt)
          
          p.value = -1
          hl = hl.test.iterate(yt, pt, 12)
          if(! is.null(hl)){
            p.value = hl$p.value
            chisq = hl$statistic
            dof = hl$parameter
            print(hl)
            
          }
          
          save.calib.curve(yt, pt, 10, "Reliability diagram(Naive Bayes)")
          
          results[nrow(results)+1, ] <-  c(df.name, length(result$y_test), "NB", "NA", "NA", "NA",
                                           unlist(measures, use.names = F), p.value, result$cv_time)
          write.csv(results, "result-NB.csv", append = TRUE)
          
 }
 
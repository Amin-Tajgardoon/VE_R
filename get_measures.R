
getMeasures = function(YTE, PTE ){
  #   GETMEASURES Summary of this function goes here
  if (sum(!is.finite(PTE))>0){
    print('there are some nan value in predictions')
  }
  
  if (sum(!is.finite(YTE))>0){
    print('there are some nan value in predictions')
  }
  
  
  idx = is.finite(YTE)&is.finite(PTE);
  YTE = YTE[idx]; PTE = PTE[idx];
  
  res = list();
  res$RMSE = getRMSE(YTE,PTE);
  res$AUC = getAUC(YTE,PTE);
  res$ACC = sum(YTE==(PTE>=0.5))/length(YTE);    
  res$MCE = getMCE(YTE,PTE); #Computing the Max Calibratio Error among all bins 
  res$ECE = getECE(YTE,PTE); #Computing Average Calinration Error on different binns
  return(res)
}

getAUC = function(Actual,Predicted){
  stopifnot( (length(unique(Actual))==2)&(max(unique(Actual))==1))
  nTarget     = sum(Actual == 1);
  nBackground = sum(Actual != 1);
  # Rank data
  R = rank(Predicted, ties.method = "average");  # % 'tiedrank' from Statistics Toolbox  
  #   Calculate AUC
  AUC = (sum(R[Actual == 1]) - (nTarget^2 + nTarget)/2) / (nTarget * nBackground);
  AUC = max(AUC,1-AUC);
  return(AUC)
}

getMCE = function ( Y, P ){
  predictions = P;
  labels = Y;
  sortObj = sort.int(predictions, index.return=TRUE)
  predictions = sortObj$x
  labels = labels[sortObj$ix]
  ordered = cbind(predictions,labels)
  N = length(predictions);
  rest = N%%10;
  S=rep(0,10)
  for (i in 1:10){
    if (i <= rest){
      startIdx = as.integer((i-1) * ceiling(N / 10) + 1)
      endIdx = as.integer(i * ceiling(N / 10))
    }else{
      startIdx = as.integer(rest + (i-1)*floor(N/10)+1)
      endIdx = as.integer(rest + i*floor(N/10))    
    }
    group = ordered[startIdx:endIdx,];
    
    n = dim(group)[1];
    observed = mean(group[,2]);
    expected = mean(group[,1]);
    S[i] = abs(expected-observed);
  }
  res = max(S);
  return(res)
}

getECE = function ( Y, P ){
  predictions = P;
  labels = Y;
  sortObj = sort.int(predictions, index.return=TRUE)
  predictions = sortObj$x
  labels = labels[sortObj$ix]
  ordered = cbind(predictions,labels)
  N = length(predictions);
  rest = N%%10;
  S=rep(0,10)
  W=rep(0,10)
  for (i in 1:10){
    if (i <= rest){
      startIdx = as.integer((i-1) * ceiling(N / 10) + 1)
      endIdx = as.integer(i * ceiling(N / 10))
    }else{
      startIdx = as.integer(rest + (i-1)*floor(N/10)+1)
      endIdx = as.integer(rest + i*floor(N/10))    
    }
    group = ordered[startIdx:endIdx,];
    
    n = dim(group)[1];
    observed = mean(group[,2]);
    expected = mean(group[,1]);
    S[i] = abs(expected-observed);
    W[i] = n/N;
  }
  res = sum(S*W);
  return(res)
}

getRMSE = function( Y, P ){
  res = sqrt(sum((Y-P)*(Y-P))/length(Y));
}
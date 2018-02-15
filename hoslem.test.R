## function for hosmer-lemeshow goodness of fit test (did not report in paper)

hoslem.test <-
  function(x, y, g=10) {
    DNAME <- paste(deparse(substitute(x)), deparse(substitute(y)), sep=", ")
    METHOD <- "Hosmer and Lemeshow goodness of fit (GOF) test"
    yhat <- y
    y <- x
    qq <- unique(quantile(yhat, probs=seq(0, 1, 1/g), na.rm = TRUE))
    if(length(qq) < 2)
      return()
    cutyhat <- cut(yhat,
                   breaks = qq, include.lowest = TRUE)
    observed <- xtabs(cbind("y0"=1 - y, "y1"=y) ~ cutyhat)
    expected <- xtabs(cbind("yhat0"=1 - yhat, "yhat1"=yhat) ~ cutyhat)
    chisq <- sum((observed - expected)^2 / expected)
    PVAL = 1 - pchisq(chisq, g - 2)
    PARAMETER <- g - 2
    names(chisq) <- "X-squared"
    names(PARAMETER) <- "df"
    structure(list(statistic = chisq, parameter = PARAMETER, 
                   p.value = PVAL, method = METHOD, data.name = DNAME, observed = observed, 
                   expected = expected), class = "htest")
  }

hl.test.iterate = function(yt, pt, g){
  p.value = NaN
  while(is.nan(p.value)){
    if(g > 50){
      print("MAXIMUM number of groups reached in HL test")
      break;
    }
    hl=hoslem.test(yt,pt, g);
    g = g+1
    if(is.null(hl))
      next
    p.value = hl$p.value
  }
  return(hl)
}
## model-averaging using booststrap sampling
## assumes training data (df) is loaded in memory
## saves the final model properties

set.seed(1234)
b.s = boot.strength(df, algorithm = "tabu",
                    algorithm.args = list(blacklist = bl, score = "bic", maxp = 8, optimized=TRUE)
                    ,R = 10, m = nrow(df))

set.seed(1234)
b.s2 = boot.strength(df, algorithm = "gs"
                     , algorithm.args = list(blacklist = bl, test = "mi", undirected = FALSE, optimized = TRUE)
                    ,R = 10, m = nrow(df))

b.s[b.s$to == "ANTIGEN_MATCHED" & b.s$strength > 0,]

b.s2[b.s2$to == "ANTIGEN_MATCHED" & b.s2$strength > 0,]

head(b.s2[b.s2$to == "ANTIGEN_MATCHED",])


a.n = averaged.network(strength = b.s, threshold = .5)


df2 = df
smbin = smbinning(df = df , x = "AGE_wdemog", y = "ANTIGEN_MATCHED", p = 0.05)
cuts = smbin$bands #c(min(df$AGE_wdemog), smbin$cuts, max(df$AGE_wdemog))
df2$AGE_wdemog = cut(x = df2$AGE_wdemog, breaks = cuts, include.lowest = TRUE)
df2$ANTIGEN_MATCHED = as.factor(df2$ANTIGEN_MATCHED)

b.s3 = boot.strength(df2,
                     algorithm = "tabu",
                     algorithm.args = 
                       list(blacklist = bl, score = "bde",
                            maxp = 8, optimized=TRUE)
                    ,R = 100, m = nrow(df2))


b.s3[b.s3$to == "ANTIGEN_MATCHED" & b.s3$strength > 0,]

bs.s3.ord = b.s3[order(b.s3$strength, decreasing = TRUE),]

bs.s3.ord[bs.s3.ord$to == "ANTIGEN_MATCHED" & bs.s3.ord$strength > 0,]

b.s3[(b.s3$to == "AGE_wdemog" | b.s3$from == "AGE_wdemog")
     & b.s3$strength > 0,]

set.seed(1234)
b.s2 = boot.strength(df2, algorithm = "fast.iamb"
                     , algorithm.args = list(blacklist = bl, test = "mi", undirected = FALSE, optimized = TRUE)
                     ,R = 100, m = nrow(df2))

bs.s2.ord = b.s2[order(b.s2$strength, decreasing = TRUE),]
bs.s2.ord[bs.s2.ord$to == "ANTIGEN_MATCHED" & bs.s2.ord$strength > 0,]


algo.const <- c( "fast.iamb", "si.hiton.pc", "mmpc") #c( "iamb", "fast.iamb", "si.hiton.pc", "mmpc", "gs")
tests <-  c( "x2-adf","mi") #, "mc-mi", "smc-mi") #c( "mc-mi", "x2-adf")

for(alg in algo.const){
  for(t in tests){
    print(paste(alg,t))
    b.s = boot.strength(df2, algorithm = alg
                         , algorithm.args = list(blacklist = bl, test = t, undirected = FALSE, optimized = TRUE)
                         ,R = 100, m = nrow(df2))
    
    bs.s.ord = b.s[order(b.s$strength, decreasing = TRUE),]
    print(bs.s.ord[bs.s.ord$to == "ANTIGEN_MATCHED" & bs.s.ord$strength > 0,])
    print("-----------------------------------------")
    print("-----------------------------------------")
    
    
  }
}


algo.score <-  c("tabu","hc") 
scores <- c("bic",  "k2", "bde") #,"aic","loglik")

set.seed(1234)

for(alg in algo.score){
  for(s in scores){
    print(paste(alg,s))
    b.s = boot.strength(df2, algorithm = alg
                        , algorithm.args = 
                          list(blacklist = bl, score = s, maxp = 8, optimized = TRUE)
                        ,R = 100, m = nrow(df2))
    
    bs.s.ord = b.s[order(b.s$strength, decreasing = TRUE),]
    res = bs.s.ord[bs.s.ord$to == "ANTIGEN_MATCHED" & bs.s.ord$strength > 0,]
    res$alg = rep(paste(alg,s,sep = "&"),nrow(res))
    res = res[,c(5,1,2,3,4)]
    write.table(res,file = "sig-parents-score-based.csv", append = TRUE, col.names = FALSE)
    
  }
}
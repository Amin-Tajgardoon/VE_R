library(bnlearn)
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/cv.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/get_measures.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/data_preprocess.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/bivariate_analysis.R")


dataDir <- 'C:/Users/mot16/projects/master/data/final/'
outDir <- 'C:/Users/mot16/projects/master/output/'
#fileName <- "gsk_108134_NoMedicine_discrete_age_tempr_MATCH_ONLY_ColNameFixed_replacedNullsInR_removedConstants"

fileName <- "gsk_108134_joined_allvars_discrete_age_tempr_final_MATCH_ONLY"

dataPath <- paste(dataDir,fileName,'.csv',sep = "")

flu <-  read.csv(file = dataPath, stringsAsFactors = F)

flu <- data_preprocess(flu)

#write.csv(x = flu,row.names = F, file = paste0(dataDir,fileName,"_cleanedInR",".csv"))

target <- "ANTIGEN_MATCHED"

# forbid edge from class variable
x_names = colnames(flu)[colnames(flu) != target]
bl = data.frame(from = rep(target,length(x_names)), to = x_names)

algo <-  "hc" ## "tabu"
score <-  "mbde" ## "loglik" ## "bde" (best so far with HC) ## "k2"

result <-  cross_validation(balance = FALSE, typeBalance="ubSMOTE", percOver=200, percUnder=200,
                            kBalance=5, percBalance=50, methodBalance="percPos", wBalance=NULL, 
                            data= flu, target= target, pos_class= "MATCH", algo = algo,
                            blacklist= bl, kfold= 10, score= score, maxp= 6, estimator= "bayes")

result$y_test  <-  as.numeric(result$y_test)
result$y_test[result$y_test == 2]  <-  0

measures <- getMeasures(result$y_test, result$probs)
print(paste0(algo," performance measures: "))
print(measures)

outputPath <-  paste0(outDir,fileName,'_',algo,'_',score,'_oputput_YTE_PTE.csv')
                     
#write.csv(x = result,row.names = F, file = outputPath)



# learn structure with hillclimbing algorithm

#hc = hc(x = flu, blacklist = bl, score = "k2", maxp = 5, optimized = TRUE)

# estimate parameters
# fit <-  bn.fit(x = hc, data = flu, method = "bayes")



#cv.k2 <-  bn.cv(flu, bn = "hc", k = 10, algorithm.args = list(score = "k2", blacklist = bl, maxp = 7, optimized = T), fit = "bayes", loss = "pred", loss.args = list(target= "ANTIGEN_MATCHED") )

#OBS = unlist(lapply(cv.k2, `[[`, "observed"))
#PRED = unlist(lapply(cv.k2, `[[`, "predicted"))

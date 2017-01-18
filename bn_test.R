library(bnlearn)
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/get_measures.R")

dataDir <- 'C:/Users/mot16/projects/master/data/final/'
outDir <- 'C:/Users/mot16/projects/master/output/'
fileName <- "gsk_108134_NoMedicine_discrete_age_tempr_MATCH_ONLY_ColNameFixed_replacedNullsInR_removedConstants"

dataPath <- paste(dataDir,fileName,'.csv',sep = "")

flu <-  read.csv(file = dataPath, stringsAsFactors = F)

#  replace Nulls with 'N'
#flu[is.na(flu)] <- 'N'
#  remove columns with no variance
#flu = flu[, sapply(flu,function(x) length(unique(x))) > 1]

# convert columns to factor 
flu <- data.frame(lapply(flu, as.factor), stringsAsFactors=FALSE)

target <- "ANTIGEN_MATCHED"

# forbid edge from class variable
x_names = colnames(flu)[colnames(flu) != target]
bl = data.frame(from = rep(target,length(x_names)), to = x_names)

       
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/cv.R")

score = "k2"
result <-  cross_validation(data= flu, target= target, pos_class= "MATCH", blacklist= bl, k= 10, score= score, maxp= 6, estimator= "bayes")

result$y_test  <-  as.numeric(result$y_test)
result$y_test[result$y_test == 2]  <-  0

measures <- getMeasures(result$y_test, result$probs)
print(measures)

outputPath <-  paste(outDir,fileName,'_HC_',score,'_oputput_YTE_PTE.csv', sep = "")
                     
#write.csv(x = result,row.names = F, file = outputPath)



# learn structure with hillclimbing algorithm

#hc = hc(x = flu, blacklist = bl, score = "k2", maxp = 5, optimized = TRUE)

# estimate parameters
# fit <-  bn.fit(x = hc, data = flu, method = "bayes")



#cv.k2 <-  bn.cv(flu, bn = "hc", k = 10, algorithm.args = list(score = "k2", blacklist = bl, maxp = 7, optimized = T), fit = "bayes", loss = "pred", loss.args = list(target= "ANTIGEN_MATCHED") )

#OBS = unlist(lapply(cv.k2, `[[`, "observed"))
#PRED = unlist(lapply(cv.k2, `[[`, "predicted"))

## reads the data and constructs the results by calling cross-validation function

library(bnlearn)
debugSource("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/cv.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/cv.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/get_measures.R")
debugSource("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/data_preprocess.R")
debugSource("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/flu_data_preprocess.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/bivariate_analysis.R")

seed = 1234
set.seed(seed)

dataDir <- 'C:/Users/mot16/projects/master/data/'
outDir <- 'C:/Users/mot16/projects/master/output/'
fileName <- "gsk_108134_joined_allvars"
dataPath <- paste(dataDir,fileName,'.csv',sep = "")
target <- "ANTIGEN_MATCHED"
pos_class <- "MATCH"

## read & clean data
df <- data.preprocess(read.csv(file = dataPath, stringsAsFactors = F), seed)

results  <-  data.frame(dataset= character(), test_size = numeric(), algorithm = character(), score = character(), estimator = character(),
                        test = character(), RMSE = numeric(), AUC = numeric(), ACC = numeric(),
                        MCE = numeric(), ECE = numeric(), R2 = numeric(),
                        run_time = character(), stringsAsFactors = F)

algo.score <-  c("tabu","hc") 
scores <- c("bic",  "k2", "bde") #,"aic","loglik")
estimators <- c( "bayes") #"mle"

algo.const <- c( "fast.iamb", "si.hiton.pc", "mmpc")
tests <-  c( "x2-adf","mi") 

optimized = TRUE
# forbid edge from class variable
x_names = colnames(df)[colnames(df) != target]
bl = data.frame(from = rep(target,length(x_names)), to = x_names)

results <- cv(df = df,df.name = df.name, target = target, pos_class = pos_class, algorithm = algo.score, max.iter = Inf,
                  score = scores, estimator = estimators, blacklist = bl, results = results, optimized = optimized, maxp = 8, 
                  restart = 5, seed = seed)

#results <- cv(df = df,df.name = df.name, target = target, pos_class = pos_class, algorithm = algo.const,
#                test = tests, estimator = estimators, blacklist = bl, results = results, optimized = optimized, seed = seed)


#results <- cv.nb(df = df,df.name = df.name, target = target, pos_class = pos_class,
#              results = results, seed = seed)
  
  

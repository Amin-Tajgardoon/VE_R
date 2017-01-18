dataPath = "C:/Users/mot16/projects/master/data/final/gsk_108134_NoMedicine_discrete_age_tempr_MATCH_ONLY_ColNameFixed.csv"
flu = read.csv(file = dataPath, stringsAsFactors = T)

fit = bn.fit(x = hc, data = flu, method = "bayes")
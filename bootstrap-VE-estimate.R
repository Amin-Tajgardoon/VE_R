## infer conditional probabilities from the Bayesian network, each for 100 times
## calculate vaccine efficacy (VE)  from probabilities, plus mean and standard deviation
## calculate vaccine efficacy from dataset, report risk-ratio and 95% confidence interval
## perform fisher exact test for each VE

library(coda)
library(fmsb)

set.seed(1234)

# Discretize AGE_wdemog
df2 = df
smbin = smbinning(df = df , x = "AGE_wdemog", y = "ANTIGEN_MATCHED", p = 0.05)
cuts = smbin$bands #c(min(df$AGE_wdemog), smbin$cuts, max(df$AGE_wdemog))
df2$AGE_wdemog = cut(x = df2$AGE_wdemog, breaks = cuts, include.lowest = TRUE)
df2$ANTIGEN_MATCHED = as.factor(df2$ANTIGEN_MATCHED)

N = 101
set.seed(1234)
b.s = boot.strength(df2,
                     algorithm = "tabu",
                     algorithm.args = 
                       list(blacklist = bl, score = "bde",
                            maxp = 8, optimized=TRUE)
                     ,R = N, m = nrow(df2))

## build average network from saved bootstrapped edge-strengths
saveRDS(b.s, file = "boot-strength.rds")
#b.s2 = readRDS(file = "boot-strength.rds")

avrg.net = averaged.network(strength = b.s, threshold = .5)
if(!directed(avrg.net))
  avrg.net = cextend(x = avrg.net, strict = TRUE)
fit <-  bn.fit(x = avrg.net, data = df2, method = "bayes")
  

ARV = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                 evidence = (VACCINE_expogn == "FLUARIX VACCINE") ))

ARU = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                 evidence = (VACCINE_expogn == "PLACEBO") ) )

ARV_AGE1 = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                      evidence = (VACCINE_expogn == "FLUARIX VACCINE" & AGE_wdemog == "[18,45]") ) )
ARV_AGE2 = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                      evidence = (VACCINE_expogn == "FLUARIX VACCINE" & AGE_wdemog == "(45,66]") ) )
ARU_AGE1 = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                      evidence = (VACCINE_expogn == "PLACEBO" & AGE_wdemog == "[18,45]") ) )
ARU_AGE2 = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                      evidence = (VACCINE_expogn == "PLACEBO" & AGE_wdemog == "(45,66]") ) )

ARV_F = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                   evidence = (VACCINE_expogn == "FLUARIX VACCINE" & SEX_wdemog == 'F') ) )
ARV_M = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                   evidence = (VACCINE_expogn == "FLUARIX VACCINE" & SEX_wdemog == 'M') ) )
ARU_F = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                   evidence = (VACCINE_expogn == "PLACEBO" & SEX_wdemog == 'F') ) )
ARU_M = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                   evidence = (VACCINE_expogn == "PLACEBO" & SEX_wdemog == "M") ) )

ARV_S = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                                        evidence = (VACCINE_expogn == "FLUARIX VACCINE" & VIAL_TYP_expogn == 'S') ) )
ARV_R = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                                        evidence = (VACCINE_expogn == "FLUARIX VACCINE" & VIAL_TYP_expogn == 'R') ) )


ARU_S = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                                        evidence = (VACCINE_expogn == "PLACEBO" & VIAL_TYP_expogn == 'S') ) )
ARU_R = sapply(1:N, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
                                        evidence = (VACCINE_expogn == "PLACEBO" & VIAL_TYP_expogn == "R") ) )


hdi.arv = HPDinterval(as.mcmc(ARV),prob = 0.95)
hdi.aru = HPDinterval(as.mcmc(ARU),prob = 0.95)
arv.m = median(ARV)
aru.m = median(ARU)
1 - c(hdi.arv[1]/hdi.aru[1], arv.m/aru.m, hdi.arv[2]/hdi.aru[2])

1 - mean(ARV)/mean(ARU)
ve = 1 - (ARV/ARU)
# ve = 1 - (ARV_F/ARU_F)
# ve = 1 - (ARV_M/ARU_M)
# ve = 1 - (ARV_AGE2/ARU_AGE2)
# ve = 1 - (ARV_AGE1/ARU_AGE1)
# ve = 1 - (ARV_S/ARU_S)
# ve = 1 - (ARV_R/ARU_R)


hdi.ve = function(ar1, ar2){
  ve = 1 - (ar1/ar2)
  hdi = HPDinterval(as.mcmc(ve),prob = 0.95)
  interval = round(c(hdi[1], mean(ve), hdi[2]),2)
  return(interval)
}

hdi.ve(ARV, ARU)
hdi.ve(ARV_F, ARU_F)
hdi.ve(ARV_M, ARU_M)
hdi.ve(ARV_AGE2, ARU_AGE2)
hdi.ve(ARV_AGE1, ARU_AGE1)
hdi.ve(ARV_S, ARU_S)
hdi.ve(ARV_R, ARU_R)

round(c(mean(ve)-2*(sqrt(var(ve))), mean(ve), mean(ve)+2*(sqrt(var(ve))))
      ,2)

t = table(df2$VACCINE_expogn,  df2$ANTIGEN_MATCHED)
# t = table(df2$VACCINE_expogn[df2$SEX_wdemog == "F"],
#           df2$ANTIGEN_MATCHED[df2$SEX_wdemog == "F"])
# t = table(df2$VACCINE_expogn[df2$SEX_wdemog == "M"],
#           df2$ANTIGEN_MATCHED[df2$SEX_wdemog == "M"])
# t = table(df2$VACCINE_expogn[df2$AGE_wdemog == "[18,45]"],
#           df2$ANTIGEN_MATCHED[df2$AGE_wdemog == "[18,45]"])
# t = table(df2$VACCINE_expogn[df2$AGE_wdemog == "(45,66]"],
#           df2$ANTIGEN_MATCHED[df2$AGE_wdemog == "(45,66]"])
# t = table(df2$VACCINE_expogn[df2$VIAL_TYP_expogn == "S"],
#           df2$ANTIGEN_MATCHED[df2$VIAL_TYP_expogn == "S"])
# t = table(df2$VACCINE_expogn[df2$VIAL_TYP_expogn == "R"],
#           df2$ANTIGEN_MATCHED[df2$VIAL_TYP_expogn == "R"])

rr = epitools::riskratio(t, rev = "r")
rr
round(1- rr$measure,2)

t = table(df2$VACCINE_expogn,  df2$ANTIGEN_MATCHED)
t = table(df2$VACCINE_expogn[df2$AGE_wdemog == "(45,66]"],
      df2$ANTIGEN_MATCHED[df2$AGE_wdemog == "(45,66]"])

rr = riskratio(X = 14, Y = 12, m1 = 1917 + 14, m2 = 931 + 12, p.calc.by.independence = TRUE)
1- rr$estimate
1- rr$conf.int

ftest = fisher.test(table(df2$VACCINE_expogn[df2$AGE_wdemog == "(45,66]"],
                                                df2$ANTIGEN_MATCHED[df2$AGE_wdemog == "(45,66]"]))

ftest = fisher.test(table(df2$VACCINE_expogn,
                          df2$ANTIGEN_MATCHED))

exp(ftest$estimate)

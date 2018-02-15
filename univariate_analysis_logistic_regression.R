## performs univarite analysis using logistic regression

y = relevel(df$ANTIGEN_MATCHED, ref = "NMATCH")

f = function(i){
  lr = glm(y ~ df[,i], 
           family = binomial(link = "logit"))
  coefs = summary(lr)$coefficients
  l = coefs[coefs[-1,4] <= 0.05, 1]
  if(length(l) > 0)
    return(i)
  return(0)
}


s = sapply(1:ncol(df) , function(x) f(x))


cols = names(df)[s[s > 0]]



f2 = function(var.name){
  formula = as.formula(paste("ANTIGEN_MATCHED", "~", var.name))
  lr = glm(formula, data = df, 
           family = binomial(link = "logit"))
  coefs = summary(lr)$coefficients[-1,]
  return(coefs)
  #return(summary(lr))
}

res = sapply(cols, function(x) f2(x))

center = relevel(df$CENTER_wdemog, ref = "1")
country = relevel(df$CTY_COD_wdemog, ref = "FI")
treatmnt = relevel(df$TREATMNT_wdemog, ref = "Placebo")
age = factor(df$AGE_wdemog, ordered = FALSE)
vaccine = relevel(df$VACCINE_expogn, ref = "PLACEBO")
vial.type = relevel(df$VIAL_TYP_expogn, ref = "R")
agecat = factor(df$AGECAT_pid, ordered = FALSE)
vac.atp = relevel(df$TYPE_VAC_wnpap, ref = "0")
diag.immune = relevel(df$diag_IMMUNE.SYSTEM..INCL.ALLERGIES..AUTOIMMUNE.DISORDERS., 
                      ref = "N")
medic_984601 = df$medic_984601
medic_845801 = df$medic_845801


var.names = c("country", "treatmnt", "age", "vaccine", "vial.type", 
              "agecat", "vac.atp")
# MODEL SUMMARY
for(nm in var.names){
  formula = as.formula(paste("y", "~", nm))
  lr = glm(formula, 
           family = binomial(link = "logit"))
  print(summary(lr))
}

# CONFINT
for(nm in var.names){
  formula = as.formula(paste("y", "~", nm))
  lr = glm(formula, family = binomial(link = "logit"))
  print(confint(lr))
}

dummies = model.matrix(~treatmnt)

lr = glm(y ~ vaccine + agecat, 
         family = binomial(link = "logit"))
print(summary(lr))
confint(lr)

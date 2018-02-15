## process and transform the raw data

source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/data_preprocess.R")
source("C:/Users/mot16/Box Sync/Master Project/BN representation of VE/code and data-docs/code/VE_R/bivariate_analysis.R")

library(discretization)

## prepares FLU dfset
data.preprocess <-  function(df, seed){
  
  set.seed(seed)
  
  target <- "ANTIGEN_MATCHED"
  
  # merge with new_medic_by_day_21
  medic.new = read.csv("C:/Users/mot16/projects/master/data/new_medic_by_day21.csv", stringsAsFactors = F)
  medicCols = grep("^medic_", names(medic.new))
  medic.new = medic.new[c("PID", names(medic.new)[medicCols])]
  
  medicCols = grep("^medic_", names(df))
  df <- df[-medicCols]

  m <- merge(x = df, y = medic.new, by = "PID", all.x = TRUE)
  
  df <- m
  
  ## list all unique values in the dfframe
  #as.character(unique(unlist(flu)))
  
  ##  replace Nulls with 'N'
  df <- replace.nulls(df, 'N')
  
  ##  replace 'ND' & 'MC' with 'N'
  df <- replace.values(df, toReplace = c('ND', 'MC', 'MISSING CONFIRMED'), value = 'N')
  
  
  df <- clean.concvac(df)
  
  ##  remove columns with no variance
  df <-  remove.constant(df)

  any_dummies <- names(df)[grep("^any_", names(df))]
  df <- df[ , !(colnames(df) %in% c(any_dummies, "PID", "DROP_OUT_wconc", "LC_GC_wconc", "SYMP_VAL_wsolpre",
                                "SYMP_SIT_wsolpre", "SAE_CNT_wconc", "AE_FLAG_wconc", "OUTCOME_wpneumo",
                                "diag_MC", "AGECAT_pid"))]
  
  df <-  prepare.target(df, target)
  
  # df <-  discretize.age(df,seed)
  # df$AGE_wdemog  <-  as.ordered(df$AGE_wdemog)
  
  #df$AGE_wdemog = cut(x = df$AGE_wdemog, breaks = 12, ordered_result = T, include.lowest = T, right = F)
  
  ## convert columns to factor 
  df <- to.factor(df)
  
  df$AGE_wdemog = as.numeric(levels(df$AGE_wdemog))[df$AGE_wdemog]
  df$ANTIGEN_MATCHED = as.numeric(df$ANTIGEN_MATCHED)
  df$ANTIGEN_MATCHED[df$ANTIGEN_MATCHED == 2] = 0
  
  #df$AGECAT_pid  <-  as.ordered(df$AGECAT_pid)

  return (df)
  
}

discretize.age <- function(df, seed){
  set.seed(seed)
  dis <- chiM(df[, c('AGE_wdemog', 'ANTIGEN_MATCHED')])
  print(paste("age cutpoits:", dis$cutp))
  
  df$AGE_wdemog <- dis$Disc.data$AGE_wdemog
  
  return (df)
  
}

prepare.target <-  function(df , target){
  
  df <-  df[ , !(colnames(df) %in% c("RAWRES", "viral_res_binary", "rt_PCR"))]
  
  df[target] <- "NMATCH"
  
  strains <- c("H1", "H3", "B")
  
  df[apply(df[,strains], MARGIN = 1, function(x) any(x == 'MATCH')), target] <- "MATCH"
  
  df <-  df[, !(colnames(df) %in% strains)]
  
  return (df)
}


clean.concvac <-  function(df){
  
  conc_vac_cols <- c('conc_vac_1','conc_vac_2','conc_vac_3','conc_vac_4','conc_vac_5',
                     'conc_vac_6','conc_vac_7','conc_vac_8','conc_vac_9','conc_vac_10',
                     'conc_vac_11','conc_vac_12','conc_vac_13','conc_vac_14','conc_vac_15')
  
  df[conc_vac_cols] <- "N"
  
  concvac_1 <- c("concvac_A.HEPATITIS", "concvac_A.HEPATITIS.VACCINE", "concvac_EPAXAL", "concvac_HAVRIX",
                 "concvac_HAVRIX.I", "concvac_HAVRIX.SECOND.DOSE", "concvac_HAWRIX", "concvac_HEPATITIS.A")
  
  concvac_2 <- c("concvac_B.HAPATITE.VACCIN.", "concvac_ENGERIX", "concvac_ENGERIX.B", "concvac_ENGERIX.B.1",
                 "concvac_HBV", "concvac_HEPATIT.B", "concvac_HEPATITIS.B", "concvac_HEPATITIS.B.VACCINE")
  
  concvac_3 <- c("concvac_DIFTERIA.TETANUS", "concvac_DIPHTERIA.TETANUS", "concvac_DIPHTHERIA.TETANUS", 
                 "concvac_DITEBOOSTER", "concvac_DT", "concvac_TD", "concvac_TETANUS...DIFTERIA", 
                 "concvac_TETANUS.AND.DIPHTERIA", "concvac_TETANUS.D", "concvac_TETANUS.DIFTERIA",
                 "concvac_TETANUS.D.1", "concvac_TETANUS.DIFTERIA.1",  "concvac_TETANUS.DIPHTERIA", 
                 "concvac_TETANUS.DIFTERIA.2", "concvac_TETENUS.D.VACCINATION" )
  
  concvac_4 <- c("concvac_HEPATIT.A.B", "concvac_TVINRIX", "concvac_TWINRIX", "concvac_TWINRIX.I",
                 "concvac_TWINRIX.II", "concvac_TWINRIX.III" )
  
  concvac_5 <- c("concvac_IPV", "concvac_POLIO", "concvac_POLIO.VACCINE")
  
  concvac_6 <- c("concvac_ENCEPUR", "concvac_ENCEPUR.II", "concvac_TBE", "concvac_TBE.VACCINATION",
                 "concvac_TBE.VACCINE", "concvac_TICK.VACCINATION", "concvac_TICOVAC", "concvac_TICOVAX")
  
  concvac_7 <- c("concvac_ALTEANA", "concvac_TETANUS", "concvac_TETANUS.VACCINE", "concvac_TETANUS.VACCINE.1")
  
  concvac_8 <- c("concvac_TETANUS.D.AND.POLIO")
  
  concvac_9 <-  c( "concvac_STAMARIL","concvac_YELLOW.FEAVER..VACCINE", "concvac_YELLOW.FEVER.VACCINE")
  
  concvac_10 <- c("concvac_HERPES.ZOSTER.VACCINE")
  
  concvac_11 <- c("concvac_IMOVAX")
  
  concvac_12 <- c("concvac_MENINGOVAX.A.C")
  
  concvac_13 <- c("concvac_MPR")
  
  concvac_14 <- c("concvac_TYPHERIX")
  
  concvac_15 <- c( "concvac_VARILIX", "concvac_VARILRIX")
  
  df$conc_vac_1[apply(df[concvac_1], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_2[apply(df[concvac_2], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_3[apply(df[concvac_3], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_4[apply(df[concvac_4], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_5[apply(df[concvac_5], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_6[apply(df[concvac_6], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_7[apply(df[concvac_7], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_8[apply(df[concvac_8], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_9[apply(df[concvac_9], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_10[apply(df[concvac_10], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_11[apply(df[concvac_11], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_12[apply(df[concvac_12], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_13[apply(df[concvac_13], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_14[apply(df[concvac_14], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  df$conc_vac_15[apply(df[concvac_15], MARGIN = 1, function(x) any(x == "Y"))] <- "Y"
  
  to.drop <- colnames(df)[grep("concvac_", colnames(df))]
  
  df <- df[ , !(colnames(df) %in% to.drop)]
  
  return (df)
  
}


list.of.dataset <- function(flu, target){
  
  ## remove medic columns
  flu.no.medic  <- flu[, ! grepl("(medic_|any_medicine)" , colnames(flu))]
  
  ## remove medic and concvac columns
  flu.no.medic.concvac  <- flu[, ! grepl("(medic_|any_medicine|conc_vac_|any_conc_vac)" , colnames(flu))]
  
  ## remove medic in correlated vars
  pval <- chi.square.test(flu, target)
  corr <- get.correlated(pval)
  corr.no_med.names <- names(corr[grep("(medic_|any_medicine)", names(corr))])
  flu.no.medic.in.corr <- flu[, !(colnames(flu) %in% corr.no_med.names)]
  
  ## keep only correlated vars. remove medicine
  pval <- chi.square.test(flu, target)
  corr <- get.correlated(pval)
  corr.no_med.names <- names(corr[grep("(medic_|any_medicine)", names(corr))])
  flu.only.corr <- flu[, colnames(flu) %in% names(corr) ]
  flu.only.corr <- flu.only.corr[, !(colnames(flu.only.corr) %in% corr.no_med.names) ]
  flu.only.corr[target] <- flu[target]
  
  return (list(flu.only.corr = flu.only.corr, flu.no.medic = flu.no.medic, flu.no.medic.concvac = flu.no.medic.concvac,
               flu.no.medic.in.corr = flu.no.medic.in.corr))
}
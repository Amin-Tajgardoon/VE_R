## performs exact inference (only approximate inference was reported in the final paper)

#install.packages("gRain")
#source("https://bioconductor.org/biocLite.R")
#biocLite("RBGL")

library("gRain")

jtree = compile(as.grain(fit))
jprop = setEvidence(jtree, nodes = c("VACCINE_expogn","AGE_wdemog"),
                   states = c("FLUARIX VACCINE", "[18,45]"))

querygrain(jprop, nodes = "ANTIGEN_MATCHED")$ANTIGEN_MATCHED


ls = sapply(1:100, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
        evidence = (VACCINE_expogn == "FLUARIX VACCINE" &
                      AGE_wdemog == "[18,45]") ) )


lw = sapply(1:100, function(i) cpquery(fitted = fit, event = (ANTIGEN_MATCHED == "1"),
        evidence = list(VACCINE_expogn = "FLUARIX VACCINE",
                      AGE_wdemog = "[18,45]"),method = "lw" ) )

mean(ls)
summary(ls)
mean(lw)
summary(lw)

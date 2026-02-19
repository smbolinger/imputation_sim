
library(brglm2)
# realDat <- read.csv("dat_complete.csv")
realDat <- readRDS("dat_complete.rds")
allMods <- readLines("modList.txt")
modList <- allMods[c(1,8,16)]

print(str(realDat))

# realDat$species <- relevel(as.factor(realDat$species), ref = "LETE")
realDat$species  <- as.factor(realDat$species)
realDat$cam_fate <- relevel(as.factor(realDat$cam_fate), ref="H")
#fitReal <- glm(is_u ~ nest_age * species + obs_int + cam_fate + fdate, family=binomial, data=realDat, method=brglm2::brglmFit)

fitReal <- function(resp, dataMod, modlist, iter=500){
  mods  <- list()
  for (m in seq_along(modlist)){
    mods[[m]] <- glm( as.formula( paste( resp, modlist[m], sep=" " )  ),
          data=dataMod,
          family=binomial,
          method=brglm2::brglm_fit,
          control=brglmControl(maxit=iter)
        )
  }
  names(mods) <- modlist
  return(mods)
  # glm(as.formulais_u ~ nest_age * species + obs_int + cam_fate + fdate, family=binomial, data=realDat, method=brglm2::brglmFit)
}

# fitIsU <- fitReal("is_u", realDat, modlist=mods4sim)
fitIsU <- fitReal("is_u", realDat, modlist=modList)
print(lapply(fitIsU, summary))
fitHM <- fitReal("HF_mis", realDat, modlist=modList)
# print(summary(fitHM))
print(lapply(fitHM, summary))
# fits <- c(fitIsU, fitHM)

betasU <- lapply(fitIsU, function(x) coef(x))
betasM <- lapply(fitHM, function(x) coef(x))
cat("\nbetas, CONI as ref level:\n")
# print(betasU)
# print(betasM)


betas <- list("is_u"=betasU, "hfm"=betasM)
print(betas)

# saveRDS(betas, "beta_list.rds") # all mods
saveRDS(betas, "betas2.rds") # all mods



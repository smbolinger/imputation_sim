
realDat <- read.csv("dat_complete.csv")
modList <- readLines("modList.txt")

realDat$species <- relevel(as.factor(realDat$species), ref = "LETE")
realDat$cam_fate <- relevel(as.factor(realDat$cam_fate), ref="H")
#fitReal <- glm(is_u ~ nest_age * species + obs_int + cam_fate + fdate, family=binomial, data=realDat, method=brglm2::brglmFit)

fitReal <- function(resp, dataMod, modlist, iter=500){
  mods <- list()
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

fitIsU <- fitReal("is_u", realDat, modlist=mods4sim)
fitHM <- fitReal("HF_mis", realDat, modlist=mods4sim)
# fits <- c(fitIsU, fitHM)

betasU <- lapply(fitIsU, function(x) coef(x))
betasM <- lapply(fitHM, function(x) coef(x))

betas <- list("is_u"=betasU, "hfm"=betasM)
saveRDS(betas, "beta_list.rds")

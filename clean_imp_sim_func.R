
###### IMPORT FUNCTIONS ################################################################
########################################################################################
options(include.full.call.stack = FALSE)               # or configure it globally
# source("missingness_tab.R")
debugging <- FALSE # for uickly setting values when working in the file with the functions

######  RANDOM FUNCTIONS ################################################################
########################################################################################

num_cols <- c("nest_age", "obs_int", "fdate")
# means <- function(realDat) sapply(realDat[,c(8:10)], mean)
means <- function(realDat) sapply(realDat[num_cols], mean)
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
  return(mods)
  # glm(as.formulais_u ~ nest_age * species + obs_int + cam_fate + fdate, family=binomial, data=realDat, method=brglm2::brglmFit)
}

###### FACTORS --> DUMMIES #############################################################
########################################################################################

if(FALSE){
  dat <- dat4sim
}
add_dummy <- function(dat, vb=0){
    dummy_vars <- model.matrix(~ cam_fate -1, data = as.data.table(dat))
    dummy_vars <- dummy_vars[,-1]
    dat$speciesCONI <- ifelse(dat$species=="CONI", 1, 0)
    dat<- dat[, c( "speciesCONI", "obs_int","nest_age", "fdate", "HF_mis", "is_u")] # if LETE is the reference level, you have speciesCONI as the dummy variable
    dat <- cbind(dat, dummy_vars)
    return(dat)
}

###### DUMMIES --> FACTORS #############################################################
########################################################################################

if(FALSE){
  dat <-amp_out_wt$amp
}
add_fact <- function(dat, facToNum=FALSE, vb=0){ 
    #cat("\n NOTE: can speed up conversion back to factor with tidytable\n")
    # I already know what the dummy var names are in this case; not a general purpose function
    dat <- as.data.frame(dat) %>% 
        mutate(cam_fate = case_when(
              cam_fateA==1 ~ "A",
                        # cam_fateH==1 ~ "H",
              cam_fateHu==1 ~ "Hu",
              cam_fateF==1 ~ "F",
              cam_fateS==1 ~ "S",
              cam_fateD==1 ~ "D",
              .default="H"), # so the NAs were becoming "H"
              species = if_else(speciesCONI==1, "CONI", "LETE")) %>%
        mutate(cam_fate = if_else(is.na(HF_mis), NA, cam_fate))

    if(facToNum){
        dat <- dat %>%
        rename(cfate = cam_fate,
             spp   = species) %>%
        mutate( cam_fate = case_match(as.character(cfate),
                        "H"  ~ 0,
                        "F"  ~ 1,
                        "D"  ~ 2,
                        "S"  ~ 3,
                        "A"  ~ 4,
                        "Hu" ~ 5,
                        "Ca" ~ 6,
                        "U"  ~ 7,
                        NA   ~ NA
        ),
        species = ifelse(spp=="LETE", 0, 1),
        isu=as.character(is_u),
        HFmis = as.character(HF_mis)
        )
    } else {
        dat <- dat %>%
        mutate(across(c(cam_fate, species), as.factor))
        dat$cam_fate <- relevel(dat$cam_fate, "H")
        dat$species  <- relevel(dat$species, "LETE")
    }
    rem_col <- na.omit(str_extract(string = names(dat), pattern = "cam_fate\\w+|species\\w+"))
    dat <- dat[,-which(names(dat) %in% rem_col)]
    return(dat)
}

###### MAKE SIMULATED DATA ###################################################################
########################################################################################

#mkResp <- function(sDat, betas, form, debug=FALSE, xdebug=FALSE){
#mkResp <- function(resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE){
#mkResp <- function(seed, resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE, mindebug=FALSE){
mkResp <- function(seed, resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, strat=TRUE, vbose=0){
    set.seed(seed=seed)
    success <- FALSE # success becomes true when dims match (allowing matrix multiplication)

    # 1. test whether dimensions match (for matrix multiplication) 
    while(!success){
        sDat <- mkSim(s_size, cMat, mList, betas, fprob, sprob, stratify=strat, vb=vbose)
        form <- as.formula(mod)
        dummySim <- model.matrix(form, data = sDat)
        success <- dim(dummySim)[2]==length(betas[[1]]) & dim(dummySim)[2]==length(betas[[2]])
    }
    # 2. pass betas to respMatMul for each response variable separately:
    sDat[,is_u := respMatMul("is_u",dummySim, betas[[1]], vb=vbose)]
    sDat[,HF_mis := respMatMul("HF_mis",dummySim, betas[[2]], vb=vbose)]
    # make response variables into factors:
    sDat[,is_u := as.factor(sDat[,is_u])]
    sDat[,HF_mis := as.factor(sDat[,HF_mis])]
    # 3. get the fit for the simulated data (to compare to true coefs):
    fitSim1 <- glm(as.formula(paste0("is_u",mod)),
                   family=binomial,
                   data=sDat,
                   method=brglm2::brglmFit)
    fitSim2 <- glm(as.formula(paste0("HF_mis",mod)),
                   family=binomial,
                   data=sDat,
                   method=brglm2::brglmFit)
    return(sDat)
}

respMatMul <- function(resp, dummySim, betas, vb=0){
     eta <- dummySim %*% betas 
     y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
     return(y)
}

mkSim <- function( s_size, cMat, mList, betas,fprob, sprob, stratify=TRUE, vb=0){
    fates <- names(fprob) # create variable for camera fates
    spp   <- names(sprob) # create variable for species names
    simDat <- list()
    if(stratify){
        for (s in seq_along(spp)){
            sp <- spp[s]
            simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=round(s_size*sprob[s]), mu=mList[[sp]], Sigma=cMat[[sp]]))
            simDat[[sp]][,cam_fate := sample(fates, size=nrow(simDat[[sp]]), replace=T, prob=fprob)]
        }
        names(simDat) <- c("CONI", "LETE")
        sDat <- data.table::rbindlist(simDat, idcol="species", use.names=T, fill=T)
    } 
    else {
        cMat <- readRDS("cmat_all.rds")
        mList <- readRDS("mlist_all.rds")
        sDat <- as.data.table(MASS::mvrnorm(n=s_size, mu=mList, Sigma=cMat))
        sDat[,species := sample(spp, size=nrow(sDat), replace=T, prob=sprob)]
        sDat[,cam_fate := sample(fates, size=nrow(sDat), replace=T, prob=fprob)]
    }
    sDat[,cam_fate := relevel(as.factor(sDat[,cam_fate]), ref="H")]
    sDat[,species  := relevel(as.factor(sDat[,species] ), ref="LETE")]
    id <- seq(1,s_size)
    sDat <- cbind(id, sDat)
    return(sDat)
}

##### AMPUTE SIMULATED DATA ###############################################################
########################################################################################

#mkSimDat <- function(seeed, nd, mpatt, wts, new_prop=0.2, patt_freq=c(0.45,0.45,0.1),wt=TRUE, test=FALSE, convFact=FALSE, facToNum=FALSE,vbose=0){
mkAmpDat <- function(seeed, nd, mpatt, wts, new_prop=0.2, patt_freq=c(0.45,0.45,0.1),wt=TRUE, test=FALSE, convFact=FALSE, facToNum=FALSE,vbose=0){
    dat4amp <- add_dummy(nd, vb=vbose)
    no_miss <- c("obs_int", "fdate", "is_u", "speciesCONI")
    is_miss <- colnames(mpatt)[!colnames(mpatt) %in% no_miss]
    new_order <- c(is_miss, no_miss)
    dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix
    suppressWarnings(amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = mpatt, freq = patt_freq,weights = wts))
    
    if(wt) {datList <- amp_out_wt} else {datList <- amp_out} ### *~*~*~*~* #######
    datList$amp
    if (convFact) datList$amp <- add_fact(dat = datList$amp, facToNum=facToNum, vb=vbose) # could probably reference the global debug instead...
    ampd <- datList$amp
    return(datList)
}


##### IMPUTE/FIT/POOL SIM DATA ##########################################################
########################################################################################
### Nest the loop inside the if statement so you aren't running the if check every loop?
# mkImpSim <- function(fullDat, ampDat,pr_list, resp_list, mod, vars, met, form_list, met_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
mkImpSim <- function(fullDat, aDat, resp_list, modd, vars, met, outFile,form_list, met_list, m=20, fam=binomial, regMet="brglm_fit", iter=500,impplot=FALSE, vbose=0){
    pr_list <-  c("species", "cam_fate", "obs_int", "nest_age", "fdate")
    ret    <- array(NA, dim=c(length(vars), 3, length(resp_list)))
    dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"), resp_list)
    if (met == "cc"){
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
            dat1 <- ampDat[complete.cases(ampDat),]
            fit = glm(as.formula(paste0(resp, modd)), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
            vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            rownames(vals) <- names(coef(fit))[-1]
            colnames(vals) <-  c("estimate", "2.5 %", "97.5 %") # this doesn't work if not a df
            
            #ret[vlist,,r]  <- as.matrix(exp(vals))# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
            ret[vlist,,r]  <- as.matrix(vals)# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
        }
        return(ret)
    } else {
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            cfates <- paste0("cam_fate", c("H", "A", "D", "F", "Hu", "S"))
            ampDat <- aDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?
            if (met=="cf_cc") ampDat <- ampDat  %>% filter(!is.na(cam_fate))
            metList <- met_list[resp,,met]
            inters <- sapply(modd,  function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))[[1]]
            inter <- paste(inters, collapse=".")
            if(met=="passive" & length(inters)!=0)  ampDat[inter] <- NA
            names(metList)[6] <- resp
            names(metList)[7] <- inter
            metList <- metList[!is.na(metList)]
            if(all(is.na(metList))) metList=NULL

            frmla <- lapply(form_list[[resp]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
            impCall <- case_when(met=="stratify"~ "impStrat(ampDat, met=metList,formulas=frmla, col_sel=cols)",
                       met=="passive" ~ "mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)",
                       .default="mice::mice(ampDat, method=metList, m=m, formulas=frmla,print=FALSE)"
            )
            imp <- eval(parse(text=impCall))
            impd <- sapply(imp$imp, function(x) any(!is.na(x)))

            if(length(imp$loggedEvents > 0)) { # ~*~*~*
                cat(sprintf("\n\n*** LOGGED EVENTS FOR MODEL %s & METHOD %s & RESP %s", modd, met, resp), file=outFile, append=TRUE)
                cat("\n\n it im dep \t meth \t\tout\n", file=outFile, append=TRUE)
                cat(paste(imp$loggedEvents, collapse=" "), file=outFile, append=TRUE)
                # cat(sprintf("\n\t*** LOGGED EVENTS FOR METHOD %s*******************************\n", met))
                # qvcalc::indentPrint(imp$loggedEvents, indent=6)
            }
            if(impplot) visImp(imp)  # ~*~*~*
            fit = with(imp,
                 glm( as.formula( paste0(resp, modd) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
                 ))
            pool = summary(mice::pool(fit), "all", conf.int=TRUE)
            rownames(pool) = pool[,"term"]
            pool <- pool %>% filter(term %in% vars)
            ret[vlist,,r] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
        }
    return(ret[vars,,])
    }
}


#printImp <- function(pool, sim, actual){
printImp <- function(impDat, simDat, trueDat){
    printImp <- array(NA, dim=c(length(vars),3, length(resp_list)))
    printImp[,1,] <- impDat[,1,]
    printImp[,2,1] <- simDat[,2,"is_u"]
    printImp[,2,2] <- simDat[,2,"HF_mis"]
    printImp[,3,] <- 
    #pool <- impDat[,1,]
    #printImp <- cbind(pool, sim, actual)
    names(printImp) <- c("POOL", "SIM", "ACTUAL")
    cat("\n\n oOoOoOo Comparison oOoOoOo\n")
    qvcalc::indentPrint(printImp, indent=6)
}

###### VISUALIZE IMPUTATIONS ###################################################################
########################################################################################

visImp <- function(imp){
  
    cat("\n\n ************* PLOTS *********************************** \n\n")
    # library(patchwork)
    # strPl <- mice::stripplot(imp, layout=c(3,1))
    # strPl <- mice::stripplot(imp)
    strPl <- mice::stripplot(imp, nest_age + cam_fate ~ .imp)
    # print(strPl)

    # denPl <- mice::densityplot(imp, layout=c(3,1))
    ### there's only one continuous predictor:
    denPl <- mice::densityplot(imp)
    # print(denPl)


    imp40 <- mice.mids(imp, maxit=35, print=F)
    # plot(imp40)
    pl40 <- mice:::plot.mids(imp40, layout=c(4,2))
    print(pl40)

    # (strPl | denPl) / pl40 # patchwork may not work with trellis plots?
    cplot <- cowplot::plot_grid(strPl, denPl, rel_widths=c(1.5, 1))
    print(cplot)


}



###### CALCULATE BIAS ###################################################################
########################################################################################

# calcBias <- function(impDat, mods4sim, vars, mets, rlist, vbose=0){
calcBias <- function(impDat, mod, vars, mets, rlist, vbose=0){
    cat("\n\t\t >-> CALCULATE BIAS VALUES <-<")
    bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(rlist) ) )
    dimnames(bias) <- list(sort(as.character(vars)),
                           # c("pmm", "rf", "cart"),
                           as.character(mets),
                           as.character(biasVals),
                           # names(mods4sim),
                           # "m1",
                           rlist)
    if(vbose>=2) cat("\n\t >>->> empty array to store bias values:" )
    if(vbose>=2) qvcalc::indentPrint(str(bias))
    if(vbose>=2) qvcalc::indentPrint(dimnames(bias))

    for(resp in rlist){
        if(vbose>=1) cat(sprintf("\n<> <> <> <> <> MODEL FORMULA: %s %s <> <> <> <> <> \n", resp, mod))
        fitReal <- glm(as.formula(paste0(resp, mod)),
                       data=trueDat,
                       family=fam,
                       method=regMet,
                       control=brglmControl(maxit=iter) )
        trueVals <- coef(fitReal)[vars] # the coefs have names associated with them

        ## Loop through the predictor variables and store the bias values to the matrix
            for(v in vars){
                simVals <- impDat[v,"sim", , "estimate", resp] # sim vals for this var & resp var, all runs
                avg <- apply(impDat[v, , , ,resp],MARGIN=c(1,3),FUN = mean, na.rm=TRUE) #print(impDat[v,,1:5,1,z])
                true <- trueVals[v] # true2 <- avg["true",] # simVal <- avg["sim",]
                avg    <- avg[mets,] # remove true and sim, which are not in "mets"
                sdev <- apply(impDat[v, mets, , ,resp],MARGIN=c(1,3),FUN = sd, na.rm=TRUE) # calculate std dev after removing 'sim' and 'true':
                true <- simVals
                bias[v, ,"value",resp] <- avg[,"estimate"]
                bias[v,,"SD",resp] <- sdev[,"estimate"]
                bias[v,mets,"bias",resp] <- rowMeans(impDat[v,mets,,"estimate",resp] - true, na.rm=params$naRm)
                bias[v,mets,"pctBias",resp] <- apply(impDat[v,mets,,"estimate",resp], 1, function(x) 100 *abs(mean((x-true)/true)))
                bias[v,mets,"covRate",resp] <- rowMeans(impDat[v,mets,,"2.5 %",resp] < true & true < impDat[v,mets,,"97.5 %",resp],na.rm=params$naRm)
                bias[v,mets,"avgWidth",resp] <- rowMeans(impDat[v,mets,,"97.5 %",resp] - impDat[v,mets,,"2.5 %",resp], na.rm=params$naRm)
                bias[v,mets,"RMSE",resp] <- 100 * sqrt(rowMeans((impDat[v,mets,,"estimate",resp] - true) ^2 ,na.rm=params$naRm))
               }

           # if(vbose >= 1) cat("\n>>> ALL BIAS VALUES: \n")
           cat("\n>>> ALL BIAS VALUES: \n")
           qvcalc::indentPrint(bias)
           # biasfile <-  sprintf("%s/bias_vals_%s_m%s_%s.rds",now_dir,resp, names(mods4sim)[z], suffix)
           biasfile <-  sprintf("%s/bias_vals_%s_m%s_%s.rds",now_dir,resp, 1, suffix)
           cat(sprintf("\n>>> saving to file: %s ", biasfile))
           saveRDS(bias, biasfile)
           biasfile1 <-  sprintf("%s/bias_vals_%s_%s_%s.csv",now_dir,resp, 1, suffix)
           if(vbose>=1) cat(sprintf("\t & to: %s \n", biasfile1))
           biasdf <- as.data.frame(bias)
           names(trueVals)[ is.na(names(trueVals)) ] <- setdiff(vars, names(trueVals))
           biasdf <- cbind(trueVals, biasdf)
           write.csv(biasdf, file = biasfile1)# write to csv in case script aborts

    }
    # }
}


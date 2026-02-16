########################################################################################
###### IMPORT FUNCTIONS ################################################################
########################################################################################

## NOTE: if you update this document, need to update the debug_imp_sim_func.R doc as well!

options(include.full.call.stack = FALSE)               # or configure it globally
# source("missingness_tab.R")
debugging <- FALSE # for uickly setting values when working in the file with the functions

########################################################################################
######  RANDOM FUNCTIONS ################################################################
########################################################################################

means <- function(realDat) sapply(realDat[,c(8:10)], mean)
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

########################################################################################
###### FACTORS --> DUMMIES #############################################################
########################################################################################

if(FALSE){
  dat <- dat4sim
}
add_dummy <- function(dat, vb=0){
#    dummy_vars <- model.matrix(~ cam_fate -1, data = dat)
    dummy_vars <- model.matrix(~ cam_fate -1, data = as.data.table(dat))
    dummy_vars <- dummy_vars[,-1]
    dummy_vars
    # dat$speciesLETE <- ifelse(dat$species=="LETE", 1, 0)
    dat$speciesCONI <- ifelse(dat$species=="CONI", 1, 0)

    # why did I get the columns in this way? I guess to reorder them?
    # dat<- dat[, c("obs_int", "nest_age", "fdate", "HF_mis", "is_u", "speciesLETE", "speciesCONI")]
    dat<- dat[, c( "speciesCONI", "obs_int","nest_age", "fdate", "HF_mis", "is_u")]
    # if LETE is the reference level, you have speciesCONI as the dummy variable
    # but here it seems like we keep all the dummy
    # dat<- dat[, c("obs_int", "nest_age", "fdate", "HF_mis", "is_u", "speciesCONI")]
    # dat4amp <- cbind(dat4amp, dummy_vars, speciesLETE, speciesCONI)
    dat <- cbind(dat, dummy_vars)
    #print(class(dat)) # ~*~*~*
    #if (debug) cat("\n\n>> add dummy variables and remove factors. new columns:\n", names(dat))
    return(dat)
}

########################################################################################
###### DUMMIES --> FACTORS #############################################################
########################################################################################

if(FALSE){
  dat <-amp_out_wt$amp
}
add_fact <- function(dat, facToNum=FALSE, vb=0){ 
    #cat("\n NOTE: can speed up conversion back to factor with tidytable\n")
    # I already know what the dummy var names are in this case; not a general purpose function
    # dat <- amp_out$amp
    # if(is.list(dat)) dat <- dat$amp
    # if(any(dim(dat))) dat <- dat$amp
    # dim(p) >0

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

    #cat("\n\t>> check new factor variable:\n\t", table(dat$cam_fate, useNA="ifany"))

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
    } else{
    dat <- dat %>%
    mutate(across(c(cam_fate, species), as.factor))
    dat$cam_fate <- relevel(dat$cam_fate, "H")
    dat$species  <- relevel(dat$species, "LETE")
    }

    # rem_col <- c(3:8, 12,13)
    rem_col <- na.omit(str_extract(string = names(dat), pattern = "cam_fate\\w+|species\\w+"))
    # names(dat)[ c(3:8, 12,13)]
    # names(dat)
    #if (debug) cat("\n>> converted dummies to factors. all columns:\n", names(dat)) 
    # if (debug) cat("\n>> columns to remove:\n", names(dat)[rem_col])
    # if (debug) cat("\n>> columns to remove:\n", rem_col)

    # dat <- 
    # dat[-c("species")]
    # dat[rem_col]
    # dat <- dat[,-rem_col]
    dat <- dat[,-which(names(dat) %in% rem_col)]
    #return(as.data.table(dat))
    return(dat)
    # there MUST be a shorter way to do this, using regex?
}

########################################################################################
##### MAKE SIMULATED DATA ###############################################################
########################################################################################

mkSimDat <- function(seeed, nd, mpatt, wts, new_prop=0.2, patt_freq=c(0.45,0.45,0.1),wt=TRUE, test=FALSE, convFact=FALSE, facToNum=FALSE,vbose=0){
  #cat("mkSimDat seed=", seeed, class(seeed))
    dat4amp <- add_dummy(nd, vb=0)
    no_miss <- c("obs_int", "fdate", "is_u", "speciesCONI")
    is_miss <- colnames(mpatt)[!colnames(mpatt) %in% no_miss]
    new_order <- c(is_miss, no_miss)
    dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix
    suppressWarnings(amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = mpatt, freq = patt_freq,weights = wts))
    #if(debug) cat("\n\nCreate more new missing values, with weighted probabilities:\n")
    # if(debug)  print(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE))
    # missing_tab("amp_out_wt", prVars)
    
    if(wt) {datList <- amp_out_wt} else {datList <- amp_out}
    ### *~*~*~*~* #######
    datList$amp
    if (convFact) datList$amp <- add_fact(dat = datList$amp, facToNum=facToNum, vb=vbose) # could probably reference the global debug instead...
    return(datList)
}

########################################################################################
##### IMPUTE/FIT/POOL SIM DATA ##########################################################
########################################################################################
if(debugging){
    met_list <- metLists[,,,mod]
    form_list <- formulas[[names(mods4sim)[mod]]]
    dat4sim <- mkSim(resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, debug=params$deb)
    convFact <- TRUE
    datNA <- mkSimDat( seeed = run+seed, nd = dat4sim, mpatt=mpatt, wts=ampwt, xdebug=params$xdeb, debug = params$deb, convFact=convFact)
    aDat <- datNA$amp
    modd <- mods4sim[2]
    fam=binomial
    regMet="brglm_fit"
    iter=500
    r <- 1
    m <- 5
    met <- "default"
}
### Nest the loop inside the if statement so you aren't running the if check every loop?
# mkImpSim <- function(fullDat, ampDat,pr_list, resp_list, mod, vars, met, form_list, met_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
mkImpSim <- function(fullDat, aDat, resp_list, modd, vars, met, outFile,form_list, met_list, m=20, fam=binomial, regMet="brglm_fit", iter=500,impplot=FALSE, vbose=0){
    pr_list <-  c("species", "cam_fate", "obs_int", "nest_age", "fdate")
    ret    <- array(NA, dim=c(length(vars), 3, length(resp_list)))
    dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"), resp_list)
    # cat("\n\n >>>> empty matrix to store output:\n") # ~*~*~*
    # print(ret)
    # dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %", "fmi"), names(mods))
    #dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"))
    if (met == "cc"){
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            #if(debug)cat("\nvars in this df:", vlist,"\t& cols for imputation:", cols) # ~*~*~*
            ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
            dat1 <- ampDat[complete.cases(ampDat),]
            fit = glm(as.formula(paste0(resp, modd)), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
            vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            rownames(vals) <- names(coef(fit))[-1]
            colnames(vals) <-  c("estimate", "2.5 %", "97.5 %") # this doesn't work if not a df
            
            ret[vlist,,r]  <- as.matrix(vals)# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
            #if(debug) cat("\nret, FILLED IN :") # ~*~*~*
            #if(debug) print(ret[,,r])
        }
        return(ret)

    } else {
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            #if (debug) cat("\n===== method:", met, "\t resp:", resp," =====\n")
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            #if(debug)cat("\nvars in this df:", vlist,"\t& cols for imputation:", cols) # ~*~*~*
            cfates <- paste0("cam_fate", c("H", "A", "D", "F", "Hu", "S"))
            # collapse pkg might be faster - https://stackoverflow.com/questions/77679416/drop-rows-with-na-in-multiple-columns
            ampDat <- aDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?
            if (met=="cf_cc") ampDat <- ampDat  %>% filter(!is.na(cam_fate))
            # if (met=="cf_cc") ampDat <- ampDat %>% filter(if_any(cfates, ~ !is.na(.)))
              #for(y in seq_along(mods)){
                # needs to be a named list; but if you keep the names, it doesn't drop NA...
            metList <- met_list[resp,,met]
            inters <- sapply(modd,  function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))[[1]]
            inter <- paste(inters, collapse=".")
            #cat("\n\n####### inter=", inter,"\n") # ~*~*~*
            if(met=="passive" & length(inters)!=0)  ampDat[inter] <- NA
            names(metList)[6] <- resp
            names(metList)[7] <- inter
            metList <- metList[!is.na(metList)]
            if(all(is.na(metList))) metList=NULL

            ## *~*~*~*~*
            #cat("\nMETHOD LIST for",met,"-") # ~*~*~*
            #print(metList)
            #cat("\n")

            frmla <- lapply(form_list[[resp]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
                ## *~*~*~*~*
            #if(debug) cat("\nmice formulas:\n") # ~*~*~*
            #if(debug) print(frmla)
            #if(debug) cat(paste(frmla, collapse="\t\t"))

            impCall <- case_when(met=="stratify"~ "impStrat(ampDat, met=metList,formulas=frmla, col_sel=cols)",
                       met=="passive" ~ "mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)",
                       #.default="mice::mice(ampDat, method=metList, m=m, print=FALSE)"
                       ### turn off print=FALSE to get output from imputations as they happen
                       #.default="mice::mice(ampDat, method=metList, m=m,print=FALSE)"
                       .default="mice::mice(ampDat, method=metList, m=m, formulas=frmla,print=FALSE)"
            )
            imp <- eval(parse(text=impCall))

            if(length(imp$loggedEvents > 0)) { # ~*~*~*
                #cat("\n")
                cat(sprintf("\n\n*** LOGGED EVENTS FOR MODEL %s & METHOD %s & RESP %sn", modd, met, resp), file=outFile, append=TRUE)
                cat("\n\n it im dep \t meth \t\tout\n", file=outFile, append=TRUE)
                cat(paste(imp$loggedEvents, collapse=" "), file=outFile, append=TRUE)
                #print(imp$loggedEvents)
            }
            # if(impplot) visImp(imp)  # ~*~*~*
            #if(debug) cat("\n>> fitting model", mods[y])
            fit = with(imp,
                 glm( as.formula( paste0(resp, modd) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
                 ))
            #if(debug) cat("\n>> pooling model fit\n") # ~*~*~*
            pool = summary(mice::pool(fit), "all", conf.int=TRUE, exponentiate=TRUE)
            # print(pool)
            rownames(pool) = pool[,"term"]
            pool <- pool %>% filter(term %in% vars)
            #cat("\nOUTPUT:\n") # ~*~*~*
            #print(pool[,c("estimate", "2.5 %", "97.5 %")])
            # cat("\nPUT IT HERE:")
            # print(ret[vlist,,r])
            ret[vlist,,r] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
            #if(debug) cat("\nret, FILLED IN:\n")
            #if(debug) print(ret[,,r])
        }
    # if(xdebug) cat("\n\nret:\n") # ~*~*~*
    #    cat("\n\nret:\n")
    #    print(head(ret))
    # if (xdebug) print(str(ret))
    return(ret[vars,,])
    }
}

########################################################################################
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

########################################################################################
###### MAKE SIMULATED DATA ###################################################################
########################################################################################



#mkResp <- function(sDat, betas, form, debug=FALSE, xdebug=FALSE){
#mkResp <- function(resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE){
#mkResp <- function(seed, resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE, mindebug=FALSE){
mkResp <- function(seed, resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, vbose=0){
    #dummySim <- as.data.table(model.matrix(form, data = sDat))
    #tryCatch(
    
    set.seed(seed=seed)
    success <- FALSE
    # 1. test whether dimensions match (for matrix multiplication) 
    while(!success){
        #sDat <- mkSim(s_size, cMat, mList, betas, fprob, sprob, debug=debug, xdebug=xdebug)
        sDat <- mkSim(s_size, cMat, mList, betas, fprob, sprob, vb=vbose)
        form <- as.formula(mod)
        dummySim <- model.matrix(form, data = sDat)
        #success <- dim(dummySim)[2]==dim(betas)[1]
        success <- dim(dummySim)[2]==length(betas[[1]]) & dim(dummySim)[2]==length(betas[[2]])
        #cat("success?", success)
    }
    # 2. pass betas to respMatMul for each response variable separately:
    sDat[,is_u := respMatMul(dummySim, betas[[1]])]
    sDat[,HF_mis := respMatMul(dummySim, betas[[2]])]
    # make response variables into factors:
    sDat[,is_u := as.factor(sDat[,is_u])]
    sDat[,HF_mis := as.factor(sDat[,HF_mis])]
    #cat("\n>> created response vars:")
    #print(sDat[,HF_mis])
    #print(sDat[,is_u])
    return(sDat)
}

respMatMul <- function( dummySim, betas, vb=0){
    #if(vb>=1) cat("\n***MATRIX MULTIPLICATION***", dim(dummySim), dim(betas))
    eta <- dummySim %*% betas 
    y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
    return(y)
}

#mkSim <- function( s_size, cMat, mList, betas,fprob, sprob, stratify=TRUE, debug=FALSE, xdebug=FALSE, mindebug=FALSE){
mkSim <- function( s_size, cMat, mList, betas,fprob, sprob, stratify=TRUE, vb=0){
    fates <- names(fprob)
    spp   <- names(sprob)
    #if(debug) cat("\nfates:", fates, "& species:", spp) # ~*~*~*
    # if(debug){ # ~*~*~*
    #     cat("\ncoefficients, class:\n", class(betas))
    #     print(betas)
    # }
    simDat <- list()
    if(stratify){
        for (s in seq_along(spp)){
            sp <- spp[s]
            # cat("\n>>> species:", sp, "\t& means:", mList[[sp]], "\t& correlation matrix:\n") # ~*~*~*
            # print(cMat[[sp]])
            simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=round(s_size*sprob[s]), mu=mList[[sp]], Sigma=cMat[[sp]]))
            #cat("\n>> adding camera fates\n") # ~*~*~*
            simDat[[sp]][,cam_fate := sample(fates, size=nrow(simDat[[sp]]), replace=T, prob=fprob)]
        }
        #cat("\n>>> created sim data") # ~*~*~*
        names(simDat) <- c("CONI", "LETE")
        sDat <- data.table::rbindlist(simDat, idcol="species", use.names=T, fill=T)
        #cat(" - merged sim data for the two species") # ~*~*~*
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
    return(sDat)
    #cat(" - releved factors") # ~*~*~*
    #form <- as.formula(mod)
    # if(xdebug) cat("\n>>>>> simulated data frame w/o dummy variables:\n") # ~*~*~*
    # if(xdebug) print(sDat)

    #dummySim <- model.matrix(form, data = sDat)
    #if(dim(dummySim)[2]==dim(betas)[1]) {
    #    sDat[,is_u := mkResp(sDat, (betas[[1]]), form, debug=debug )]
    #    sDat[,HF_mis := mkResp(sDat, (betas[[2]]), form, debug=debug )]
    #}else{
    #}

    #cat(' - added response variables\n') # ~*~*~*
    #for (r in seq_along(resp_list)){
    #    #resp <- mkResp(sDat, betas, form, debug=debug)
    #    sDat[,resp_list[r] := mkResp(sDat, betas[[r]], form, debug=debug)]
    #}

    #dummySim <- model.matrix(form, data = sDat)
    #if(debug) cat("\ndesign matrix, class:\n", str(dummySim))
    #if(debug) print(dummySim)
    #eta <- dummySim %*% betas
    #if(debug) cat("\nmultiplied by beta values:\n")
    #print(eta)
    #y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
    #sDat[,is_u := y]
    #sDat[,is_u := resp]
    # if(xdebug){ # ~*~*~*
    #     fitReal2 <- summary(glm(is_u ~ nest_age * species + obs_int + cam_fate + fdate,
    #                        family=binomial,
    #                        data=sDat,
    #                        method=brglm2::brglmFit))
    #     cat("\n>>>> model fit for simulated data:\n")
    #     print(fitReal2)
    # }
    #return(sDat)
}

########################################################################################
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



########################################################################################
###### CALCULATE BIAS ###################################################################
########################################################################################

calcBias <- function(impDat, mods4sim, vars, mets, rlist){
    biasVals <- c("value","bias","pctBias","covRate","avgWidth")
    #bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(mods4sim), length(rlist) ) )
    #dimnames(bias) <- list(sort(as.character(vars)),
   #                        # c("pmm", "rf", "cart"),
    #                       as.character(mets),
    #                       as.character(biasVals),
    #                       names(mods4sim),
    #                       rlist)
    #if(params$debug) cat("\n    empty array to store bias values:" )
    #if(params$debug) qvcalc::indentPrint(str(bias))
    #if(params$debug) qvcalc::indentPrint(dimnames(bias))


    for(z in seq_along(mods4sim)){
        cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
        cat("\n\t\t\t MODEL:", mods4sim[z])
        cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
        mod     <- mods4sim[z]
        cat("\n\t\t\t DATA:")
        qvcalc::indentPrint(str(impDat))
        cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
            #res <- array(NA, dim = c(length(vars), params$nrun,
            #             3, length(mods4sim), length(resp_list)))
            #nrun <- nruns[]
            #cat("nrun for this model:", nrun, "correct number:", nrun==params$nrun)
            #allSeeds <- seq(seed,seed+params$nrun)
            #cat("all seeds used for simulating data:", allSeeds)

            # no mod nums in coeefe files name
            #seedMods <- str_extract(seedCoef, "(?<=mod)\\d+")
            #readFile <- seedCoef[seedMods==mod]
            #cat("\n>>> attempt to merge:\n")
            #impDat <- abind::abind(flist[mods==mod], along=3)
            ## double bracket doesn't work'
            #impDat <- abind::abind(flist[[mods==mod]], along=3)
            ## Separate by response variable and get the real model output
        for(resp in rlist){
            cat(sprintf("\n<> <> <> <> <> MODEL: %s %s <> <> <> <> <> \n", resp, mod))
                #sim_val <- simVals[,,,z,resp]
                #sim_val <- impDat[,"sim",]
                #cat("\n    >> sim values for this model:\n")
                #qvcalc::indentPrint(sim_val, indent=8)
                #fitReal <- glm(as.formula(paste0(resp, mod)),
                #               data=fullDat,
                #               family=fam,
                #               method=regMet,
                #               control=brglmControl(maxit=iter) )
                # don't want coef - want coefficients from summary.glm - need exp(coef())?
                #trueVals <- exp(coef(fitReal)[vars]) # the coefs have names associated with them
                #if(params$debug) qvcalc::indentPrint(summary(fitReal)) 
                #cat("\n    true values:\n")
                #qvcalc::indentPrint(trueVals)
            bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(mods4sim), length(rlist) ) )
            dimnames(bias) <- list(sort(as.character(vars)),
                                   # c("pmm", "rf", "cart"),
                                   as.character(mets),
                                   as.character(biasVals),
                                   names(mods4sim),
                                   rlist)
            if(params$debug) cat("\n    empty array to store bias values:" )
            if(params$debug) qvcalc::indentPrint(str(bias))
            if(params$debug) qvcalc::indentPrint(dimnames(bias))

            #for(r in seq(1, params$nrun)){
            #    sim <- exp(sim_val[,r,1]) # 1 = estimate
            #    cat(sprintf("\n    sim vals (estimate) for run %s (exp):\n", r))
            #    qvcalc::indentPrint(sim)

            ## Loop through the predictor variables and store the bias values to the matrix
                for(v in vars){
                    cat(sprintf("\n    ----- VARIABLE: %s ------------------------------------------------\n", v))
                    if(params$debug) cat(sprintf("\n    >> output for variable %s:\n", v))
                    # now I'm doing this in the other script (debug_whatever.R)
                    #impDat[v,"cc",,,z,resp] <- exp(impDat[v,"cc",,,z,resp]) 
                    if(params$debug) qvcalc::indentPrint(str(impDat[v,,r,,z,resp]), indent=8)
                    #avg <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
                    avg <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
                    sdev <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = sd, na.rm=TRUE)
                    #print(impDat[v,,1:5,1,z])
                    #true <- trueVals[v]
                    #true <- exp(sim_val[v,,"Estimate"])
                    true <- impDat[v,"sim",,]
                    cat("\n    vals for comparison (exp)=", true,"\n    & dimensions:", dim(true), dim(impDat[v,,,"estimate",z,resp]))
                    if(params$debug){
                          cat(sprintf("\n    average of estimate for mod %s:\n",z))
                          qvcalc::indentPrint(apply(impDat[v,,,,z,resp], FUN=mean, MARGIN=c(1,3)))
                          cat("\n    dimensions:\n")
                          qvcalc::indentPrint(dim(impDat[v,,,,z,resp]))
                          cat("\n    avg:\n")
                          qvcalc::indentPrint(avg)
                    }
                    cat("\n    >>> storing vals to matrix\n")
                    bias[v,,"value",z,resp] <- avg[,"estimate"]
                    #bias[v,,"bias",z,resp] <- avg[,"estimate"] - true
                    bias[v,,"bias",z,resp] <- rowMeans(impDat[v,,,"estimate",z,resp] - true)
                    #bias[v,,"bias",z,resp] <- rowMeans(impDat[v,,,"estimate",z,resp] - true, dims=2)
                    #bias[v,,"pctBias",z,resp] <- 100 * abs((avg[,"estimate"] - true) / true )
                    bias[v,,"pctBias",z,resp] <- 100 * abs(rowMeans((impDat[v,,,"estimate",z,resp] - true) / true ))
                    bias[v,,"covRate",z,resp] <- rowMeans(impDat[v,,,"2.5 %",z,resp] < true & true < impDat[v,,,"97.5 %",z,resp])
                    bias[v,,"avgWidth",z,resp] <- rowMeans(impDat[v,,,"97.5 %",z,resp] - impDat[v,,,"2.5 %",z,resp])
                    #bias[v,,"RMSE",z,resp] <- sqrt((avg[,"estimate"] - true)^2)
                    #bias[v,,"RMSE",z,resp] <- 100 * sqrt((rowMeans(impDat[v,,,"estimate",]) - true) ^2 )
                    bias[v,,"RMSE",z,resp] <- 100 * sqrt(rowMeans(impDat[v,,,"estimate",z,resp] - true) ^2 )
                    bias[v,,"SD",z,resp] <- sdev[,"estimate"]
                    #if(parrams$debug) cat(sprintf("\nbias values for %s and model %s %s:\n\n", v, resp, mods4sim[z]))
                    cat(sprintf("\n>>> bias values for model %s %s for variable %s\n",
                              resp, mods4sim[z], v))
                    qvcalc::indentPrint(bias[v,,,z,resp], indent=8)
                }

            }
            biasfile <-  sprintf("%s/bias_vals_%s_m%s_%s.rds",now_dir,resp, names(mods4sim)[z], suffix)
            cat(sprintf("\n>>> saving to file: %s ", biasfile))
            saveRDS(bias, biasfile)
            biasfile1 <-  sprintf("%s/bias_vals_%s_%s_%s.csv",now_dir,resp, names(mods4sim)[z], suffix)
            cat(sprintf("\t & to: %s \n", biasfile1))
            biasdf <- as.data.frame(bias)
            names(trueVals)[ is.na(names(trueVals)) ] <- setdiff(vars, names(trueVals))
            biasdf <- cbind(trueVals, biasdf)
            write.csv(biasdf, file = biasfile1)# write to csv in case script aborts

        #}
    }

}


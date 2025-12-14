



########################################################################################
###### IMPORT FUNCTIONS ################################################################
########################################################################################
options(include.full.call.stack = FALSE)               # or configure it globally
# source("missingness_tab.R")
debugging <- FALSE # for uickly setting values when working in the file with the functions

########################################################################################
###### FACTORS --> DUMMIES #############################################################
########################################################################################

if(FALSE){
  dat <- dat4sim
}
add_dummy <- function(dat, debug=FALSE){
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
    if (debug) cat("\n\t>> add dummy variables and remove factors. new columns:\n\t\t", names(dat))
    return(dat)
}

########################################################################################
###### DUMMIES --> FACTORS #############################################################
########################################################################################

if(FALSE){
  dat <-amp_out_wt$amp
}
add_fact <- function(dat, facToNum=FALSE, debug=FALSE){ 
    #cat("\n NOTE: can speed up conversion back to factor with tidytable\n")
    # I already know what the dummy var names are in this case; not a general purpose function
    #cat("\n\t>> cam fates before conversion:\n\t", colSums(dat))
    #if(debug) cat("\n\t>> get positions of NAs:\n\t")
    #qvcalc::indentPrint(head(dat[!complete.cases(dat)]))
    #naPos <- which(is.na(dat$))
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

    if(debug) cat("\n\t>> check new factor variable:\n\t", table(dat$cam_fate, useNA="ifany"))

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
    if (debug) cat("\n\t>> converted dummies to factors. all columns:\n\t\t", names(dat)) 
    #if (debug) cat("\n>> columns to remove:\n", names(dat)[rem_col])
    if (debug) cat("\n\t>> columns to remove:\n\t\t", rem_col)

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
##### AMPUTE SIMULATED DATA ###############################################################
########################################################################################

mkSimDat <- function(seeed, nd, mpatt, wts, new_prop=0.2, patt_freq=c(0.45,0.45,0.1),wt=TRUE, test=FALSE, debug=FALSE, xdebug=FALSE, convFact=FALSE, facToNum=FALSE){
    if(debug) cat("\n<><><><><><><><><><><><><><><> MAKE AMPUTED DATA: <><><><><><><><><><><><><><><>  ")
    if(debug) cat("\n\t>> mkSimDat seed=", seeed, class(seeed))
  # if(method=="amp"){
    dat4amp <- add_dummy(nd, debug=debug)
    set.seed(seed=seeed)
    # no_miss <- c("obs_int", "fdate", "is_u", "speciesLETE", "speciesCONI")
    no_miss <- c("obs_int", "fdate", "is_u", "speciesCONI")
    is_miss <- colnames(mpatt)[!colnames(mpatt) %in% no_miss]
    new_order <- c(is_miss, no_miss)
    ### *~*~*~*~* #######
    if(debug) cat("\n\t>> reorder columns:\n\t\t", new_order)
    dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix
    
    suppressWarnings(amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = mpatt, freq = patt_freq,weights = wts))
    
    if(debug) cat("\n\t>> Create more new missing values, with weighted probabilities:\n")
    if(xdebug)  qvcalc::indentPrint(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE), indent=8)
    #if(test)  (mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE), indent=8) #figuree out how to save to file
    # missing_tab("amp_out_wt", prVars)
    
    if(wt) {datList <- amp_out_wt} else {datList <- amp_out}
    ### *~*~*~*~* #######
    datList$amp
    if(xdebug) cat("\n\t>> the amputed data before converting to factors:\n")
    if(xdebug) qvcalc::indentPrint(str(datList$amp), indent=8)
    if(xdebug) qvcalc::indentPrint(colSums(is.na(datList$amp)), indent=8)
    if(xdebug) qvcalc::indentPrint(names(datList$amp), indent=8)
    # if (debug) print(class(datList$amp))
    if (convFact) datList$amp <- add_fact(dat = datList$amp, facToNum=facToNum, debug=debug) # could probably reference the global debug instead...
    if(debug) cat("\n\n<><><><><><><><><><<><><><><<><><><><> AMPUTED DATA: <><><><><><><><><><<><><><><<><><><><>>>\n ")
    if(debug) qvcalc::indentPrint(str(datList$amp), indent=8)
    if(debug) qvcalc::indentPrint(colSums(is.na(datList$amp)), indent=8)
    ampd <- datList$amp
    if (test) vmiss <- naniar::vis_miss(ampd)
    if (test) now <- format(Sys.time(), "%d%b%H%M%S")
    if (test) fname <- sprintf("figs/ampDat_%s.svg",now)
    if (test) ggsave(fname, vmiss, device="svg")

    #if(params$deb) cat("\n>> species for this run:")
    #if(debug) qvcalc::indentPrint(table(datList$amp$speciesCONI), indent=8)
   # if(debug) qvcalc::indentPrint(table(ampd$species), indent=8)
    #if(debug) cat("\n>> cam fate vars for this run:")
    #if(debug) qvcalc::indentPrint(table(datList$amp$cam_fate), indent=8)
    #if(debug) qvcalc::indentPrint(table(ampd$cam_fate), indent=8)
    if(debug) cat("\n\t data summary:\n")
    qvcalc::indentPrint(summary(datList$amp), indent=8)
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
mkImpSim <- function(fullDat, aDat, resp_list, modd, vars, met, outFile,form_list, met_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
    #cat("\n\t>>>>>> prepare for imputation")# ~*~*~
    #print(str(ampDat))
    #if(xdebug) cat("- amputed data:")
    #cat("\n")
    #if(xdebug) bprint(aDat)
    #if(xdebug) qvcalc::indentPrint(str(aDat))
    #if(xdebug) qvcalc::indentPrint(rowSums(is.na(aDat)))
    pr_list <-  c("species", "cam_fate", "obs_int", "nest_age", "fdate")
    ret    <- array(NA, dim=c(length(vars), 3, length(resp_list)))
    #dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %","n"), resp_list)
    dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"), resp_list)
    #cat("\n\t>>>> make empty matrix to store output") # ~*~*~*
    if (met == "cc"){
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            if (debug) cat("\n===========================================================================================")
            if(debug) cat("\n==================== complete-case analysis \t resp:", resp," =====================")# ~*~*~*
            if (debug) cat("\n===========================================================================================")
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            if(debug)cat("\n\t>> vars in this df:", vlist,"\n\t\t& cols for analysis:", cols) # ~*~*~*
            ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
            dat1 <- ampDat[complete.cases(ampDat),]
            fit = glm(as.formula(paste0(resp, modd)), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
            vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            if(xdebug) cat("\n\t>> output vals:\n") # ~*~*~*
            if(xdebug) qvcalc::indentPrint(vals, indent=8)# ~*~*~
            rownames(vals) <- names(coef(fit))[-1]
            colnames(vals) <-  c("estimate", "2.5 %", "97.5 %") # this doesn't work if not a df
            
            ret[vlist,,r]  <- as.matrix(vals)# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
            if(xdebug) cat("\n\tret, FILLED IN:\n") # ~*~*~*
            if(xdebug) qvcalc::indentPrint(ret[,,r], indent=8)
        }
        return(ret)

    } else {
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            if (debug) cat("\n===========================================================================================")
            if (debug) cat("\n======================= method:", met, "\t resp:", resp," =============================== ")
            if (debug) cat("\n===========================================================================================")
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            if(debug)cat("\n\t>> vars in this df:", vlist,"\n\t\t& cols for imputation:", cols) # ~*~*~*
            cfates <- paste0("cam_fate", c("H", "A", "D", "F", "Hu", "S"))
            # collapse pkg might be faster - https://stackoverflow.com/questions/77679416/drop-rows-with-na-in-multiple-columns
            if(debug) cat("\n\t>>>> vars with empty levels:", names(which(sapply(aDat, function(x) length(setdiff(levels(x),unique(x))))>0)))
            if(debug) cat("\n\t>>>> levels of factors:", levels(aDat$cam_fate),"\t", levels(aDat$HF_mis),"\t", levels(aDat$species))
            ampDat <- aDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?
            if(debug) cat("\n\t>>>> levels of factors:", levels(aDat$cam_fate),"\t", levels(aDat$HF_mis),"\t", levels(aDat$species))
            if (met=="cf_cc") ampDat <- ampDat  %>% filter(!is.na(cam_fate))
            if(debug) cat("\n\t>>>> vars that are missing values:", names(which(colSums(is.na(ampDat))>0)))
            metList <- met_list[resp,,met]
            inters <- sapply(modd,  function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))[[1]]
            inter <- paste(inters, collapse=".")
            if(xdebug) cat("\n\t\t >>>>> inter=", inter) # ~*~*~*
            if(met=="passive" & length(inters)!=0)  ampDat[inter] <- NA
            names(metList)[6] <- resp
            names(metList)[7] <- inter
            metList <- metList[!is.na(metList)]
            if(all(is.na(metList))) metList=NULL

            if(debug) cat("\n<><><><><><><><><><><><><><><> imputation with MICE <><><><><><><><><><><><><><><>  ")
            if(xdebug){ ## *~*~*~*~*
                cat("\n\tMETHOD LIST for",met,":\n") # ~*~*~*
                qvcalc::indentPrint(metList, indent=8)
                #cat("\n")
            }

            frmla <- lapply(form_list[[resp]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
            if(xdebug) cat("\n\tmice formulas:\n") # ~*~*~*
            if(xdebug) cat("\t\t",paste(frmla, collapse="\n\t\t"))

            impCall <- case_when(met=="stratify"~ "impStrat(ampDat, met=metList,formulas=frmla, col_sel=cols)",
                       #met=="passive" ~ "mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)",
                       met=="passive" ~ "mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)",
                       #.default="mice::mice(ampDat, method=metList, m=m, print=FALSE)"
                       ### turn off print=FALSE to get output from imputations as they happen
                       #.default="mice::mice(ampDat, method=metList, m=m, formulas=frmla,print=FALSE)"
                       #.default="mice::mice(ampDat, method=metList, m=m,print=FALSE)"
                       .default="mice::mice(ampDat, method=metList, m=m, formulas=frmla,print=FALSE)"
            )
            if(xdebug) cat("\n\n\t>> call:", impCall,"\n")
            imp <- eval(parse(text=impCall))
            if(debug) cat("\n<><><><><><><><><><><><><><><> MICE results <><><><><><><><><><><><><><><> ")
            if(xdebug) cat("\n\t>> the imputed values:\n")
            #if(xdebug) qvcalc::indentPrint(str(imp$imp))
            #impd <- sapply(imp$imp, function(x) any(rowSums(x)>0))
            #impd <- sapply(imp$imp, function(x) any(x>0))
            impd <- sapply(imp$imp, function(x) any(!is.na(x)))
            if(xdebug) qvcalc::indentPrint(impd, indent=8)
            if(debug) cat(">>>>> vars that were imputed:", names(impd)[impd])

            if(length(imp$loggedEvents > 0)) { # ~*~*~*
                cat(sprintf("\n\n*** LOGGED EVENTS FOR MODEL %s & METHOD %s & RESP %s", modd, met, resp), file=outFile, append=TRUE)
                cat("\n\n it im dep \t meth \t\tout\n", file=outFile, append=TRUE)
                cat(paste(imp$loggedEvents, collapse=" "), file=outFile, append=TRUE)
                cat(sprintf("\n\t*** LOGGED EVENTS FOR METHOD %s*******************************\n", met))
                qvcalc::indentPrint(imp$loggedEvents, indent=8)
            }
            if(impplot) visImp(imp)  # ~*~*~*
            if(debug) cat("\n<><><><><><><><><><><><><><><> Analysis <><><><><><><><><><><><><><><> ")
            if(debug) cat("\n\t>> fitting model", modd)
            fit = with(imp,
                 glm( as.formula( paste0(resp, modd) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
                 ))
            if(debug) cat("\n\t>> pooling model fit") # ~*~*~*
            pool = summary(mice::pool(fit), "all", conf.int=TRUE, exponentiate=TRUE)
            rownames(pool) = pool[,"term"]
            pool <- pool %>% filter(term %in% vars)
            if(xdebug){
                cat("\n\tpool() OUTPUT:\n") # ~*~*~*
                qvcalc::indentPrint(pool[,c("estimate", "2.5 %", "97.5 %")], indent=8)
                cat("\n\tPUT IT HERE:\n")
                qvcalc::indentPrint(ret[vlist,,r], indent=8)
            }
            ret[vlist,,r] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
            if(xdebug) cat("\n\tret, FILLED IN:\n")
            if(xdebug) qvcalc::indentPrint(ret[,,r], indent=8)
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
    qvcalc::indentPrint(printImp, indent=8)
}

########################################################################################
###### MAKE SIMULATED DATA ###################################################################
########################################################################################

#mkResp <- function(sDat, betas, form, debug=FALSE, xdebug=FALSE){
#mkResp <- function(resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE){
mkResp <- function(seed, resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE){
    #dummySim <- as.data.table(model.matrix(form, data = sDat))
    #tryCatch(
    set.seed(seed=seed)
    success <- FALSE
    while(!success){
        sDat <- mkSim(s_size, cMat, mList, betas, fprob, sprob, debug=debug, xdebug=xdebug)
        form <- as.formula(mod)
        if(debug) cat("\n\t>>>> making response variables with model ", mod)
        dummySim <- model.matrix(form, data = sDat)
        if(debug) cat("\n\t>> do dims match?", dim(dummySim), length(betas[[1]]), length(betas[[2]]))
        success <- dim(dummySim)[2]==length(betas[[1]]) & dim(dummySim)[2]==length(betas[[2]])
        if (debug) cat("\n\t>>> ***success?***", success)
    }
    sDat[,is_u := respMatMul(dummySim, betas[[1]])]
    sDat[,HF_mis := respMatMul(dummySim, betas[[2]])]
    sDat[,is_u := as.factor(sDat[,is_u])]
    sDat[,HF_mis := as.factor(sDat[,HF_mis])]
    #cat("\n>> added response variables, which should be factors:", sapply(sDat[,c(HF_mis, is_u)], is.factor))
    # DOESN'T work with data.table:
    #cat("\n>> added response variables, which should be factors:", sapply(as.data.frame(sDat)[c(HF_mis, is_u)], is.factor))
    #fitReal2 <- summary(glm(is_u ~ species * nest_age + obs_int + cam_fate + fdate,
    fits <- readRDS("fits.rds")
    #modnum <- as.numericnames(mod)
    modnum<- gsub("\\D+", "", names(mod))
    #mnums <- as.numeric(gsub("\\D+", "", names(fits)))
    #idx <- gsub("\\D+", "", names(fits))
    mnums <- rep(c("1", "8", "16"))
    if(debug) cat("\n\tmodel number:", modnum, "and model nums from list:", paste(mnums, collapse=","))
    if(debug) cat("\n\t>> Actual model fit:\n")
    #qvcalc::indentPrint(sapply(fits[[mnums==names(m)]], summary))
    # it's already a summary
    qvcalc::indentPrint(fits[mnums==modnum], indent=8)
    fitSim1 <- summary(glm(as.formula(paste0("is_u",mod)),
                       family=binomial,
                       data=sDat,
                       method=brglm2::brglmFit))
    fitSim2 <- summary(glm(as.formula(paste0("HF_mis",mod)),
                       family=binomial,
                       data=sDat,
                       method=brglm2::brglmFit))
    if(debug){
        cat("\n\n>>>> model coefficients for simulated data, model = ",mod,"\n")
        qvcalc::indentPrint(fitSim1[[1]], indent=8)
        qvcalc::indentPrint(fitSim2[[1]], indent=8)
        return(list(sDat, fitSim1, fitSim2))
    }
    #sDat[,HF_mis := respMatMul(dummySim, betas[[2]])]
    #sDat[,cam_fate := relevel(as.factor(sDat[,cam_fate]), ref="H")]
    #dummySim <- model.matrix(form, data = sDat)
    #dummySim <- as.data.table(model.matrix(form, data = sDat))
    # if(xdebug) cat("\n>> and with dummy variables:\n") # ~*~*~*
    # print(class(dummySim))
    #if(debug) aprint(dummySim)
    #if(debug) print.data.table(dummySim)
    # if(debug) bprint(dummySim)
    #if(debug) print(dummySim, topn=5, trunc.cols=T, )
     #if(debug) qvcalc::indentPrint(head(eta), indent=8)
     #if(debug) qvcalc::indentPrint(tail(eta), indent=8)
    #y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
    #return(y)
    return(sDat)
}

respMatMul <- function( dummySim, betas, debug=FALSE){
    if(debug) cat("\n***MATRIX MULTIPLICATION***", dim(dummySim), dim(betas))
    eta <- dummySim %*% betas 
    y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
     #if(debug) cat("\n\t>>>> beta values, class:", class(betas), "\n\n") # ~*~*~*
     if(debug) cat("\n\t>>>> beta values multiplied by design matrix = eta:\n\n") # ~*~*~*
     #if (debug) qvcalc::indentPrint(cbind(betas, dummySim, eta))
     if (debug) qvcalc::indentPrint(paste(betas, dummySim, eta), indent=8)
     #if(debug) qvcalc::indentPrint(rbind(head(betas), tail(betas)), indent=8)
     #if(debug) qvcalc::indentPrint(paste(betas,collapse=" "), indent=8)
     #if(debug) qvcalc::indentPrint(tail(betas), indent=8)
    # print(head(betas))
     #if(debug) cat("\n\t>>> multiplied by design matrix, class:", class(dummySim), "\n")
     #if(debug) qvcalc::indentPrint(rbind(head(dummySim), tail(dummySim)), indent=8)
     #if(debug) qvcalc::indentPrint(str(dummySim), indent=8)
    #if(debug) aprint(dummySim)
    #eta <- dummySim %*% betas
     #if(debug) cat("\n\t>>> equals eta:\n") # ~*~*~*
     #if(debug) qvcalc::indentPrint(rbind(head(eta), tail(eta)), indent=8)
     return(y)
}

mkSim <- function( s_size, cMat, mList, betas,fprob, sprob, stratify=TRUE, debug=FALSE, xdebug=FALSE){
#mkSim <- function(resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE){
    fates <- names(fprob)
    spp   <- names(sprob)
    #if(debug) cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>  ")
    #if(debug) cat(sprintf( "\n<><><><><><><><><><> REP:%s <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>  ",  ))
    if(debug) cat("\n\n + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ")
    if(debug) cat("\n<><><><><><><><><> NEXT REP <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>  ")
    if(debug) cat("\n + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ")
    #if(debug) cat("\n\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>  ")
    if(debug) cat("\n\n<><><><><><><><><><><><><><><> MAKE SIMULATED DATA <><><><><><><><><><><><><><><> ")
    if(debug) cat("\n\tfates:", fates, "& species:", spp) # ~*~*~*
    if(debug){ # ~*~*~*
        #cat("\ncoefficients, class:\n", class(betas))
        cat("\n\tcoefficients:\n")
        #qvcalc::indentPrint(paste(betas, collapse=" "), indent=8)
        qvcalc::indentPrint(betas, indent=8)
    }
    simDat <- list()
    if (stratify) {
        for (s in seq_along(spp)){
            sp <- spp[s]
            #simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=s_size*sprob[s], mu=mList[[sp]], Sigma=cMat[[sp]]))
            if(debug) cat("\n\t>>> species:", sp, "\t& means:", mList[[sp]], "\n\t\t& correlation matrix:\n") # ~*~*~*
            qvcalc::indentPrint(cMat[[sp]], indent=12)
            simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=round(s_size*sprob[s]), mu=mList[[sp]], Sigma=cMat[[sp]]))
            # cat("\n>> adding camera fates\n") # ~*~*~*
            simDat[[sp]][,cam_fate := sample(fates, size=nrow(simDat[[sp]]), replace=T, prob=fprob)]
        }
        if(debug) cat("\n\t>>> created sim data") # ~*~*~*
        names(simDat) <- c("CONI", "LETE")
        sDat <- data.table::rbindlist(simDat, idcol="species", use.names=T, fill=T)
        #cat(" - merged sim data for the two species") # ~*~*~*
    } 
    else {
        if(debug) cat( "\n\t>> means:\n")
        qvcalc::indentPrint(mList, indent=8)
        if(debug) cat("\n\t& correlation matrix:\n") # ~*~*~*
        qvcalc::indentPrint(cMat, indent=8)
        sDat <- as.data.table(MASS::mvrnorm(n=s_size, mu=mList, Sigma=cMat))
        sDat[,species := sample(spp, size=nrow(sDat), replace=T, prob=sprob)]
        sDat[,cam_fate := sample(fates, size=nrow(sDat), replace=T, prob=fprob)]
    }
#mkSim <- function(  s_size, cMat, mList, betas,fprob, sprob,  debug=FALSE, xdebug=FALSE){
    sDat[,cam_fate := relevel(as.factor(sDat[,cam_fate]), ref="H")]
    sDat[,species  := relevel(as.factor(sDat[,species] ), ref="LETE")]
    if(debug) cat(" - releved factors") # ~*~*~*
    #dummySim <- model.matrix(mod, data = sDat)
    id <- seq(1,s_size)
    sDat <- cbind(id, sDat)
    if(xdebug) cat("\n\t>>>>> simulated data frame w/o dummy variables:\n") # ~*~*~*
    if(xdebug) qvcalc::indentPrint(sDat, indent=8)
    #if(debug) print(str(sDat))
    return(sDat)
    #sDat[,is_u := mkResp(sDat, as.matrix(betas[[1]]), form, debug=debug )]
    #sDat[,HF_mis := mkResp(sDat, as.matrix(betas[[2]]), form, debug=debug )]
    #sDat[,is_u := mkResp(sDat, (betas[[1]]), form, debug=debug )]
    #sDat[,HF_mis := mkResp(sDat, (betas[[2]]), form, debug=debug )]
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

#m <- mkSim(resp="is_u", mod=mods4sim[[1]], s_size=500, cMat=cMat, mList=mList, betas=betas, fprob=fprob, sprob=sprob, prList=prList, debug=TRUE)


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


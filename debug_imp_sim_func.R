



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
    if (vb>=1) cat("\n\t>> add dummy variables and remove factors. new columns:\n\t\t", names(dat))
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

    if(vb>=2) cat("\n\t>> check new factor variable & NA distribution after converting:\n")
    if(vb>=2) qvcalc::indentPrint( table(dat$cam_fate, useNA="ifany"), indent=12)
    cat("\n")
    if(vb>=2) qvcalc::indentPrint( colSums(is.na(dat)), indent=12)

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
    if (vb>=1) cat("\n\t>> converted dummies to factors. all columns:\n\t\t", names(dat)) 
    #if (debug) cat("\n>> columns to remove:\n", names(dat)[rem_col])
    if (vb>=2) cat("\n\t>> columns to remove:\n\t\t", rem_col)

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

mkSimDat <- function(seeed, nd, mpatt, wts, new_prop=0.2, patt_freq=c(0.45,0.45,0.1),wt=TRUE, test=FALSE, convFact=FALSE, facToNum=FALSE,vbose=0){
    if(vbose>=1) cat("\n\n<><><><><><><><><><><><><><><> MAKE AMPUTED DATA: <><><><><><><><><><><><><><><> \n ")
    #if(vbose>=3) cat("\n\t>> mkSimDat seed=", seeed, class(seeed))
    dat4amp <- add_dummy(nd, vb=0)
    no_miss <- c("obs_int", "fdate", "is_u", "speciesCONI")
    is_miss <- colnames(mpatt)[!colnames(mpatt) %in% no_miss]
    new_order <- c(is_miss, no_miss)
    if(vbose>=3) cat("\n\t>> reorder columns:\n\t\t", new_order) ### *~*~*~*~* #######
    dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix
    suppressWarnings(amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = mpatt, freq = patt_freq,weights = wts))
    if(vbose>=2) cat("\n\t>> Create more new missing values, with weighted probabilities:\n")
    if(vbose>=3)  qvcalc::indentPrint(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE), indent=12)
    #if(test)  (mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE), indent=8) #figuree out how to save to file
    # missing_tab("amp_out_wt", prVars)
    
    if(wt) {datList <- amp_out_wt} else {datList <- amp_out} ### *~*~*~*~* #######
    datList$amp
    if(vbose>=3) cat("\n\t>> the amputed data before converting to factors:\n")
    if(vbose>=3) qvcalc::indentPrint(summary(datList$amp), indent=12)
    if(vbose>=3) cat("\n\t>> & the NAs:\n")
    if(vbose>=3) qvcalc::indentPrint(colSums(is.na(datList$amp)), indent=12)
    if (convFact) datList$amp <- add_fact(dat = datList$amp, facToNum=facToNum, vb=vbose) # could probably reference the global debug instead...
    if(vbose>=2) cat("\n\n<><><><><><><><><><<><><><><<><><><><> AMPUTED DATA: <><><><><><><><><><<><><><><<><><><><>>>\n ")
    if(vbose>=1) qvcalc::indentPrint(str(datList$amp), indent=12)
    if(vbose>=1) cat("\n\tNA summary:\n")
    if(vbose>=1) qvcalc::indentPrint(colSums(is.na(datList$amp)), indent=12)
    ampd <- datList$amp
    if (test) vmiss <- naniar::vis_miss(ampd)
    if (test) now <- format(Sys.time(), "%d%b%H%M%S")
    if (test) fname <- sprintf("figs/ampDat_%s.svg",now)
    if (test) ggsave(fname, vmiss, device="svg")
    if(test) cat("\n    > missing value figure saved to figs directory\n")


    #if(params$deb) cat("\n>> species for this run:")
    #if(debug) qvcalc::indentPrint(table(datList$amp$speciesCONI), indent=8)
   # if(debug) qvcalc::indentPrint(table(ampd$species), indent=8)
    #if(debug) cat("\n>> cam fate vars for this run:")
    #if(debug) qvcalc::indentPrint(table(datList$amp$cam_fate), indent=8)
    #if(debug) qvcalc::indentPrint(table(ampd$cam_fate), indent=8)
    if(vbose>=1) cat("\n\t AMPUTED DATA-summary:\n")
    if(vbose>=1) qvcalc::indentPrint(summary(datList$amp), indent=12)
    if(vbose>=3) cat("\n\tALL THE AMPUTED DATA:\n")
    if(vbose>=3) qvcalc::indentPrint(datList$amp, indent=12)
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
    if(vbose>=2) sprintf(cat("\n\n<><><><><><><><><><><><><><> MULTIPLE IMPUTATION: <><><><><><><><><><><><><><><><><> m %s <><><><><><><><><>>>\n ", modd))
    if(vbose>=4) cat("\t\n*.*.*.VARS:",vars,".*.*.*\n")
    pr_list <-  c("species", "cam_fate", "obs_int", "nest_age", "fdate")
    ret    <- array(NA, dim=c(length(vars), 3, length(resp_list)))
    #dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %","n"), resp_list)
    dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"), resp_list)
    if(vbose>=4) cat("\n\t>>>> make empty matrix to store output") # ~*~*~*
    if(vbose>=4) qvcalc::indentPrint(ret, indent=12)
    if (met == "cc"){
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            #if(vbose==2) cat("\n    ===========================================================================================")
            if(vbose>=2) cat("\n    ==================== complete-case analysis \t resp:", resp," =====================")# ~*~*~*
            if(vbose>=2) cat("\n    ===========================================================================================\n")
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            if(vbose>=2) cat("\n\t>> vars in this df:", vlist,"\n\t\t& cols for analysis:", cols) # ~*~*~*
            ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
            if(vbose>=3) cat("\n\tcheck incoming data:\n")
            if(vbose>=3) qvcalc::indentPrint(summary(ampDat), indent=12)
            dat1 <- ampDat[complete.cases(ampDat),]
            fit = glm(as.formula(paste0(resp, modd)), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
            vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            if(vbose>=2) cat("\n\t>> output vals:\n") # ~*~*~*
            if(vbose>=2) qvcalc::indentPrint(vals, indent=12)# ~*~*~
            if(vbose>=2) cat("\n\t>> model fit summary:\n") # ~*~*~*
            if(vbose>=3) qvcalc::indentPrint(summary(fit), indent=12)# ~*~*~
            rownames(vals) <- names(coef(fit))[-1]
            colnames(vals) <-  c("estimate", "2.5 %", "97.5 %") # this doesn't work if not a df
            
            ret[vlist,,r]  <- as.matrix(exp(vals))# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
            if(vbose>=3) cat("\n\tret, FILLED IN (exponentiated):\n") # ~*~*~*
            if(vbose>=3) qvcalc::indentPrint(ret[,,r], indent=12)
            #if(vbose==2) cat("\n=== complete cases: ===\n")
            #if(vbose==2) qvcalc::indentPrint(ret[,,r], indent=16)
        }
        return(ret)

    } else {
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            #if (debug) cat("\n    ===========================================================================================")
            if(vbose>=2) cat("\n    ======================= method:", met, "\t resp:", resp," =============================== ")
            if(vbose>=2) cat("\n    ===========================================================================================\n")
            vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
            cols <- c(resp, pr_list)
            if(vbose>=2)cat("\n\t>> vars in this df:", vlist,"\n\t\t& cols for imputation:", cols) # ~*~*~*
            cfates <- paste0("cam_fate", c("H", "A", "D", "F", "Hu", "S"))
            # collapse pkg might be faster - https://stackoverflow.com/questions/77679416/drop-rows-with-na-in-multiple-columns
            if(vbose>=3) cat("\n\t>>>> vars with empty levels:", names(which(sapply(aDat, function(x) length(setdiff(levels(x),unique(x))))>0)))
            if(vbose>=3) cat("\n\t>>>> levels of factors:", levels(aDat$cam_fate),"\t", levels(aDat$HF_mis),"\t", levels(aDat$species))
            ampDat <- aDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?
            if(vbose>=3) cat("\n\t>>>> levels of factors after droplevels():", levels(aDat$cam_fate),"\t", levels(aDat$HF_mis),"\t", levels(aDat$species))
            if (met=="cf_cc") ampDat <- ampDat  %>% filter(!is.na(cam_fate))
            if(vbose>=3) cat("\n\tcheck incoming data:\n")
            if(vbose>=3) qvcalc::indentPrint(summary(ampDat), indent=12)
            if(vbose>=2) cat("\n\t>>>> vars that are missing values:", names(which(colSums(is.na(ampDat))>0)))
            metList <- met_list[resp,,met]
            inters <- sapply(modd,  function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))[[1]]
            inter <- paste(inters, collapse=".")
            if(vbose>=2) cat("\n\t\t >>>>> inter=", inter) # ~*~*~*
            if(met=="passive" & length(inters)!=0)  ampDat[inter] <- NA
            names(metList)[6] <- resp
            names(metList)[7] <- inter
            metList <- metList[!is.na(metList)]
            if(all(is.na(metList))) metList=NULL

            if(vbose>=2) cat("\n\n\t<><><><><><><><><><><><><><><> imputation with MICE <><><><><><><><><><><><><><><> \n ")
            if(vbose>=2){ ## *~*~*~*~*
                cat("\n\tMETHOD LIST for",met,":\n") # ~*~*~*
                qvcalc::indentPrint(metList, indent=12)
                #cat("\n")
            }

            frmla <- lapply(form_list[[resp]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
            if(vbose>=3) cat("\n\tmice formulas:\n") # ~*~*~*
            if(vbose>=3) cat("\t\t",paste(frmla, collapse="\n\t\t"))

            impCall <- case_when(met=="stratify"~ "impStrat(ampDat, met=metList,formulas=frmla, col_sel=cols)",
                       #met=="passive" ~ "mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)",
                       met=="passive" ~ "mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)",
                       #.default="mice::mice(ampDat, method=metList, m=m, print=FALSE)"
                       ### turn off print=FALSE to get output from imputations as they happen
                       #.default="mice::mice(ampDat, method=metList, m=m, formulas=frmla,print=FALSE)"
                       #.default="mice::mice(ampDat, method=metList, m=m,print=FALSE)"
                       .default="mice::mice(ampDat, method=metList, m=m, formulas=frmla,print=FALSE)"
            )
            if(vbose>=2) cat("\n\n\t>> call:", impCall,"\n")
            imp <- eval(parse(text=impCall))
            #if(vbose>=2) cat("\n\n    <><><><><><><><><><><><><><><> MICE results <><><><><><><><><><><><><><><> \n")
            if(vbose>=2) cat("\n\n\t<><><><><><><><><><><><><><><> MICE results <><><><><><><><><><><><><><><> \n")
            if(vbose>=3) cat("\n\t>> the imputed values:\n")
            #if(xdebug) qvcalc::indentPrint(str(imp$imp))
            #impd <- sapply(imp$imp, function(x) any(rowSums(x)>0))
            #impd <- sapply(imp$imp, function(x) any(x>0))
            impd <- sapply(imp$imp, function(x) any(!is.na(x)))
            if(vbose>=3) qvcalc::indentPrint(impd, indent=12)
            if(vbose>=2) cat("\n\t>>>>> vars that were imputed:", names(impd)[impd])

            if(length(imp$loggedEvents > 0)) { # ~*~*~*
                cat(sprintf("\n\n*** LOGGED EVENTS FOR MODEL %s & METHOD %s & RESP %s", modd, met, resp), file=outFile, append=TRUE)
                cat("\n\n it im dep \t meth \t\tout\n", file=outFile, append=TRUE)
                cat(paste(imp$loggedEvents, collapse=" "), file=outFile, append=TRUE)
                cat(sprintf("\n\t*** LOGGED EVENTS FOR METHOD %s*******************************\n", met))
                qvcalc::indentPrint(imp$loggedEvents, indent=12)
            }
            if(impplot) visImp(imp)  # ~*~*~*
            if(vbose>=2) cat("\n\n\t<><><><><><><><><><><><><><><> Analysis <><><><><><><><><><><><><><><>\n ")
            if(vbose>=2) cat("\n\t>> fitting model", modd)
            fit = with(imp,
                 glm( as.formula( paste0(resp, modd) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
                 ))
            if(vbose>=2) cat("\t>> pooling model fit") # ~*~*~*
            if(vbose>=2) cat("\n\t    pool() OUTPUT:\n") # ~*~*~*
            pool = summary(mice::pool(fit), "all", conf.int=TRUE, exponentiate=TRUE)
            if(vbose>=2) qvcalc::indentPrint(pool,indent=16)
            rownames(pool) = pool[,"term"]
            pool <- pool %>% filter(term %in% vars)
            #cat("\n")

            #print(pool)
            if(vbose>=2){
                cat("\n\t    pool() OUTPUT with vars selected:\n") # ~*~*~*
                qvcalc::indentPrint(pool[,c("estimate", "2.5 %", "97.5 %")], indent=16)
                #outp <- pool[,c("estimate", "2.5 %", "97.5 %")]
                #comp <- cbind()
            }
            if(vbose>=3){
                cat("\n\t    PUT IT HERE:\n")
                qvcalc::indentPrint(ret[vlist,,r], indent=12)
            }
            ret[vlist,,r] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
            if(vbose>=3) cat("\n\t    ret, FILLED IN:\n")
            if(vbose>=3) qvcalc::indentPrint(ret[,,r], indent=16)
            #if(mindebug) cat(sprintf("\n=== imputation results for %s & %s: ===\n", met, resp))
            #if(mindebug) qvcalc::indentPrint(ret[,,r], indent=16)
        }
    #if(vbose>=2) print(ret[,1,])
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
    qvcalc::indentPrint(printImp, indent=12)
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
    while(!success){
        #sDat <- mkSim(s_size, cMat, mList, betas, fprob, sprob, debug=debug, xdebug=xdebug)
        sDat <- mkSim(s_size, cMat, mList, betas, fprob, sprob, vb=vbose)
        form <- as.formula(mod)
        if(vbose>=2) cat("\n\t>>>> making response variables with model ", mod, "- formula:", form)
        dummySim <- model.matrix(form, data = sDat)
        if(vbose>=2) cat("\n\t>> do dims match?", dim(dummySim), length(betas[[1]]), length(betas[[2]]))
        success <- dim(dummySim)[2]==length(betas[[1]]) & dim(dummySim)[2]==length(betas[[2]])
        if(vbose>=2) cat("\n\t>>> ***success?***", success)
    }
    if(vbose>=2) cat("\n\t>>> is_u betas:", betas[[1]])
    if(vbose>=2) cat("\n\t>>> HF_mis betas:", betas[[2]])
    sDat[,is_u := respMatMul(dummySim, betas[[1]], vb=vbose)]
    sDat[,HF_mis := respMatMul(dummySim, betas[[2]], vb=vbose)]
    sDat[,is_u := as.factor(sDat[,is_u])]
    sDat[,HF_mis := as.factor(sDat[,HF_mis])]
    if(vbose>=1){
        cat("\n\t>> simulated data with response variables:\n ")
        qvcalc::indentPrint(sDat, indent=12)
        cat("\n")
        qvcalc::indentPrint(summary(sDat), indent=12)
        cat("\n\t    >> camera fate levels (check for small cell sizes):\n ")
        qvcalc::indentPrint(table(sDat$cam_fate), indent=8)
    }
    #cat("\n>> added response variables, which should be factors:", sapply(sDat[,c(HF_mis, is_u)], is.factor))
    # DOESN'T work with data.table:
    #cat("\n>> added response variables, which should be factors:", sapply(as.data.frame(sDat)[c(HF_mis, is_u)], is.factor))
    #fitReal2 <- summary(glm(is_u ~ species * nest_age + obs_int + cam_fate + fdate,
    #fits <- readRDS("fits.rds")
    #if(vbose>=3) cat("\nmodel fits imported from file:\n")
    #if(vbose>=3) qvcalc::indentPrint(fits, indent=12)
    #coefss <- sapply(fits, "[[", 12) 
    #modnum <- as.numericnames(mod)
    #modnum<- gsub("\\D+", "", names(mod))
    #mnums <- as.numeric(gsub("\\D+", "", names(fits)))
    #idx <- gsub("\\D+", "", names(fits))
    #mnums <- rep(c("1", "8", "16"))
    #if(vbose>=2) cat("\n\tmodel number:", modnum, "and model nums from list:", paste(mnums, collapse=","))
    #if(debug) cat("\n\t>> Actual model fit:\n")
    #if(vbose>=3) cat("\n\t>> Actual model coefs (from file):\n")
    #qvcalc::indentPrint(sapply(fits[[mnums==names(m)]], summary))
    # it's already a summary
    #qvcalc::indentPrint(fits[mnums==modnum], indent=16)
    #if(vbose>=3) qvcalc::indentPrint(coefss[mnums==modnum], indent=12)
    #fitSim1 <- summary(glm(as.formula(paste0("is_u",mod)),
    fitSim1 <- glm(as.formula(paste0("is_u",mod)),
                   family=binomial,
                   data=sDat,
                   method=brglm2::brglmFit)
    fitSim2 <- glm(as.formula(paste0("HF_mis",mod)),
                   family=binomial,
                   data=sDat,
                   method=brglm2::brglmFit)
    if(vbose>=1){
        cat("\n\t>>>> model coefficients for simulated data, model = ",mod,"\n\n")
        #qvcalc::indentPrint(coef(fitSim1[[1]]), indent=16)
        #qvcalc::indentPrint(coef(fitSim1[[1]]), indent=12)
        cat("\nis_u & exponentiated is_u:  ")
        qvcalc::indentPrint(coef(fitSim1), indent=12)
        cat("\n")
        qvcalc::indentPrint(exp(coef(fitSim1)), indent=12)
        cat("\nHF_mis & exponentiated HF_mis:  ")
        qvcalc::indentPrint(coef(fitSim2), indent=12)
        cat("\n")
        qvcalc::indentPrint(exp(coef(fitSim2)), indent=12)
    #    return(list(sDat, fitSim1, fitSim2))
    }
    return(sDat)
}

respMatMul <- function( dummySim, betas, vb=0){
    if(vb>=3) cat("\n\t***MATRIX MULTIPLICATION***", dim(dummySim), dim(betas))
    eta <- dummySim %*% betas 
    y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
     #if(debug) cat("\n\t>>>> beta values, class:", class(betas), "\n\n") # ~*~*~*
     if(vb>=2) cat("\n\t>>>> beta values multiplied by design matrix = eta:\n\n") # ~*~*~*
     if(vb>=2) qvcalc::indentPrint(matlib::printMatEqn(betas, "*", dummySim, "=", eta), indent=8)
     if(vb>=2) cat("\n\n    >> y values are chosen from a binomial distribution with inverse link eta:\n")
     if(vb>=2) qvcalc::indentPrint(y, indent=8)
     #if (debug) qvcalc::indentPrint(cbind(betas, dummySim, eta))
     #cat("\ndummySim:")
     #qvcalc::indentPrint(str(dummySim), indent=8)
     #cat( "\neta:")
     #qvcalc::indentPrint(str(eta), indent=8)
     return(y)
}

#mkSim <- function( s_size, cMat, mList, betas,fprob, sprob, stratify=TRUE, debug=FALSE, xdebug=FALSE, mindebug=FALSE){
mkSim <- function( s_size, cMat, mList, betas,fprob, sprob, stratify=TRUE, vb=0){
#mkSim <- function(resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE){
    fates <- names(fprob)
    spp   <- names(sprob)
    if(vb>=1) cat("\n\n + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ")
    if(vb>=1) cat("\n + + + + + + + + + NEXT REP  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  ")
    if(vb>=1) cat("\n + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ")
    
    if(vb>=1) cat("\n\n<><><><><><><><><><><><><><><> MAKE SIMULATED DATA <><><><><><><><><><><><><><><>\n ")
    if(vb>=2) cat("\n\tfates:", fates, "& species:", spp) # ~*~*~*
    if(vb>=1){ # ~*~*~*
        #cat("\ncoefficients, class:\n", class(betas))
        cat("\n    > coefficients:\n")
        qvcalc::indentPrint(betas, indent=8)
        cat("\n    > exponentiated coefficients:\n")
        qvcalc::indentPrint(exp(betas), indent=8)
    }
    simDat <- list()
    if (stratify) {
        for (s in seq_along(spp)){
            sp <- spp[s]
            #simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=s_size*sprob[s], mu=#mList[[sp]], Sigma=cMat[[sp]]))
            if(vb>=1) cat("\n    >>> species:", sp, "\t& means:", mList[[sp]], "\n\t\t& correlation matrix:\n") # ~*~*~*
            qvcalc::indentPrint(cMat[[sp]], indent=8)
            simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=round(s_size*sprob[s]), mu=mList[[sp]], Sigma=cMat[[sp]]))
            # cat("\n>> adding camera fates\n") # ~*~*~*
            simDat[[sp]][,cam_fate := sample(fates, size=nrow(simDat[[sp]]), replace=T, prob=fprob)]
        }
        if(vb>=1) cat("\n    >>> created sim data") # ~*~*~*
        names(simDat) <- c("CONI", "LETE")
        sDat <- data.table::rbindlist(simDat, idcol="species", use.names=T, fill=T)
        #cat(" - merged sim data for the two species") # ~*~*~*
    } 
    else {
        if(vb>=2) cat( "\n\t>> means:\n")
        qvcalc::indentPrint(mList, indent=12)
        if(vb>=2) cat("\n\t& correlation matrix:\n") # ~*~*~*
        qvcalc::indentPrint(cMat, indent=12)
        sDat <- as.data.table(MASS::mvrnorm(n=s_size, mu=mList, Sigma=cMat))
        sDat[,species := sample(spp, size=nrow(sDat), replace=T, prob=sprob)]
        sDat[,cam_fate := sample(fates, size=nrow(sDat), replace=T, prob=fprob)]
    }
    #mkSim <- function(  s_size, cMat, mList, betas,fprob, sprob,  debug=FALSE, xdebug=FALSE){
    sDat[,cam_fate := relevel(as.factor(sDat[,cam_fate]), ref="H")]
    sDat[,species  := relevel(as.factor(sDat[,species] ), ref="LETE")]
    if(vb>=1) cat("    >>> releved factors") # ~*~*~*
    #dummySim <- model.matrix(mod, data = sDat)
    id <- seq(1,s_size)
    sDat <- cbind(id, sDat)
    if(vb>=3) cat("\n    >>>>> simulated data frame w/o dummy variables:\n") # ~*~*~*
    if(vb>=3) qvcalc::indentPrint(sDat, indent=12)
    cat("\n")
    if(vb>=3) qvcalc::indentPrint(summary(sDat), indent=12)
    #if(debug) print(str(sDat))
    return(sDat)
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

calcBias <- funtion(impDat, mods4sim, vars, mets, rlist){
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


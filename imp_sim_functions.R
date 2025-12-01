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
    print(class(dat))
    # if (debug) cat("\n\n>> add dummy variables and remove factors. new columns:\n", names(dat))
    return(dat)
}

########################################################################################
###### DUMMIES --> FACTORS #############################################################
########################################################################################

if(FALSE){
  dat <-amp_out_wt$amp
}
add_fact <- function(dat, facToNum=FALSE, debug=FALSE){ 
    cat("\n NOTE: can speed up conversion back to factor with tidytable\n")
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
      .default="H"),
      species = if_else(speciesCONI==1, "CONI", "LETE"))

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
                    "U"  ~ 7
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
    # if (debug) cat("\n>> converted dummies to factors. all columns:\n", names(dat)) 
    # if (debug) cat("\n>> columns to remove:\n", names(dat)[rem_col])
    # if (debug) cat("\n>> columns to remove:\n", rem_col)

    # dat <- 
    # dat[-c("species")]
    # dat[rem_col]
    # dat <- dat[,-rem_col]
    dat <- dat[,-which(names(dat) %in% rem_col)]
    return(as.data.table(dat))
    # there MUST be a shorter way to do this, using regex?
    # fate_letter <- str_extract_all(dat$)
    # dat <- dat %>% mutate(cam_fate = case_when())
    # dat$cam_fate <- case
}

########################################################################################
##### MAKE SIMULATED DATA ###############################################################
########################################################################################

mkSimDat <- function(seeed, nd, mpatt, wts, new_prop=0.2, patt_freq=c(0.45,0.45,0.1),wt=TRUE, debug=FALSE, xdebug=FALSE, convFact=FALSE, facToNum=FALSE){
  cat("mkSimDat seed=", seeed, class(seeed))
  # if(method=="amp"){
    dat4amp <- add_dummy(nd, debug=debug)
    set.seed(seed=seeed)
    # no_miss <- c("obs_int", "fdate", "is_u", "speciesLETE", "speciesCONI")
    no_miss <- c("obs_int", "fdate", "is_u", "speciesCONI")
    is_miss <- colnames(mpatt)[!colnames(mpatt) %in% no_miss]
    new_order <- c(is_miss, no_miss)
    ### *~*~*~*~* #######
    if(xdebug) cat("\n\n>> reorder columns:", new_order)
    dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix
    
    suppressWarnings(amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = mpatt, freq = patt_freq,weights = wts))
    
    if(FALSE){
      out <- add_fact(amp_out_wt$amp,facToNum = T, debug=T)
      # levels(out$cam_fate)
      table(out$cam_fate)
      table(out$HF_mis)
      table(out$is_u)
      xtabs(formula = ~ cam_fate + species, data = out)
    }
    ### *~*~*~*~* #######
    # if(debug) cat("\n\nCreate more new missing values, with weighted probabilities:\n")
    # if(debug)  print(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE))
    # missing_tab("amp_out_wt", prVars)
    
    if(wt) {datList <- amp_out_wt} else {datList <- amp_out}
    ### *~*~*~*~* #######
    datList$amp
    # if(debug) print(str(datList$amp))
    # if (debug) print(class(datList$amp))
    if (convFact) datList$amp <- add_fact(dat = datList$amp, facToNum=facToNum, debug=debug) # could probably reference the global debug instead...
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
# mkImpSim <- function(fullDat, ampDat, cols, resp, mods, vars, met, m=20, fam=binomial, regMet="brglm_fit", iter=500, passive="both", debug=FALSE, xdebug=FALSE, impplot=FALSE){
#mkImpSim <- function(fullDat, ampDat, cols, resp, mods, vars, met, form_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
# mkImpSim <- function(fullDat, ampDat,pr_list, resp_list, mod, vars, met, form_list, met_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
mkImpSim <- function(fullDat, aDat, resp_list, modd, vars, met, form_list, met_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
    # cat("data for imputation:\n")# ~*~*~
   # print(str(ampDat))
    # bprint(aDat)
    # vlist <- colnames(ampDat)
    vlist <- colnames(model.matrix(as.formula(paste0(resp, modd)),data = aDat))[-1]
    # if(debug) cat("\n*** vars in the df:  ", vlist)
    # colnames(vlist)
    # ampDat <- add_fact(ampDat) # convert dummies back to factor
    pr_list <-  c("species", "cam_fate", "obs_int", "nest_age", "fdate")
    # cat("\nvars in this df:", vlist)
    #ret    <- array(NA, dim=c(length(vars), 3))
    ret    <- array(NA, dim=c(length(vars), 3, length(resp_list)))
    dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"), resp_list)
    # cat("\n\n >>>> empty matrix to store output:\n") # ~*~*~*
    # print(ret)
    # dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %", "fmi"), names(mods))
    #dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"))
    if (met == "cc"){
        # ### *~*~*~*~* #######
        if(debug) cat("\n===== complete-case analysis \t resp:", resp," =====\n")# ~*~*~*
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            cols <- c(resp, pr_list)
            cat("\nvars in this df:", vlist)
            cat("\t& cols for imputation:", cols) # ~*~*~*
            ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
            dat1 <- ampDat[complete.cases(ampDat),]
            # ### *~*~*~*~* #######
            # if(debug) cat("\n\n>> data:\n") # ~*~*~*
            # if(debug) print(str(dat1))
            fit = glm(as.formula(paste0(resp, modd)), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
            vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            cat("\n>> vals:")
            print(vals)# ~*~*~
            # cat("\n>> ret to fill:")
            # print(ret[vlist,,r])# ~*~*~
            rownames(vals) <- names(coef(fit))[-1]
            colnames(vals) <-  c("estimate", "2.5 %", "97.5 %") # this doesn't work if not a df
            
            #    vmatch <- match(rownames(vals), rownames(ret)) # col 1 of vals is the row names
            ret[vlist,,r]  <- as.matrix(vals)# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
            if(debug) cat("\nret, FILLED IN :")
            if(debug) print(ret[,,r])
        }
        return(ret)

    } else {
        # if (debug) cat("\n===== method:", met)
        # for(r in resp_list){
        for(r in seq_along(resp_list)){
            resp <- resp_list[r]
            if (debug) cat("\n===== method:", met)
            if(debug) cat("\t resp:", resp," =====\n")
            cols <- c(resp, pr_list)
            if(debug)cat("\nvars in this df:", vlist)
            if(debug) cat("\t& cols for imputation:", cols) # ~*~*~*
            cfates <- paste0("cam_fate", c("H", "A", "D", "F", "Hu", "S"))
            # if (met=="cf_cc") ampDat <- ampDat %>% select(all_of(cols)) %>% filter(!is.na(cam_fate))
            # collapse pkg might be faster - https://stackoverflow.com/questions/77679416/drop-rows-with-na-in-multiple-columns
            ampDat <- aDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?
            if (met=="cf_cc") ampDat <- ampDat  %>% filter(!is.na(cam_fate))
            # if (met=="cf_cc") ampDat <- ampDat %>% filter(if_any(cfates, ~ !is.na(.)))
              #for(y in seq_along(mods)){
                # needs to be a named list; but if you keep the names, it doesn't drop NA...
            # print(met_list)
            # cat("\n<><><><><>")
            # print(dimnames(met_list))
            # metList <- met_list[resp,,met,mod]
            metList <- met_list[resp,,met]
            # cat("\ngot method list")
            # ampDat <- ampDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?
            inters <- sapply(modd,  function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))[[1]]
            inter <- paste(inters, collapse=".")
            cat("inter=", inter)
            # if(met=="passive" & length(inters)!=0)  ampDat[inter] <- NA
            if(met=="passive" & length(inters)!=0)  ampDat[, ..inter := NA]
            # if(met=="passive" & length(inters)!=0)  ampDat[, (inter) := NA]
            # if(met=="passive" & length(inters)!=0) set(ampDat, "inter", NA)
            # if(met=="passive" & length(inters)!=0) set(ampDat, as.character(inter), NA)
            
            print(ampDat[,inter])
            names(metList)[6] <- resp
            names(metList)[7] <- inter
            metList <- metList[!is.na(metList)]
            if(all(is.na(metList))) metList=NULL

            ## *~*~*~*~*
            cat("\nMETHOD LIST for",met,"-") # ~*~*~*
            print(metList)
            cat("\n")

            # frmla <- lapply(form_list[[mod]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
            # cat("formula list:")
            # print(form_list)
            frmla <- lapply(form_list[[resp]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
                ## *~*~*~*~*
            if(debug) cat("\nmice formulas:\n") # ~*~*~*
            if(debug) print(frmla)

                # 3: In all(metList) : coercing argument of type 'character' to logical
                # if(is.na(all(metList))) metList=NULL # this is incorrect I think
                # mlt <- metList[!is.na(metList)]
                # dimnames(metList)
            impCall <- case_when(met=="stratify"~ "impStrat(ampDat, met=metList,formulas=frmla, col_sel=cols)",
                       met=="passive" ~ "mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)",
                       .default="mice::mice(ampDat, method=metList, m=m, print=FALSE)"
            )
            imp <- eval(parse(text=impCall))

            if(length(imp$loggedEvents > 0)) { # ~*~*~*
                cat(sprintf("\n*** LOGGED EVENTS FOR METHOD %s*******************************\n", met))
                print(imp$loggedEvents)
            }
            # if(impplot) visImp(imp)  # ~*~*~*
            # if(debug) cat("\n>> fitting model", mods[y])
            fit = with(imp,
                 glm( as.formula( paste0(resp, modd) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
                 ))
            # if(debug) cat("\n>> pooling model fit\n") # ~*~*~*
            pool = summary(mice::pool(fit), "all", conf.int=TRUE, exponentiate=TRUE)
            # print(pool)
            rownames(pool) = pool[,"term"]
            pool <- pool %>% filter(term %in% vars)
            cat("\nOUTPUT:") # ~*~*~*
            print(pool[,c("estimate", "2.5 %", "97.5 %")])
            # cat("\nPUT IT HERE:")
            # print(ret[vlist,,r])
            #vmatch <- match(pool[,1], rownames(ret)) # col 1 of vals is the row names
            # ret[vmatch,,y] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %", "fmi")])
            #ret[vmatch,,y] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
            ret[vlist,,r] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
            if(debug) cat("\nret, FILLED IN :")
            if(debug) print(ret[,,r])
            # if(xdebug) { # ~*~*~*
            #    cat("\n\n>> pooled model fit:\n")
            #str(pool)
            #    print(head(pool[,c("term", "estimate", "2.5 %", "97.5 %")]))
            # }
            #}
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


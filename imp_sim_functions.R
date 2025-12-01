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
    dummy_vars

    dat$speciesLETE <- ifelse(dat$species=="LETE", 1, 0)
    dat$speciesCONI <- ifelse(dat$species=="CONI", 1, 0)

    # why did I get the columns in this way? I guess to reorder them?
    dat<- dat[, c("obs_int", "nest_age", "fdate", "HF_mis", "is_u", "speciesLETE", "speciesCONI")]
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
    # I already know what the dummy var names are in this case; not a general purpose function
    # dat <- amp_out$amp
    # if(is.list(dat)) dat <- dat$amp
    # if(any(dim(dat))) dat <- dat$amp
    # dim(p) >0

    dat <- dat %>% 
    mutate(cam_fate = case_when(cam_fateA==1 ~ "A",
                cam_fateH==1 ~ "H",
                cam_fateHu==1 ~ "Hu",
                cam_fateF==1 ~ "F",
                cam_fateS==1 ~ "S",
                cam_fateD==1 ~ "D"),
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
    return(dat)
    # there MUST be a shorter way to do this, using regex?
    # fate_letter <- str_extract_all(dat$)
    # dat <- dat %>% mutate(cam_fate = case_when())
    # dat$cam_fate <- case
}

########################################################################################
####### MAKE METHOD LIST #############################################################
########################################################################################

if(FALSE){
    fateMet <- "pmm"
    ageMet <- "pmm"
    misMet <- "pmm"
    resp <- "HF_mis"
    col_sel[6] <- "HF_mis"
    col_sel <- c(col_sel, "speciesLETE:is_u")
    col_sel <- names(ampDat)
    c<-col_sel[2]
    c <- inter
    met <- "rf"
    met <- "caliber"
    met <- "passive"
    met <- "stratify"
    dat <- ampDat
    int <- NULL
    int <- paste(inters, collapse="*")
}
# mkMetList <- function(met, dat, col_sel, int=NULL, debug=FALSE){
mkMetList <- function(met, dat, int=NULL, debug=FALSE){
    col_sel <- names(dat)
    metList <- rep("",length(col_sel))
    catList <- c("HF_mis", "cam_fate", "species", "is_u")
    names(metList) <- col_sel

    colNA <- names(dat)[colSums(is.na(dat)) > 0]
    is_num <- names(dat)[sapply(dat, is.numeric)]
    is_bin <- names(dat)[sapply(dat, function(x) length(levels(x)))==2]
    is_cat <- names(dat)[sapply(dat, function(x) length(levels(x)))>2]


    for(c in col_sel){
    if(met=="caliber"){
    if(c %in% catList ) met1 = "rfcat" else met1 = "rfcont"
    }else if(met %in% c("default.int", "pmm.int")){
    met1 <- str_extract(met, pattern = "\\w+(?=.)") # everything up until the period
    # }else if(met == "passive.int"){
    }else if(met == "passive" | met=="stratify" | met=="cf_cc"){
    # met1 <- ""
    # met1 <- ifelse(is.numeric())
    met1 <- case_when(c %in% is_num ~ "pmm",
            c %in% is_bin ~ "logreg",
            # c %in% is_cat ~ "polyreg")
            c %in% is_cat ~ "pmm") # because I have many cells w/ <10 obs
    # met2 <- paste("~I(", int,")")
    } else {met1 <- met}

    interTrue <- ifelse((grepl(".",c,fixed=T)), TRUE, FALSE)
    # make NAs so that those columns can be dropped later
    metList[c] <- case_when(
    # interTrue & met=="passive" ~ paste("~I(",int,")"),
    met=="cf_cc" & c=="cam_fate" ~ "",
    met =="default" ~ NA,
    met == "passive" ~ NA,
    met=="stratify" & c=="species" ~ NA, # this should still exist, but take species out of formula?
    # met =="default" ~ "",
    interTrue & met !="passive" ~ NA,
    c=="inter" ~ NA,
    !interTrue & c%in%colNA ~ met1,
    .default = ""
    )

    }
    return(metList)
}

########################################################################################
##### MAKE SIMULATED DATA ###############################################################
########################################################################################

if (FALSE){
    # nd <- ndGLM_scl_cc
    nd <- dat4sim
    wt <- TRUE
    convFact <- TRUE
    debug <- params$deb
    debug <- TRUE
    # facToNum <- TRUE
    facToNum <-FALSE 
    vars <- var_list
    seeed <- 614
}

mkSimDat <- function(seeed, nd, vars, facToNum=FALSE, method="amp", wt=TRUE, debug=FALSE, xdebug=FALSE, convFact=FALSE){
    if(method=="amp"){
        cat("\n>>> making amputed data\n") # ~*~*~*
        #nd <- as.matrix(nd)
        dat4amp <- add_dummy(nd, debug=debug)
        set.seed(seed=seeed)
        bprint(dat4amp)# ~*~*~*
        #suppressWarnings(amp_out1 <- mice::ampute(dat4amp))
        suppressMessages(suppressWarnings(amp_out1 <- mice::ampute(dat4amp)))
        bprint(amp_out1$amp)# ~*~*~*

        if(FALSE){
            out <- add_fact(amp_out1$amp,facToNum = T, debug=T)
            # levels(out$cam_fate)
            table(out$cam_fate)
            xtabs(formula = ~ cam_fate + species, data = out)
            xtabs(formula = ~ cam_fate + species, data = nd)
        }
        new_patt <- amp_out1$patterns
        no_miss <- c("obs_int", "fdate", "is_u", "speciesLETE", "speciesCONI")
        # The problem with doing this with indices is that (as just happened), I forget how
        # I did things and then add new variables w/o accounting for them     # Easier to reference names?
        is_miss <- names(new_patt)[!names(new_patt) %in% no_miss]
        is_miss

        # the order: missing age only; missing fate/HF_mis only; missing both
        miss_pat <- list(c(0, rep(1,7)), c(1,rep(0,7)), c(rep(0,8)) )
        miss_pat

        m <- do.call(rbind,miss_pat) # create the matrix of the missingness patterns 
        p <- do.call(rbind, rep(list(c(rep(1,5))), 3)) # create additional columns for the vars not missing values

        new_prop <- 0.16
        miss_patt_mat <- cbind(m,p)
        colnames(miss_patt_mat) <- c(is_miss, no_miss)
        patt_freq <- c(0.45,0.45,0.1) # missing: age only, fate only, both
        ### *~*~*~*~* #######
        #if(debug & seeed==1){
        rr=0
        if(debug & rr<1){
            cat("\n\n>> missingness matrix:\n")
            print(miss_patt_mat)
            cat("\n>> proportion missing:", new_prop, "\n\n>> frequency of each pattern:", patt_freq)
            rr = rr +1
        }
        new_order <- c(is_miss, no_miss)
        ### *~*~*~*~* #######
        if(debug) cat("\n\n>> reorder columns:", new_order) # ~*~*~*
        dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix

        suppressWarnings(amp_out <- mice::ampute(dat4amp, prop = new_prop, patterns = miss_patt_mat, freq = patt_freq))
        if(FALSE){
            out <- add_fact(amp_out$amp,facToNum = T, debug=T)
            # levels(out$cam_fate)
            table(out$cam_fate)
            xtabs(formula = ~ cam_fate + species, data = out)
        }
        ### *~*~*~*~* #######
        # if(xdebug) cat("\n\nCreate new missing values:\n") # ~*~*~*
        # if(xdebug) print(mice::md.pattern(amp_out$amp, rotate.names = TRUE)) # ~*~*~*

        # 1 = complete; 0 = has missings
        wts <- amp_out$weights 
        wts[1,] <- c(0, 0, 0, 0.8, 0.2, rep(0, 8)) # missing age only - what vars contribute
        # wts[3,] <- c(0, 0, 0, 0.8, 0.2, rep(0, 8))
        ### *~*~*~*~* #######
         if(debug){ # ~*~*~*
           cat("\n>> new weights:\n")
           print(wts)
           cat("\n>> to go with the patterns:\n")
           print(miss_patt_mat)
         }
        # what exactly are the weights doing?
        suppressWarnings(amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = miss_patt_mat, freq = patt_freq,weights = wts))

        if(FALSE){
            out <- add_fact(amp_out_wt$amp,facToNum = T, debug=T)
            # levels(out$cam_fate)
            table(out$cam_fate)
            table(out$HF_mis)
            table(out$is_u)
            xtabs(formula = ~ cam_fate + species, data = out)
        }
        ### *~*~*~*~* #######
        # if(debug) cat("\n\nCreate more new missing values, with weighted probabilities:\n") # ~*~*~*
        # if(debug)  print(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE)) # ~*~*~*
        # missing_tab("amp_out_wt", prVars)

        # datList <- ifelse(wt, amp_out_wt, amp_out) # ifelse is not what I need here
        if(wt) {datList <- amp_out_wt} else {datList <- amp_out}
        ### *~*~*~*~* #######
        if(debug) print(names(datList)) # ~*~*~*
        if(debug) bprint(amp_out_wt)
        #if(debug) print(str(datList$amp)) # ~*~*~*
        # if (debug) print(class(datList$amp)) # ~*~*~*
        if (convFact) datList$amp <- add_fact(dat = datList$amp, facToNum=facToNum, debug=debug) # could probably reference the global debug instead...
        # datList$amp <- datList$amp %>%
        # levels(datList$amp$cam_fate)
        # str(datList)
        # datList
        # return(as.matrix(dat))
        return(datList)

    } else {

    ageNA <- sample(x=nd$nest, size = 10, replace = FALSE)
    nd$nest_age[!is.na(match(nd$nest,ageNA))] = NA
    fateNA <- sample(x=nd$nest, size = 12, replace = FALSE)
    nd$cam_fate[!is.na(match(nd$nest,fateNA))] = NA

    # ### *~*~*~*~* #######
    # if (debug){ # ~*~*~*
    #   cat("new number missing nest_age:", sum(is.na(nd$nest_age)),
    #       "new percent missing nest_age:", sum(is.na(nd$nest_age)) / nrow(nd),
    #       "original percent missing:", sum(is.na(ndGLM_scl_all$nest_age)) / nrow(ndGLM_scl_all) )
    #   cat("\nnew number missing cam_fate:", sum(is.na(nd$cam_fate)),
    #       "new percent missing cam_fate:", sum(is.na(nd$cam_fate)) / nrow(nd),
    #       "original percent missing:", sum(is.na(ndGLM_scl_all$cam_fate)) / nrow(ndGLM_scl_all) )
    #   
    # }
    return(nd)
  }
}

########################################################################################
##### IMPUTE/FIT/POOL SIM DATA ##########################################################
########################################################################################

if(debugging){
    # dat=mkSimDat(ndGLM_scl_cc, convFact = TRUE)$amp
    # ampDat=mkSimDat(ndGLM_scl_cc, convFact = TRUE)$amp
    s = 614
    s=618
    ampDat = mkSimDat(seeed=s, nd=dat4sim, vars=var_list, convFact=TRUE)$amp
    xtabs(~ cam_fate + species, data=ampDat)
    # ampDat <- datList$amp
    # ampDat <- datNA
    # aDat = sim_dat$amp
    # aDat = sim_dat$amp[1:20,]
    resp="is_u"
    # resp = "HF_mis"
    # resp="isu"
    # resp = "HFmis"
    mod=modList[1]
    mods = mods4sim
    m=20
    met="pmm"
    met="rf"
    met="default"
    # met="cc"
    # met = "caliber"
    met="passive"
    met = "stratify"
    mets <- met_list
    # seed=61389

    fam=binomial
    regMet="brglm_fit"
    iter=500
    cols <- col_list
    debug = TRUE
    xdebug=TRUE
    y=1
    # why only these vars?
    # vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
    vars <- var_list
    impplot = TRUE
    fcToNum <- TRUE
    # mod = modList[1]
    # ampDat <- dat
    # ampDat <- datList$amp
}

### Nest the loop inside the if statement so you aren't running the if check every loop?
# mkImpSim <- function(fullDat, ampDat, cols, resp, mods, vars, met, m=20, fam=binomial, regMet="brglm_fit", iter=500, passive="both", debug=FALSE, xdebug=FALSE, impplot=FALSE){
#mkImpSim <- function(fullDat, ampDat, cols, resp, mods, vars, met, form_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
mkImpSim <- function(fullDat, ampDat,pr_list, resp_list, mod, vars, met, form_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
    cat("data for imputation:\n")# ~*~*~
   # print(str(ampDat))
    bprint(ampDat)
    #ret    <- array(NA, dim=c(length(vars), 3))
    ret    <- array(NA, dim=c(length(vars), 3, length(resp_list)))
    # dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %", "fmi"), names(mods))
    #dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"))
    dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"), resp_list)
    if (met == "cc"){
        # ### *~*~*~*~* #######
        # if(debug) cat(" complete-case analysis:\n") # ~*~*~*
        for(r in resp_list){
            cols <- c(r, pr_list)
            cat("\n\t& cols for imputation:", cols) # ~*~*~*
            ampDat <- ampDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
            dat1 <- ampDat[complete.cases(ampDat),]
            # ### *~*~*~*~* #######
            # if(debug) cat("\n\n>> data:\n") # ~*~*~*
            # if(debug) print(str(dat1))
            fit = glm(as.formula(paste0(resp, mod)), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
            vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            print(vals)# ~*~*~
            print(ret)# ~*~*~
            rownames(vals) <- names(coef(fit))[-1]
            colnames(vals) <-  c("estimate", "2.5 %", "97.5 %") # this doesn't work if not a df
            #    vmatch <- match(rownames(vals), rownames(ret)) # col 1 of vals is the row names
            ret[,,r]  <- as.matrix(vals)# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
        }
        return(ret)

    } else {
        for(r in resp_list){
            cols <- c(r, pr_list)
            if (met=="cf_cc") ampDat <- ampDat %>% select(all_of(cols)) %>% filter(!is.na(cam_fate))
            #for(y in seq_along(mods)){
                # needs to be a named list; but if you keep the names, it doesn't drop NA...
            metList <- metLists[resp,,met,mod]
            ampDat <- ampDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?
            inters <- sapply(mod,  function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))[[1]]
            inter <- paste(inters, collapse=".")
            if(met=="passive" & length(inters)!=0)  ampDat[inter] <- NA
            names(metList)[6] <- resp
            names(metList)[7] <- inter
            metList <- metList[!is.na(metList)]

            ## *~*~*~*~*
            cat("\nMETHOD:",met,"-") # ~*~*~*
            print(metList)
            cat("\n")

            frmla <- lapply(form_list[[mod]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
                ## *~*~*~*~*
            if(debug) cat("\nmice formulas:\n") # ~*~*~*
            if(debug) print(frmla)

                # 3: In all(metList) : coercing argument of type 'character' to logical
                # if(is.na(all(metList))) metList=NULL # this is incorrect I think
            if(all(is.na(metList))) metList=NULL
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
                 glm( as.formula( paste0(resp, mod) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
                 ))
            # if(debug) cat("\n>> pooling model fit\n") # ~*~*~*
            pool = summary(mice::pool(fit), "all", conf.int=TRUE, exponentiate=TRUE)
            pool <- pool %>% filter(term %in% vars)
            print(pool)
            print(ret)
            #vmatch <- match(pool[,1], rownames(ret)) # col 1 of vals is the row names
            # ret[vmatch,,y] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %", "fmi")])
            #ret[vmatch,,y] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
            ret[,,r] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
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
    return(ret)
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


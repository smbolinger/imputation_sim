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
    dummy_vars <- model.matrix(~ cam_fate -1, data = dat)
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
        dat4amp <- add_dummy(nd, debug=debug)
        set.seed(seed=seeed)
        #suppressWarnings(amp_out1 <- mice::ampute(dat4amp))
        suppressMessages(suppressWarnings(amp_out1 <- mice::ampute(dat4amp)))

        if(FALSE){
            out <- add_fact(amp_out1$amp,facToNum = T, debug=T)
            # levels(out$cam_fate)
            table(out$cam_fate)
            xtabs(formula = ~ cam_fate + species, data = out)
            xtabs(formula = ~ cam_fate + species, data = nd)
        }
        new_patt <- amp_out1$patterns
        new_patt
        # is_miss <- c("")
        no_miss <- c("obs_int", "fdate", "is_u", "speciesLETE", "speciesCONI")
        # The problem with doing this with indices is that (as just happened), I forget how
        # I did things and then add new variables w/o accounting for them
        # Easier to reference names?
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
        # if(debug & seeed==1){
        #   cat("\n\n>> missingness matrix:\n")
        #   print(miss_patt_mat)
        #   cat("\n>> proportion missing:", new_prop, "\n\n>> frequency of each pattern:", patt_freq)
        # }

        new_order <- c(is_miss, no_miss)
        ### *~*~*~*~* #######
        # if(xdebug) cat("\n\n>> reorder columns:", new_order)
        dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix

        suppressWarnings(amp_out <- mice::ampute(dat4amp, prop = new_prop, patterns = miss_patt_mat, freq = patt_freq))
        if(FALSE){
        out <- add_fact(amp_out$amp,facToNum = T, debug=T)
        # levels(out$cam_fate)
        table(out$cam_fate)
        xtabs(formula = ~ cam_fate + species, data = out)
        }
        ### *~*~*~*~* #######
        # if(xdebug) cat("\n\nCreate new missing values:\n")
        # if(xdebug) print(mice::md.pattern(amp_out$amp, rotate.names = TRUE))

        # 1 = complete; 0 = has missings

        wts <- amp_out$weights 
        wts[1,] <- c(0, 0, 0, 0.8, 0.2, rep(0, 8)) # missing age only - what vars contribute
        # wts[3,] <- c(0, 0, 0, 0.8, 0.2, rep(0, 8))
        ### *~*~*~*~* #######
        # if(xdebug){
        #   cat("\n>> new weights:\n")
        #   print(wts)
        #   cat("\n>> to go with the patterns:\n")
        #   print(miss_patt_mat)
        # }
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
        # if(debug) cat("\n\nCreate more new missing values, with weighted probabilities:\n")
        # if(debug)  print(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE))
        # missing_tab("amp_out_wt", prVars)

        # datList <- ifelse(wt, amp_out_wt, amp_out) # ifelse is not what I need here
        if(wt) {datList <- amp_out_wt} else {datList <- amp_out}

        # str(amp_out_wt)
        # amp_out_wt
        ### *~*~*~*~* #######
        # if(debug) print(str(datList$amp))
        # if (debug) print(class(datList$amp))
        # class(datList)
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
    # if (debug){
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

# full versions with all comments (old code) in fate_GLM_imp_simulation.Rmd
### Nest the loop inside the if statement so you aren't running the if check every loop?
### 
# mkImpSim <- function(fullDat, ampDat, cols, resp, mods, vars, met, m=20, fam=binomial, regMet="brglm_fit", iter=500, passive="both", debug=FALSE, xdebug=FALSE, impplot=FALSE){
mkImpSim <- function(fullDat, ampDat, cols, resp, mods, vars, met, form_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
    ret    <- array(NA, dim=c(length(vars), 3, length(mods)))
    # metLists <- readRDS("met_lists.rds")

    # dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %", "fmi"), names(mods))
    dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %"), names(mods))
    # ret
    # browser()
    if (met == "cc"){
        # ### *~*~*~*~* #######
        # if(debug) cat(" complete-case analysis:\n")
        ampDat <- ampDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
        dat1 <- ampDat[complete.cases(ampDat),]

        # ### *~*~*~*~* #######
        # if(debug) cat("\n\n>> data:\n")
        # if(debug) print(str(dat1))

        for(y in seq_along(mods)){
            fit = glm(as.formula(paste0(resp, mods[y])), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
            # confint(fit)[-1,]

            vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            # vals <- cbind(vals, rep(-1.0, nrow(vals))) # fmi for complete-case analysis
            # also so that numeric stays numeric when you add names?
            rownames(vals) <- names(coef(fit))[-1]
            colnames(vals) <-  c("estimate", "2.5 %", "97.5 %") # this doesn't work if not a df
            vmatch <- match(rownames(vals), rownames(ret)) # col 1 of vals is the row names
            ret[vmatch, , y]  <- as.matrix(vals)# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
            ##################
            # if(FALSE){
            #   ret
            #   vals
            #   rownames(ret)
            #   rownames(vals)
            #   vmatch
            #   y=1
            #   vals
            #   ret[vmatch,,y]
            #   ret
            #   y=y+1
            #   names(coef(fit))
            # }
            ##################
        }

    } else {
        if (met=="cf_cc") ampDat <- ampDat %>% select(all_of(cols)) %>% filter(!is.na(cam_fate))
        for(y in seq_along(mods)){
            # needs to be a named list; but if you keep the names, it doesn't drop NA...
            metList <- metLists[resp,,met,y]
            ampDat <- ampDat %>% select(all_of(cols)) %>% droplevels() # need to select the cols that are relevant for mod?

            # metList <- metList[!is.na(metList), drop=F] 
            # dimnames(metList)
            # if(length(metList)==0) metList=NULL
            # if(metList=='') metList=NULL
            # if(nchar(metList)==0) metList=NULL
            # length(metList)

            # if(xdebug & y==1) cat("\n\n>> data:\n")
            # if(xdebug & y==1) print(str(ampDat))
            # noMet <- rep(NA,)
            # mett <- metL[,met]
            # metList <- if_else(met=="default", c(rep(NA,6)), na.omit(metL[,met]))
            # cat(met,"-", metList,";")
            # frmla=NULL
            # rm(metL)
            # don't need the if statement?
            inters <- sapply(mods[y],  function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))[[1]]
            inter <- paste(inters, collapse=".")
            if(met=="passive" & length(inters)!=0)  ampDat[inter] <- NA

            names(metList)[6] <- resp
            # dimnames(metList)[[2]][7] <- paste(inters, collapse=".")
            names(metList)[7] <- inter
            metList <- metList[!is.na(metList)]

            ## *~*~*~*~*
            # cat("\n",met,"-")
            # print(metList)
            # cat("\n")

            frmla <- lapply(form_list[[y]][[met]], function(x) as.formula(paste(x[[1]], "~", x[[2]])))
            ## *~*~*~*~*
            # if(xdebug) cat("\nmice formulas:\n")
            # if(xdebug) print(frmla)

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

            if(length(imp$loggedEvents > 0)) {
                cat(sprintf("\n*** LOGGED EVENTS FOR METHOD %s*******************************\n", met))
                print(imp$loggedEvents)
            }
            ## *~*~*~*~*
            # if(impplot) visImp(imp) 
            ## *~*~*~*~*
            # if(debug) cat("\n>> fitting model", mods[y])
            fit = with(imp,
                 glm( as.formula( paste0(resp, mods[y]) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
                 ))
            ## *~*~*~*~*
            # if(debug) cat("\n>> pooling model fit\n")
            pool = summary(mice::pool(fit), "all", conf.int=TRUE, exponentiate=TRUE)
            pool <- pool %>% filter(term %in% vars)
            vmatch <- match(pool[,1], rownames(ret)) # col 1 of vals is the row names
            # ret[vmatch,,y] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %", "fmi")])
            ret[vmatch,,y] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %")])
            ## *~*~*~*~*
            # if(xdebug) {
            #    cat("\n\n>> pooled model fit:\n")
            #str(pool)
            #    print(head(pool[,c("term", "estimate", "2.5 %", "97.5 %")]))
            # }
        }
    }
    ## *~*~*~*~*
    # if(xdebug) cat("\n\nret:\n")
    #    cat("\n\nret:\n")
    #    print(head(ret))
    # if (xdebug) print(str(ret))
    return(ret)
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
###### RUN SIMULATION ###################################################################
########################################################################################

##################
# if(FALSE){
#   # nruns <- params$nrun
#   # debug <- params$deb
#   # xdebug <- params$xdeb
#   # ipl <- params$ipl
#   resp <- r
#   fullDat <- dat4sim
#   run <- 1
#   datNA <- ampDat
#   mods <- mods4sim
#   forms <- form_list
#   col_sel <- col_list
# # (nd = fullDat, seeed = run, vars=vars, method = "amp", wt = TRUE, xdebug=xdebug, debug = debug, convFact = TRUE)
#   # nruns=10
#   x <- "default"
#   mod=modList[1]
#   ndat=ndGLM_scl_cc
#   mets <- c("default","pmm", "rf", "cart", "caliber","cc")
#   mets <- met_list
#   # does not include the reference levels:
#   # vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
#   vars <- var_list
#   mLists <- metLists
#   # datNA <- dat
#   # datNA <- ndGLM_scl_cc
#   # debug=TRUE
#   # datNA <- ampDat
#   # seed=82985
# 
#     # run=1
#     # m=1
#     # x=1
#     # x = "rf"
#     # x = met
#     # vals <- imp_sim
#     # run <- 9
# par <- list(nrun=5, 
#             hdir="/home/wodehouse/projects/fate_glm/",
#             deb=TRUE,
#             xdeb=TRUE, # for obsessively checking things - not useful for normal debugging & lots of text output 
#             ipl=FALSE,
#             lin=FALSE,
#             m=5)
# # run =  run + 1
# }

##################
# for some reason, this can't find the global vars in fate_GLM1 when I knit that document.
# so, pass them all as arguments? but they already are...
# arguments: function to make sim data, real data, response, predictors, model, nruns
runSim <- function(fullDat, col_sel, resp, vars, mods,forms, mLists,par,fcToNum=FALSE){
  m=par$m
  nruns=par$nrun
  debug = par$deb
  xdebug=par$xdeb
  # ipl   = par$impPlot
  ipl   = par$ipl
  mets <- unlist(dimnames(mLists)[3])
  
  # res <- array(NA, dim = c(length(vars), length(mets), nruns, 4, length(mods)))
  res <- array(NA, dim = c(length(vars), length(mets), nruns, 3, length(mods)))
  # dimnames(res) <- list(c("pmm", "rf"),
  dimnames(res) <- list(sort(as.character(vars)),
                        # c("pmm", "rf", "cart"),
                        as.character(mets),
                        as.character(1:nruns),
                        # c("estimate", "2.5 %","97.5 %","fmi"),
                        c("estimate", "2.5 %","97.5 %"),
                        names(mods)
  )
  # if(debug) cat("res matrix is formatted as follows:",str(res))
  ## *~*~*~*~*
  # if(xdebug) cat("\n\n>>>>> res matrix is formatted as follows:\n")
  # if(xdebug) print(str(res))
  # datNA <- datNA %>% select(all_of(col_sel))
  ### does this need to be here?
  # fullDat <- fullDat %>% select(all_of(col_sel))
  cat(sprintf(">>> running simulation %s times\n", nruns))
  for(run in 1:nruns){
    cat(run)
    ## *~*~*~*~*
    # if(debug){
    #   cat("\n\n********************************************************************************************")
    #   cat(sprintf("\n//////// RUN %s ////////////////////////////////////////////////////////////////////////////\n", run))
    #   cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
    #   cat("\n\n********************************************************************************************")
    # }
    # ndat = nDat
    # if(missing(datNA)) datNA <- mkSimDat(ndat)
    # if the complete data was passed as an argument, add the NAs; otherwise, should be data with NA
    # if(all(colSums(is.na(datNA)) == 0)) datNA <- mkSimDat(datNA, debug = TRUE, convFact = TRUE)$amp
    # MAKE DATA WITHIN THE RUN_SIM FUNCTION?
    # datNA <- mkSimDat(datNA, debug = debug, convFact = TRUE, fcToNum=fcToNum)$amp
    datNA <- mkSimDat(nd = fullDat, seeed = run, vars=vars, method = "amp", wt = TRUE, xdebug=xdebug, debug = debug, convFact = TRUE)
    ##################
    # if(FALSE ){
    #   run=run+1
    # }
    ##################
    datNA <- datNA$amp
    # datNA <- datNA %>% select(all_of(col_sel))
    # datNA1 <- datNA
    # for(x in seq_along(mets)){
    # for(x in mets){ # this makes it match by name, not just by index
    for(x in seq_along(mets)){ # does matching by index help the trycatch statement?
      ## *~*~*~*~*
      # if (xdebug) cat("\n\n>>>> method:", x)
      skiptoNext <- FALSE
      
      tryCatch(
        expr = {
          vals <- mkImpSim(ampDat=datNA,cols=col_sel,resp=resp, form_list =forms, vars=vars, mods=mods, met=mets[x], debug = debug,m=m, xdebug=xdebug, impplot=ipl)
          # imp <- eval(parse(text=impCall))
          # skiptoNext <- FALSE
          # if(length(imp$loggedEvents > 0)) print(imp$loggedEvents)
        },
        error = function(e){
          cat("\nERROR:", conditionMessage(e), "\n")
          next
          # return(NULL)
          # skiptoNext <- TRUE
          skiptoNext <<- TRUE  # superassignment operator- not sure if necessary
          # imp <- list(imp=NA)
          # ret[,,y]
          # continue()
        }
      )
      if(skiptoNext) next
      # vals <- mkImpSim(ampDat=datNA,cols=col_sel,resp=resp, form_list =forms, vars=vars, mods=mods, met=mets[x], debug = debug,m=m, xdebug=xdebug, impplot=ipl)
      ##################
      # if(FALSE){
      #   vals
      #   rownames(res)
      #   print(vals) # before, had the var names...
      #   vals
      #   vals[,,"m16"]
      #   dim(vals)
      #   dimnames(vals)
      #   dimnames(res)
      #   vals[1,,]
      #   res["nest_age", x, , "estimate","m1"] 
      #   res["nest_age", x, , "estimate","m1"] <- c(-2, -2, -2)
      #   res["nest_age", x, , "2.5 %","m1"] 
      # }
      ## *~*~*~*~*
      # if(xdebug){ 
      #   cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
      #   cat(sprintf("\n>> output of mkImpSim for all models for run %s and method %s:\n", run, x))
      #   str(vals)
      #   }
      # vmatch <- match(vals[,1], rownames(res)) # col 1 of vals is the row names
      ##################
      vmatch <- match(rownames(vals), rownames(res)) # col 1 of vals is the row names
      # vals <- as.matrix(vals[,-1]) # remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
      # dim(res[vmatch, x, r,,])  
      res[vmatch, mets[x], run,,]  <- vals
      ## *~*~*~*~*
      # if(xdebug){ 
      #   cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
      #   cat("\n>> res matrix filled in:\n")
      #   # print(res[vmatch,x,r,])
      #   # print(res[,x,r,])
      #   print(res[,,run,,])
      # }
    }
    ## *~*~*~*~*
    # if(debug & run==nruns){
    #   cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
    #   cat("\n>>>> full res matrix:\n")
    #   print(res[,,,,])
    # }
    # used x bc m already exists, so for testing it was confusing. doesn't matter once fxn works.
    # res[vmatch ,x,,]
    # res[,x,r,]
    # 
    if(run %% 100 == 0){
      # if(run %% 4 == 0){
      begn <- run-100
      endd <- run-0
      nowtime <- format(Sys.time(), "%d%b%H%M")
      fname <- paste(sprintf("out/runs%sto%s_%s.rds", begn, endd, nowtime))
      saveRDS(res[,,begn:endd,,], fname)
    }
    ##################
    # if(FALSE){
    #   run <- run + 1
    #   res[,,c(begn,endd),,]
    #   resss <- readRDS(fname)
    #   resss
    # }
    ##################
    ### print memory usage every 10 runs:
    ## *~*~*~*~*
    # if(run %% 10 == 0) print(gc()) # modulo operator %%
    
    # res[, 1, r,] <- as.matrix(mkImpSim(dat=datNA, resp=resp, vars=vars, mod=mod, met="pmm"))
    # res[, 2, r,] <- as.matrix(mkImpSim(dat=datNA, resp=resp, vars=vars, mod=mod, met="rf"))
    # res[, 3, r,] <- as.matrix(mkImpSim(dat=datNA, resp=resp, vars=vars, mod=mod, met="cart"))
    # res
  }
  return(res)
}

if(FALSE){
  # res1 <- runSim(datNA = dat, resp = resp, vars = vars, mod = mod, mets = mets, nruns = 2, debug=TRUE)
  res1 <- runSim(datNA = dat, resp = resp, vars = vars, mod = mod, mets = mets, nruns = 5)
}

########################################################################################
###### GET AVERAGE BIAS ################################################################
########################################################################################

if(FALSE){
  # dat=simDatNA
  dat=imp_sim
  fitReal <- glm(as.formula(paste0(resp, mod)), data=ndGLM_scl_cc, family=fam, method=regMet, control=brglmControl(maxit=iter) )
  coef(fitReal)
  trueVals <- coef(fitReal)[vars]
  vars
  v="nest_age"
  vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
  v=3
  vars[v]
  dat[vars[v],,,]
  dat[vars[v],,,"2.5 %"]
  dat[vars[v],,,"2.5 %"] < trueVals[vars[v]]
  dat[vars[v],,,"97.5 %"] > trueVals[vars[v]]
  rowMeans(dat[vars[v],,,"2.5 %"] < trueVals[vars[v]])
  biasVals <- c("value","bias", "pctBias", "covRate", "avgWidth", "RMSE","SD")
  # mets <- c("default","pmm", "rf", "cart", "caliber","cc")
  # trueVals <- data.frame(vars=vars, value=)
  # fullDat <- ndGLM_scl_cc
  # impDat  <- res1
  # trueVals
  v <- "cam_fateA"
  # impDat <- res
  resp="is_u"
  resp = "HF_mis"
  # mod=modList[1]
  
  
  #OLD met list:
  mets <- c("default","pmm", "rf", "cart", "caliber","passive", "stratify","cc")# don't need full here?
  # if working from fateGLM_impsim.R:
  resp = r
  fullDat <- dat4sim 
  impDat <- imp_sim
  z=1
  v="nest_age"
  # mod=mods4sim[z]
  mods <- mods4sim
  fam=binomial
  regMet="brglm_fit"
  iter=500
  cols <- col_list
  # why only these vars?
  # vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
  biasVals = bias_names
  mets = met_list
  vars <- var_list
  debug=TRUE
  hdir <- params$hdir
}

parAvg <- function(fullDat, impDat, hdir, resp, vars, mods, regMet="brglm_fit", fam=binomial, iter=500, mets, biasVals, debug=FALSE, xdebug=FALSE){
  # trueVals is optional
  # if(missing(trueVals)){
  # bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(mods) ) )
  for(z in seq_along(mods)){ # loop through mods, save output of each to file
    # mod     <- mods4sim[modnum]
    mod     <- mods4sim[z]
    fitReal <- glm(as.formula(paste0(resp, mod)),
                   # data=ndGLM_scl_cc,
                   data=fullDat,
                   family=fam,
                   method=regMet,
                   control=brglmControl(maxit=iter) )
    # saveRDS(fitReal, sprintf("out/fitReal_%s_m%s.rds", resp, modnum ))
    now_dir = paste0(hdir, "out/", format(Sys.time(), "%d%b"))
    # now_dir
    # if(!dir.exists(paste0(hdir, "out/now_"))) dir.create(paste0(hdir, "out/now_"))
    if(!dir.exists(now_dir)) dir.create(now_dir)
    cat("\n>>> creating directory:", now_dir, "exists?", exists(now_dir))
    # saveRDS(fitReal, sprintf("%s/fitReal_%s_m%s.rds", now_dir, resp, z))
    trueVals <- coef(fitReal)[vars] # the coefs have names associated with them
    # but this ^ creates NAs
    # trueVals <- coef(fitReal)[-1] # the coefs have names associated with them
    # trueVals
    # names(trueVals) <-  names(fitReal[[1]])[vars]
    
    # trueVals
    # vars
    # }
    # mList <- modList[c(1, 8, 16)]
    # for (v in vars){ # v = var names, not indices
    # avg <- list()
    # bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals) ) )
    # bias <- array(NA, dim = c(length(vars), length(mets) +1, length(biasVals) ) )
    # bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals) ) )
    bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(mods) ) )
    
    dimnames(bias) <- list( sort(as.character(vars)),
                            # c("pmm", "rf", "cart"),
                            # c("real",as.character(mets)),
                            as.character(mets),
                            as.character(biasVals),
                            names(mods)
    )
    ## *~*~*~*~*
    # if (xdebug){
    #   cat("\n>> empty bias matrix:\n")
    #   print(bias) 
    #   cat("\n>> impDat:\n")
    #   str(impDat)
    # }
    # for (v in seq_along(vars)){ # v = var names, not indices
    for (v in vars){ # v = var names, not indices
      # avg[v] <- apply(dat[vars[v], , ,], c(1,3), mean, na.rm=TRUE)
      # vars[v]
      # avg <- apply(dat[vars[v], , ,], c(1,3), mean, na.rm=TRUE)
      # avg <- apply(impDat[v, , ,], c(1,3), mean, na.rm=TRUE)
      # which(dimnames(impDat) %in% vars)
      # dim(impDat)[dimnames(impDat)==vars]
      # dimnames(impDat)[[1]] %in% vars
      
      # ddd <- dimnames(impDat)
      # dname <- lapply(dimnames(impDat), function(x) x%in%c(vars,seq(1:nrun)))
      # dname <- lapply(dimnames(impDat), function(x) x%in%biasVals)
      # # which(dimnames(impDat) == vars)
      # # dd <- dim(impDat)
      # # which(dim(impDat))
      # marg <- c(which(dimnames(impDat) %in% vars), which(dimnames(impDat) %in% biasVals))
      # impDat[v, , , ,] 
      avg <- apply(impDat[v, , , ,z],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
      # avg
      sdev <- apply(impDat[v, , , ,],MARGIN=c(1,3),FUN = sd, na.rm=TRUE)
      if(FALSE){
        impDat["nest_age", "default", , "estimate","m1"] 
        impDat["nest_age", , , , ]
        mean(impDat["nest_age", "default", , "estimate","m1"] )
        
        impDat["nest_age",, , "estimate","m1"] 
        impDat["nest_age", "default",1 , "estimate","m1"] <-  -2
        impDat["nest_age", "pmm",2 , "estimate","m1"] <-  -2
        # impDat["nest_age", "default", , "estimate","m1"] <- c(-2, -2, -2)
        impDat[v,,,,]
        avg
        dimnames(impDat)
        sdev
        bias[v,,"value",]
        true
        impDat[v, , , , z]
        impDat[v, "pmm", , , z]
        apply(impDat[v, "pmm", , , z], MARGIN=c(3,4), FUN=mean, na.rm=TRUE)
      }
      # impDat[v,1,,1]
      # impDat[]
      # mean(impDat[v,1,,1]) # this is the mean you are taking, but applied over all methods (dim 2) & parameters (dim 4)
      # for a given variable, these are dim 1 and dim 3.
      # mean(impDat[v,"pmm",,"estimate"]) # This is the same
      
      # mean(impDat[v,"pmm",,"2.5 %"]) # Second value
      # mean(impDat[v,"pmm",,2]) # This is the same
      # only the dimension being averaged over is not specified
      # avg
      # trueVals
      # true <- trueVals[vars[v]]
      true <- trueVals[v] # using v instead of vars[v] matches it by name
      # if(debug) cat("true param value:",trueVals[v])
      # true
      # bias[vars[v], , "bias"] <- avg[,"estimate"] - trueVals[vars[v]]
      # dimnames(avg)
      bias[v, , "value", ] <- avg[,"estimate"]
      bias[v, , "bias", ] <- avg[,"estimate"] - true
      bias
      # bias[v, , "bias"] <- rowMeans(impDat[v,,,"estimate"]) - true #should be equivalent to the line above
      bias[v, , "pctBias",] <- 100 * abs((avg[,"estimate"] - true) / true )
      # bias[vars[v], , "covRate"] <- avg[,"2.5 %"] < true & true > avg[,"97.5 %"]
      # bias[vars[v], , "avgWidth"] <- avg[,"97.5 %"] - avg[,"2.5 %"] 
      bias[v, , "covRate",] <- rowMeans(impDat[v,,,"2.5 %",] < true & true < impDat[v,,,"97.5 %",])
      # bias[v, , "avgWidth"] <- avg[,"97.5 %"] - avg[,"2.5 %"]
      bias[v, , "avgWidth",] <- rowMeans(impDat[v,,,"97.5 %",] - impDat[v,,,"2.5 %",]) # gives the same output
      bias[v, , "RMSE",] <- sqrt((avg[,"estimate"] - true)^2)
      bias[v, , "SD",] <- sdev[,"estimate"]
      
      ## *~*~*~*~*
      # if (xdebug) cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
      # if (xdebug) cat("\n\n>> true param values:\n")
      # if (xdebug) print(trueVals)
      # if (xdebug) cat("\n\n>> bias values:\n")
      # if (xdebug) print(bias)
      # biasfile <- paste0(params$hdir, sprintf("%sout/bias_vals_%s_%s_%s_.rds", hdir,r, names(mods4sim)[z], suffix))
      # biasfile <-  sprintf("%sout/%s/bias_vals_%s_%s_%s_.rds", hdir,now_,r, names(mods4sim)[z], suffix)
      cat(sprintf("\n>>> saving bias values for model %s with response %s for variable %s", 
                  names(mods4sim)[z], r, v))
      biasfile <-  sprintf("%s/bias_vals_%s_%s_%s.rds",now_dir,r, names(mods4sim)[z], suffix)
      # biasfile
      # saveRDS(bias_out, sprintf("out/bias_vals_%s_%s.rds",r, names(mods4sim)[z]))
      saveRDS(bias, biasfile)
      # biasfile1 <- paste0(params$hdir, sprintf("%sout/bias_vals_%s_%s_%s_.csv", hdir, r, names(mods4sim)[z], suffix))
      # biasfile1 <-  sprintf("%sout/%s/bias_vals_%s_%s_%s_.csv", hdir, now_, r, names(mods4sim)[z], suffix)
      biasfile1 <-  sprintf("%s/bias_vals_%s_%s_%s.csv",now_dir,r, names(mods4sim)[z], suffix)
      
      biasdf <- as.data.frame(bias)
      names(trueVals)[ is.na(names(trueVals)) ] <- setdiff(vars, names(trueVals))
      # trueVals
      biasdf <- cbind(trueVals, biasdf)
      # biasdf <- merge
      # trVals <- (trueVals[!is.na(trueVals)])
      
      write.csv(biasdf, file = biasfile1)# write to csv in case script aborts 
      if(FALSE){
        true
        # avg[,"2.5 %"]
        # avg[,"97.5 %"]
        # is value within avg confidence interval VS. how many times is it within the confidence interval out of all the runs
        impDat[v,,,"2.5 %",]
        impDat[v,,,"97.5 %",]
        bias[v,,,]
      }
      # bias
    }
  }
  # return(bias)
  
}

########################################################################################
########################################################################################
test_avg <- function(simDatNA, fitReal){
  # colMeans(simDatNA["cam_fateA", , , ])
  cat(">> all 5 replicates for cam_fateA:\n\n")
  print(simDatNA["cam_fateA", , , ])
  cat(">> & colnames:\n\n")
  print(colnames(simDatNA["cam_fateA", , , ]))
  
  cat("\n\n>> all replicates of estimate for cam_fateA & colnames & rownames:\n\n")
  print(simDatNA["cam_fateA", , ,"estimate"])
  print(colnames(simDatNA["cam_fateA", , , "estimate"]))
  print(rownames(simDatNA["cam_fateA", , , "estimate"]))
  
  cat("\n\n>> averages of estimate from cam_fateA:\n\n")
  print(rowMeans(simDatNA["cam_fateA", , , "estimate"]))
  
  # simDatNA["cam_fateA", , ,]
  cat("\n\n>> all averages from cam_fateA:\n\n")
  simDatAvg <- apply(simDatNA["cam_fateA", , ,], c(1,3), mean, na.rm=TRUE)
  print(simDatAvg)
  
  cat("\n\n>> bias of estimate from cam_fateA:\n\n")
  trueVal <- coef(fitReal)["cam_fateA"]
  print(simDatAvg[,"estimate"] - trueVal)
}

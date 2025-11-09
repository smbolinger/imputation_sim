########################################################################################
###### IMPORT FUNCTIONS ################################################################
########################################################################################

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
  if (debug) cat(">> add dummy variables and remove factors. new columns:\n", names(dat))
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
  if (debug) cat("\n>> converted dummies to factors. all columns:\n", names(dat)) 
  # if (debug) cat("\n>> columns to remove:\n", names(dat)[rem_col])
  if (debug) cat("\n>> columns to remove:\n", rem_col)
  
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
# mkMetList <- function(ageMet = "pmm", fateMet = "pmm", misMet = "pmm", resp=resp){
# mkMetList <- function(met="pmm", resp=resp){
# mkMetList <- function(mType, cols){
# mkMetList <- function(metL, col_sel, debug=FALSE){
mkMetList <- function(met, dat, col_sel, int=NULL, debug=FALSE){
  
  # if (mType=="default"){
  #   met=NULL
  # } else if (mType=="cc"){
  #   
  # }
  # names(dat4imp[])
  # u<-names(colSums(is.na(dat4imp)))
  # metList <- rep("", 6)
  metList <- rep("",length(col_sel))
  catList <- c("HF_mis", "cam_fate", "species", "is_u")
  # metList
  names(metList) <- col_sel
  # metList["cam_fate"] <- fateMet
  # metList["nest_age"] <- ageMet
  # metList["cam_fate"] <- met
  # metList["nest_age"] <- met
  
  colNA <- names(dat)[colSums(is.na(dat)) > 0]
  is_num <- names(dat)[sapply(dat, is.numeric)]
  # is_bin <- names(dat)[length(sapply(dat, levels))==2]
  is_bin <- names(dat)[sapply(dat, function(x) length(levels(x)))==2]
  is_cat <- names(dat)[sapply(dat, function(x) length(levels(x)))>2]
  # is_num <- names(dat)[where(is.numeric(dat))]
  # which(is.numeric(
  
  # if (resp %in% colNA) metList[resp] <- misMet
  # metList[resp] <- ifelse(resp %in% colNA, misMet, "")
  # metList[resp] <- ifelse(resp %in% colNA, met, "")

  # if(met=="caliber"){
  #   for(c in col_sel){
  #     catList <- c("HF_mis", "cam_fate", "species", "is_u")
  #     # if(c %in% catList ) metList[c] = "rfcat" 
  #     if(c %in% catList ) metList[c] = "rfcat"
  #   }
  # }
  # met2 <- ""
  for(c in col_sel){
    if(met=="caliber"){
      # catList <- c("HF_mis", "cam_fate", "species", "is_u")
      if(c %in% catList ) met1 = "rfcat" else met1 = "rfcont"
    }else if(met %in% c("default.int", "pmm.int")){
      met1 <- str_extract(met, pattern = "\\w+(?=.)") # everything up until the period
    # }else if(met == "passive.int"){
    }else if(met == "passive" | met=="stratify"){
      # met1 <- ""
      # met1 <- ifelse(is.numeric())
      met1 <- case_when(c %in% is_num ~ "pmm",
                         c %in% is_bin ~ "logreg",
                         c %in% is_cat ~ "polyreg")
      # met2 <- paste("~I(", int,")")
    } else {met1 <- met}
    
    # print(met)
  # for(c in cols){
    # metList[c]<- ifelse(c %in% colNA, metL[c], "")
    if(grepl(".", c, fixed=TRUE)){
      metList[c] = paste("~I(", int,")") 
    } else {
      metList[c]<- ifelse(c %in% colNA, met1, "")
    }
    # metList[c] <- ifelse(grepl(".", c, fixed=TRUE), paste("~I(", int,")"), metList[c])
    # metList[c]<- case_when(c %in% colNA  met1, "")
    # metList[c]<- ifelse(grepl(":", c,fixed = TRUE), met2, "")
    # if(met=="passive" & grepl(":", c, fixed=TRUE)) metList[c] = met2 
    # if(grepl(":", c, fixed=TRUE)) metList[c] = met2 
    # grepl(".", c, fixed=TRUE)
    # this is not being evaluated for some reason??
    
  }
  if(debug){
    cat(sprintf("methods for each variable (%s total):\n", length(metList)))
    print( metList)
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
}

# mkSimDat <- function(nd,col_sel, method="amp", wt=TRUE, debug=FALSE, convFact=FALSE){
mkSimDat <- function(nd, vars, facToNum=FALSE, method="amp", wt=TRUE, debug=FALSE, convFact=FALSE){
  if(method=="amp"){
    dat4amp <- add_dummy(nd, debug=debug)
    # dat4amp <- as.numeric(dat4amp)
    # dat4amp <- apply(dat4amp,MARGIN=2, FUN=function(x) as.numeric(x))
    
    # dat4amp <-nd %>% select(all_of(col_sel)) %>% add_dummy(debug=TRUE)
    # amp_out <- mice::ampute(dat4amp, run = FALSE)
    suppressWarnings(amp_out1 <- mice::ampute(dat4amp))
    
    new_patt <- amp_out1$patterns
    new_patt
    # is_miss <- c("")
    no_miss <- c("obs_int", "fdate", "is_u", "speciesLETE", "speciesCONI")
    # The problem with doing this with indices is that (as just happened), I forget how
    # I did things and then add new variables w/o accounting for them
    # Easier to reference names?
    # no_mis_c
    # is_mis <- vars[-c(no_mis)]
    # is_miss <- vars[!vars %in% no_mis]
    is_miss <- names(new_patt)[!names(new_patt) %in% no_miss]
    is_miss
    # no_miss <- names(new_patt)[c(1, 3, 5:7)]
    # is_miss <- names(new_patt)[-c(1,3,5:7)] # Reve Coffee Lab in BR
    
    # the order: missing age only; missing fate/HF_mis only; missing both
    miss_pat <- list(c(0, rep(1,7)), c(1,rep(0,7)), c(rep(0,8)) )
    
    miss_pat
    
    m <- do.call(rbind,miss_pat) # create the matrix of the missingness patterns 
    # rep(c(rep(1,5)),3)
    p <- do.call(rbind, rep(list(c(rep(1,5))), 3)) # create additional columns for the vars not missing values
    
    new_prop <- 0.2
    miss_patt_mat <- cbind(m,p)
    colnames(miss_patt_mat) <- c(is_miss, no_miss)
    # rownames(miss_patt_mat) <- c("nest_age", )
    patt_freq <- c(0.45,0.45,0.1) # missing: age only, fate only, both
    if(debug){
      cat("\n\n>> missingness matrix:\n")
      print(miss_patt_mat)
      cat("\n>> proportion missing:", new_prop, "\n\n>> frequency of each pattern:", patt_freq)
    }
    
    new_order <- c(is_miss, no_miss)
    if(debug) cat("\n\n>> reorder columns:", new_order)
    dat4amp <- dat4amp %>% select(all_of(new_order)) # reorder the columns to match the matrix
    
    suppressWarnings(amp_out <- mice::ampute(dat4amp, prop = new_prop, patterns = miss_patt_mat, freq = patt_freq))
    if(debug) cat("\n\nCreate new missing values:\n")
    if(debug) print(mice::md.pattern(amp_out$amp, rotate.names = TRUE))
    # print(missing_tab("amp_out", prVars))
    # missing_tab("amp_out", prVars)
    # 1 = complete; 0 = has missings
    
    wts <- amp_out$weights 
    wts[1,] <- c(0, 0, 0, 0.8, 0.2, rep(0, 8)) # missing age only - what vars contribute
    # wts[3,] <- c(0, 0, 0, 0.8, 0.2, rep(0, 8))
    if(debug){
      cat("\n>> new weights:\n")
      print(wts)
      cat("\n>> to go with the patterns:\n")
      print(miss_patt_mat)
    }
    # what exactly are the weights doing?
    suppressWarnings(amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = miss_patt_mat, freq = patt_freq,weights = wts))
    
    if(debug) cat("\n\nCreate more new missing values, with weighted probabilities:\n")
    if(debug)  print(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE))
    # missing_tab("amp_out_wt", prVars)
    
    # datList <- ifelse(wt, amp_out_wt, amp_out) # ifelse is not what I need here
    if(wt) {datList <- amp_out_wt} else {datList <- amp_out}
    
    # str(amp_out_wt)
    # amp_out_wt
    if(debug) print(str(datList$amp))
    if (debug) print(class(datList$amp))
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
    
    if (debug){
      cat("new number missing nest_age:", sum(is.na(nd$nest_age)),
          "new percent missing nest_age:", sum(is.na(nd$nest_age)) / nrow(nd),
          "original percent missing:", sum(is.na(ndGLM_scl_all$nest_age)) / nrow(ndGLM_scl_all) )
      cat("\nnew number missing cam_fate:", sum(is.na(nd$cam_fate)),
          "new percent missing cam_fate:", sum(is.na(nd$cam_fate)) / nrow(nd),
          "original percent missing:", sum(is.na(ndGLM_scl_all$cam_fate)) / nrow(ndGLM_scl_all) )
      
      return(nd)
    }
  }
}
########################################################################################
##### IMPUTE/FIT/POOL SIM DATA ##########################################################
########################################################################################
if(debugging){
  # dat=mkSimDat(ndGLM_scl_cc, convFact = TRUE)$amp
  # ampDat=mkSimDat(ndGLM_scl_cc, convFact = TRUE)$amp
  aDat = sim_dat$amp
  # aDat = sim_dat$amp[1:20,]
  # resp="is_u"
  resp = "HF_mis"
  # resp="isu"
  # resp = "HFmis"
  mod=modList[1]
  mods = mods4sim
  m=20
  # met="rf"
  met="default"
  met="cc"
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

# mkImpSim <- function(ampDat, resp, mod, vars, m=30, met="rf", fam=binomial, regMet="brglm_fit", iter=500, seed=NULL, passive="both", debug=FALSE){
# mkImpSim <- function(ampDat, resp, mod, vars, m=30, metList=rep("", 6),  fam=binomial, regMet="brglm_fit", iter=500, seed=NULL, passive="both", debug=FALSE){
# mkImpSim <- function(ampDat, cols, resp, mod, vars, met, m=30, fam=binomial, regMet="brglm_fit", iter=500, seed=NULL, passive="both", debug=FALSE){
# mkImpSim <- function(fullDat, ampDat, cols, resp, mod, vars, met, m=20, fam=binomial, regMet="brglm_fit", iter=500, passive="both", debug=FALSE){
# mkImpSim <- function(fullDat, ampDat, cols, resp, mods, vars, met, fcToNum=FALSE, m=20, fam=binomial, regMet="brglm_fit", iter=500, passive="both", debug=FALSE){
# mkImpSim <- function(fullDat, aDat, cols, resp, mods, vars, met, m=20, fam=binomial, regMet="brglm_fit", iter=500, passive="both", debug=FALSE, xdebug=FALSE){

### Nest the loop inside the if statement so you aren't running the if check every loop?
### 
mkImpSim <- function(fullDat, aDat, cols, resp, mods, vars, met, m=20, fam=binomial, regMet="brglm_fit", iter=500, passive="both", debug=FALSE, xdebug=FALSE, impplot=FALSE){
  ret    <- array(NA, dim=c(length(vars), 4, length(mods)))
  
  dimnames(ret) <- list(vars, c( "estimate", "2.5 %", "97.5 %", "fmi"), names(mods))
  # ret
  # browser()
  if (met == "cc"){
    
    # dat <- ampDat[!is.na(ampDat),]
    if(debug) cat(" complete-case analysis:\n")
    # dat1 <- ifelse(met=="cc",ampDat[complete.cases(ampDat),],fullDat) # I forget why this makes a list
    ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
    dat1 <- ampDat[complete.cases(ampDat),]
    
    if(debug) cat("\n\n>> data:\n")
    if(debug) print(str(dat1))
    for(y in seq_along(mods)){
      fit = glm(as.formula(paste0(resp, mods[y])), data=dat1, family=fam, method=regMet, control=brglmControl(maxit=iter))
      # confint(fit)[-1,]

      vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
      vals <- cbind(vals, rep(-1.0, nrow(vals)))
      # vals <- data.frame(vals) # why does it need to be a dataframe? oh, maybe bc of filter?
      # also so that numeric stays numeric when you add names?
      rownames(vals) <- names(coef(fit))[-1]
      colnames(vals) <-  c("estimate", "2.5 %", "97.5 %", "fmi") # this doesn't work if not a df
      vmatch <- match(rownames(vals), rownames(ret)) # col 1 of vals is the row names
      # vmatch <- na.omit(match(rownames(vals), rownames(ret))) # col 1 of vals is the row names
      ret[vmatch, , y]  <- as.matrix(vals)# remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
      if(FALSE){
        ret
        vals
        rownames(ret)
        rownames(vals)
        vmatch
        y=1
        vals
        ret[vmatch,,y]
        ret
        y=y+1
        names(coef(fit))
      }
    }
  } else {
    # ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
    for(y in seq_along(mods)){
      #ampDat resets each time:
      ampDat <- aDat %>% select(all_of(cols)) # need to select the cols that are relevant for mod?
    if(xdebug & y==1) cat("\n\n>> data:\n")
    if(xdebug & y==1) print(str(ampDat))
      if(grepl("*", mods[y], fixed=TRUE)){ # the defaults are for this model
        # don't know if as.vector is necessary
        inters <- sapply(mods[y], 
                         # function(x) as.vector(unlist( str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))))
                         # function(x) unlist( str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+")))
                         function(x)  str_extract_all(x, "\\w+(?=\\s\\*)|(?<=\\*\\s)\\w+"))
        
        inters <- inters[[1]]
        inter <- paste(inters, collapse=".")
        inter2 <- paste(inters, collapse=":")
        # inters2 <- str_replace(inters, "species", "speciesCONI")
        # metList 
      }
      if (met=="default"){
        # if(debug) cat("\n\n>> default method\n")
        if(debug & y==1) cat(" default method:\n")
        metList <- NULL
        if(debug) cat("\n>>metList=", metList)
        
      } else if(met=="passive"){
          if(debug & y==1) cat(" passive imputation:\n")
          addInt <- paste(" +", inter2, sep=" ")
          frmla <- vector("list", ncol(ampDat))
          names(frmla) <- names(ampDat)
          # n=1
          ### I'll just make my own formulas instead of trying to update those....
          ### updating them is such a headache
          # for(n in seq_along(names(ampDat))){
          if(xdebug) cat("\nmice formulas:\n")
          for(n in seq_along(names(ampDat))){
          # for(n in names(ampDat)){
            nm <- names(ampDat)[n]
            rhside <- paste(names(ampDat)[-c(n)], collapse="+")
            frmla[[n]] <- as.formula(paste(nm, rhside, sep="~"))
            if(xdebug) print(frmla[[n]])
          }
          # # for(form in seq_along(inImp$formulas)){
          # #   frmla[[form]] <- ifelse(names(inImp$formulas)[form] %in% inters,
          # #                           as.character(inImp$formulas[form]),
          # #                           paste(as.character(inImp$formulas[form]), addInt, sep=" "))
          # for(form in seq_along(ff)){
          #   frmla[[form]] <- ifelse(names(ff)[form] %in% inters,
          #                           ff[form],
          #                           paste(as.character(ff[form]), addInt, sep=" "))
          # }
          # names(frmla) <- names(ff)
          # if(xdebug) cat("\nmice formulas:\n")
          # if(xdebug) print(frmla)
          # inImp$formulas <- sapply(inImp$formulas, paste(resp, mod))
          
          # pre <- inImp$pred
          # # pre[c(inters), inter]
          # pre[c(inters),inter] <- 0
          metList <- mkMetList(met=met, dat=ampDat, int=paste(inters, collapse="*"),col_sel = names(ampDat), debug=debug)
          # imp <- mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)
          
      } else {
        if (debug & y==1) cat(sprintf(" %s method:\n", met))
        metList <- mkMetList(met=met, dat=ampDat, col_sel=names(ampDat), debug=debug)
        if(FALSE){
          imp <- mice::mice(ampDat, method=NULL, m=m, print=FALSE, seed=13)
          imp_comp <- mice::complete(imp,action="long") 
          imp2 <- mice::mice(ampDat, m=m, print=FALSE, seed=13)
          imp_comp2 <- mice::complete(imp,action="long") 
        }
      }
#######################################################################################
      # if(debug) print(str(imp$imp))
      # if(debug) print(plot(imp))
      # create one long dataframe of all imputed datasets stacked, so you have an "impVar" for ppsfmi
      # imp_comp <- mice::complete(imp,action="long") 
      # if(debug){
      #   imp40 <- mice.mids(imp, maxit=35, print=F)
      #   print(plot(imp40))
      # }
      # imp$blocks
    # }
    
        # imp_comp <- imp_comp %>% 
        #   mutate( cam_fate = case_match(as.character(cam_fate),
        #                                              "H"  ~ 0,
        #                                              "F"  ~ 1,
        #                                              "D"  ~ 2,
        #                                              "S"  ~ 3,
        #                                              "A"  ~ 4,
        #                                              "Hu" ~ 5,
        #                                              "Ca" ~ 6,
        #                                              "U"  ~ 7
        #   ),
        #   species = ifelse(species=="LETE", 0, 1)
        #   )
        # levels(imp_comp$cam_fate)
        # levels(inters[2])
        # psdat <- psfmi::lbpmilr
        
        # formula here requires categorical vars to be wrapped in factor([var])
        # pool <- psfmi::psfmi_lr(data=imp_comp, nimp=m, impvar=".imp", formula=mods[y])
        # this function thinks cam_fate has 2 levels, for some asinine reason??
        # mod_form
      #   pool <- psfmi::psfmi_lr( data=imp_comp, 
      #                           # formula=mod_form, 
      #                           nimp=m, 
      #                           impvar=".imp", 
      #                           Outcome=resp,
      #                           # predictors=c("nest_age", "obs_int", "fdate"), 
      #                           # predictors=c("species","nest_age", "obs_int", "fdate", "cam_fate"),
      #                           predictors=c("species","nest_age", "obs_int", "fdate"),
      #                           # predictors=c("spcs","nest_age", "obs_int", "fdate"),
      #                           cat.predictors = c("cam_fate"),
      #                           # cat.predictors = c("cfate"),
      #                           int.predictors = inter,
      #                           method="D1")
      # if(FALSE){
      #   pool$multiparm
      # }
#######################################################################################
      if(debug) cat(sprintf("\n>> making %s imputed datasets", m))
      if(met=="stratify"){
        # imp <- impStrat("ampDat", met = metList, col_sel = cols)
        cat("stratified")
        imp <- impStrat(ampDat, met = metList, col_sel = cols)
        imp$formulas
        # I don't think just passing the name really works (from the command line or 
        # another environment). can't find actual object
      } else if (met=="passive"){
        # imp <- mice::mice(ampDat, method=metList, m=m, pred=pre, maxit=10, print=FALSE)
        cat("passive")
        imp <- mice::mice(ampDat, method=metList, m=m, formulas=frmla, maxit=10, print=FALSE)
      } else {
        cat(met)
        imp <- mice::mice(ampDat, method=metList, m=m, print=FALSE)
      }
      if(!is.null(imp$loggedEvents)) print(imp$loggedEvents)
      # if(xdebug){
      if(impplot){
        imp40 <- mice.mids(imp, maxit=35, print=F)
        # plot(imp40)
        mice:::plot.mids(imp40)
      }
      # if(debug) print(plot(imp))
      imp_comp <- mice::complete(imp,action="long") 
      if(debug) cat("\n>> fitting model:", mods[y])
      fit = with(imp,
                 glm( as.formula( paste0(resp, mods[y]) ),
                      family=fam,
                      method = regMet,
                      control=brglmControl(maxit=iter)
      ))
      if(debug) cat("\n>> pooling model fit\n")
      # other pooling & related functions worth looking into:
      # psfmi::psfmi_lr
      # Hmisc:aRegImpute
      pool = summary(mice::pool(fit), "all", conf.int=TRUE)
      pool <- pool %>% filter(term %in% vars)
      vmatch <- match(pool[,1], rownames(ret)) # col 1 of vals is the row names
      ret[vmatch,,y] <- as.matrix(pool[, c("estimate", "2.5 %", "97.5 %", "fmi")])
      if(FALSE){
        rownames(ret)
        pool[,1]
        fit$analyses[[1]]
        pool
        y=1
        ret
        y=y+1
      }
      # if(debug) {
      #   cat("\n\n>> pooled model fit:\n")
      #   str(pool)
      # }
      # ret[m,] = pool[, c("term", "estimate", "2.5 %", "97.5 %", "fmi")] # maybe do need term?
      # for(y in seq_along(mods)){
      #   
      # }
      
      
      
    }
      }
  
  # fits <- list()
  # if (passive=="both" | passive=="no"){
  #   fit = with(imp,
  #              glm( as.formula( paste0(resp, mod) ),
  #                   family=fam,
  #                   method = regMet,
  #                   control=brglmControl(maxit=iter)
  #   ))
  #   append(fits, fit)
  # }
  # if (passive=="both" | passive=="yes"){
  #   fit = with(imp,
  #              glm( as.formula( paste0(resp, mod) ),
  #                   family=fam,
  #                   method = regMet,
  #                   control=brglmControl(maxit=iter)
  #   ))
  #   append(fits, fit)
  # }
  # if (debug) cat("length of fits:", length(fits))
  # pool[[m]] = summary(mice::pool(fit), "all", conf.int=TRUE)
  # pool = summary(mice::pool(fit), "all", conf.int=TRUE)
  # why did I make impV?
  # impV = as.character(pool$term[pool$term %in% vars]) # don't need the levels
  # impV
  # pool <- pool[order(as.character(pool$term)),]
  # if(debug) {
  #   cat("pooled model fit:\n")
  #   print(pool)
  #   str(pool)
  # }
  # pool <- pool %>% filter(term %in% impV)
  # pool <- pool %>% filter(term %in% vars)
  # pool
  # pool2
  # pool = pool[, c("term", "estimate", "2.5 %", "97.5 %", "fmi")] # maybe do need term?
  # don't need term? - term is the name in the matrix
  # pool = pool[, c("estimate", "2.5 %", "97.5 %", "fmi")]
  # pool
  if(debug) cat("\n\nret:\n")
  if (debug) print(str(ret))
  return(ret)
}

########################################################################################
###### RUN SIMULATION ###################################################################
########################################################################################
if(FALSE){
  nruns=10
  mod=modList[1]
  ndat=ndGLM_scl_cc
  mets <- c("default","pmm", "rf", "cart", "caliber","cc")
  mets <- met_list
  # does not include the reference levels:
  # vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
  vars <- var_list
  # datNA <- dat
  # datNA <- ndGLM_scl_cc
  # debug=TRUE
  datNA <- ampDat
  # seed=82985
  nruns <- params$nrun
  debug <- params$deb
  mods <- mods4sim

    r=1
    # m=1
    # x=1
    x = "rf"
    x = met
    vals <- imp_sim
}

# for some reason, this can't find the global vars in fate_GLM1 when I knit that document.
# so, pass them all as arguments? but they already are...
# arguments: function to make sim data, real data, response, predictors, model, nruns
# runSim <- function(datNA, col_sel, resp, vars, mod, mets, nruns=100, seed=NULL, passive="both", debug=FALSE){
# runSim <- function(dat, datNA, col_sel, resp, vars, mod, mets, nruns=100, passive="both", debug=FALSE){
# runSim <- function(datNA, col_sel, resp, vars, mods, mets,fcToNum=FALSE,m=25, nruns=100, passive="both", debug=FALSE, xdebug=FALSE){
runSim <- function(datNA, col_sel, resp, vars, mods, mets,par,fcToNum=FALSE){
   m=par$m
   nruns=par$nrun
   debug = par$deb
   xdebug=par$xdeb
   # ipl   = par$impPlot
   ipl   = par$ipl
  
  res <- array(NA, dim = c(length(vars), length(mets), nruns, 4, length(mods)))
  # dimnames(res) <- list(c("pmm", "rf"),
  dimnames(res) <- list(sort(as.character(vars)),
                        # c("pmm", "rf", "cart"),
                        as.character(mets),
                        as.character(1:nruns),
                        c("estimate", "2.5 %","97.5 %","fmi"),
                        names(mods)
                        )
  # if(debug) cat("res matrix is formatted as follows:",str(res))
  if(debug) cat("\n\n>>>> res matrix is formatted as follows:\n")
  if(debug) print(str(res))
  datNA <- datNA %>% select(all_of(col_sel))
  cat(sprintf(">>> running simulation %s times", nruns))
  for(r in 1:nruns){
    # ndat = nDat
    # if(missing(datNA)) datNA <- mkSimDat(ndat)
    # if the complete data was passed as an argument, add the NAs; otherwise, should be data with NA
    # if(all(colSums(is.na(datNA)) == 0)) datNA <- mkSimDat(datNA, debug = TRUE, convFact = TRUE)$amp
    if(all(colSums(is.na(datNA)) == 0)) datNA <- mkSimDat(datNA, debug = debug, convFact = TRUE, fcToNum=fcToNum)$amp
    # datNA1 <- datNA
    # for(x in seq_along(mets)){
    for(x in mets){ # this makes it match by name, not just by index
      if (xdebug) cat("\n\n>>>> method:", x)
      # could try using map to make sure the vars match up?
      # vals <- as.matrix(mkImpSim(dat=datNA, resp=resp, vars=vars, mod=mod, met=mets[m]))
    # if you include term in the mkImpSim output, as.matrix coerces all columns to character. 
    # need all num, so remove "term" before converting
      # vals <- as.matrix(mkImpSim(ampDat=datNA, resp=resp, vars=vars, mod=mod, met=x))
      # if (x == "default"){
      #   vals <- mkImpSim(ampDat=datNA, resp=resp, vars=vars, mod=mod, metList=c("none"), seed=seed, debug = debug)
      #   
      # }
      # metL <- mkMetList(met = x, cols = col_sel)
      # vals <- mkImpSim(ampDat=datNA, resp=resp, vars=vars, mod=mod, metList=metL, seed=seed, debug = debug)
      # vals <- mkImpSim(ampDat=datNA, cols=col_sel,resp=resp, vars=vars, mod=mod, met=x, seed=seed, debug = debug)
      # vals <- mkImpSim(fullDat=dat,ampDat=datNA, cols=col_sel,resp=resp, vars=vars, mod=mod, met=x, debug = debug)
      # vals <- mkImpSim(ampDat=datNA,cols=col_sel,resp=resp, vars=vars, mod=mod, met=x, debug = debug)
      # vals <- mkImpSim(ampDat=datNA,cols=col_sel,resp=resp, vars=vars, mods=mods, met=x, debug = debug,m = 75)
      vals <- mkImpSim(aDat=datNA,cols=col_sel,resp=resp, vars=vars, mods=mods, met=x, debug = debug,m=m, xdebug=xdebug, impplot=ipl)
      if(FALSE){
        vals
        rownames(res)
        print(vals) # before, had the var names...
        vals
        vals[,,"m16"]
        dim(vals)
        dimnames(vals)
        dimnames(res)
        vals[1,,]
        res["nest_age", x, , "estimate","m1"] 
        res["nest_age", x, , "estimate","m1"] <- c(-2, -2, -2)
        res["nest_age", x, , "2.5 %","m1"] 
      }
      if(xdebug){ 
        cat(sprintf("\n>> output of mkImpSim for all models for run %s and method %s:\n", r, x))
        str(vals)
        }
      # vmatch <- match(vals[,1], rownames(res)) # col 1 of vals is the row names
      vmatch <- match(rownames(vals), rownames(res)) # col 1 of vals is the row names
      # vals <- as.matrix(vals[,-1]) # remove chr column AFTER match so others aren't coerced to chr when you convert to matrix
      # dim(res[vmatch, x, r,,])  
      res[vmatch, x, r,,]  <- vals
      if(xdebug){ 
        cat("\n>> res matrix filled in:\n")
        # print(res[vmatch,x,r,])
        # print(res[,x,r,])
        print(res[,,r,,])
        }
      # used x bc m already exists, so for testing it was confusing. doesn't matter once fxn works.
      # res[vmatch ,x,,]
      # res[,x,r,]
    }
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
  cols <- col_sel
  # why only these vars?
  # vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
  biasVals = bias_names
  mets = met_list
  vars <- var_list
  debug=TRUE
  hdir <- params$hdir
}

# parAvg <- function(dat, vars, mets, biasVals, trueVals){
# parAvg <- function(dat, resp, vars, mod, mets, biasVals, trueVals){
# parAvg <- function(fullDat, impDat, resp, vars, mod, regMet="brglm_fit", fam=binomial, iter=500, mets, biasVals, debug=FALSE){
# parAvg <- function(fullDat, impDat, resp, vars, modnum, regMet="brglm_fit", fam=binomial, iter=500, mets, biasVals, debug=FALSE){
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
    saveRDS(fitReal, sprintf("%s/fitReal_%s_m%s.rds", now_dir, resp, z))
    trueVals <- coef(fitReal)[vars] # the coefs have names associated with them
    # trueVals
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
    if (xdebug){
      cat("\n>> empty bias matrix:\n")
      print(bias) 
      cat("\n>> impDat:\n")
      str(impDat)
    }
    # for (v in seq_along(vars)){ # v = var names, not indices
    for (v in vars){ # v = var names, not indices
      # avg[v] <- apply(dat[vars[v], , ,], c(1,3), mean, na.rm=TRUE)
      # vars[v]
      # avg <- apply(dat[vars[v], , ,], c(1,3), mean, na.rm=TRUE)
      # avg <- apply(impDat[v, , ,], c(1,3), mean, na.rm=TRUE)
      avg <- apply(impDat[v, , , ,],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
      # avg
      sdev <- apply(impDat[v, , , ,],MARGIN=c(1,3),FUN = sd, na.rm=TRUE)
      if(FALSE){
        impDat["nest_age", "default", , "estimate","m1"] 
        impDat["nest_age",, , "estimate","m1"] 
        impDat["nest_age", "default",1 , "estimate","m1"] <-  -2
        impDat["nest_age", "pmm",2 , "estimate","m1"] <-  -2
        # impDat["nest_age", "default", , "estimate","m1"] <- c(-2, -2, -2)
        impDat[v,,,,]
        avg
        sdev
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
      bias[v, , "value", ] <- avg[,"estimate"]
      bias[v, , "bias", ] <- avg[,"estimate"] - true
      # bias
      # bias[v, , "bias"] <- rowMeans(impDat[v,,,"estimate"]) - true #should be equivalent to the line above
      bias[v, , "pctBias",] <- 100 * abs((avg[,"estimate"] - true) / true )
      # bias[vars[v], , "covRate"] <- avg[,"2.5 %"] < true & true > avg[,"97.5 %"]
      # bias[vars[v], , "avgWidth"] <- avg[,"97.5 %"] - avg[,"2.5 %"] 
      bias[v, , "covRate",] <- rowMeans(impDat[v,,,"2.5 %",] < true & true < impDat[v,,,"97.5 %",])
      # bias[v, , "avgWidth"] <- avg[,"97.5 %"] - avg[,"2.5 %"]
      bias[v, , "avgWidth",] <- rowMeans(impDat[v,,,"97.5 %",] - impDat[v,,,"2.5 %",]) # gives the same output
      bias[v, , "RMSE",] <- sqrt((avg[,"estimate"] - true)^2)
      bias[v, , "SD",] <- sdev[,"estimate"]
      if(debug) cat("true param values:")
      if(debug) print(trueVals)
      if(debug) cat("\n\nbias values:\n")
      if (debug) print(bias)
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
      write.csv(bias, file = biasfile1)# write to csv in case script aborts 
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

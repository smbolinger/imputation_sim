########################################################################################
###### IMPORT FUNCTIONS ################################################################
########################################################################################

# source("missingness_tab.R")

########################################################################################
###### FACTORS --> DUMMIES #############################################################
########################################################################################
add_dummy <- function(dat, debug=FALSE){
  dummy_vars <- model.matrix(~ cam_fate -1, data = dat)
  dummy_vars
  
  dat$speciesLETE <- ifelse(dat$species=="LETE", 1, 0)
  dat$speciesCONI <- ifelse(dat$species=="CONI", 1, 0)
  
  dat<- dat[, c("obs_int", "nest_age", "fdate", "HF_mis", "is_u", "speciesLETE", "speciesCONI")]
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
add_fact <- function(dat, debug=FALSE){ 
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
           species = if_else(speciesCONI==1, "CONI", "LETE")) %>%
    mutate(across(c(cam_fate, species), as.factor))
  
  dat$cam_fate <- relevel(dat$cam_fate, "H")
  
  rem_col <- c(3:8, 12,13)
  # names(dat)[, c(3:8, 12,13)]
  if (debug) cat("\n>> converted dummies to factors. all columns:\n", names(dat)) 
  if (debug) cat("\n>> columns to remove:\n", names(dat)[rem_col])
  
  
  dat <- dat[,-rem_col]
  return(dat)
  # there MUST be a shorter way to do this, using regex?
  # fate_letter <- str_extract_all(dat$)
  # dat <- dat %>% mutate(cam_fate = case_when())
  # dat$cam_fate <- case
}

########################################################################################
##### MAKE SIMULATED DATA ###############################################################
########################################################################################
if (FALSE){
  nd <- ndGLM_scl_cc
  wt <- TRUE
  convFact <- TRUE
}

mkSimDat <- function(nd, method="amp", wt=TRUE, debug=FALSE, convFact=FALSE){
  if(method=="amp"){
    dat4amp <- add_dummy(nd, debug=TRUE)
    # amp_out <- mice::ampute(dat4amp, run = FALSE)
    amp_out1 <- mice::ampute(dat4amp)
    
    new_patt <- amp_out1$patterns
    no_miss <- names(new_patt)[c(1, 3, 5:7)]
    is_miss <- names(new_patt)[-c(1,3,5:7)]
    
    # the order: missing age only; missing fate/HF_mis only; missing both
    miss_pat <- list(c(0, rep(1,7)), c(1,rep(0,7)), c(rep(0,8)) )
    
    m <- do.call(rbind,miss_pat) # create the matrix of the missingness patterns 
    # rep(c(rep(1,5)),3)
    p <- do.call(rbind, rep(list(c(rep(1,5))), 3)) # create additional columns for the vars not missing values
    
    new_prop <- 0.2
    miss_patt_mat <- cbind(m,p)
    patt_freq <- c(0.45,0.45,0.1) # missing: age only, fate only, both
    if(debug){
      cat("\n\n>> missingness matrix:\n")
      print(miss_patt_mat)
      cat("\n>> proportion missing:", new_prop, "\n\n>> frequency of each pattern:", patt_freq)
    }
    
    new_order <- c(is_miss, no_miss)
    if(debug) cat("\n\n>> reorder columns:", new_order)
    dat4amp <- dat4amp %>% select(new_order) # reorder the columns to match the matrix
    
    amp_out <- mice::ampute(dat4amp, prop = new_prop, patterns = miss_patt_mat, freq = patt_freq)
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
    amp_out_wt <- mice::ampute(dat4amp, prop = new_prop, patterns = miss_patt_mat, freq = patt_freq,weights = wts)
    
    if(debug) cat("\n\nCreate more new missing values, with weighted probabilities:\n")
    if(debug)  print(mice::md.pattern(amp_out_wt$amp, rotate.names = TRUE))
    # missing_tab("amp_out_wt", prVars)
    
    # datList <- ifelse(wt, amp_out_wt, amp_out) # ifelse is not what I need here
    if(wt) {datList <- amp_out_wt} else {datList <- amp_out}
    
    # str(amp_out_wt)
    # amp_out_wt
    str(datList$amp)
    class(datList$amp)
    class(datList)
    if (convFact) datList$amp <- add_fact(dat = datList$amp, debug=TRUE)
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
if(FALSE){
  dat=mkSimDat(ndGLM_scl_cc, convFact = TRUE)$amp
  resp="is_u"
  mod=modList[1]
  m=30
  met="rf"
  fam=binomial
  regMet="brglm_fit"
  iter=500
  vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
  mod = modList[1]
  ampDat <- datList$amp
}

# full versions with all comments (old code) in fate_GLM_imp_simulation.Rmd

mkImpSim <- function(ampDat, resp, mod, vars, m=30, met="rf", fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE){
  imp <- mice::mice(ampDat, method=met, m=m, print=FALSE)
  fit = with(imp,
             glm( as.formula( paste0(resp, mod) ),
                  family=fam,
                  method = regMet,
                  control=brglmControl(maxit=iter)
  ))
  
  # pool[[m]] = summary(mice::pool(fit), "all", conf.int=TRUE)
  pool = summary(mice::pool(fit), "all", conf.int=TRUE)
  # why did I make impV?
  # impV = as.character(pool$term[pool$term %in% vars]) # don't need the levels
  # impV
  # pool <- pool[order(as.character(pool$term)),]
  if(debug) {
    cat("pooled model fit:\n")
    print(pool)
    str(pool)
  }
  # pool <- pool %>% filter(term %in% impV)
  pool <- pool %>% filter(term %in% vars)
  # pool
  # pool2
  pool = pool[, c("term", "estimate", "2.5 %", "97.5 %", "fmi")] # maybe do need term?
  # don't need term? - term is the name in the matrix
  # pool = pool[, c("estimate", "2.5 %", "97.5 %", "fmi")]
  # pool
  # return(pool)
}

########################################################################################
###### RUN SIMULATION ###################################################################
########################################################################################
if(FALSE){
  nruns=10
  mod=modList[1]
  ndat=ndGLM_scl_cc
  mets <- c("pmm", "rf")
  # does not include the reference levels:
  vars= c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesLETE", "speciesLETE:nest_age")
  datNA <- dat
  datNA <- ndGLM_scl_cc
    r=1
    # m=1
    # x=1
    x = "rf"
}

# for some reason, this can't find the global vars in fate_GLM1 when I knit that document.
# so, pass them all as arguments? but they already are...
# arguments: function to make sim data, real data, response, predictors, model, nruns
runSim <- function(datNA, resp, vars, mod, mets, nruns=100, debug=FALSE){
  
  res <- array(NA, dim = c(length(vars), length(mets), nruns, 4))
  # dimnames(res) <- list(c("pmm", "rf"),
  dimnames(res) <- list(sort(as.character(vars)),
                        # c("pmm", "rf", "cart"),
                        as.character(mets),
                        as.character(1:nruns),
                        c("estimate", "2.5 %","97.5 %","fmi")
                        )
  if(debug) cat("res matrix is formatted as follows:",str(res))
  for(r in 1:nruns){
    # ndat = nDat
    # if(missing(datNA)) datNA <- mkSimDat(ndat)
    # if the complete data was passed as an argument, add the NAs; otherwise, should be data with NA
    if(all(colSums(is.na(datNA)) == 0)) datNA <- mkSimDat(datNA, debug = TRUE, convFact = TRUE)$amp
    # datNA1 <- datNA
    # for(x in seq_along(mets)){
    for(x in mets){ # this makes it match by name, not just by index
      # could try using map to make sure the vars match up?
      # vals <- as.matrix(mkImpSim(dat=datNA, resp=resp, vars=vars, mod=mod, met=mets[m]))
    # if you include term in the mkImpSim output, as.matrix coerces all columns to character. 
    # need all num, so remove "term" before converting
      # vals <- as.matrix(mkImpSim(ampDat=datNA, resp=resp, vars=vars, mod=mod, met=x))
      vals <- mkImpSim(ampDat=datNA, resp=resp, vars=vars, mod=mod, met=x, debug = debug)
      vmatch <- match(vals[,1], rownames(res)) # col 1 of vals is the row names
      vals <- as.matrix(vals[,-1]) # remove chr column AFTER match so others aren't coerced to chr
      # for (v in vmatch){
      # }
      # res[vmatch, x, r,]  <- vals[,-1]
      res[vmatch, x, r,]  <- vals
      if(debug){ 
        cat(sprintf("\noutput of mkImpSim for run %s and method %s:\n", r, x))
        print(vals)
        str(vals)
        cat("\nres matrix filled in:\n")
        print(res[vmatch,x,r,])
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
  dat=simDatNA
  dat
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
  biasVals <- c("bias", "pctBias", "covRate", "avgWidth", "RMSE")
  # trueVals <- data.frame(vars=vars, value=)
  fullDat <- ndGLM_scl_cc
  impDat  <- res1
  impDat
  trueVals
}

# parAvg <- function(dat, vars, mets, biasVals, trueVals){
# parAvg <- function(dat, resp, vars, mod, mets, biasVals, trueVals){
parAvg <- function(fullDat, impDat, resp, vars, mod, regMet="brglm_fit", fam=binomial, iter=500, mets, biasVals, debug=FALSE){
  # trueVals is optional
  # if(missing(trueVals)){
  fitReal <- glm(as.formula(paste0(resp, mod)),
                 # data=ndGLM_scl_cc,
                 data=fullDat,
                 family=fam,
                 method=regMet,
                 control=brglmControl(maxit=iter) )
  trueVals <- coef(fitReal)[vars] # the coefs have names associated with them
  # trueVals
  # }
# mList <- modList[c(1, 8, 16)]
  # for (v in vars){ # v = var names, not indices
  # avg <- list()
  bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals) ) )
  
  dimnames(bias) <- list( sort(as.character(vars)),
                          # c("pmm", "rf", "cart"),
                          as.character(mets),
                          as.character(biasVals)
  )
  if (debug){
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
    avg <- apply(impDat[v, , ,],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
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
    trueVals[v]
    # true
    # bias[vars[v], , "bias"] <- avg[,"estimate"] - trueVals[vars[v]]
    avg[,"estimate"]
    bias[v, , "bias"] <- avg[,"estimate"] - true
    bias[v, , "pctBias"] <- 100 * abs((avg[,"estimate"] - true) / true )
    # bias[vars[v], , "covRate"] <- avg[,"2.5 %"] < true & true > avg[,"97.5 %"]
    # bias[vars[v], , "avgWidth"] <- avg[,"97.5 %"] - avg[,"2.5 %"] 
    bias[v, , "covRate"] <- rowMeans(impDat[v,,,"2.5 %"] < true & true < impDat[v,,,"97.5 %"])
    bias[v, , "avgWidth"] <- avg[,"97.5 %"] - avg[,"2.5 %"] 
    bias[v, , "RMSE"] <- sqrt((avg[,"estimate"] - true)^2)
    bias
  }
  return(bias)
  
}
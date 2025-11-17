
  # col_sel <- c(prVars,r) # columns to select, as strings
# source("fateGLM_impsim.R")
library(brglm2) # penalized logistic regression
library(CALIBERrfimpute) # could access using ::
library(gt)
library(gtsummary)
library(tidyverse)
source("mice_functions.R")
source("imp_sim_functions.R")
source("missing_data.R")
source("other_imp.R")

#########################################################################################

params <- list(nrun=100, hdir="/home/wodehouse/projects/fate_glm/",
               deb=FALSE,
               xdeb=FALSE, # for obsessively checking things - not useful for normal debugging & lots of text output
               ipl=FALSE,
               lin=FALSE,
               seed = NULL,
               resp = NULL,
               j = 50,
               m=20)

arg <- commandArgs(trailingOnly=TRUE)
if(length(arg)==0){
  # stop("at minimum, needs argument for 'lin' (T or F)")
  cat("\n\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
  cat("\n** NOTE ** no arguments provided - using default of linux = FALSE\n\n")
  cat("/////////////////////////////////////////////////////////////////////////////////////////////\n\n\n")
} else if(length(arg) > 0){
  cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n\n")
  cat("\t\targ =")
  print(arg)
  cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n\n")
  for(a in arg){
    # if(is.numeric(a)) params$nrun <- a
    if(grepl("^\\d+$", a)) params$nrun  <- a
    else if(a=="lin") params$lin        <- TRUE
    else if(a=="deb") params$deb        <- TRUE
    else if(a=="xdeb") params$xdeb        <- TRUE
    else if(a=="ipl") params$ipl       <- TRUE
    else if(a=="is_u" | a == "HF_mis") params$resp       <- a
    else if(grepl("j\\d+", a)) params$j <- as.numeric(str_extract(a, "\\d+"))
    else if(grepl("s\\d+", a)) params$seeds <- c(as.numeric(str_extract(a, "\\d+")))
    else if(grepl("m\\d+", a)) params$m <- as.numeric(str_extract(a, "\\d+"))
    # else if()
  }
}
#########################################################################################

debugging <- FALSE # for uickly setting values when working in the file with the functions
suffix <- sprintf("%sruns", params$nrun)
# cat("\n\n>> home directory:", params$hdir, "\t> & number of runs:", params$nrun, "\t> & output suffix:", suffix)
# cat("\n\n>> home directory:", params$hdir)
# "\t>> & output suffix:", suffix)
if(params$j > params$nrun) stop("j must be less than or equal to nrun!")
if(is.null(params$resp)) stop("no response variable specified!")

modList <- readLines("modList.txt")
mods4sim <- modList[c(1,8,16) ]
names(mods4sim) <- c("m1", "m8", "m16")
formulas <- readRDS("form_lists.rds")
metLists <- readRDS("met_lists.rds")

dat4sim <- read.csv("dat_complete.csv", stringsAsFactors = TRUE)
dat4sim$cam_fate <- relevel(dat4sim$cam_fate, ref="H") # make 'H' the reference category
dat4sim$species <- relevel(dat4sim$species, ref="LETE")
levels(dat4sim$HF_mis) <- c(0,1)
levels(dat4sim$is_u)   <- c(0,1)
 
prVars <- c("species", "cam_fate", "obs_int", "nest_age", "fdate")
vars <-  c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") # all vars
# mets <- c("default","pmm", "rf", "cart", "caliber","passive", "stratify","cc")# don't need full here?
mets <- c("default","pmm", "rf", "cart", "caliber","passive", "stratify","cf_cc","cc")# don't need full here?

# resp <- "is_u"
col_list<- c(prVars,params$resp )# columns to select, as strings
form_list <- formulas[[params$resp]]
if(is.null(params$seeds)){
  params$seeds <- c( 11153, 71358, 102891, 82985, 61389)
  # params$seeds <- c( 102891, 82985, 61389, )
}

#########################################################################################

# cat("\n\n>> number of imputations:", params$m, class(params$m))
cat("\n********************************************************************************************")
cat("\n>> methods:", mets)
# cat("\t\t>> & number of imputations:", params$m, class(params$m))
cat("\t>> & no. imputations:", params$m)
# cat("\n\n>> bias to be calculated:", bias_names, "\n")
# cat("\n\n>>>> date & time:", format(Sys.time(), "%d-%b %H:%M"))

cat("\n\n********************************************************************************************")
cat("\n>> response:", params$resp,"\n\t& columns for imputation:", col_list)
cat("\n\n********************************************************************************************")
cat("\n>> output will be saved every", params$j, "runs to home directory:", params$hdir)
cat("\n\n********************************************************************************************")
cat("\n>>>> date & time:", format(Sys.time(), "%d-%b %H:%M"))

res <- array(NA, dim = c(length(vars), length(mets), params$nrun, 3, length(mods4sim)))
# dimnames(res) <- list(c("pmm", "rf"),
dimnames(res) <- list(sort(as.character(vars)),
                      # c("pmm", "rf", "cart"),
                      as.character(mets),
                      as.character(1:params$nrun),
                      # c("estimate", "2.5 %","97.5 %","fmi"),
                      c("estimate", "2.5 %","97.5 %"),
                      names(mods4sim)
                      )

# cat(sprintf("\t\t>>> running simulation %s times. seed = %s\n\n", params$nrun, params$seed))
# cat("\n********************************************************************************************")

for (seed in params$seeds){
  cat(sprintf("\t>>>>>> running simulation %s times. seed = %s >>>>>>>>\n\n", params$nrun, seed))
  # cat("seed=",params$seed, " - ")
  for(run in 1:params$nrun){
    cat(run)
    # repeat this until you get a dataset w/o missing levels??
    datNA <- mkSimDat(nd = dat4sim, seeed = run+seed, vars=vars, method = "amp", wt = TRUE, xdebug=params$xdeb, debug = params$deb, convFact = TRUE)
    datNA <- datNA$amp
    
    for(x in seq_along(mets)){ # does matching by index help the trycatch statement?
      ## *~*~*~*~*
      # if (xdebug) cat("\n\n>>>> method:", x)
      skiptoNext <- FALSE
      
      # the error seems to come from levels with zero observations
      # so maybe ampute() is removing all the camera fates of category X (F in this case)
      tryCatch(
        expr = {
          vals <- mkImpSim(ampDat=datNA,
                           cols=col_list,
                           resp=params$resp, 
                           form_list =form_list,
                           vars=vars,
                           mods=mods4sim,
                           met=mets[x],
                           debug = params$debug,
                           m=params$m, 
                           xdebug=params$xdebug,
                           impplot=params$ipl
                           )
          # imp <- eval(parse(text=impCall))
          # skiptoNext <- FALSE
          # if(length(imp$loggedEvents > 0)) print(imp$loggedEvents)
        },
        error = function(e){
          cat("\nERROR:", conditionMessage(e), "\n")
          # next
          # return(NULL)
          # skiptoNext <- TRUE
          skiptoNext <<- TRUE  # superassignment operator- not sure if necessary
          # cat(levels)
          # imp <- list(imp=NA)
          # ret[,,y]
          # continue()
        }
      )
      if(skiptoNext) next
      
      vmatch <- match(rownames(vals), rownames(res)) # col 1 of vals is the row names
      res[vmatch, mets[x], run,,]  <- vals
    }
    
    # j <- 50
    # j <-2 
    if(run %% params$j == 0){
    # if(run %% 4 == 0){
      begn <- run-params$j
      endd <- run-0
      nowtime <- format(Sys.time(), "%d%b%H%M")
      fname <- paste(sprintf("out/runs%sto%s_resp%s_seed%s_%s.rds", begn, endd, params$resp, seed, nowtime))
      saveRDS(res[,,begn:endd,,], fname)
      cat(sprintf("\n>>>>>> saved runs %s to %s to file!\n", begn, endd))
      # cat("\n********************************************************************************************")
    }
    # return(res)
  }
}


# }
# imp_sim <- runSim(fullDat=dat4sim,col_sel = col_list,mets = met_list,forms=form_list, resp = r, vars = var_list, mods = mods4sim,par=params) # don't want to set seed
# if(debug) Rprof()
# imp_sim_m <- runSim(fullDat=dat4sim,col_sel = col_list,mLists = metLists,forms=form_list, resp = r, vars = var_list, mods = mods4sim,par=params) # don't want to set seed
# if(debug) Rprof(NULL)
# if(params$deb){
#   cat("\n\n********************************************************************************************\n")
#   cat(">>>>> BIAS VALUES: \n")
#   cat("********************************************************************************************\n")
# }
# # bias_out <- parAvg(fullDat = dat4sim, impDat = imp_sim,resp = r, vars = var_list, mod = mods4sim[z], mets = met_list, biasVals = bias_names, debug = params$deb)
# # bias_out <- parAvg(fullDat = dat4sim, impDat = imp_sim,resp = r, vars = var_list, modnum = z, mets = met_list, biasVals = bias_names, debug = params$deb)
# bias_out <- parAvg(fullDat = dat4sim, impDat = imp_sim,hdir = params$hdir,resp = r, vars = var_list, mods=mods4sim, mets = met_list, biasVals = bias_names, debug = params$deb, xdebug=params$xdeb)
# }
# }

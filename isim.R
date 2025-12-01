

  # col_sel <- c(prVars,r) # columns to select, as strings
# source("fateGLM_impsim.R")

library(brglm2) # penalized logistic regression
library(CALIBERrfimpute) # could access using ::
library(gt)
library(gtsummary)
suppressMessages(library(mice))
suppressMessages(library(tidyverse))
# look into tidytable for tidyverse syntax w/ data.table
source("mice_functions.R")
source("imp_sim_functions.R")
source("missing_data.R")
source("other_imp.R")
source("gen_sim.R") ### new data simulation function

#########################################################################################

params <- list(nrun=100,
	       hdir="/home/wodehouse/Projects/fate_glm/",
	       #outdir = paste0(format(Sys.time(), "%d%b"),"/"),
	       outdir = format(Sys.time(), "%d%b"),
               deb=FALSE,
               xdeb=FALSE, # for obsessively checking things - not useful for normal debugging & lots of text output
               ipl=FALSE,
               win=FALSE,
               test=FALSE,
               seeds = NULL,
               resp = NULL,
               j = 50,
               m=20)

arg <- commandArgs(trailingOnly=TRUE)
if(length(arg)==0){
	cat("\n\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
	cat("** NOTE ** no arguments provided - using default of windows = FALSE\n")
	cat("/////////////////////////////////////////////////////////////////////////////////////////////\n\n")
} else if(length(arg) > 0){
  #cat("\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
  cat("\n/////////////////")
  cat("  arg =  ", paste(unlist(arg), collapse=" ; "))
  #print(arg)
  cat("  /////////////////\n")
  #cat("/////////////////////////////////////////////////////////////////////////////////////////////\n")
  for(a in arg){
	  if(grepl("r\\d+$", a)) params$nrun  <- as.numeric(str_extract(a, "\\d+"))
	  else if(a=="lin") params$win        <- TRUE
	  else if(a=="deb") params$deb        <- TRUE
	  else if(a=="xdeb") params$xdeb        <- TRUE
	  else if(a=="ipl") params$ipl       <- TRUE
	  else if(a=="is_u" | a == "HF_mis") params$resp       <- a
          else if(a=="test") params$test <- TRUE
	  else if(grepl("j\\d+", a)) params$j <- as.numeric(str_extract(a, "\\d+"))
	  else if(grepl("s\\d+", a)) params$seeds <- c(as.numeric(str_extract(a, "\\d+")))
	  else if(grepl("m\\d+", a)) params$m <- as.numeric(str_extract(a, "\\d+"))
  }
}
#########################################################################################

if(params$win) cat("\n>>>> date & time:", format(Sys.time(), "%d-%b %H:%M"))
debugging <- FALSE # for uickly setting values when working in the file with the functions
suffix <- sprintf("%sruns", params$nrun)
if(params$j > params$nrun) stop("j must be less than or equal to nrun!")
if(is.null(params$resp)) stop("no response variable specified!")

modList <- readLines("modList.txt")
mods4sim <- modList[c(1,8,16) ]
names(mods4sim) <- c("m1", "m8", "m16")
formulas <- readRDS("form_lists.rds")
metLists <- readRDS("met_lists.rds")
betas    <- readRDS("betas.rds")

cat("\n>> methods for individual variables:\n")
print(metLists)

#dat4sim <- read.csv("dat_complete.csv", stringsAsFactors = TRUE)
#dat4sim$cam_fate <- relevel(dat4sim$cam_fate, ref="H") # make 'H' the reference category
#dat4sim$species <- relevel(dat4sim$species, ref="LETE")
#levels(dat4sim$HF_mis) <- c(0,1)
#levels(dat4sim$is_u)   <- c(0,1)
#dat4sim <- dat4sim %>% mutate(nest_age = if_else(nest %in% c(10240, 20301, 30903, 30929), NA, nest_age))
#cat("\n>> new nest ages:\n", table(dat4sim$nest_age, useNA="ifany"))
#cat("\n>> new missing nest age:\n", dat4sim$nest[is.na(dat4sim$nest_age)])
#dat4sim <- dat4sim[complete.cases(dat4sim),]
#dat4sim <- dat4sim[!is.na(dat4sim$nest_age),]
#dat4sim <- dat4sim %>% filter(nest_age == !is.na(nest_age))
#cat("\n>>>> remove NAs:\n", str(dat4sim))

#dat4fates <- dat4sim %>% 
#  mutate(cam_fate = if_else(cam_fate=="S"|cam_fate=="Hu", "F", cam_fate)) %>% 
#  mutate(cam_fate= relevel(as.factor(cam_fate), ref="H"))
# dat3fate
 
prVars <- c("species", "cam_fate", "obs_int", "nest_age", "fdate")

cat("sort vars or not???")
#vars <- sort( c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") )# all vars
vars <-  c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") # all vars
#mets <- c("default","pmm", "rf", "cart", "caliber","passive", "stratify","cf_cc","cc")# don't need full here?
mets <- ifelse(params$test,c("default", "passive", "stratify"), c("default", "cart", "caliber","passive", "stratify","cf_cc","cc"))# don't need full here?

# resp <- "is_u"
#col_list<- c(prVars,params$resp )# columns to select, as strings
form_list <- formulas[[params$resp]]
seed_out <- FALSE
now_dir <- paste(params$hdir, params$outdir, sep="out/")
#cat("\n>>> save directory:", now_dir,"\n")
if(!dir.exists(now_dir)) dir.create(now_dir)
if(is.null(params$seeds)){
	seed_out <- TRUE
	params$seeds <- c(71358, 102891, 82985, 61389, 11153)
	last_seed <- as.numeric(grep(readLines("seed.flag"), params$seeds))
	cat("last seed:", last_seed)
	#new_order <- params$seeds + last_seed
	new_order <- seq(1,length(params$seeds)) + last_seed
	cat("new order:", new_order)
	new_order[new_order > 5] = new_order[new_order > 5] - 5
	cat("new order:", new_order)
	params$seeds <- params$seeds[new_order]
	cat("\n****************************************************************************")
	cat("\n>>> using default seed list, in new order:", params$seeds, "\t")
}

#########################################################################################

# cat("\n\n>> number of imputations:", params$m, class(params$m))
#cat("\n****************************************************************************")
cat("\n")
cat("\n>> methods:", mets)
# cat("\t\t>> & number of imputations:", params$m, class(params$m))
#cat("  >> & no. imp.:", params$m)
# cat("\n\n>> bias to be calculated:", bias_names, "\n")
# cat("\n\n>>>> date & time:", format(Sys.time(), "%d-%b %H:%M"))

#cat("\n****************************************************************************")
cat("\n")
cat(sprintf("\n>> total reps: %s x %s", length(params$seeds), params$nrun)) 
cat("\t>> response:", params$resp,"\n\t& cols for imputation:", col_list)
cat("\n\n****************************************************************************")
#cat("\n>> output will be saved every", params$j, "runs to dir:", paste0(params$hdir,"out/",params$outdir))
cat("\n>> output will be saved every", params$j, "runs to dir:", now_dir)
#cat("\n****************************************************************************")
cat("\n")
#cat("\n>>>> date & time:", format(Sys.time(), "%d-%b %H:%M"))

# cat(sprintf("\t\t>>> running simulation %s times. seed = %s\n\n", params$nrun, params$seed))
# cat("\n********************************************************************************************")

for (seed in params$seeds){
    res <- array(NA, dim = c(length(vars), length(mets), params$nrun, 3, length(mods4sim)))
    dimnames(res) <- list(sort(as.character(vars)),
                          # c("pmm", "rf", "cart"),
                          as.character(mets),
                          as.character(1:params$nrun),
                          # c("estimate", "2.5 %","97.5 %","fmi"),
                          c("estimate", "2.5 %","97.5 %"),
                          names(mods4sim)
                          )

    cat(sprintf(">>> running simulation %s times \t>>> seed = %s ", params$nrun, seed))
    cat("\t>> & no. imp.:", params$m, "\n\n")
    for(mod in mods4sim){
        for(run in 1:params$nrun){
                cat(run)
            # repeat this until you get a dataset w/o missing levels??
            #mkSim <- function(resp, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug)
            dat4sim <- mkSim(params$resp, mod, 150, cMat, mList, betas, fprob, sprob, prList, debug=params$deb)
            datNA <- mkSimDat(nd = dat4sim, seeed = run+seed, vars=vars, method = "amp", wt = TRUE, xdebug=params$xdeb, debug = params$deb, convFact = TRUE)
            datNA <- datNA$amp
            
            # prv<- profvis::profvis(
            for(x in seq_along(mets)){ # does matching by index help the trycatch statement?
              ## *~*~*~*~*
              # if (xdebug) cat("\n\n>>>> method:", x)
              skiptoNext <- FALSE
              
              # the error seems to come from levels with zero observations
              # so maybe ampute() is removing all the camera fates of category X (F in this case)
              tryCatch(
                expr = {
#mkImpSim <- function(fullDat, ampDat,pr_list, resp_list, mod, vars, met, form_list, m=20, fam=binomial, regMet="brglm_fit", iter=500, debug=FALSE, xdebug=FALSE, impplot=FALSE){
                  vals <- mkImpSim(ampDat=datNA,
                                   #cols=col_list,
                                   pr_list=prVars,
                                   resp=params$resp, 
                                   form_list =form_list,
                                   vars=vars,
                                   #mods=mods4sim,
                                   mod=mod,
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
            
          #  vmatch <- match(rownames(vals), rownames(res)) # col 1 of vals is the row names
            #print(head(res))
            print(res[,mets[x],run,,mod])
            print(vals)
            #res[vmatch, mets[x], run,,]  <- vals
            res[, mets[x], run,,mod]  <- vals
          }
            
            if(run %% params$j == 0){
                    begn <- run-params$j
                    endd <- run-0
                    #cat("\n>>> creating directory:", now_dir, "exists?", exists(now_dir))
                    nowtime <- format(Sys.time(), "%H_%M")
                    fname <- paste0(now_dir,sprintf("runs%sto%s_%s_seed%s_%s.rds", begn, endd, params$resp, seed, nowtime))
                    saveRDS(res[,,begn:endd,,], fname)
                    cat(sprintf("\n>>>>>> %s - saved runs %s to %s to file: %s\n\n",nowtime, begn, endd, fname))
            }
      }
	if (seed_out){
		writeLines(as.character(seed), "seed.flag")
		cat("\n>> wrote seed value to file\n")
	}
    }
}



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

### PROFILING THE CODE
# if(FALSE){
#   for(i in 1:2){
#     # resp <- "is_u"
#     for(resp in resp_list){
#       col_list<- c(prVars,resp )# columns to select, as strings
#       form_list <- formulas[[resp]]
#       filename <- sprintf("out/rprof_%s_%s.out", resp, i)
#       Rprof(filename)
#       runSim(fullDat=dat4sim, 
#              col_sel=col_list,
#              mLists=metLists,
#              forms=form_list,
#              resp=resp,
#              vars=var_list,
#              mods=mods4sim,
#              par=params
#              )
#       Rprof(NULL)
#     }
#     # filename <- paste0("rprof_",i,".out")
#     
#   }
# }
if(FALSE){
    for(resp in resp_list){
      col_list<- c(prVars,resp )# columns to select, as strings
      form_list <- formulas[[resp]]
      filename <- sprintf("out/prvis_%s.out", resp)
      prv<- profvis::profvis(
        runSim(fullDat=dat4sim, 
               col_sel=col_list,
               mLists=metLists,
               forms=form_list,
               resp=resp,
               vars=var_list,
               mods=mods4sim,
               par=params
               ),
        prof_output = filename
        )
      htmlwidgets::saveWidget(prv, sprintf("profile_%s.html", resp))
    }
}

if(FALSE){
resp_list <- c("is_u", "HF_mis")
    for(resp in resp_list){
      col_list<- c(prVars,resp )# columns to select, as strings
      form_list <- formulas[[resp]]
      filename <- sprintf("out/prvis_%s.out", resp)
      aDat <- mkSimDat(seeed=13, nd=dat4sim, vars=vars, convFact=TRUE)$amp
      for(met in mets){
        prv<- profvis::profvis(
          mkImpSim(fullDat=dat4sim,
                   ampDat=aDat,
                   cols=col_list,
                   resp=resp,
                   mods=mods4sim,
                   vars=vars,
                   met=met,
                   form_list=form_list
                   
                   ),
          prof_output = filename
          )
        htmlwidgets::saveWidget(prv, sprintf("prof/profile2_%s_%s.html", resp, met))
      }
    }
}

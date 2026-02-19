

  # col_sel <- c(prVars,r) # columns to select, as strings
# source("fateGLM_impsim.R")

library(brglm2) # penalized logistic regression
library(CALIBERrfimpute) # could access using ::
library(gt)
library(gtsummary)
library(data.table)
suppressMessages(library(mice))
suppressMessages(library(tidyverse))
# look into tidytable for tidyverse syntax w/ data.table
options(width=99)
#if(params$deb|params$xdeb) func = "debug_imp_sim_func.R" else func = "imp_sim_functions.R"
#source(func)
#source("mice_functions.R")
#source("imp_sim_functions.R")
#source("missing_data.R")
#source("other_imp.R")
#source("gen_sim.R") ### new data simulation function

#########################################################################################

params <- list(nrun=100,
	       hdir="/home/wodehouse/Projects/fate_glm/",
	       #outdir = paste0(format(Sys.time(), "%d%b"),"/"),
               suffix = "",
	       outdir = format(Sys.time(), "%d%b"),
               ampWt=1,
               deb=FALSE,
               xdeb=FALSE, # for obsessively checking things - not useful for normal debugging & lots of text output
               ipl=FALSE,
               win=FALSE,
               test=FALSE,
               seeds = NULL,
               # resp = NULL,
               j = 50,
               m=20)

arg <- commandArgs(trailingOnly=TRUE)
# Rscript isim.R win test r3 s613 m5 j1
if(length(arg)==0){
    cat("\n\n/////////////////////////////////////////////////////////////////////////////////////////////\n")
    cat("** NOTE ** no arguments provided - using default of windows = FALSE\n")
    cat("/////////////////////////////////////////////////////////////////////////////////////////////\n\n")
}else if(length(arg) > 0){
    cat("\n\n/////////////////////")
    cat("  arg =  ", paste(unlist(arg), collapse=" ; "))
    cat("  ////////////////////////////////\n")
  for(a in arg){
	  if(grepl("r\\d+$", a)) params$nrun  <- as.numeric(str_extract(a, "\\d+"))
	  else if(a=="win") params$win        <- TRUE
	  else if(a=="deb") params$deb        <- TRUE
	  else if(a=="xdeb") params$xdeb        <- TRUE
	  else if(a=="ipl") params$ipl       <- TRUE
	  # else if(a=="is_u" | a == "HF_mis") params$resp       <- a
          else if(a=="test") params$test <- TRUE
	  else if(grepl("suff_\\w+", a)) params$suffix <- as.numeric(str_extract(a, "(?<=suff_)\\w+"))
	  else if(grepl("aw\\d", a)) params$ampWt <- as.numeric(str_extract(a, "\\d"))
	  else if(grepl("j\\d+", a)) params$j <- as.numeric(str_extract(a, "\\d+"))
	  else if(grepl("s\\d+", a)) params$seeds <- c(as.numeric(str_extract(a, "\\d+")))
	  else if(grepl("m\\d+", a)) params$m <- as.numeric(str_extract(a, "\\d+"))
  }
}

#########################################################################################

s_files <- c("mice_functions.R", "missing_data.R", "other_imp.R", "imp_sim_functions.R")
if(params$deb|params$xdeb) s_files[4] = "debug_imp_sim_func.R" 
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
#cat("\n###################################################################")
cat("\n>>> sourcing files:", paste(s_files,collapse=" ; "))
#lapply(s_files, function(x) invisible(source(x)))
invisible(lapply(s_files, function(x) source(x)))

#########################################################################################

bprint <- function(x) print(rbind(head(x, 3), tail(x,3)))
debugging <- FALSE # for uickly setting values when working in the file with the functions
if(params$win){
    cat("\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n")
    cat("\n>>>> date & time:", format(Sys.time(), "%d-%b %H:%M")) 
    cat("\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n") 
    params$hdir <- "C:/Users/sarah/Dropbox/Models/fate_glm/"
}
if(params$j > params$nrun) stop("j must be less than or equal to nrun!")
# if(is.null(params$resp)) stop("no response variable specified!")
if(params$test){
    #mets <- c("cc", "cart", "caliber", "passive")
    params$outdir <- paste("t",params$outdir,sep="-")
}

suffix <- paste(sprintf("%sruns", params$nrun), params$suffix, sep="_")
now_dir <- paste(params$hdir, params$outdir, sep="out/")
if(!dir.exists(now_dir)) dir.create(now_dir)

seed_out <- FALSE
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


modList <- readLines("modList.txt")
mods4sim <- modList[c(1,8,16) ]
names(mods4sim) <- c("m1", "m8", "m16")

#### for the imputation: ####################################################################
formulas <- readRDS("form_lists.rds")
metLists <- readRDS("met_lists.rds")
#cat("\n>> methods for individual variables:\n")
#print(metLists)

#### for creating the simulated data: #######################################################
betas <- readRDS("betas.rds") # this is a list of vectors or something
mList <- readRDS('means.rds')
cMat  <- readRDS('cormat.rds')
fprob <- c( 'H'=0.53,'A'=0.1, 'D'=0.13,'F'=0.06, 'Hu'=0.13,'S'=0.06 )
sprob <- c('CONI'=0.32, 'LETE'=0.68)
nnest <- ifelse(params$test==T, 200, 200) # number of simulated nests
mpatt  <- readRDS("misPatt.rds")
awFile <- case_when(params$ampWt==1 ~ "ampWts.rds",
                    params$ampWt==2 ~ "ampWts2.rds",
                    params$ampWt==3 ~ "ampWts3.rds")
ampwt <- readRDS(awFile)
cat("\n\n<><><><><><><><><><> Missingness pattern & variable weights: <><><><><><><><><><><><><><><><><><><><>\n")
print(mpatt)
print(ampwt)

#########################################################################################

#cat("sort vars or not???")
vars <-  c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") # all vars
prVars <- c("species", "cam_fate", "obs_int", "nest_age", "fdate")
mets <- c("default", "cart", "caliber","passive", "stratify","cf_cc","cc")# don't need full here?
resp_list <- c("is_u", "HF_mis")

#form_list <- formulas[[params$resp]]
#########################################################################################

cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
#cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
cat("\n>> METHODS:", mets)
cat(sprintf("\t>> TOTAL REPS: %s x %s WITH %s NESTS EACH", length(params$seeds), params$nrun, nnest)) 
if(params$deb) cat("\t- debug ON")
if(params$xdeb) cat(" + EXTRA")
#cat("\t>> response:", params$resp,"\n\t& cols for imputation:", col_list)
#cat("\n\n****************************************************************************")
cat("\n>> output will be saved every", params$j, "runs to dir:", now_dir,"\n\n")

for (seed in params$seeds){
    outf <- paste0(now_dir,sprintf("/%s_loggedEvents.out",seed ))
    #if(params$test) outf <- paste0(now_dir,sprintf("/%s_loggedEvents-TEST.out",seed ))
    res <- array(NA, dim = c(length(vars), length(mets), params$nrun, 3, length(mods4sim), length(resp_list)))
    # dimnames(res) <- list(sort(as.character(vars)),
    dimnames(res) <- list(as.character(vars),# need to be in same order as vals
                          as.character(mets),
                          as.character(1:params$nrun),
                          c("estimate", "2.5 %","97.5 %"),
                          names(mods4sim),
                          resp_list
                          )
    camFateVars <- vars[grepl(pattern="cam_fate", x=vars, fixed=TRUE)]
    # keep a count of sample size for each category
    varInfo <- array(NA, dim=c(length(camFateVars), params$nrun, length(mods4sim)))
    dimnames(varInfo) <- list( camFateVars, seq(1,params$nrun), mods4sim)

#cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
    cat(sprintf(">>>>>>>>>>>>>> running simulation %s times \t>>> seed = %s ", params$nrun, seed))
    cat("\t>> & no. imp.:", params$m, "\n\n")
    for(mod in seq_along(mods4sim)){
        # select the correct list of formulas and list of beta values!
        #cat("\n()()()()()() MODEL", mods4sim[mod], " ()()()()()()()()()()()()()()()() \n\n") ## *~*~*~*~*
        cat("\n\n::::::::::::::::::: MODEL", mods4sim[mod], ":::::::::::::::::::::::::::::::\n\n") ## *~*~*~*~*
        beta_list <- betas[[mod]]
        form_list <- formulas[[names(mods4sim)[mod]]]
        cat("\n>> formula list, beta list:") ## *~*~*~*~*
        print(names(form_list))
        print(str(beta_list))
        for(run in 1:params$nrun){
            cat(run)
            # repeat this until you get a dataset w/o missing levels??
            #dat4sim <- mkSim(resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, debug=params$deb)
            dat4sim <- mkResp(seed=run+seed,resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, debug=params$deb)
            datNA <- mkSimDat( seeed = run+seed, nd = dat4sim, mpatt=mpatt, wts=ampwt, xdebug=params$xdeb, debug = params$deb, convFact = TRUE)
            #if(params$deb) cat("\n*** datNA:\n")
            #if(params$deb) cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
            #if(params$deb) cat("\n\n<><><><><><><><><><<><><><><<><><><><> AMPUTED DATA: <><><><><><><><><><<><><><><<><><><><>>>\n ")
            #if(params$deb) qvcalc::indentPrint(str(datNA$amp), indent=8)
            #if(params$deb) qvcalc::indentPrint(colSums(is.na(datNA$amp)), indent=8)
            datNA <- datNA$amp
            #if(params$deb) qvcalc::indentPrint(str(datNA), indent=8)
            #if(params$debug) cat("\n>> cam fate vars for this run:")
            #cat("created simulated nest data")
            #if(params$debug) qvcalc::indentPrint(datNA$cam_fate, indent=8)
            #if(params$deb) cat("\n>> cam fate vars for this run:")
            #if(params$deb) qvcalc::indentPrint(table(datNA$cam_fate), indent=8)
            
            for(x in seq_along(mets)){ # does matching by index help the trycatch statement?
                skiptoNext <- FALSE
                ## *~*~*~*~*
                # if (params$deb) cat("\n\n<><><><><><><><> method:", mets[x], "\t- ")
                #if (params$deb) cat("\n\n<><><><><><><><> method:", mets[x], "<><><><><><><><><><><><><><><><><><><><>") ## *~*~*~*~*

                # the error seems to come from levels with zero observations
                # so maybe ampute() is removing all the camera fates of category X (F in this case)
                tryCatch(
                expr = {
                    vals <- mkImpSim(aDat=datNA,
                                   #cols=col_list,
                                   # pr_list=prVars,
                                   # resp=params$resp, 
                                   resp_list=resp_list,
                                   form_list =form_list,
                                   met_list=metLists[,,,mod],
                                   vars=vars,
                                   outFile=outf,
                                   #mods=mods4sim,
                                   modd=mods4sim[mod],
                                   met=mets[x],
                                   debug = params$deb,
                                   m=params$m, 
                                   xdebug=params$xdeb,
                                   impplot=params$ipl
                                   )
                },
                error = function(e){
                    cat("\nERROR:", conditionMessage(e), "\n")
                    skiptoNext <<- TRUE  # superassignment operator- not sure if necessary
                }
                )
                if(skiptoNext) next
                res[, mets[x], run,,mod,]  <- vals
            }
            
            if(params$test) cat(sprintf("\n\n::::::::::::::::::::::::::::::: ALL OUTPUT - RUN %s, MODEL %s:::::::::::::::::::::::::::::::\n\n",run,mod ))## *~*~*~*~*
            if(params$test) qvcalc::indentPrint(res[,, run,,mod,])
            #if(run %% params$j == 0){
                    #begn <- run-params$j
                    #endd <- run-0
                    ##cat("\n>>> creating directory:", now_dir, "exists?", exists(now_dir))
                    #nowtime <- format(Sys.time(), "%H_%M")
                    ## fname <- paste0(now_dir,sprintf("runs%sto%s_%s_seed%s_%s.rds", begn, endd, params$resp, seed, nowtime))
                    #fname <- paste0(now_dir,sprintf("/runs%sto%s_mod%s_seed%s_%s_%s.rds", begn, endd,mod, seed, nowtime, suffix))
                    #saveRDS(res[,,begn:endd,,,], fname)
                    #cat(sprintf("\n>>>>>> %s - saved runs %s to %s to file: %s\n\n",nowtime, begn, endd, fname))
            #}
            #varInfo <- array(NA, dim=c(length(camFateVars), params$nrun, length(mods4sim), length(resp_list)))
            #if(params$deb) cat("\n>> cam fate vars for this run:", table(datNA$cam_fate))
            #if(params$deb) cat("\n>> species for this run:", table(datNA$speciesCONI))
            #if(params$deb) cat("\n>> add to matrix:", varInfo[,run,mod]) 
            varInfo[,run,mod] <- table(datNA$cam_fate)[-1]
            #if(params$deb) cat("\n>> in matrix:", varInfo[,run,mod])
        }
	if (seed_out){
            writeLines(as.character(seed), "seed.flag")
            cat("\n>> wrote seed value to file\n")
	}
        #if(params$test) cat("\n >> cam fate vars for this model:")
        if(params$test) qvcalc::indentPrint(varInfo[,,mod])
        nowtime <- format(Sys.time(), "%H_%M")
        fname <- paste0(now_dir, sprintf("/vals_mod%s_seed%s_%s_%s.rds", mod, seed, nowtime, suffix))
        cat("\n ~ ~ ~ saving output to file:", fname)
        saveRDS(res, fname)

        fnameV <- paste0(now_dir, sprintf("/camFate_mod%s_seed%s_%s_%s.rds", mod, seed, nowtime, suffix ))
        cat("\n ~ ~ ~ saving cam fate vars to file:", fnameV)
        saveRDS(varInfo, fnameV)

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
#if(FALSE){
    #for(resp in resp_list){
    #  col_list<- c(prVars,resp )# columns to select, as strings
    #  form_list <- formulas[[resp]]
    #  filename <- sprintf("out/prvis_%s.out", resp)
    #  prv<- profvis::profvis(
    #    runSim(fullDat=dat4sim, 
    #           col_sel=col_list,
    #           mLists=metLists,
    #           forms=form_list,
    #           resp=resp,
    #           vars=var_list,
    #           mods=mods4sim,
    #           par=params
    #           ),
    #    prof_output = filename
    #    )
    #  htmlwidgets::saveWidget(prv, sprintf("profile_%s.html", resp))
    #}
#}

#if(FALSE){
#resp_list <- c("is_u", "HF_mis")
#    for(resp in resp_list){
#      col_list<- c(prVars,resp )# columns to select, as strings
#      form_list <- formulas[[resp]]
#      filename <- sprintf("out/prvis_%s.out", resp)
#      aDat <- mkSimDat(seeed=13, nd=dat4sim, vars=vars, convFact=TRUE)$amp
#      for(met in mets){
#        prv<- profvis::profvis(
#          mkImpSim(fullDat=dat4sim,
#                   ampDat=aDat,
#                   cols=col_list,
#                   resp=resp,
#                   mods=mods4sim,
#                   vars=vars,
#                   met=met,
#                   form_list=form_list
#                   
#                   ),
#          prof_output = filename
#          )
#        htmlwidgets::saveWidget(prv, sprintf("prof/profile2_%s_%s.html", resp, met))
#      }
#    }
#
#}
#

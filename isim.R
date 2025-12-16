
library(brglm2) # penalized logistic regression
library(CALIBERrfimpute) # could access using ::
library(gt)
library(gtsummary)
library(data.table)
suppressMessages(library(mice))
suppressMessages(library(tidyverse)) # look into tidytable for tidyverse syntax w/ data.table
options(width=199)
###### PARAMS & ARGS ###########################################################################
params <- list(nrun=100,
	       hdir="/home/wodehouse/Projects/fate_glm/",
	       #outdir = paste0(format(Sys.time(), "%d%b"),"/"),
               suffix = "",
	       outdir = format(Sys.time(), "%d%b"),
               ampWt=1,
               #deb=FALSE,
               #xdeb=FALSE, # for obsessively checking things - not useful for normal debugging & lots of text output
               #mdeb=FALSE,
               vbose=0,
               ipl=FALSE,
               win=FALSE,
               test=FALSE,
               seeds = NULL,
               # resp = NULL,
               j = 50,
               m=20)
arg <- commandArgs(trailingOnly=TRUE)
if(length(arg)==0){
    cat("\n\n/////////////////////// ** NOTE ** no arguments provided - using defaults of nrun=1, m=5, verbose=FALSE //////////////////////////\n\n")
    params$nrun=1
    params$j=1
    params$m=5
    # isn't the same as "test" - this set of params is for when I source the script within R and want little output
    #params$deb=TRUE # willl also use functions from debug script
}else if(length(arg) > 0){
    cat("\n\n/////////////////////")
    cat("  arg =  ", paste(unlist(arg), collapse=" ; "))
    cat("  ////////////////////////////////\n")
    for(a in arg){
	  if(grepl("r\\d+$", a)) params$nrun  <- as.numeric(str_extract(a, "\\d+"))
	  else if(a=="win") params$win        <- TRUE
	  #else if(a=="deb"|a=="debug") params$deb        <- TRUE
	  #else if(a=="xdeb"|a=="xdebug") params$xdeb        <- TRUE
	  else if(a=="v1"|a=="vbose1"|a=="vb1") params$vbose <- 1
	  else if(a=="v2"|a=="vbose2"|a=="vb2") params$vbose <- 2
	  else if(a=="v3"|a=="vbose3"|a=="vb3") params$vbose <- 3
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
if(params$j > params$nrun) stop("j must be less than or equal to nrun!")
# if(is.null(params$resp)) stop("no response variable specified!")
if(params$win){
    cat("\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n")
    cat("\n>>>> date & time:", format(Sys.time(), "%d-%b %H:%M")) 
    cat("\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n") 
    params$hdir <- "C:/Users/sarah/Dropbox/Models/fate_glm/"
}
if(params$test){
    #mets <- c("cc", "cart", "caliber", "passive")
    params$outdir <- paste("t",params$outdir,sep="-")
    #params$mdeb <- TRUE
}
####### source files: ##########################################################################
s_files <- c("mice_functions.R", "missing_data.R", "other_imp.R", "imp_sim_functions.R")
#if(params$deb|params$xdeb|params$mdeb) s_files[4] = "debug_imp_sim_func.R" 
if(params$vbose>=1) s_files[4] = "debug_imp_sim_func.R" 
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
cat("\n>>> sourcing files:", paste(s_files,collapse=" ; "))
invisible(lapply(s_files, function(x) source(x)))
#########################################################################################
bprint <- function(x) print(rbind(head(x, 3), tail(x,3)))
debugging <- FALSE # for uickly setting values when working in the file with the functions
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
formulas <- readRDS("form_lists.rds")
metLists <- readRDS("met_lists.rds")
####### imports for making sim data: #######################################################
if(TRUE){
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
}
#########################################################################################
#cat("sort vars or not???")
vars <-  c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") # all vars
prVars <- c("species", "cam_fate", "obs_int", "nest_age", "fdate")
mets <- c("default", "cart", "caliber","passive", "stratify","cf_cc","cc")# don't need full here?
resp_list <- c("is_u", "HF_mis")
#########################################################################################
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
cat("\n>> METHODS:", mets)
cat(sprintf("\t>> TOTAL REPS: %s x %s WITH %s NESTS EACH", length(params$seeds), params$nrun, nnest)) 
if(params$vbose >= 1) cat("\t *** verbosity = ", params$vbose)
#cat(sprintf("\n >>>> ALL VARIABLES: %s ; PREDICTOR VARIABLES: %s", vars, prVars))
#if(params$deb) cat("\t- debug ON")
#if(params$xdeb) cat(" + EXTRA")
#if(params$mdeb) cat("\t ** debug = minimum")
cat("\n>> output will be saved every", params$j, "runs to dir:", now_dir,"\n\n")

for (seed in params$seeds){
    outf <- paste0(now_dir,sprintf("/%s_loggedEvents.out",seed ))
    res <- array(NA, dim = c(length(vars), length(mets), params$nrun, 3, length(mods4sim), length(resp_list)))
    dimnames(res) <- list(as.character(vars),# need to be in same order as vals
                          as.character(mets),
                          as.character(1:params$nrun),
                          c("estimate", "2.5 %","97.5 %"),
                          names(mods4sim),
                          resp_list
                          )
    camFateVars <- vars[grepl(pattern="cam_fate", x=vars, fixed=TRUE)]
    varInfo <- array(NA, dim=c(length(camFateVars), params$nrun, length(mods4sim))) # keep a count of sample size for each category
    dimnames(varInfo) <- list( camFateVars, seq(1,params$nrun), mods4sim)
    cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
    cat(sprintf(">>>>>>>>>>>>>> running simulation %s times \t>>> seed = %s ", params$nrun, seed))
    cat("\t>> & no. imp.:", params$m, "\n")
    for(mod in seq_along(mods4sim)){
        cat("\n\n::::::::::::::::::: MODEL", mods4sim[mod], ":::::::::::::::::::::::::::::::\n\n") ## *~*~*~*~*
        beta_list <- betas[[mod]]
        form_list <- formulas[[names(mods4sim)[mod]]]
        for(run in 1:params$nrun){
            cat(run)
            # repeat this until you get a dataset w/o missing levels??
            #dat4sim <- mkSim(resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, debug=params$deb)
            #dat4sim <- mkResp(seed=run+seed,resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList,  debug=params$deb, mindebug=params$mdeb)
            dat4sim <- mkResp(seed=run+seed,resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, vbose=params$vbose)
            #dat4sim <- dat4sim[[1]] # comment out except when debugging - mkResp returns list(sDat, fitSim1, fitSim2)
            #datNA <- mkSimDat( seeed = run+seed, nd = dat4sim, mpatt=mpatt, wts=ampwt, xdebug=params$xdeb, debug = params$deb, convFact = TRUE, mindebug=params$mdeb)
            datNA <- mkSimDat( seeed = run+seed, nd = dat4sim, mpatt=mpatt, wts=ampwt,convFact = TRUE, vbose=params$vbose)
            #if(params$deb) cat("\n*** datNA:\n")
            #if(params$deb) cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
            #if(params$deb) cat("\n\n<><><><><><><><><><<><><><><<><><><><> AMPUTED DATA: <><><><><><><><><><<><><><><<><><><><>>>\n ")
            #if(params$deb) qvcalc::indentPrint(str(datNA$amp), indent=16)
            #if(params$deb) qvcalc::indentPrint(colSums(is.na(datNA$amp)), indent=8)
            datNA <- datNA$amp
            for(x in seq_along(mets)){ # does matching by index help the trycatch statement?
                skiptoNext <- FALSE
                ## *~*~*~*~*
                #if (params$deb) cat("\n\n<><><><><><><><> method:", mets[x], "<><><><><><><><><><><><><><><><><><><><>") ## *~*~*~*~*
                # the error seems to come from levels with zero observations
                # so maybe ampute() is removing all the camera fates of category X (F in this case)
                tryCatch(
                expr = {
                    vals <- mkImpSim(aDat=datNA,
                                     fullDat=dat4sim,
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
                                   #debug = params$deb,
                                   #mindebug=params$mdeb,
                                   m=params$m, 
                                   vbose=params$vbose,
                                   #xdebug=params$xdeb,
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
            
            #impDat[v,"cc",,,z,resp] <- exp(impDat[v,"cc",,,z,resp]) 
            #cat("\nbefore exp:\n")
            #print(res[,"cc",run,,z,])
            #res[,"cc",run,,z,] <- exp(res[,"cc",run,,z,]) 
            #cat("\nafter exp:\n")
            #print(res[,"cc",run,,z,])
            if(params$test) cat(sprintf("\n\n::::::::::::::::::::::::::::::: ALL OUTPUT - RUN %s, MODEL %s:::::::::::::::::::::::::::::::\n\n",run,mod ))## *~*~*~*~*
            if(params$test) qvcalc::indentPrint(res[,, run,,mod,])
            if(params$test){
                regMet <- brglm2::brglmFit
                iter <- 500
                coefFull <- list()
                #cat("\n\timporting true & sim values")
                coefTrue <- list()
                #mvars <- c()
                for(r in seq_along(resp_list)){
                    resp <- resp_list[r]
                    #fitFull <- glm(as.formula(paste(resp, mods4sim[mod])), data=fullDat, family=binomial, method=regMet, control=brglmControl(maxit=iter))
                    fitFull <- glm(as.formula(paste(resp, mods4sim[mod])), data=dat4sim, family=binomial, method=regMet, control=brglmControl(maxit=iter))
                    if(params$vbose>=3) cat("\n\t get sim coef values (exponentiated)")
                    #coefFull[[r]] <- coef(fitFull)[-1]
                    coefFull[[r]] <- exp(coef(fitFull))[-1]
                    trueDat <- readRDS("dat_complete.rds")
                    #trueDat <- read.csv("dat_complete.csv")
                    trueDat$cam_fate <- relevel(as.factor(trueDat$cam_fate), ref="H")
                    trueDat$species <- relevel(as.factor(trueDat$species), ref="LETE")
                    if(params$vbose>=3) cat("\n\treleveled trueDat")
                    fitTrue <- glm(as.formula(paste(resp,mods4sim[mod])), data=trueDat, family=binomial, method=regMet, control=brglmControl(maxit=iter))
                    if(params$vbose>=3) cat("\n\t get true coef values (exponentiated)")
                    coefTrue[[r]] <- exp(coef(fitTrue))[-1]
                    }
                cat(sprintf("\n\n::::::::::::::::::::::::::::::: MICE OUTPUT COMPARED TO TRUE/SIM COEFS :::::::::::::::::::::::::::::::\n\n",run,mod ))## *~*~*~*~*
               # qvcalc::indentPrint(res[,, run,,mod,])
                mvars <- names(coefTrue[[1]])
                coef1 <- cbind(coefTrue[[1]], coefFull[[1]])
                #out2 <- cbind(coefTrue[[2]][,1], coefFull[[2]][,1], ret[,1,2])
                # this one doesn't have the CIs anyway
                #coef2 <- cbind(coefTrue[[2]], coefFull[[2]])
                coef2 <- cbind(coefFull[[2]], coefTrue[[2]])
                colnames(coef1) <- c("sim", "true")
                colnames(coef2) <- c("sim", "true")
                if(params$vbose >=3) cat("\n\tcheck match:\n")
                if(params$vbose >=3) print(coefTrue[[1]])
                ##print(coef1[mvars,])
                if(params$vbose >=3) print(coef1)
                #print(res[mvars,,run,1,mod,1])
                if(params$vbose >=3) print(res[,,run,1,mod,1])
                #out1 <- sapply(res)
                # vars should be in same order:
                # make sur the vars are in th same order:
                out1 <- cbind(res[mvars,,run,1,mod,1], coef1[mvars,]) # 1=estimate, 1=is_u
                out2 <- cbind(res[mvars,,run,1,mod,2], coef2[mvars,]) # 1=estimate, 2=HF_mis
                cat("\n\tMICE output compared to true and sim values:\n")
                cat("\n")
                qvcalc::indentPrint(out1, indent=8)
                cat("\n")
                qvcalc::indentPrint(out2, indent=8)
            }
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


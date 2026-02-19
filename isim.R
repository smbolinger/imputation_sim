

library(brglm2) # penalized logistic regression
library(CALIBERrfimpute) # could access using ::
library(gt)
library(gtsummary)
library(data.table)
# suppressMessages(library(mice))
suppressPackageStartupMessages(library(mice))
suppressMessages(library(tidyverse)) # look into tidytable for tidyverse syntax w/ data.table
options(width=180)
options(scipen=999) # discourage R from using scientific notation?
source("/home/wodehouse/.local/bin/r_func.R")
source("config.R")

###### PARAMS & ARGS ###########################################################################
arg <- commandArgs(trailingOnly=TRUE)
# print(str(arg))
arg <- unlist(strsplit(commandArgs(trailingOnly=TRUE), split=" "))
# print(str(arg))
if(length(arg)==0){
    cat("\n\n/////////////////////// ** NOTE ** no arguments provided - using defaults of nrun=1, m=5, verbose=FALSE //////////////////////////\n\n")
    # isn't the same as "test" - this set of params is for when I source the script within R and want little output
    params$nrun=1
    params$j=1
    params$m=5
    params$seeds = c(71358)
    params$nNest = 50
    params$vbose = 1
    #params$deb=TRUE # willl also use functions from debug script
}else if(length(arg) > 0){
    cat("\n\n/////////////////////")
    #cat("  arg =  ", paste(unlist(arg), collapse=" ; "))
    cat("  arg =  ", paste(unlist(arg), collapse=" ; "))
    cat("  ////////////////////////////////\n")
    for(a in arg){
        cat("arg =",a)
        if(grepl("r\\d+$", a)) params$nrun  <- as.numeric(str_extract(a, "\\d+"))
        else if(a=="win") params$win        <- TRUE
        else if(a=="ipl") params$ipl       <- TRUE
        else if(a=="strat"|a=="str") params$strat_sim <- TRUE
        else if(a=="nostrat"|a=="nostr") params$strat_sim <- FALSE
        else if(a=="test") params$test <- TRUE
        else if(grepl("suff_\\w+", a)) params$suffix <- str_extract(a, "(?<=suff_)\\w+")
        else if(grepl("aw\\d", a)) params$ampWt <- as.numeric(str_extract(a, "\\d"))
        else if(grepl("nn\\d", a)) params$nNest <- as.numeric(str_extract(a, "\\d+"))
        else if(grepl("j\\d+", a)) params$j <- as.numeric(str_extract(a, "\\d+"))
        else if(grepl("s\\d+", a)) params$seeds <- c(as.numeric(str_extract(a, "\\d+")))
        else if(grepl("m\\d+", a)) params$m <- as.numeric(str_extract(a, "\\d+"))
        else if(grepl("v\\d+", a)) params$vbose <- as.numeric(str_extract(a, "\\d+"))
        #else if(grepl("suff_\\w+", a)) params$suffix <- as.numeric(str_extract(a, "(?<=suff_)\\w+"))
        # else if(a=="is_u" | a == "HF_mis") params$resp       <- a
        #else if(a=="deb"|a=="debug") params$deb        <- TRUE
        #else if(a=="xdeb"|a=="xdebug") params$xdeb        <- TRUE
        #else if(a=="v1"|a=="vbose1"|a=="vb1") params$vbose <- 1
        #else if(a=="v2"|a=="vbose2"|a=="vb2") params$vbose <- 2
        #else if(a=="v3"|a=="vbose3"|a=="vb3") params$vbose <- 3
    }
    #print(params)
    cat("\n\nPARAMS:",paste(names(params), params, sep="=", collapse="    |    "))
}
#cat("params$j & params$nrun:", params$j, params$nrun)
# if(params$j > params$nrun) stop(sprintf("j (%s) must be less than or equal to nrun (%s)!", params$j, params$nrun))

####### SOURCE FILES: ##########################################################################
s_files <- c("mice_functions.R", "missing_data.R", "other_imp.R")
# if(params$test) s_files[4] = "debug_imp_sim_func.R" 
s_paths <- paste0("functions/", s_files)
if(params$test) s_paths[4] = "debug_imp_sim_func.R" else s_paths[4] = "clean_imp_sim_func.R"
# s_paths[4] = "debug_imp_sim_func.R" 
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
# cat("\n>>> sourcing files:", paste(s_files,collapse=" ; "))
cat("\n>>> sourcing files:", paste(s_paths,collapse=" ; "))
invisible(lapply(s_paths, function(x) source(x)))

#########################################################################################
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
cat("\n>> METHODS:", mets)
cat(sprintf("\t>> TOTAL REPS: %s x %s WITH %s NESTS EACH; \n", length(params$seeds), params$nrun, nnest)) 
if(params$test){
    cat("\t.:.:.:. TEST MODE .:.:.:. \t")
    #mets <- c("cc", "cart", "caliber", "passive")
    params$outdir <- paste("t",params$outdir,sep="-")
    #params$mdeb <- TRUE
}
cat(sprintf("\t *** stratify = %s ***\t\t *** verbosity = %s ***", params$strat_sim,params$vbose))
suffix <- paste(sprintf("%sruns", params$nrun), params$suffix, sep="_")
now_dir <- paste(params$hdir, params$outdir, sep="out/")
outf <- paste0(now_dir,sprintf("/%s_loggedEvents.out",seed ))
if(!dir.exists(now_dir)) dir.create(now_dir)
cat("\n>> output will be saved every", params$j, "runs to dir:", now_dir,"\n\n")
# writeLines(now_dir, "last_out.txt")

########## MAKE MATRIX ###################################################################
res <- array(NA, dim=c(length(vars), length(mets)+2, params$nrun, 3, length(resp_list)))
dimnames(res) <- list(as.character(vars), 
                      c(as.character(mets), "sim", "true"),
                      as.character(1:params$nrun),
                      c("estimate", "2.5 %", "97.5 %"),
                      as.character(resp_list)
)

########## LOOP THROUGH RUNS ###################################################################
for(run in 1:params$nrun){
    if(params$vbose>=1) cat("\n\n+ + + + + + + + + + + + + + + + + + + RUN: ")
    cat(run)
    if(params$vbose>=1) cat(" + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +")
    ########## MAKE SIM DATA ###################################################################
    # resp  <- resp_list[1]
    # dat4sim <- mkResp(seed=run+seed,resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, strat=params$strat_sim, vbose=params$vbose)
    dat4sim <- mkResp(seed=run+seed,resp_list, oneMod, nnest, cMat, mList, beta_list, fprob, sprob, prList, strat=params$strat_sim, vbose=vb$mk_sim)
    datNA <- mkAmpDat( seeed = run+seed, nd = dat4sim, mpatt=mpatt, test=params$test, wts=ampwt,convFact = TRUE, vbose=vb$mk_amp) #if(params$deb) cat("\n*** datNA:\n")
    datNA <- datNA$amp

    ########## IMPUTATION ###################################################################
    for(x in seq_along(mets)){ # does matching by index help the trycatch statement?
        skiptoNext <- FALSE
        tryCatch(
        expr = {
            vals <- mkImpSim(aDat=datNA,
                             fullDat=dat4sim,
                           #cols=col_list,
                           # pr_list=prVars,
                           # resp=params$resp, 
                           resp_list=resp_list,
                           # resp_list=c("is_u"),
                           form_list =form_list,
                           met_list=metLists[,,,mod],
                           vars=vars,
                           outFile=outf,
                           #mods=mods4sim,
                           # modd=mods4sim[mod],
                           modd=oneMod,
                           met=mets[x],
                           #debug = params$deb,
                           #mindebug=params$mdeb,
                           m=params$m, 
                           # vbose=params$vbose,
                           vbose=vb$mk_imp,
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
        # res[, mets[x], run,,mod,]  <- vals
        res[, mets[x], run,,]  <- vals
    }

    if(TRUE){
        coefFull <- list()
        coefTrue <- list()
        for(r in seq_along(resp_list)){
            #if(params$vbose >=2) cat("\n\n    <><><><><><><><><> TRUE / SIM COEFS <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n\n")
            resp <- resp_list[r]
            fitTrue <- glm(as.formula(paste(resp,oneMod[mod])),
                           data=trueDat,
                           family=binomial,
                           method=brglm2::brglmFit,
                           control=brglmControl(maxit=500))
            fitFull <- glm(as.formula(paste(resp, oneMod[mod])),
                           data=dat4sim,
                           family=binomial,
                           method= brglm2::brglmFit,
                           control=brglmControl(maxit=500))
            # if(params$vbose>=3){
            #     cat("\n\t trueDat & fit:\n")
            #     qvcalc::indentPrint(names(trueDat),indent=12)
            #     qvcalc::indentPrint(summary(fitTrue))
            #     cat("\n\t dat4sim & fit:\n")
            #     qvcalc::indentPrint(names(dat4sim))
            #     qvcalc::indentPrint(summary(fitFull),indent=12)
            # }
            # if(params$vbose>=3) cat("\n\t get sim coef values ")
            #coefFull[[r]] <- exp(coef(fitFull))[-1]
            coefFull[[resp]] <- cbind(coef(fitFull)[-1], confint(fitFull)[-1,])
            #coefFull[[resp]] <- cbind(exp(coef(fitFull))[-1], exp(confint(fitFull))[-1,])
            # does having exp here prevent the rownames from carrying over?
            #vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
            # if(params$vbose>=3) cat("\t get true coef values ")
            coefTrue[[resp]] <- cbind(coef(fitTrue)[-1], confint(fitTrue)[-1,])
            }
        mvars <- rownames(coefTrue[[1]])
        # if(params$vbose >=3){
        #     cat("\n    >> print sim/true coef vals with CIs (not exp):\n")
        #     qvcalc::indentPrint(coefFull, indent=16)
        #     qvcalc::indentPrint(coefTrue, indent=16)
        #     cat("\n\tcheck match for is_u & true:\n")
        #     qvcalc::indentPrint(coefTrue[[1]][mvars,], indent=8)
        #     qvcalc::indentPrint(res[mvars,"true",run,1,1], indent=8)
        # }
        #cat(sprintf("\n\n::::::::::::::::::::::::::::::: MICE OUTPUT COMPARED TO TRUE/SIM COEFS :::::::::::::::::::::::::::::::\n\n",run,mod ))## *~*~*~*~*
       # qvcalc::indentPrint(res[,, run,,mod,])
        ##print(coef1[mvars,])
        #if(params$vbose >=3) print(coef1)
        #print(res[mvars,,run,1,mod,1])
        # if(params$vbose >=3) qvcalc::indentPrint(res[mvars,"true",run,1,mod,1], indent=8)
        # vars should be in same order:
        res[mvars,"sim",run,,"is_u"] <-  coefFull[[1]][mvars,]
        #res[mvars,"sim",run,,mod,"is_u"] <-  exp(coefFull[[1]][mvars,])
        res[mvars,"sim",run,,"HF_mis"] <-  coefFull[[2]][mvars,]
        #res[mvars,"sim",run,,mod,"HF_mis"] <-  exp(coefFull[[2]][mvars,])
        #res[mvars,"true",run,,mod,"is_u"] <-  exp(coefTrue[[1]][mvars,])
        #res[mvars,"true",run,,mod,"HF_mis"] <-  exp(coefTrue[[2]][mvars,])
        res[mvars,"true",run,,"is_u"] <-  coefTrue[[1]][mvars,]
        res[mvars,"true",run,,"HF_mis"] <-  coefTrue[[2]][mvars,]

        # if(params$vbose >=3) {
        #     cat("\n\texponentiated values:\n")
        #     qvcalc::indentPrint(res[mvars,,run,,])
        # }
    }

    ###### PRINT ##########################################################################
    # if(params$test) cat(sprintf("::::::::::::::::::::::::::::::: ALL OUTPUT - RUN %s, MODEL %s:::::::::::::::::::::::::::::::\n",run,mod ))## *~*~*~*~*
    if(params$test) {
        cat(sprintf("::::::::::::::::::::::::::::::: ALL OUTPUT - RUN %s, MODEL: %s :::::::::::::::::::::::::::::::\n",run,mod))## *~*~*~*~*
        qvcalc::indentPrint(round(res[,, run,,],3))
        cat("\n    >> cam fate vars for this run:\n")
        qvcalc::indentPrint(table(datNA$cam_fate), indent=8)
        cat("\n >>> 'estimate' for entire res:\n")
        qvcalc::indentPrint(round(res[,,,"estimate",],3))
    }
    # if(params$test) qvcalc::indentPrint(round(res[,, run,,,],2))
    # fname <- paste0(now_dir, sprintf("/vals_mod%s_seed%s_%s_%s.rds", mod, seed, nowtime, suffix))

}

fname <- paste0(now_dir, sprintf("/vals_mod%s_seed%s_%s_%s.rds", mod, seed, nowtime, suffix))
writeLines(fname, "last_fname.txt")
# if(params$vbose>=3) cat("\n >>> 'estimate' for entire res:\n")
# if(params$vbose>=3) qvcalc::indentPrint(res[,,,"estimate",])
cat("\n ~ ~ ~ saving output to file:", fname)
saveRDS(round(res,3), fname)

########## CALCULATE BIAS #####################################################################

biasVals <- calcBias(res, mods4sim[1], vars, mets, resp_list, vbose=vb$calc_bias)
# bfname <- paste0(params$hdir, sprintf("out/bias/bias_vals_seed%s_%s.rds", seed, nowtime))
# saveRDS(round(biasVals, 3), bfname)


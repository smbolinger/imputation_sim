

library(brglm2) # penalized logistic regression
library(CALIBERrfimpute) # could access using ::
library(gt)
library(gtsummary)
library(data.table)
# suppressMessages(library(mice))
suppressPackageStartupMessages(library(mice))
suppressMessages(library(tidyverse)) # look into tidytable for tidyverse syntax w/ data.table
options(width=130)
options(scipen=999) # discourage R from using scientific notation?

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
               #strat_sim=FALSE,
               strat_sim=TRUE,
               nNest=50,
               j = 50,
               m=20
)
arg <- commandArgs(trailingOnly=TRUE)
print(str(arg))
arg <- unlist(strsplit(commandArgs(trailingOnly=TRUE), split=" "))
print(str(arg))

####### SOURCE FILES: ##########################################################################
s_files <- c("mice_functions.R", "missing_data.R", "other_imp.R", "imp_sim_functions.R")
if(params$test) s_files[4] = "debug_imp_sim_func.R" 
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
cat("\n>>> sourcing files:", paste(s_files,collapse=" ; "))
invisible(lapply(s_files, function(x) source(x)))

######## MODELS/METHODS #####################################################################
modList <- readLines("modList.txt")
mods4sim <- modList[c(1,8,16) ]
names(mods4sim) <- c("m1", "m8", "m16")
formulas <- readRDS("form_lists.rds")
metLists <- readRDS("met_lists.rds")

######## NEST DATA #####################################################################
if(TRUE){
    trueDat <- readRDS("dat_complete.rds")
    if(params$vbose>=2) cat("\nrelevel cam_fate\n")
    trueDat$cam_fate <- relevel(as.factor(trueDat$cam_fate), ref="H")
    # if(params$vbose>=2) cat("\n    relevel species\n")
    # trueDat$species <- relevel(as.factor(trueDat$species), ref="LETE")
    if(params$vbose>=2) cat("\nmake categorical vars factors:\n")
    trueDat$species <- as.factor(trueDat$species)
    trueDat$is_u <- as.factor(trueDat$is_u)
    trueDat$HF_mis <- as.factor(trueDat$HF_mis)
    if(params$vbose>=2) qvcalc::indentPrint(summary(trueDat), indent=4)
    if(params$vbose>=2) qvcalc::indentPrint(str(trueDat), indent=4)
}

######## SIMULATION INPUT #####################################################################
if(TRUE){
    betas <- readRDS("betas.rds") # this is a list of vectors or something
    mList <- readRDS('means.rds')
    cMat  <- readRDS('cormat.rds')
    if(params$vbose>=2) cat("\nbetas, means, & matrix:\n")
    if(params$vbose>=2) qvcalc::indentPrint(betas, indent=4)
    if(params$vbose>=2) qvcalc::indentPrint(mList, indent=4)
    if(params$vbose>=2) qvcalc::indentPrint(cMat, indent=4)
    fprob <- c( 'H'=0.53,'A'=0.1, 'D'=0.13,'F'=0.06, 'Hu'=0.13,'S'=0.06 )
    sprob <- c('CONI'=0.32, 'LETE'=0.68)
    #nnest <- ifelse(params$test==T, 200, 200) # number of simulated nests
    nnest <- params$nNest
    mpatt  <- readRDS("misPatt.rds")
    awFile <- case_when(params$ampWt==1 ~ "ampWts.rds",
                        params$ampWt==2 ~ "ampWts2.rds",
                        params$ampWt==3 ~ "ampWts3.rds")
    ampwt <- readRDS(awFile)
    cat("\n\n<><><><><><><><><><> Missingness pattern & variable weights: <><><><><><><><><><><><><><><><><><><><>\n\n")
    print(mpatt)
    print(ampwt)
}

######## VAR/METHOD LISTS ##################################################################### cat("sort vars or not???")
vars <-  c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") # all vars
prVars <- c("species", "cam_fate", "obs_int", "nest_age", "fdate")
# mets <- c("default", "cart", "caliber","passive", "stratify","cf_cc","cc")# don't need full here?
mets <- c("cart","passive", "cc")# don't need full here?
resp_list <- c("is_u", "HF_mis")

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
if(!dir.exists(now_dir)) dir.create(now_dir)
cat("\n>> output will be saved every", params$j, "runs to dir:", now_dir,"\n\n")
# writeLines(now_dir, "last_out.txt")

########## MAKE SIM DATA ###################################################################
resp  <- resp_list[1]
dat4sim <- mkResp(seed=run+seed,resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, strat=params$strat_sim, vbose=params$vbose)
datNA <- mkAmpDat( seeed = run+seed, nd = dat4sim, mpatt=mpatt, test=params$test, wts=ampwt,convFact = TRUE, vbose=params$vbose) #if(params$deb) cat("\n*** datNA:\n")
datNA <- datNA$amp

########## MAKE MATRIX ###################################################################
res <- array(NA, dim=c(length(vars), length(mets), params$nrun, 3))
dimnames(res) <- list(as.character(vars), 
                      as.character(mets),
                      as.character(1:params$nrun),
                      c("estimate", "2.5 %", "97.5 %"))

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

if(TRUE){
    coefFull <- list()
    coefTrue <- list()
    for(r in seq_along(resp_list)){
        if(params$vbose >=2) cat("\n\n    <><><><><><><><><> TRUE / SIM COEFS <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n\n")
        resp <- resp_list[r]
        fitTrue <- glm(as.formula(paste(resp,mods4sim[mod])),
                       data=trueDat,
                       family=binomial,
                       method=brglm2::brglmFit,
                       control=brglmControl(maxit=500))
        fitFull <- glm(as.formula(paste(resp, mods4sim[mod])),
                       data=dat4sim,
                       family=binomial,
                       method= brglm2::brglmFit,
                       control=brglmControl(maxit=500))
        if(params$vbose>=3) cat("\n\t get sim coef values ")
        #coefFull[[r]] <- exp(coef(fitFull))[-1]
        coefFull[[resp]] <- cbind(coef(fitFull)[-1], confint(fitFull)[-1,])
        #coefFull[[resp]] <- cbind(exp(coef(fitFull))[-1], exp(confint(fitFull))[-1,])
        # does having exp here prevent the rownames from carrying over?
        #vals <- cbind(coef(fit)[-1], confint(fit)[-1,]) # confint has 2 columns, so need comma
        if(params$vbose>=3) cat("\t get true coef values ")
        coefTrue[[resp]] <- cbind(coef(fitTrue)[-1], confint(fitTrue)[-1,])
        }
    if(params$vbose >=3) cat("\n    >> print sim/true coef vals with CIs (not exp):\n")
    if(params$vbose >=3) qvcalc::indentPrint(coefFull, indent=16)
    if(params$vbose >=3) qvcalc::indentPrint(coefTrue, indent=16)
    #cat(sprintf("\n\n::::::::::::::::::::::::::::::: MICE OUTPUT COMPARED TO TRUE/SIM COEFS :::::::::::::::::::::::::::::::\n\n",run,mod ))## *~*~*~*~*
   # qvcalc::indentPrint(res[,, run,,mod,])
    mvars <- rownames(coefTrue[[1]])
    if(params$vbose >=3) cat("\n\tcheck match for is_u & true:\n")
    if(params$vbose >=3) qvcalc::indentPrint(coefTrue[[1]][mvars,], indent=8)
    ##print(coef1[mvars,])
    #if(params$vbose >=3) print(coef1)
    #print(res[mvars,,run,1,mod,1])
    if(params$vbose >=3) qvcalc::indentPrint(res[mvars,"true",run,1,mod,1], indent=8)
    # vars should be in same order:
    res[mvars,"sim",run,,mod,"is_u"] <-  coefFull[[1]][mvars,]
    #res[mvars,"sim",run,,mod,"is_u"] <-  exp(coefFull[[1]][mvars,])
    res[mvars,"sim",run,,mod,"HF_mis"] <-  coefFull[[2]][mvars,]
    #res[mvars,"sim",run,,mod,"HF_mis"] <-  exp(coefFull[[2]][mvars,])
    #res[mvars,"true",run,,mod,"is_u"] <-  exp(coefTrue[[1]][mvars,])
    #res[mvars,"true",run,,mod,"HF_mis"] <-  exp(coefTrue[[2]][mvars,])
    res[mvars,"true",run,,mod,"is_u"] <-  coefTrue[[1]][mvars,]
    res[mvars,"true",run,,mod,"HF_mis"] <-  coefTrue[[2]][mvars,]

    if(params$vbose >=3) {
        cat("\n\texponentiated values:\n")
        qvcalc::indentPrint(res[mvars,,run,,mod,])

    }
}

###### PRINT ##########################################################################
# if(params$test) cat(sprintf("::::::::::::::::::::::::::::::: ALL OUTPUT - RUN %s, MODEL %s:::::::::::::::::::::::::::::::\n",run,mod ))## *~*~*~*~*
if(params$test) cat(sprintf("::::::::::::::::::::::::::::::: ALL OUTPUT - RUN %s, MODEL: %s %s :::::::::::::::::::::::::::::::\n",run,mod ))## *~*~*~*~*
if(params$test) qvcalc::indentPrint(round(res[,, run,,mod,],2))
if(params$vbose>=1) cat("\n    >> cam fate vars for this run:\n")
if(params$vbose>=1) qvcalc::indentPrint(table(datNA$cam_fate), indent=8)

fname <- paste0(now_dir, sprintf("/vals_mod%s_seed%s_%s_%s.rds", mod, seed, nowtime, suffix))
cat("\n ~ ~ ~ saving output to file:", fname)
saveRDS(round(res,3), fname)


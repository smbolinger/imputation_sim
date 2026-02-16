
library(brglm2)
library(tidyverse)
library(gt)
naRm <- TRUE
unexp <- FALSE
options(width=190)
source("imp_sim_functions.R")
arg <- commandArgs(trailingOnly=TRUE)
#if(length(arg)==0) stop("needs name of output directory within out directory, and nrun value to use; add debug as final argument to get print statements")
if(length(arg)==0) stop("needs name of output directory within out directory, and nrun value to use; (optional) add verbosity level prepended by 'v'")
#if(length(arg)==0) stop("needs name of output directory within out directory")
#cat("\n===========================================================================")
cat("\n\n /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ /|||\\ ")
cat("\n++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++\n")
cat("\n\t>> date & time:",format(Sys.time(), "%d-%b %H:%M"))
params <- list(hdir = "/home/wodehouse/Projects/fate_glm/", 
               #outdir = arg[1],
	       dir_ext = paste0("out/",arg[1]),
               #nrun    = 100,
               # nrun    = arg[2],
               nrun    = as.numeric(str_extract(arg[2], "\\d+")),
	       #debug   = FALSE)
	       vbose   = 1)

#if(length(arg)>2) params$debug=TRUE
#if(length(arg)>2) params$vbose=arg[3]
if(length(arg)>2) params$vbose = str_extract(arg[3], "\\d+")

#cat("\n===========================================================================")
cat("\t           *** CALCULATE BIAS VALUES - calc_bias.R ***")
#cat("\n===========================================================================")
cat(sprintf("\n\n\t>> *** verbosity = %s *** \t>> model output from: %s \t w/ nrun = %s\n", params$vbose, params$dir_ext, params$nrun))
cat("\n++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++\n")
cat(" \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \\|||/ \n\n")
#cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
suffix <- ""

rlist <- c("is_u", "HF_mis")
vars <-  c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") # all vars
#mets <- c("default","pmm", "rf", "cart", "caliber","passive","stratify","cf_cc","cc")# don't need full here?
mets <- c("default", "cart", "caliber","passive","stratify","cf_cc","cc")# don't need full here?
biasVals <- c("value","bias", "pctBias", "covRate", "avgWidth", "RMSE", "SD")

modList <- readLines("modList.txt")
mods4sim <- modList[c(1,8,16) ]
names(mods4sim) <- c("m1", "m8", "m16")
cat("\n>> models:\n")
qvcalc::indentPrint(mods4sim, indent=4)

#patt <- "100runs|250runs"
patt <- "vals"
#fnlist <- list.files(path = paste0(params$hdir, params$dir_ext), pattern = "^runs0to50.*|^runs50to100.*", full.names =T)
fnlist <- list.files(path = paste0(params$hdir, params$dir_ext), pattern = patt, full.names =T)
#if(params$debug) cat("\n\t>> filename list:")
cat("\n>> filename list:\n")
qvcalc::indentPrint(fnlist, indent=4)
flist <- lapply(fnlist, FUN=function(x) readRDS(x))
seeds <- unlist(lapply(fnlist, function(x) str_extract(x, "(?<=seed)\\d+")))
modnums <- unlist(lapply(fnlist, function(x) str_extract(x, "(?<=mod)\\d+")))
runs <- unlist(lapply(fnlist, function(x) str_extract(x, "\\d+(?=runs)")))
#mods <- names(mods4sim)[as.numeric(modnums)]
mods <- mods4sim[as.numeric(modnums)]
cat("\n>> seeds from the files:\n")
qvcalc::indentPrint(seeds)
qvcalc::indentPrint(class(seeds))
cat("\n>> models from the files:\n")
#qvcalc::indentPrint(mods, indent=4)
qvcalc::indentPrint(paste(names(mods), mods, sep=": "), indent=4)
cat("\n>> nruns from the files:\n")
qvcalc::indentPrint(runs)
#print(str(flist))
#print(dim(flist))
#fl1 <- flist[[1]]
#str(fl1)
if(params$vbose > 1) cat(sprintf("\n\n    >>files with %s runs: \n", params$nrun))
fnlist2 <- fnlist[runs==params$nrun]
if(params$vbose > 1) qvcalc::indentPrint(fnlist2, indent=8)
flist2 <- lapply(fnlist2, FUN=function(x) readRDS(x))
if(params$vbose > 1) qvcalc::indentPrint(str(flist2), indent=8)

#dat4sim <- read.csv("dat_complete.csv", stringsAsFactors = TRUE)
dat4sim <- readRDS("dat_complete.rds")
if(params$vbose>=2) cat("\nrelevel cam_fate:\n")
dat4sim$cam_fate <- relevel(as.factor(dat4sim$cam_fate), ref="H")
# if(params$vbose>=2) cat("\nrelevel species:\n")
# dat4sim$species <- relevel(as.factor(dat4sim$species), ref="LETE")
if(params$vbose>=2) cat("\nmake categorical vars factors:\n")
dat4sim$species <- as.factor(dat4sim$species)
dat4sim$is_u <- as.factor(dat4sim$is_u)
dat4sim$HF_mis <- as.factor(dat4sim$HF_mis)
# dat4sim$cam_fate <- relevel(dat4sim$cam_fate, ref="H") # make 'H' the reference category
# dat4sim$species <- relevel(dat4sim$species, ref="LETE")
# levels(dat4sim$HF_mis) <- c(0,1)
# levels(dat4sim$is_u)   <- c(0,1)
fullDat <- dat4sim
cat("\n>>>> releveled categorical variables\n")
qvcalc::indentPrint(summary(dat4sim), indent=8)
cat("\n>>> species vs other variables:\n")
qvcalc::indentPrint(xtabs(~ is_u + species, data=fullDat))
qvcalc::indentPrint(xtabs(~ HF_mis + species, data=fullDat))
fam <- binomial
regMet <- brglm2::brglm_fit
iter <- 500

now_dir = paste0(params$hdir, "out/bias/", format(Sys.time(), "%d%b"))
if(!dir.exists(now_dir)) dir.create(now_dir)
cat("\n\n>>>> output directory:", now_dir, "exists?", dir.exists(now_dir))
cat("\n")

#resp_list <- unlist(lapply(fnlist, function(x) str_extract(x, paste(rlist,collapse="|"))))
#print(resp_list)

#for(resp in seq_along(resp_list)){
#for(resp in unique(resp_list)){

#for(resp in rlist)
#for(seed in seeds){
    #cat(sprintf("\n<> <> <> <> <> SEED: %s <> <> <> <> <> \n", seed))
    #cat("\n>> matrices to bind:\n", seeds==seed)
    #impDat <- abind::abind(flist[seeds=seed], along=3)
    #cat("\n>> matrices combined:\n")
    #print(str(impDat))

## Group the different seed output together by model 
for(z in seq_along(mods4sim)){
    mod     <- mods4sim[z]
    if(params$vbose > 1){
        cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
        cat("\n\t\t MODEL:", mods4sim[z])
        cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
        cat("\n\n    >> files for this model:", mods==mod)
    }
    if(params$vbose > 1) cat("\n\n    >> number of files for this model:", length(fnlist[mods==mod]), "\n")
    if(params$vbose > 1) cat(sprintf("\n\n    >> number of files for this model for %s runs: %s\n", params$nrun, length(fnlist2[mods==mod])))
    if(length(flist2[mods==mod])<1) stop(sprintf("No files matching model %s. Exiting script.", mod))
    if(params$vbose > 2) qvcalc::indentPrint(flist2[mods==mod], indent=8)
    if(params$vbose > 2) cat("\n    >>> attempt to merge:\n")
    impDat <- abind::abind(flist2[mods==mod], along=3)
    ## double bracket doesn't work'
    #impDat <- abind::abind(flist[[mods==mod]], along=3)
    if(params$vbose >= 2) cat("\n\n    >> matrices combined:\n")
    # if(params$vbose >= 2) qvcalc::indentPrint(head(impDat), indent=8)
    if(params$vbose > 2) qvcalc::indentPrint(str(impDat), indent=8)
    if(params$vbose >= 2) qvcalc::indentPrint(dimnames(impDat), indent=8)

    impDatfile <-  sprintf("%s/impdat_vals.rds",now_dir)
    cat(sprintf("\n>>> saving to file: %s ", impDatfile))
    saveRDS(impDat, impDatfile)


    if(unexp==TRUE){
        cat("\n***NOTE*** for impDat, using transformed values from file.\n")
        impDat <- readRDS("ln_out.rds")
    }

    ## Separate by response variable and get the real model output
    for(resp in rlist){
        cat(sprintf("\n<> <> <> <> <> MODEL: %s %s <> <> <> <> <> \n", resp, mod))
        #cat(sprintf("\n\n============== MODEL: %s %s ==================\n",resp,mod))
        fitReal <- glm(as.formula(paste0(resp, mod)),
                       data=fullDat,
                       family=fam,
                       method=regMet,
                       control=brglmControl(maxit=iter) )

        # don't want coef - want coefficients from summary.glm - need exp(coef())?
        #trueVals <- exp(coef(fitReal)[vars]) # the coefs have names associated with them
        trueVals <- coef(fitReal)[vars] # the coefs have names associated with them
        if(params$vbose > 1) qvcalc::indentPrint(summary(fitReal), indent=8) 
        cat("\n    true values (not exp):\n")
        qvcalc::indentPrint(trueVals, indent=8)
        #cat("\n    avg sim values:\n")
        #qvcalc::indentPrint(rowMeans(impDat[,"sim", , "estimate", z, resp]), indent=8)

        bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(mods4sim), length(rlist) ) )
        dimnames(bias) <- list(sort(as.character(vars)),
                               # c("pmm", "rf", "cart"),
                               as.character(mets),
                               as.character(biasVals),
                               names(mods4sim),
                               rlist)
        if(params$vbose > 1) cat("\n    >>> empty array to store bias values:" )
        if(params$vbose > 2) qvcalc::indentPrint(str(bias), indent=8)
        if(params$vbose > 1) qvcalc::indentPrint(dimnames(bias), indent=8)

        ## Loop through the predictor variables and store the bias values to the matrix
        for(v in vars){
            cat(sprintf("\n    ----- VARIABLE: %s ------------------------------------------------\n", v))
            cat("\n    >--> sim values:\n")
            simVals <- impDat[v,"sim", , "estimate", z, resp]
            qvcalc::indentPrint(simVals, indent=8)
            if(params$vbose > 2) cat(sprintf("\n    >> output for variable %s:\n", v))
            #if(params$debug) print(impDat[v,,,,z,resp])
            #if(params$debug) print(str(impDat))
            #if
            # now I'm doing this in the other script (debug_whatever.R)
            #impDat[v,"cc",,,z,resp] <- exp(impDat[v,"cc",,,z,resp]) 
            #if(params$vbose > 2) qvcalc::indentPrint(str(impDat[v,,,,z,resp]), indent=12)
            if(params$vbose > 2) qvcalc::indentPrint(impDat[v,,,"estimate",z,resp], indent=8)
            avg <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
            #print(impDat[v,,1:5,1,z])
            true <- trueVals[v]
            true2 <- avg["true",]
            simVal <- avg["sim",]
            if(params$vbose > 1){
                cat("\n\taverages - true, true, sim:\n")
                qvcalc::indentPrint(true, indent=12)
                qvcalc::indentPrint(true2, indent=12)
                #qvcalc::indentPrint(simVal, indent=12)
                qvcalc::indentPrint(simVals, indent=12)
            }
            avg    <- avg[mets,] # remove true and sim, which are not in "mets"
            #impDat <- impDat[,mets,,,,]
            # calculate std dev after removing 'sim' and 'true':
            #cat("mets = ", mets, dim(mets))
            sdev <- apply(impDat[v, mets, , ,z,resp],MARGIN=c(1,3),FUN = sd, na.rm=TRUE)
            #cat(sprintf("\nbias for %s, value, and %s:\n", v, z))
            #print(bias[v,,"value",z])
                  #cat("\nfirst 5 reps:\n")
                  #print(impDat[v,,1:5,,z])
            true <- simVals
            #if(params$vbose > 1) cat("\tset 'true' to the sim values (one for each run):", true[1,])
            if(params$vbose > 1) cat("\n\tset 'true' to the sim values (one for each run):", true)
            if(params$vbose > 1){
                  cat("\n\n    >--> estimated values:")
                  cat("\n\n\tdimensions of impDat[v,,,,z,resp]:\n")
                  qvcalc::indentPrint(dim(impDat[v,,,,z,resp]), indent=12)

                  cat(sprintf("\n\taverage of estimate for mod %s:\n",z))
                  #print(apply(impDat[v,,,,z], FUN=mean, MARGIN=c(1,3)))
                  #print(impDat[v,,,z])
                  qvcalc::indentPrint(apply(impDat[v,,,,z,resp], FUN=mean, MARGIN=c(1,3)), indent=12)

                  cat("\n\tavg:\n")
                  qvcalc::indentPrint(avg, indent=12)
            }
            #print(str(avg))
            #if(params$vbose > 1) qvcalc::indentPrint(bias[v,,"value",z,resp])
            #if(params$vbose > 1) qvcalc::indentPrint(avg[,"estimate"])
            if(FALSE){ # compare the matrix dimensions to dimensions of vals to store
                mat1 = bias[v,,"bias",z,resp]
                #mat2 = rowMeans(impDat[v,,,"estimate",z,resp] - true, na.rm=naRm)
                mat2 = sdev[,"estimate"]
                if(params$vbose > 1) qvcalc::indentPrint(mat1)
                if(params$vbose > 1) qvcalc::indentPrint(mat2)
            }
            if(params$vbose > 2) cat("\n\t>>> storing vals to matrix\n")
            bias[v, ,"value",z,resp] <- avg[,"estimate"]
            bias[v,,"SD",z,resp] <- sdev[,"estimate"]
            qvcalc::indentPrint(impDat[v,mets,,"estimate",z,resp], indent=8)
            #bias[v,,"bias",z,resp] <- avg[,"estimate"] - true
            #bias[v,met,"bias",z,resp] <- rowMeans(impDat[v,mets,,"estimate",z,resp] - true, na.rm=naRm) # exclude sim and true by matching to 'mets'
            #bias[v,met,"bias",z,resp] <- rowMeans(impDat[v,met,,"estimate",z,resp] - true, na.rm=naRm)
            qvcalc::indentPrint(rowMeans(impDat[v,mets,,"estimate",z,resp] - true, na.rm=naRm), indent=12)
            # could also apply() over all the vars
            ## as in out <- apply(m1new, 2, function(x) rowMeans(x-simVals))
            bias[v,mets,"bias",z,resp] <- rowMeans(impDat[v,mets,,"estimate",z,resp] - true, na.rm=naRm)
            #bias <- apply(impDat[,v,mets,,"estimate",z,resp])
            #bias[v,,"pctBias",z,resp] <- 100 * abs((avg[,"estimate"] - true) / true )
            #bias[v,,"pctBias",z,resp] <- 100 * abs(mean((impDat[v,met,,"estimate",z,resp] - true) / true , na.rm=naRm))
            #PB <-                          100 * abs((rowMeans(res[,, "estimate"]) - true)/ true)
            #bias[v,mets,"pctBias",z,resp] <- 100 * abs((rowMeans(impDat[v,mets,,"estimate",z,resp]) - true) / true )
            #bias[v,mets,"pctBias",z,resp] <- apply(impDat[v,mets,,"estimate",z,resp], 1, function(x) 100 *abs((mean(x)-true)/true))
            bias[v,mets,"pctBias",z,resp] <- apply(impDat[v,mets,,"estimate",z,resp], 1, function(x) 100 *abs(mean((x-true)/true)))
            #bias[v,,"pctBias",z,resp] <- 100 * abs(rowMeans((impDat[v,,,"estimate",] - true) / true ))
            #CR <-                          rowMeans(res[,, "2.5 %"] < true & true < res[,, "97.5 %"])
            bias[v,mets,"covRate",z,resp] <- rowMeans(impDat[v,mets,,"2.5 %",z,resp] < true & true < impDat[v,mets,,"97.5 %",z,resp],na.rm=naRm)
            bias[v,mets,"avgWidth",z,resp] <- rowMeans(impDat[v,mets,,"97.5 %",z,resp] - impDat[v,mets,,"2.5 %",z,resp], na.rm=naRm)
            #bias[v,,"RMSE",z,resp] <- sqrt((avg[,"estimate"] - true)^2)
            bias[v,mets,"RMSE",z,resp] <- 100 * sqrt(rowMeans((impDat[v,mets,,"estimate",z,resp] - true) ^2 ,na.rm=naRm))
            #biasByMet <- array(NA, dim=) # store individual methods' bias vals?
            #for(met in mets){
            #    cat(sprintf("\n    estimates for %s & %s & %s & %s:\n", v, met, resp, z))
            #    qvcalc::indentPrint(impDat[v,met,,"estimate",z,resp], indent=8)
            #    #bias[v,,"bias",z,resp] <- avg[,"estimate"] - true
            #    #bias[v,met,"bias",z,resp] <- rowMeans(impDat[v,mets,,"estimate",z,resp] - true, na.rm=naRm) # exclude sim and true by matching to 'mets'
            #    #bias[v,met,"bias",z,resp] <- rowMeans(impDat[v,met,,"estimate",z,resp] - true, na.rm=naRm)
            #    bias[v,met,"bias",z,resp] <- mean(impDat[v,met,,"estimate",z,resp] - true, na.rm=naRm)
            #    #bias[v,,"pctBias",z,resp] <- 100 * abs((avg[,"estimate"] - true) / true )
            #    bias[v,,"pctBias",z,resp] <- 100 * abs(mean((impDat[v,met,,"estimate",z,resp] - true) / true , na.rm=naRm))
            #    #bias[v,,"pctBias",z,resp] <- 100 * abs(rowMeans((impDat[v,,,"estimate",] - true) / true ))
            #    bias[v,,"covRate",z,resp] <- mean(impDat[v,met,,"2.5 %",z,resp] < true & true < impDat[v,met,,"97.5 %",z,resp],na.rm=naRm)
            #    bias[v,,"avgWidth",z,resp] <- mean(impDat[v,met,,"97.5 %",z,resp] - impDat[v,met,,"2.5 %",z,resp], na.rm=naRm)
            #    #bias[v,,"RMSE",z,resp] <- sqrt((avg[,"estimate"] - true)^2)
            #    bias[v,,"RMSE",z,resp] <- 100 * sqrt(mean((impDat[v,met,,"estimate",z,resp] - true) ^2 ,na.rm=naRm))
            #}

            #if(params$debug) cat(sprintf("\nbias values for %s and model %s %s:\n\n", v, resp, mods4sim[z]))
            if(params$vbose > 1)cat(sprintf("\n>>> bias values for model %s %s for variable %s\n",
                      resp, mods4sim[z], v))

            if(params$vbose > 1) qvcalc::indentPrint(bias[v,,,z,resp], indent=4)




        }

        biasfile <-  sprintf("%s/bias_vals_%s_m%s_%s.rds",now_dir,resp, names(mods4sim)[z], suffix)
        cat(sprintf("\n>>> saving to file: %s ", biasfile))
        saveRDS(bias, biasfile)

        biasfile1 <-  sprintf("%s/bias_vals_%s_%s_%s.csv",now_dir,resp, names(mods4sim)[z], suffix)
        cat(sprintf("\n\t & to: %s \n", biasfile1))
        biasdf <- as.data.frame(bias)
        names(trueVals)[ is.na(names(trueVals)) ] <- setdiff(vars, names(trueVals))
        biasdf <- cbind(trueVals, biasdf)
        write.csv(biasdf, file = biasfile1)# write to csv in case script aborts

        if(params$vbose > 1) cat(sprintf("\n\n>>>>>> BIAS VALUES FOR MODEL %s %s:\n", resp, mods4sim[z]))
        if(params$vbose > 1) qvcalc::indentPrint(bias[,,,z,resp], indent=4)

    }
}

#datt <- readRDS()

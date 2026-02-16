
library(brglm2)
library(tidyverse)
library(gt)

options(width=100)
source("imp_sim_functions.R")
params <- list(hdir = "/home/wodehouse/Projects/fate_glm/", 
#              #outdir = arg[1],
               outdir = NULL,
	       #dir_ext = paste0("out/",arg[1]),
               coefdir = "out/16Dec",
               nrun    = 100,
	       debug   = FALSE)

arg <- commandArgs(trailingOnly=TRUE)
#if(length(arg)==0) stop("needs name of output directory within out directory")
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
if(length(arg)==0){
    #outd <- readLines("last_out.txt")
    #cat("\n*** NO DIR NAME PROVIDED - defaulting to latest output directory:", outd)
    outd <- paste0(params$hdir,"out/12Dec")
    cat("\n*** NO DIR NAME PROVIDED - defaulting to:", outd)
} else {
    outd <- paste0(params$hdir,"out/", arg[1])
    #cat(sprintf("\tgetting model output from: out/%s\n", arg[1]))
    cat("\n\t>> getting model output from:",outd) 
}
cat("\n>> getting coef vals from dir:", params$coefdir)
cat("\n\n",format(Sys.time(), "%d-%b %H:%M"))
params$outdir = outd
if(length(arg)>1) params$debug=TRUE
cat("\n\n\t\tCALCULATING BIAS VALUES - calc_bias.R")
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
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
qvcalc::indentPrint(mods4sim)

#patt <- "100runs|250runs"
patt <- "vals"
#fnlist <- list.files(path = paste0(params$hdir, params$dir_ext), pattern = "^runs0to50.*|^runs50to100.*", full.names =T)
#fnlist <- list.files(path = paste0(params$hdir, params$dir_ext), pattern = patt, full.names =T)
fnlist <- list.files(path = params$outdir, pattern = patt, full.names =T)
#if(params$debug) cat("\n\t>> filename list:")
cat("\n>> filename list:\n")
qvcalc::indentPrint(fnlist, indent=8)
flist <- lapply(fnlist, FUN=function(x) readRDS(x))
nrun_list <- unlist(lapply(fnlist, function(x) str_extract(x, "\\d+(?=runs)")))
# search for spcific nrun valu instead of all:
#nrun_list <- unlist(lapply(fnlist, function(x) str_extract(x, "100(?=runs)")))
seeds <- unlist(lapply(fnlist, function(x) str_extract(x, "(?<=seed)\\d+")))
modnums <- unlist(lapply(fnlist, function(x) str_extract(x, "(?<=mod)\\d+")))
#mods <- names(mods4sim)[as.numeric(modnums)]
mods <- mods4sim[as.numeric(modnums)]
#cat("\n>> number of runs in each file:", length(nrun_list),"\n")
#print(nrun_list)
#cat("\n>> seeds from the files:", length(seeds),"\n")
#print(seeds)
#print(class(seeds))
#cat("\n>> models from the files:", length(mods),"\n")
#print(mods)
cat("\nall combos & unique combos:\n")
combb <- paste(nrun_list, seeds, modnums)
qvcalc::indentPrint(combb)
cat("\t\n")
qvcalc::indentPrint(unique(combb), indent=8)

#biasFiles <- list.files(path = params$outdir, pattern="sim_bias", full.names=T)
#cat("\nBIAS FILES:\n")
#print(biasFiles)
coefFiles <- list.files(path = params$coefdir, pattern="sim_coefs", full.names=T)
cat("\nCOEF FILES\n:")
qvcalc::indentPrint(coefFiles)

nrun_idx <- nrun_list==params$nrun
#print(nrun_idx)
nrun_list <- nrun_list[nrun_idx]
seeds <- seeds[nrun_idx]
mods  <- mods[nrun_idx]
flist <- flist[nrun_idx]
cat("\n>>>> num files:", length(flist))
cat("    & filenamees (correct nrun):\n")
qvcalc::indentPrint(fnlist[nrun_idx], indent=8)
fnlist <- fnlist[nrun_idx]

cat("\nnew combos & unique combos w/ only correct nrun numbers:\n")
combb <- paste(nrun_list, seeds, modnums)
qvcalc::indentPrint(combb)
cat("\t\n")
qvcalc::indentPrint(unique(combb), indent=8)
#combb <- paste(nrun_list, seeds, modnums)
#print(combb)
#cat("\n")
#print(unique(combb))
#print(str(flist))
#print(dim(flist))
#fl1 <- flist[[1]]
#str(fl1)

dat4sim <- read.csv("dat_complete.csv", stringsAsFactors = TRUE)
dat4sim$cam_fate <- relevel(dat4sim$cam_fate, ref="H") # make 'H' the reference category
dat4sim$species <- relevel(dat4sim$species, ref="LETE")
levels(dat4sim$HF_mis) <- c(0,1)
levels(dat4sim$is_u)   <- c(0,1)
fullDat <- dat4sim
cat("\n>>>> releveled categorical variables")

fam <- binomial
regMet <- brglm2::brglm_fit
iter <- 500

now_dir = paste0(params$hdir, "out/bias/", format(Sys.time(), "%d%b"))
if(!dir.exists(now_dir)) dir.create(now_dir)
cat("\n>>>> output directory:", now_dir, "exists?", dir.exists(now_dir))

    #resp_list <- unlist(lapply(fnlist, function(x) str_extract(x, paste(rlist,collapse="|"))))
    #print(resp_list)
    
    #for(resp in seq_along(resp_list)){
    #for(resp in unique(resp_list)){
    
    #for(resp in rlist)
    #for(seed in seeds){
    #    cat(sprintf("\n<> <> <> <> <> SEED: %s <> <> <> <> <> \n", seed))
    #    cat("\n>> matrices to bind:\n", seeds==seed)
    #    #cat("\n>> first, ")
    #    seedDat <- abind::abind(flist[seeds=seed], along=3)
    #    cat("\n>> matrices combined:\n")
    #    print(str(seedDat))
    #}
    ## Group the different seed output together by model 
for(seed in seeds){
    cat(sprintf("\n\n<> <> <> <> <> SEED: %s <> <> <> <> <> \n", seed))
    idx <- seeds==seed
    #fbind <- flist[seeds==seed]
    fbind <- flist[idx]
    cat("\n\t>>>> num files:", length(fbind))
    cat("    & filenamees (correct nrun):\n")
    qvcalc::indentPrint(fnlist[idx], indent=8)

    #nruns <- nrun_list[seeds==seed]
    nruns <- nrun_list[idx]
    modss <- mods[idx]
    cat("\n    mod & nrun for this seed:\n")
    qvcalc::indentPrint(modss, indent=8)
    qvcalc::indentPrint(nruns, indent=8)
    cat(sprintf("\n    >> matrices to bind - length %s:\n", length(fbind)))
    #qvcalc::indentPrint(fbind, indent=8)
    qvcalc::indentPrint(str(fbind))
    #cat("\n>> first, ")
    #seedDat <- abind::abind(flist[seeds=seed], along=3)
    seedDat <- abind::abind(fbind, along=3)
    cat("\n    >> matrices combined:\n")
    qvcalc::indentPrint(str(seedDat), indent=8)
#    fnlist <- list.files(path = params$outdir, pattern = patt, full.names =T)
    #seedBias <- str_extract(biasFiles, seed)
    seedCoef <- grep(seed, coefFiles, fixed=T, value=T)
    #cat("\nseed coef files:\n")
    qvcalc::indentPrint(seedCoef)
    if(length(seedCoef)<1) {
        cat("\nNo files for this seed, skip to next.\n")
        next
    }
    fname <- paste0(params$hdir,seedCoef)
    cat("\nfile to read for sim vals for this seed (should be 1):", fname)
    simVals <- readRDS(fname)
    cat("\n simulated vals read from file:\n")
    #qvcalc::indentPrint(simVals)
    qvcalc::indentPrint(str(simVals))

    for(z in seq_along(mods4sim)){
        cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
        cat("\n\t\t\t MODEL:", mods4sim[z])
        cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
        mod     <- mods4sim[z]
        cat("\n\n    >> all files for this model:\n", fnlist[mods==mod],length(flist[mods==mod]))
        cat("\n\n    >> seed files for this model (should be 1; others will be duplicates):\n", length(fbind[modss==mod]), "\n")
        #if(length(flist[mods==mod])<1) stop(sprintf("No files matching model %s. Exiting script.", mod))
        if(length(fbind[modss==mod])<1) cat("\nNo files for model, skip to next.\n")
        if(length(fbind[modss==mod])<1) next
        impDat <- fbind[modss==mod][[1]]
        qvcalc::indentPrint(str(impDat))
        #res <- array(NA, dim = c(length(vars), params$nrun,
        #             3, length(mods4sim), length(resp_list)))
        #nrun <- nruns[]
        #cat("nrun for this model:", nrun, "correct number:", nrun==params$nrun)
        #allSeeds <- seq(seed,seed+params$nrun)
        #cat("all seeds used for simulating data:", allSeeds)

        # no mod nums in coeefe files name
        #seedMods <- str_extract(seedCoef, "(?<=mod)\\d+")
        #readFile <- seedCoef[seedMods==mod]
        #cat("\n>>> attempt to merge:\n")
        #impDat <- abind::abind(flist[mods==mod], along=3)
        ## double bracket doesn't work'
        #impDat <- abind::abind(flist[[mods==mod]], along=3)
        ## Separate by response variable and get the real model output
        for(resp in rlist){
            cat(sprintf("\n<> <> <> <> <> MODEL: %s %s <> <> <> <> <> \n", resp, mod))
            sim_val <- simVals[,,,z,resp]
            cat("\n    >> sim values for this model:\n")
            qvcalc::indentPrint(sim_val, indent=8)
            fitReal <- glm(as.formula(paste0(resp, mod)),
                           data=fullDat,
                           family=fam,
                           method=regMet,
                           control=brglmControl(maxit=iter) )
            # don't want coef - want coefficients from summary.glm - need exp(coef())?
            trueVals <- exp(coef(fitReal)[vars]) # the coefs have names associated with them
            if(params$debug) qvcalc::indentPrint(summary(fitReal)) 
            cat("\n    true values:\n")
            qvcalc::indentPrint(trueVals)
            bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(mods4sim), length(rlist) ) )
            dimnames(bias) <- list(sort(as.character(vars)),
                                   # c("pmm", "rf", "cart"),
                                   as.character(mets),
                                   as.character(biasVals),
                                   names(mods4sim),
                                   rlist)
            if(params$debug) cat("\n    empty array to store bias values:" )
            if(params$debug) qvcalc::indentPrint(str(bias))
            if(params$debug) qvcalc::indentPrint(dimnames(bias))

            #for(r in seq(1, params$nrun)){
            #    sim <- exp(sim_val[,r,1]) # 1 = estimate
            #    cat(sprintf("\n    sim vals (estimate) for run %s (exp):\n", r))
            #    qvcalc::indentPrint(sim)

            ## Loop through the predictor variables and store the bias values to the matrix
                for(v in vars){
                    cat(sprintf("\n    ----- VARIABLE: %s ------------------------------------------------\n", v))
                    if(params$debug) cat(sprintf("\n    >> output for variable %s:\n", v))
                    # now I'm doing this in the other script (debug_whatever.R)
                    #impDat[v,"cc",,,z,resp] <- exp(impDat[v,"cc",,,z,resp]) 
                    if(params$debug) qvcalc::indentPrint(str(impDat[v,,r,,z,resp]), indent=8)
                    #avg <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
                    avg <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
                    sdev <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = sd, na.rm=TRUE)
                    #print(impDat[v,,1:5,1,z])
                    #true <- trueVals[v]
                    true <- exp(sim_val[v,,"Estimate"])
                    cat("\n    vals for comparison (exp)=", true,"\n    & dimensions:", dim(true), dim(impDat[v,,,"estimate",z,resp]))
                    if(params$debug){
                          cat(sprintf("\n    average of estimate for mod %s:\n",z))
                          qvcalc::indentPrint(apply(impDat[v,,,,z,resp], FUN=mean, MARGIN=c(1,3)))
                          cat("\n    dimensions:\n")
                          qvcalc::indentPrint(dim(impDat[v,,,,z,resp]))
                          cat("\n    avg:\n")
                          qvcalc::indentPrint(avg)
                    }
                    cat("\n    >>> storing vals to matrix\n")
                    bias[v,,"value",z,resp] <- avg[,"estimate"]
                    #bias[v,,"bias",z,resp] <- avg[,"estimate"] - true
                    bias[v,,"bias",z,resp] <- rowMeans(impDat[v,,,"estimate",z,resp] - true)
                    #bias[v,,"bias",z,resp] <- rowMeans(impDat[v,,,"estimate",z,resp] - true, dims=2)
                    #bias[v,,"pctBias",z,resp] <- 100 * abs((avg[,"estimate"] - true) / true )
                    bias[v,,"pctBias",z,resp] <- 100 * abs(rowMeans((impDat[v,,,"estimate",z,resp] - true) / true ))
                    bias[v,,"covRate",z,resp] <- rowMeans(impDat[v,,,"2.5 %",z,resp] < true & true < impDat[v,,,"97.5 %",z,resp])
                    bias[v,,"avgWidth",z,resp] <- rowMeans(impDat[v,,,"97.5 %",z,resp] - impDat[v,,,"2.5 %",z,resp])
                    #bias[v,,"RMSE",z,resp] <- sqrt((avg[,"estimate"] - true)^2)
                    #bias[v,,"RMSE",z,resp] <- 100 * sqrt((rowMeans(impDat[v,,,"estimate",]) - true) ^2 )
                    bias[v,,"RMSE",z,resp] <- 100 * sqrt(rowMeans(impDat[v,,,"estimate",z,resp] - true) ^2 )
                    bias[v,,"SD",z,resp] <- sdev[,"estimate"]
                    #if(parrams$debug) cat(sprintf("\nbias values for %s and model %s %s:\n\n", v, resp, mods4sim[z]))
                    cat(sprintf("\n>>> bias values for model %s %s for variable %s\n",
                              resp, mods4sim[z], v))
                    qvcalc::indentPrint(bias[v,,,z,resp], indent=8)
                }

            }
            biasfile <-  sprintf("%s/bias_vals_%s_m%s_%s.rds",now_dir,resp, names(mods4sim)[z], suffix)
            cat(sprintf("\n>>> saving to file: %s ", biasfile))
            saveRDS(bias, biasfile)
            biasfile1 <-  sprintf("%s/bias_vals_%s_%s_%s.csv",now_dir,resp, names(mods4sim)[z], suffix)
            cat(sprintf("\t & to: %s \n", biasfile1))
            biasdf <- as.data.frame(bias)
            names(trueVals)[ is.na(names(trueVals)) ] <- setdiff(vars, names(trueVals))
            biasdf <- cbind(trueVals, biasdf)
            write.csv(biasdf, file = biasfile1)# write to csv in case script aborts

        #}
    }
}


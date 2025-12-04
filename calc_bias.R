
library(brglm2)
library(tidyverse)
library(gt)

options(width=100)
source("imp_sim_functions.R")
arg <- commandArgs(trailingOnly=TRUE)
if(length(arg)==0) stop("needs name of output directory within out directory")
#cat("\n===========================================================================")
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat(format(Sys.time(), "%d-%b %H:%M"))
cat(sprintf("\tgetting model output from: out/%s\n", arg[1]))
params <- list(hdir = "/home/wodehouse/Projects/fate_glm/", 
#               outdir = arg[1],
	       dir_ext = paste0("out/",arg[1]),
	       debug   = FALSE)

if(length(arg)>1) params$debug=TRUE

cat("\n===========================================================================")
cat("\nCALCULATING BIAS VALUES - calc_bias.R")
cat("\n===========================================================================")
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
print(mods4sim)

patt <- "100runs|250runs"
#fnlist <- list.files(path = paste0(params$hdir, params$dir_ext), pattern = "^runs0to50.*|^runs50to100.*", full.names =T)
fnlist <- list.files(path = paste0(params$hdir, params$dir_ext), pattern = patt, full.names =T)
#if(params$debug) cat("\n\t>> filename list:")
cat("\n>> filename list:\n")
print(fnlist)
flist <- lapply(fnlist, FUN=function(x) readRDS(x))
seeds <- unlist(lapply(fnlist, function(x) str_extract(x, "(?<=seed)\\d+")))
modnums <- unlist(lapply(fnlist, function(x) str_extract(x, "(?<=mod)\\d+")))
#mods <- names(mods4sim)[as.numeric(modnums)]
mods <- mods4sim[as.numeric(modnums)]
cat("\n>> seeds from the files:")
print(seeds)
print(class(seeds))
cat("\n>> models from the files:")
print(mods)
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
    #cat(sprintf("\n<> <> <> <> <> SEED: %s <> <> <> <> <> \n", seed))
    #cat("\n>> matrices to bind:\n", seeds==seed)
    #impDat <- abind::abind(flist[seeds=seed], along=3)
    #cat("\n>> matrices combined:\n")
    #print(str(impDat))

## Group the different seed output together by model 
for(z in seq_along(mods4sim)){
    cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    cat("\n\t\t\t MODEL:", mods4sim[z])
    cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    mod     <- mods4sim[z]
    cat("\n\n>> files for this model:\n", mods==mod)
    cat("\n\n>> files for this model:\n", length(flist[mods==mod]), "\n")
    if(length(flist[mods==mod])<1) stop(sprintf("No files matching model %s. Exiting script.", mod))
    impDat <- abind::abind(flist[mods==mod], along=3)
    ## double bracket doesn't work'
    #impDat <- abind::abind(flist[[mods==mod]], along=3)
    if(params$debug) cat("\n\n>> matrices combined:\n")
    if(params$debug) print(str(impDat))
    if(params$debug) print(dimnames(impDat))

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
        trueVals <- exp(coef(fitReal)[vars]) # the coefs have names associated with them
        if(params$debug) print(summary(fitReal)) 
        cat("\ntrue values:\n")
        print(trueVals)

        bias <- array(NA, dim = c(length(vars), length(mets), length(biasVals), length(mods4sim), length(rlist) ) )
        dimnames(bias) <- list(sort(as.character(vars)),
                               # c("pmm", "rf", "cart"),
                               as.character(mets),
                               as.character(biasVals),
                               names(mods4sim),
                               rlist)
        if(params$debug) cat("\nempty array to store bias values:" )
        if(params$debug) print(str(bias))
        if(params$debug) print(dimnames(bias))

        ## Loop through the predictor variables and store the bias values to the matrix
        for(v in vars){
            cat(sprintf("\n----- VARIABLE: %s ------------------------------------------------\n", v))
            if(params$debug) cat(sprintf("\n>> output for variable %s:\n", v))
            #if(params$debug) print(impDat[v,,,,z,resp])
            #if(params$debug) print(str(impDat))
            #if
            impDat[v,"cc",,,z,resp] <- exp(impDat[v,"cc",,,z,resp]) 
            if(params$debug) print(str(impDat[v,,,,z,resp]))
            avg <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = mean, na.rm=TRUE)
            sdev <- apply(impDat[v, , , ,z,resp],MARGIN=c(1,3),FUN = sd, na.rm=TRUE)
            #print(impDat[v,,1:5,1,z])
            true <- trueVals[v]
            #cat(sprintf("\nbias for %s, value, and %s:\n", v, z))
            #print(bias[v,,"value",z])
                  #cat("\nfirst 5 reps:\n")
                  #print(impDat[v,,1:5,,z])
            if(params$debug){
                  cat(sprintf("\naverage of estimate for mod %s:\n",z))
                  #print(apply(impDat[v,,,,z], FUN=mean, MARGIN=c(1,3)))
                  #print(impDat[v,,,z])
                  print(apply(impDat[v,,,,z,resp], FUN=mean, MARGIN=c(1,3)))

                  cat("\ndimensions:\n")
                  print(dim(impDat[v,,,,z,resp]))
                  cat("\navg:\n")
                  print(avg)
            }
            #print(str(avg))
            cat("\n>>> storing vals to matrix\n")
            bias[v, ,"value",z,resp] <- avg[,"estimate"]
            bias[v,,"bias",z,resp] <- avg[,"estimate"] - true
            bias[v,,"pctBias",z,resp] <- 100 * abs((avg[,"estimate"] - true) / true )
            #bias[v,,"pctBias",z,resp] <- 100 * abs(rowMeans((impDat[v,,,"estimate",] - true) / true ))
            bias[v,,"covRate",z,resp] <- rowMeans(impDat[v,,,"2.5 %",z,resp] < true & true < impDat[v,,,"97.5 %",z,resp])
            bias[v,,"avgWidth",z,resp] <- rowMeans(impDat[v,,,"97.5 %",z,resp] - impDat[v,,,"2.5 %",z,resp])
            bias[v,,"RMSE",z,resp] <- sqrt((avg[,"estimate"] - true)^2)
            #bias[v,,"RMSE",z,resp] <- 100 * sqrt((rowMeans(impDat[v,,,"estimate",]) - true) ^2 )
            bias[v,,"SD",z,resp] <- sdev[,"estimate"]

            #if(params$debug) cat(sprintf("\nbias values for %s and model %s %s:\n\n", v, resp, mods4sim[z]))
            cat(sprintf("\n>>> bias values for model %s %s for variable %s\n",
                      resp, mods4sim[z], v))

            print(bias[v,,,z,resp])




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

	}
}


library(brglm2) # penalized logistic regression
library(CALIBERrfimpute) # could access using ::
library(gt)
library(gtsummary)
library(data.table)

suppressMessages(library(mice))
suppressMessages(library(tidyverse)) # look into tidytable for tidyverse syntax w/ data.table
options(width=99)
########## ARGS & PARAMS ####################################################################
params <- list(nrun=100,
	       hdir="/home/wodehouse/Projects/fate_glm/",
	       #outdir = paste0(format(Sys.time(), "%d%b"),"/"),
               suffix = "",
	       outdir = format(Sys.time(), "%d%b"),
               ampWt=3,
               #deb=FALSE,
               #xdeb=FALSE, # for obsessively checking things - not useful for normal debugging & lots of text output
               ipl=FALSE,
               vbose=0,
               win=FALSE,
               test=FALSE,
               seeds = c(666),
               # resp = NULL,
               j = 50,
               m=20)
arg <- commandArgs(trailingOnly=TRUE)
if(length(arg) > 0){ # ex: Rscript isim.R win test r3 s613 m5 j1
    cat("\n\n///////////////////// arg =  ", paste(unlist(arg), collapse=" ; "),"  ////////////////////////////////\n")
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
} else {
    cat("\n\n//////////////// *** NOTE *** no arguments provided - using default of nrun=3, imputations=5, j=1, debug=F ////////////////\n\n")
    # mostly for source-ing the file in R console
    #params$test <- TRUE
    params$nrun <- 3
    params$j    <- 1
    params$m    <- 5
    #params$deb  <- FALSE
    params$vbose <-0
    params$suffix <- "rr"

}
########## SOURCE FILES ####################################################################
s_files <- c("visualize_data.R","mice_functions.R", "missing_data.R", "other_imp.R", "imp_sim_functions.R")
#if(params$deb|params$xdeb) s_files[4] = "debug_imp_sim_func.R" 
if(params$vbose>=1) s_files[4] = "debug_imp_sim_func.R" 
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
cat("\n>>> sourcing files:", paste(s_files,collapse=" ; "))
invisible(lapply(s_files, function(x) source(x)))
######### other vals: ####################################################################
bprint <- function(x) print(rbind(head(x, 3), tail(x,3)))
debugging <- FALSE # for uickly setting values when working in the file with the functions
suffix <- paste(sprintf("%sruns", params$nrun), params$suffix, sep="_")
#now_dir <- paste(params$hdir, params$outdir, sep="out/")
checkdir <- paste0(sprintf("check_sim_%s/", params$outdir))
now_dir <- paste(params$hdir,checkdir, sep="out/")
if(!dir.exists(now_dir)) dir.create(now_dir)
seed_out <- FALSE
########## check params: ##########################################################################
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
}
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
########## import the model specification: ####################################################################
if(TRUE){
    modList <- readLines("modList.txt")
    mods4sim <- modList[c(1,8,16) ]
    names(mods4sim) <- c("m1", "m8", "m16")
    fits <- readRDS("fits.rds")
    mnums <- rep(c("1", "8", "16"))
    cat("model nums from file:", mnums)
}
########## import the full complete-case data: ####################################################################
if(TRUE){
    fullDat <- read.csv("dat_complete.csv", stringsAsFactors = TRUE)
    fullDat$cam_fate <- relevel(fullDat$cam_fate, ref="H") # make 'H' the reference category
    fullDat$species <- relevel(fullDat$species, ref="LETE")
    levels(fullDat$HF_mis) <- c(0,1)
    levels(fullDat$is_u)   <- c(0,1)
    cat("\n>>>> releveled categorical variables")
}
########## for creating the simulated data: #######################################################
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
    cat("\n\n<><><><><><><><><><> Betas, Means & Correlation Matrices: <><><><><><><><><><><><><><><><><><><><>\n")
    qvcalc::indentPrint(head(betas))
    qvcalc::indentPrint(mList)
    qvcalc::indentPrint(cMat)
    #cat("\n\n<><><><><><><><><><> Missingness pattern & variable weights: <><><><><><><><><><><><><><><><><><><><>\n")
    #print(mpatt)
    #print(ampwt)
}
#########################################################################################
#cat("sort vars or not???")
vars <-  c("nest_age", "cam_fateD", "cam_fateA", "cam_fateF", "cam_fateHu", "cam_fateS", "speciesCONI", "speciesCONI:nest_age", "speciesCONI:obs_int", "obs_int", "fdate") # all vars
prVars <- c("species", "cam_fate", "obs_int", "nest_age", "fdate")
resp_list <- c("is_u", "HF_mis")
#########################################################################################
cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
cat(sprintf("\t>> TOTAL REPS: %s x %s WITH %s NESTS EACH", length(params$seeds), params$nrun, nnest)) 
#nif(params$deb) cat("\t- debug ON")
#if(params$xdeb) cat(" + EXTRA")
if(params$vbose >= 1) cat("\t *** verbosity = ", params$vbose)
cat("\n>> output will be saved every", params$j, "runs to dir:", now_dir,"\n\n")

for (seed in params$seeds){
    outf <- paste0(now_dir,sprintf("/%s_loggedEvents.out",seed ))
    res <- array(NA, dim = c(length(vars), params$nrun, 3, length(mods4sim), length(resp_list)))
    #res <- array(NA, dim = c(length(vars), params$nrun, length(mods4sim), length(resp_list)))
    #biasN <- c("bias", "pctBias", "cov", "width")
#  vars need to be in same order as vals
    dimnames(res) <- list(as.character(vars), as.character(1:params$nrun),c("Estimate", "2.5 %", "97.5 %"), names(mods4sim), resp_list)
    resAvg <- array(NA, dim=c(length(vars), 4, length(mods4sim), length(resp_list)))
    biasNames <-  c("avgBias", "avgPctBias", "covRate", "avgWidth")
    dimnames(resAvg) <- list(vars, biasNames, names(mods4sim), resp_list)
    #if (params$deb) print(dimnames(resAvg))
    if (params$vbose>=2) print(dimnames(resAvg))
    fTabs <- array(NA, dim=c(8, params$nrun, 2, 2)) ## dim - levels of categorical vars, nrun, 2 resp vars, T/F
    catVars <- c( "speciesLETE","speciesCONI","cam_fateH", "cam_fateA", "cam_fateD", "cam_fateF", "cam_fateHu", "cam_fateS")
    dimnames(fTabs) <- list(catVars, seq(1, params$nrun), resp_list, c("YES", "NO"))
    camFateVars <- vars[grepl(pattern="cam_fate", x=vars, fixed=TRUE)]
    varInfo <- array(NA, dim=c(length(camFateVars), params$nrun, length(mods4sim))) # keep a count of sample size for each category
    dimnames(varInfo) <- list( camFateVars, seq(1,params$nrun), mods4sim)
    cat(sprintf(">>>>>>>>>>>>>> running simulation %s times \t>>> seed = %s ", params$nrun, seed))
    cat("\t>> & no. imp.:", params$m, "\n\n")
    for(mod in seq_along(mods4sim)){
        cat("\n\n::::::::::::::::::: MODEL", mods4sim[mod], ":::::::::::::::::::::::::::::::\n\n") ## *~*~*~*~*
        beta_list <- betas[[mod]]
        #form_list <- formulas[[names(mods4sim)[mod]]]
        modnum <- gsub("\\D+", "", names(mods4sim)[mod])
        cat("\n<> <> model number:", modnum)
        fit_isU <-  glm(as.formula(paste("is_u", mods4sim[mod])), family=binomial, data=fullDat, method=brglm2::brglmFit)
        fit_mis <-  glm(as.formula(paste("HF_mis", mods4sim[mod])), family=binomial, data=fullDat, method=brglm2::brglmFit)
        cat("\n>> beta list, model fits:") ## *~*~*~*~*
        if(TRUE){
            #qvcalc::indentPrint(names(form_list))
            qvcalc::indentPrint(str(beta_list))
            qvcalc::indentPrint(summary(fit_isU))
            qvcalc::indentPrint(summary(fit_mis))
        }
        for(run in 1:params$nrun){
            cat(run)
            #dat4sim <- mkResp(seed=seed+run,resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, debug=FALSE)
            dat4sim <- mkResp(seed=run+seed,resp_list, mods4sim[mod], nnest, cMat, mList, beta_list, fprob, sprob, prList, vbose=params$vbose)
            if(params$vbose>=2) cat(mod)
            if(params$vbose>=2) cat("\n>> make freeequency tables:")
            #ftabs <- mk_fr_tab(dat4sim, vars="HF_mis", save=FALSE, debug=params$deb)
            ftabs <- mk_fr_tab(dat4sim, vars="HF_mis", save=FALSE, debug=FALSE)
            #qvcalc::indentPrint(sapply(ftabs, function(x) dim(x)))
            if(params$vbose>=2) qvcalc::indentPrint(ftabs)
                #catVars <- c("cam_fateH", "cam_fateA", "cam_fateD", "cam_fateF", "cam_fateS", "cam_fateHu", "speciesLETE","speciesCONI")
                #dimnames(fTabs) <- list(catVars, seq(1, params$nrun), resp_list, c("TRUE", "FALSE"))
            fTabs[catVars, run, "is_u",] <- ftabs[[1]]
            fTabs[catVars, run, "HF_mis",] <- ftabs[[2]]
            #ftabs <- sapply(ftabs, function(x))
            # if you do summary, you get the standard error and all that as well
            if(params$vbose>=2) cat("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
            if(params$vbose>=2) cat("\n> fitting simulated data:")
            #fitSim1 <- summary(glm(as.formula(paste0("is_u",mods4sim[mod])),
            fitSim1 <- glm(as.formula(paste0("is_u",mods4sim[mod])),
                           family=binomial,
                           data=dat4sim,
                           method=brglm2::brglmFit)
            #fitSim2 <- summary(glm(as.formula(paste0("HF_mis",mods4sim[mod])),
            fitSim2 <- glm(as.formula(paste0("HF_mis",mods4sim[mod])),
                           family=binomial,
                           data=dat4sim,
                           method=brglm2::brglmFit)


            #cat("\n > > classes: ", class(fit_real[[1]]), class(fitSim1))
            if(params$vbose>=2) cat("\n+ + + coefs from real data vs from sim data:\n")
            if(params$vbose>=2) qvcalc::indentPrint(data.frame("real"=coef(fit_isU), "sim"=coef(fitSim1), "difference"=coef(fit_isU)-coef(fitSim1)))
            if(params$vbose>=2) qvcalc::indentPrint(data.frame("real"=coef(fit_mis), "sim"=coef(fitSim2), "difference"=coef(fit_mis)-coef(fitSim2)))
            #print(coef(fit_mis)[-1])
            #print(coef(fitSim2)[-1])
            if(params$vbose>=2) cat("\n + + + calculate & fill in the estimates: + + + \n")
            #print(coef(fit_isU)[-1] - coef(fitSim1)[-1])
            vlist1 <- colnames(model.matrix(as.formula(paste0("is_u", mods4sim[mod])),data = fullDat))[-1]
            vlist2 <- colnames(model.matrix(as.formula(paste0("HF_mis", mods4sim[mod])),data = fullDat))[-1]
            if(params$vbose>=2) cat("\n>>>> var lists:", vlist1, vlist2)
            if(params$vbose>=2) cat("\n\t >> percent bias:", 100*( ( coef(fit_isU)[-1] - coef(fitSim1)[-1] ) / coef(fit_isU)[-1] ))
            if(params$vbose>=2) cat("\n\t >> coefs:")
            if(params$vbose>=2) qvcalc::indentPrint(cbind(coef(fitSim2), confint(fitSim2)))

            res[vlist1,run,,mod,"is_u"]      <-cbind( coef(fitSim1)[-1], confint(fitSim1)[-1,])
            res[vlist1,run,,mod,"HF_mis"]      <-cbind( coef(fitSim2)[-1], confint(fitSim2)[-1,])
            if(params$vbose>=2) cat("\n res so far:\n")
            if(params$vbose>=2) qvcalc::indentPrint(res[,run,,mod,])
            #res[vlist2,run,,mod,"HF_mis"]    <- coef(fitSim2)[-1]
            #res[vlist1,run,mod,"is_u"]      <- coef(fit_isU)[-1] - coef(fitSim1)[-1]
            #res[vlist2,run,mod,"HF_mis"]    <- coef(fit_mis)[-1] - coef(fitSim2)[-1]
            #res[vlist1,run,mod,"is_u"]   <- 100*( ( coef(fit_isU)[-1] - coef(fitSim1)[-1] ) / coef(fit_isU)[-1] )
            #res[vlist2,run,mod,"HF_mis"] <- 100*( ( coef(fit_mis)[-1] - coef(fitSim2)[-1] ) / coef(fit_mis)[-1] )
            #if(params$xdeb) cat("\n\t >> upper CI",confint(fitSim1)[2][-1] )
            #if(params$xdeb) cat("\n\t >> lower CI",confint(fitSim1)[1][-1] )
            #uci1 <- rep(confint(fitSim1)[2], length(vlist1)+1)
            ##lci1 <- rep(confint(fitSim1)[1], length(vlist1)+1)
            #if(params$deb) cat("\n\t >> upper CI, vector:", uci1,length(uci1), "\n\t & coefs:", coef(fit_isU), length(coef(fit_isU)))
            #res[vlist1,run,"cov",mod,"is_u"] <- ( rep(confint(fitSim1)[2], length(vlist1)+1) > coef(fit_isU) ) && ( coef(fit_isU) > rep(confint(fitSim1)[1], length(vlist1)+1) )
            #res[vlist1,run,"cov",mod,"HF_mis"] <- (rep(confint(fitSim2)[2], length(vlist2)+1) > coef(fit_mis)) && (coef(fit_mis) > rep(confint(fitSim2)[1], length(vlist2)+1))
            #for (v in vars){
            #    res[v,run,"cov",mod,"is_u"] <- ( uci1 > coef(fit_isU) ) && ( coef(fit_isU) > lci1 )
            #    res[v,run,"cov",mod,"HF_mis"] <- (rep(confint(fitSim2)[2], length(vlist2)+1) > coef(fit_mis)) && (coef(fit_mis) > rep(confint(fitSim2)[1], length(vlist2)+1))
            #    #res[vlist2,run,"cov",mod,"HF_mis"] <- confint(fit_mis)[2] > coef(fit_mis) && coef(fit_mis)  >  confint(fitSim2)[1]
            #    res[v,run,"width",mod,"is_u"] <- confint(fit_isU)[2] - confint(fitSim1)[1]
            #    res[v,run,"width",mod,"HF_mis"] <- confint(fit_isU)[2] - confint(fitSim1)[1]
            #}
        }
        if(params$vbose>=2) cat("\n>>> calculating bias & storing vals to matrix\n")
    #dimnames(res) <- list(as.character(vars), as.character(1:params$nrun),3, names(mods4sim), resp_list)
    #resAvg <- array(NA, dim=c(length(vars), 4, length(mods4sim), length(resp_list)))
    #biasNames <-  c("avgbias", "avgPctBias", "covRate", "avgWidth")
        for(resp in resp_list){
            if(resp=="is_u") fit <- fit_isU[[1]] else fit <- fit_mis[[1]]
            #fit <- ifelse(resp=="is_u", fit_isU[[1]], fit_mis[[1]])
            if(params$vbose>=2) cat(sprintf("fit for %s:\n ", resp))
            if(params$vbose>=2) qvcalc::indentPrint(fit)
            for(v in vars){
                if(params$vbose>=2) cat("\n<> <> <> <> <> <> <> <> <> <> <> <> VARIABLE = ", v, "<> <> <> <> <> <> <> <> <> <> <> <> \n")
                #true <- coef(fit)[v]
                true <- (fit)[v]
                #resAvg[v, "value",mod,resp] <- avg[,"estimate"]
                #resAvg[v,"bias",mod,resp] <- rowMeans(res[v,,mod,resp ]) - true
                if(params$vbose>=2) cat("\n>var mod resp")
                if(params$vbose>=2) cat(v, mod, resp)
                if(params$vbose>=2) cat("\n>res\n")
                if(params$vbose>=2) qvcalc::indentPrint(res[v,,"Estimate",mod,resp])
                if(params$vbose>=2) cat("\n>res - true\n")
                if(params$vbose>=2) qvcalc::indentPrint(mean(res[v,,"Estimate",mod,resp] - true))
                if(params$vbose>=2) cat("\n>resAvg\n")
                if(params$vbose>=2) qvcalc::indentPrint(resAvg[,, mod,resp])
                if(params$vbose>=2) cat("\n>store to matrix\n")
                resAvg[v,"avgBias",mod,resp] <- mean(res[v,,"Estimate",mod,resp ] - true)
                resAvg[v,"avgPctBias",mod,resp] <- 100 * abs((mean(res[v,,"Estimate",mod,resp] -true)/true))
                #bias[v,,"pctBias",z,resp] <- 100 * abs(rowMeans((impDat[v,,,"estimate",] - true) / true ))
                resAvg[v,"covRate",mod,resp] <- mean(res[v,,"2.5 %",mod,resp]< true & true < res[v,,"97.5 %",mod,resp] )
                resAvg[v,"avgWidth",mod,resp] <- mean(res[v,,"97.5 %",mod,resp]-res[v,,"2.5 %",mod,resp])
                if(params$vbose>=2) qvcalc::indentPrint(resAvg[v,,mod,resp])
            }
        }
        cat("\n bias vals (simulated vs true data):")
        qvcalc::indentPrint(resAvg)
        fname1 <- paste0(now_dir, sprintf("/sim_bias_%s_mod%s.rds", seed, mod))
        cat("\n >> saving bias vals to file:", fname1)
        saveRDS(resAvg, fname1)

        #for(b in seq_along(biasNames)){
        #    if(params$deb) cat("\navg bias vals, is_u:") 
        #    if(params$deb) qvcalc::indentPrint(rowMeans(res[,,biasN[b],mod,"is_u"]))
        #    resAvgMis[,biasNames[b],mod,"is_u"] <- rowMeans(res[,,biasN[b],mod,"is_u"])
        #    if(params$deb) cat("\navg bias vals, HF_mis:")h
        #    if(params$deb) qvcalc::indentPrint(rowMeans(res[,,biasN[b],mod,"HF_mis"]))
        #    resAvgMis[,biasNames[b],mod,"HF_mis"] <- rowMeans(res[,,biasN[b],mod,"HF_mis"])
        #}
    }
    if(params$vbose>=2) cat("\n>> res, all filled in:")
    if(params$vbose>=2) qvcalc::indentPrint(res, indent=8)
    fnamee <- paste0(now_dir, sprintf("/sim_coefs_%s.rds", seed))
    cat("\n\n >> saving res to file:", fnamee)
    saveRDS(res, fnamee)

    if(params$vbose>=2) cat("\n>> frequencies, all filled in:")
    if(params$vbose>=2) qvcalc::indentPrint(fTabs, indent=8)
    fname <- paste0(now_dir, sprintf("/fr_tabs_%s.rds", seed))
    cat("\n\n >> saving fTabs to file:", fname)
    saveRDS(fTabs, fname)
}

#saveRDS(res, "sim_bias_.rds")



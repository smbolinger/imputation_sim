
######### PARAMS ############################################################################
params <- list(hdir = "/home/wodehouse/Projects/imputation_sim/", 
               dir_ext = NULL,
               outdir = format(Sys.time(), "%d%b"),
               suffix = "",
               test=FALSE,
               seeds = NULL,
               one_mod = TRUE,
               mod_num = 6,
               # resp = NULL,
               strat_sim=TRUE,
               nrun    = 100,
               nNest=250,
               j = 50,
               m=20,
               ampWt=1,
               naRm <- TRUE,
               unexp <- FALSE,
               ipl=FALSE,
               win=FALSE,
               vbose   = 0
)
######## VERBOSITY LEVELS ##################################################################### cat("sort vars or not???")
vb <- list(
             calc_bias = 3,
             mk_sim    = 1,
             mk_amp    = 1, 
             mk_imp    = 1
)
######## VAR/METHOD LISTS ##################################################################### cat("sort vars or not???")
# mets <- c("cart","passive", "cc")# don't need full here?
mets      <- c("default", "cart", "caliber","passive", "stratify","cf_cc","cc")# don't need full here?
resp_list <- c("is_u", "HF_mis") # resp_list <- c("is_u")
biasVals  <- c("value","bias", "pctBias", "covRate", "avgWidth", "RMSE", "SD")

######## NEST DATA #####################################################################
if(TRUE){
  trueDat <- readRDS("simdata/dat_complete.rds")
  if(params$vbose>=2) cat("\n\t>=> relevel cam_fate")
  trueDat$cam_fate <- relevel(as.factor(trueDat$cam_fate), ref="H")
  if(params$vbose>=2) cat("\t>=> relevel species")
  trueDat$species <- relevel(as.factor(trueDat$species), ref="LETE")
  if(params$vbose>=2) cat("\t>=> make categorical vars factors\n")
  # trueDat$species <- as.factor(trueDat$species)
  trueDat$is_u <- as.factor(trueDat$is_u)
  trueDat$HF_mis <- as.factor(trueDat$HF_mis)
  if(params$vbose>=2) qvcalc::indentPrint(summary(trueDat), indent=4)
  if(params$vbose>=2) qvcalc::indentPrint(str(trueDat), indent=4)
}
######## MODELS/METHODS #####################################################################
formulas <- readRDS("simdata/form_lists.rds")
metLists <- readRDS("simdata/met_lists.rds")
modList  <- readLines("simdata/modList.txt")
if(one_mod){
    mods4sim <- c(modList[mod_num])
    oneMod    <- modList[params$mod_num]
    mod       <- params$mod_num
    form <- as.formula(mod)
    vars <- colnames(model.matrix(form, data=trueDat))
    prVars <- strsplit(mod, split="\\W")
    cat(sprintf("\n\t>-> using only 1 model: %s \n\t\t> vars=%s & prVars=%s", oneMod, vars, prVars))
    beta_list <- lapply(resp_list, function(x) beta_lists[[x]][[params$mod_num]]) # just for mod 1 
    form_list <- formulas[[params$mod_num]]
} else {
    mods4sim <- modList[c(1,8,16) ]
    names(mods4sim) <- c("m1", "m8", "m16")
    # this part could be applied to either:
    vars <- unique(unlist(sapply(mods4sim, function(x) colnames(model.matrix(as.formula(x), data=trueDat)))))[-1]
    prVars <- unique(unlist(sapply(mods4sim, function(x) strsplit(x, split="\\W"))))[-1]
    cat("\n\t>> models:\n")
    qvcalc::indentPrint(mods4sim, indent=4)
    cat(sprintf("\n\t\t> vars=%s & prVars=%s", vars, prVars))
}
######## SIMULATION INPUT #####################################################################
if(TRUE){
    # beta_lists <- readRDS("simdata/betas2.rds") # 2 has CONI as ref level instead of LETE, so it messes up the varnames
    beta_lists <- readRDS("simdata/betas.rds") # this is a list of vectors or something
    mList <- readRDS('simdata/means.rds')
    cMat  <- readRDS('simdata/cormat.rds')
    if(params$vbose>=2) cat("\nbetas, means, & matrix:\n")
    # if(params$vbose>=2) qvcalc::indentPrint(betas, indent=4)
    if(params$vbose>=2) qvcalc::indentPrint(beta_lists, indent=4)
    if(params$vbose>=2) qvcalc::indentPrint(mList, indent=4)
    if(params$vbose>=2) qvcalc::indentPrint(cMat, indent=4)
    fprob <- c( 'H'=0.53,'A'=0.1, 'D'=0.13,'F'=0.06, 'Hu'=0.13,'S'=0.06 )
    sprob <- c('CONI'=0.32, 'LETE'=0.68)
    #nnest <- ifelse(params$test==T, 200, 200) # number of simulated nests
    nnest <- params$nNest
    mpatt  <- readRDS("simdata/misPatt.rds")
    awFile <- case_when(params$ampWt==1 ~ "simdata/ampWts.rds",
                        params$ampWt==2 ~ "simdata/ampWts2.rds",
                        params$ampWt==3 ~ "simdata/ampWts3.rds")
    ampwt <- readRDS(awFile)
    cat("\n\n<><><><><><><><><><> Missingness pattern & variable weights: <><><><><><><><><><><><><><><><><><><><>\n\n")
    print(mpatt)
    print(ampwt)
}
######## OTHER #####################################################################
nowtime <- format(Sys.time(), "%H_%M")
fam     <- binomial
regMet  <- brglm2::brglm_fit
iter    <- 500
seed    <- params$seeds[1]



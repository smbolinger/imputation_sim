

# top_num <- aicTabImp$num[1]
# regtab <- function(aicTabImp, modList, howManyMod){
regtab <- function(impDat, aicTabImp, modList, debug=FALSE){
  # browser()
  
    modnum  <- as.numeric(rownames(aicTabImp)[m])
    # NOTE - will change based on whether there is a space in the modname
    intVar1 <- str_extract(string = modList[modnum], pattern = "\\w+(?=\\s\\*)")
    intVar2 <- str_extract(string = modList[modnum], pattern = "(?<=\\*\\s)\\w+")
    cat("\nmodel number", modnum,"=", modList[modnum], ", with these interaction vars:", intVar1, intVar2, ", and delta AIC:", aicTabImp$deltaAIC[m], "\n")
    fit <- with(impDat, glm(as.formula(paste(resp, modList[[modnum]])), family=binomial, method=brglm2::brglmFit)) 
    # pass results to tbl_regression before pooling, and it will handle pooling + tidying:
    vars=c('species', 'nest_age', 'obs_int', 'cam_fate', 'fdate')
    varNames=c('SPECIES', 'NEST_AGE', 'FINAL_INTERVAL', 'TRUE_FATE', 'END_DATE')
    vnames = as.list(varNames)    
    names(vnames) <- vars
    # tabOR <- gtsummary::tbl_regression(pool.fit, label=vnames, exponentiate=TRUE)
    summ.fit <- gtsummary::tbl_regression(fit, label=vnames, exponentiate=TRUE)
    # summ.fit
    # fn <- sprintf("_q1_reg_imp_%s_%s.rtf", datNameImp, now)
    headr <- paste("Pooled MI Regression Summary for Model:",
                   modnames(list(x = mods[[modnum]]), null=FALSE),
                   sep="\n")
    fn <- sprintf("%s_regtab_imp_%s-%s_%s.rtf", quest, datNameImp,modnum, now)
    summ.fit %>%
      gtsummary::as_gt() %>%
      gt::tab_header(title=headr) %>%
      gt::gtsave( filename=fn, path="analysis/" )
    # return(summ.fit)
    if (debug) print(summ.fit)
    return(fit)
  }

#####################################################################################
# impStrat <- function(datNameImp, met=pmm, col_sel){
impStrat <- function(datImp, formulas, met=pmm, col_sel){
  # stratify by species
# impStrat <- function(datNameImp, met){
  # col_sel <- c(prVars, resp) # columns to select, as strings
  if(FALSE){
    datNameImp <- "ampDat"
    met <- metList
  }
  # datLT <- get(datNameImp) %>% filter(species=="LETE") %>% select(all_of(col_sel))
  datLT <- datImp %>% filter(species=="LETE") %>% select(all_of(col_sel))
  impLT <- mice::mice(datLT, seed=613, m=30, formulas=formulas, method=met, printFlag=FALSE)
  
  datCO <- datImp %>% filter(species=="CONI") %>% select(all_of(col_sel))
  impCO <- mice::mice(datCO, seed=613, m=30,formulas=formulas,method=met, printFlag=FALSE)
  
  # impInt <- mice::cbind(impLT, impCO)
  # impIntPMM <- mice::rbind(impLT, impCO)
  impInt <- mice::rbind(impLT, impCO)

  # imp_plot(impInt, plType = "density")
  # imp_plot(impInt, nest_age ~ fdate, plType="xy")
}

#####################################################################################
imp_small <- function(){
  impDat_small <- imputeDat( impDat1 = imppDat, return= "") # default m value is 10
  smallFit <- impFit(impDat_small, modList)
  smallPool <- poolImp(smallFit, modList)
  # Error in `brglm2::brglm_fit`(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
  #  could not find function "brglm2::brglm_fit"
  #       MAYBE HELPFUL:
  # pkg::name returns the value of the exported variable name in namespace pkg, whereas pkg:::name returns the value of the internal variable name. 
  impDat <- impDat_small
  impMod <- smallFit
  impModPool <- smallPool
}

#####################################################################################
runImpMods <- function(
    impDat1,
    m=100,
    method="pmm",
    modList,
    fam=binomial,
    regMethod="glm.fit",
    iter=500,
    resp,
    impute = TRUE,
    pool = FALSE,
    fit = TRUE,
    howMany=FALSE
    ){

  if (howMany) {
    m = 20
    pool = FALSE
  }
  if(impute | howMany){
    impDat  <- mice::mice(impDat1, seed=613, m=m, method=method, printFlag=FALSE)
    mice::complete(impDat)
    print(plot(impDat))
    print(attributes(impDat))
  } else {
    impDat = impDat1 # if you don't need to impute, the data is already in the proper format
  }

#####################################################################################
  impMod <- vector(mode="list", length=length(modList))
  if(regMethod == "brglm2") {

    regMet = brglm2::brglm_fit

    for (m in seq_along(modList)){
      impMod[[m]] = with(impDat, glm(as.formula(
        # paste0(resp, modList[[m]])), family=fam, method=brglm2::brglm_fit))
        paste0(resp, modList[[m]])), family=fam, method = regMet,
        control=brglmControl(maxit=iter)
        ))
  }
    cat("using brglmFit method")
  } else if(regMethod == "jeffreys") {
    regMet = brglm2::brglm_fit
    type   = "MPL_Jeffreys"
    cat("using Jeffreys's prior")
    for (m in seq_along(modList)){
      impMod[[m]] = with(impDat, glm(as.formula(
        # paste0(resp, modList[[m]])), family=fam, method=brglm2::brglm_fit))
        paste0(resp, modList[[m]])), family=fam, method = regMet, type=type,
        control=brglmControl(maxit=iter)
        ))
      # num[[m]] = how_many_imputations(impMod[[m]], cv=0.05, alpha=0.05)
  }
  } else if (regMethod == "glm.fit") {
    regMet = "glm.fit"
    for (m in seq_along(modList)){
      impMod[[m]] = with(impDat, glm(as.formula(
        # paste0(resp, modList[[m]])), family=fam, method=brglm2::brglm_fit))
        paste0(resp, modList[[m]])), family=fam, method = regMet
        ))
    }
  }
  if(howMany) {
    # num = howManyImputations::how_many_imputations()

    # impDat  <- mice::mice(impDat1, seed=613, m=20, method=method, printFlag=FALSE)
    num = list()
    for(m in seq_along(modList)){
      num[[m]] = how_many_imputations(impMod[[m]], cv=0.05, alpha=0.05)
    }
    # cat("optimal number of imputations (von Hippel 2020): ", num)
  }
  if(pool){
    saveRDS(impMod, "dataframes/impMod.rds")

    impModPool <- vector(mode="list", length=length(modList))

    for(m in seq_along(modList)){
      impModPool[[m]] = mice::pool(impMod[[m]])
      summary(impModPool[[m]])
      # impModPool[[m]] %>%
      #   gt::gt() # was this new? not working 12 aug 25
    }
    cat("returning pooled data sets")

    return(impModPool)
  } else if(howMany){
    cat("returning optimal number of imputations")
    return(num)
  } else if(fit) {
    cat("returning the fitted model results")
    return(impMod)
  } else if(impute) {
     cat(sprintf("returning all %s imputed data sets for each model", m))
    return(impDat)
  } else {
    cat("make sure only one of howMany and impute is set to TRUE")
    # cat("returning the fitted model results")
  }
}

#####################################################################################
impModTab <- function(impModPool, modList, fileSuffix=""){
  for(m in seq_along(impModPool)){
    imp_mod_tab = impModPool[[m]]$pooled

    fname = imp_mod_tab %>%
      gt::gt() %>%
      gt::fmt_number(decimals = 2) %>%
      # gt::cols_label() %>%
      # nestGLM::save_tab(rtf=TRUE)
      save_tab(rtf=TRUE)
  }
}

#####################################################################################
aicImpMods <- function(
    impModPool, modList, mods, dir, filePrefix=""
    ){

  aicImpMod <- list()
  mVal <- function(x) mean(as.numeric(vals[grep(x, names(vals))]))
  # mVal("glanced.AIC")
  extrNames <- c("glanced.AIC",           # names to extract
                 "glanced.deviance",
                 "glanced.df.residual",
                 "glanced.logLik")
  for(m in seq_along(impModPool)){
    vals = unlist(x=impModPool[[m]])[-1] # also remove the first element, which is a list
    aicImpMod[[m]] = lapply(extrNames, mVal)
    names(aicImpMod[[m]]) <- extrNames
  }
  saveRDS(aicImpMod, file = "dataframes/aicImpMod.rds")
  aicImpMod <- unlist(aicImpMod) # unlist so you can go through the values easily
  impAICtab <- data.frame(num = c(1:length(mods)),
                          modname    = modnames(mods),
                          dfResid    = as.integer(aicImpMod[grep(extrNames[3], names(aicImpMod))]),
                          logLik     = aicImpMod[grep(extrNames[4], names(aicImpMod))],
                          aic        = aicImpMod[grep(extrNames[1], names(aicImpMod))],
                          deltaAIC   = numeric(length=length(modList)),
                          akaikeWt   = numeric(length=length(modList)),
                          rdev       = aicImpMod[grep(extrNames[2], names(aicImpMod))]
                          )

  # impAICtab <- data.frame(modName=modList, aic=unlist(aicImpMod), aicWt=akaikeWts)
  impAICtab$akaikeWt <- MuMIn::Weights(impAICtab$aic)
  # impAICtab$
  impAICtab <- impAICtab[order(impAICtab$aic),] # order isn't working anymore?
  # impAICtab$deltaAIC <- c(0, diff(impAICtab$aic)) $ NOT RIGHT - calculates difference from previous value, not min
  impAICtab$deltaAIC <- c(impAICtab$aic - min(impAICtab$aic))
  saveRDS(impAICtab, file="dataframes/impAICtab.rds")
  impAIC <- gt::gt(impAICtab)  %>%
    gt::fmt_number(columns=c( logLik, aic, deltaAIC, akaikeWt),
                   decimals=2) %>%
    gt::cols_label(
    # cols_label( # why does this not have gt:: in front of it?
      modname  = "Model",
      dfResid  = "Degrees of\nFreedom",
      logLik   = "Log\nLikelihood",
      aic      = "AICc",
      deltaAIC = "delta AIC",
      akaikeWt = "Akaike Weight",
      rdev     = "Residual\nDeviance"
      )
  # filename <- impAIC %>% nestGLM::save_tab(dir = dir, suffix = filePrefix, rtf = TRUE)
  filename <- impAIC %>% save_tab(dir = dir, suffix = filePrefix, rtf = TRUE)
  return(filename)
}


#####################################################################################
# create the interaction plots
# grpColor <- c(LETE="#AA4499", CONI="#44AA99")
intPlot <- function(dat, 
                    modnum, 
                    plot_type = "pred", 
                    vars=c("species", "nest_age"),
                    y_lab = "Predicted probability of\n misclassification",
                    x_lab = "Nest age - centered (days)",
                    # grid=intGrid1 # probably shouldn't have the same name as a function?
                    # grid1=intGrid1 # probably shouldn't have the same name as a function?
                    grid1
                    ){
  
    model   <- with(dat,
                    glm(as.formula(paste0(resp, modList[[modnum]])),
                    family=binomial,
                    method=regMet
                    ))
    
    #  from the documentation: 
    #  Warning: Slopes and elasticities can only be calculated for continuous numeric variables. The slopes() functions will automatically revert to comparisons() for binary or categorical variables.
    #  
    #  Why are all the predictions only a few values?
    # intPred <- marginaleffects::predictions(model,
    #                                         # variables=c("species", "nest_age"),
    #                                         variables=c("species"), # try just including X
    #                                         newdata=grid1
    #                                         # so was this always working with these 2 variables?
    #                                         # I guess it was ignoring it when it was "nest_age"
    #                                         # except when the grid also varied nest_age, 
    #                                         # which is why the species*nest_age interaction plot
    #                                         # wasn't working
    #                                         # but anyway, I don't think I need the variables arg?
    #                                         # nope, still not working for the species*nest_age interaction
    #                                         # and the CI are much narrower now for the ones that work
    #                                         # could be because I changed the interaction order?
    #                                         # why does nest_age*species work, but not the reverse?
    #                                         # and why is the other one obs_int*species, not the reverse?
    #                                         # changing the order without putting the variables
    #                                         # argument back in didn't help anything...
    #                                         # variables=vars
    #                                         )
    
    # this one gives a list as output:
    # plDat   <- ifelse(plot_type=="slope",as.data.frame(intSlope),as.data.frame(intPred))
    # plDat <- intPred
    # if(plot_type=="slope") plDat <- intSlope
    # the confidence intervals in this plot don't seem wide enough...
    # plDat <- case_when(plot_type == "pred" ~ intPred,
    #                    plot_type == "slope" ~ intSlope)
    # ggplot(intPred,
    mm <- sym(vars[1])
    xx <- sym(vars[2])
    # colores <- c("")
    # g <- ggplot(plDat,
    
    # I didn't change anything about the ggplot function call
    # or really anything about this function definition, except for making it so that
    # the variables used for predictions aren't always species & nest age 
     if(plot_type=="slope"){
      
      intSlope <- marginaleffects::slopes(model,
                                          variables=c("species"),
                                          # variables=vars,
                        # variables=c("species", "obs_int"),
                        newdata=grid1)
      pl <- ggplot(intSlope,
             aes(x=!!xx, y=estimate, ymin=conf.low, ymax=conf.high,color=!!mm)) +
        labs(x=x_lab, y=y_lab)+
        scale_color_manual(values=grpColor,
                           labels=c("Common\nNighthawk", "Least Tern"),
                           name="Species") +
        geom_ribbon(alpha=0.2) +
        geom_line() +
        theme_classic() +
        theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20), # can't do margins for x and y axis titles at once
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        axis.title.y = element_text(margin=margin(r=5)),
        axis.title.x = element_text(margin=margin(t=5))) 
        
    } else {  # ggplot(plDat,
      intPred <- marginaleffects::predictions(model,
                                              # variables=c("species", "nest_age"),
                                              variables=c("species"), # try just including X
                                              newdata=grid1
      )  
      pl <- ggplot(intPred,
             # aes(x=obs_int, y=estimate, ymin=conf.low, ymax=conf.high, color=species)) +
             # aes(x=as.numeric(obs_int), y=estimate, ymin=conf.low, ymax=conf.high, color=species)) +
             aes(x=!!xx, y=estimate, color=!!mm)) +
        # geom_point() +
        labs(x=x_lab, y=y_lab)+
        scale_color_manual(values=grpColor,
                           labels=c("Common\nNighthawk", "Least Tern"),
                           name="Species") +
        geom_smooth(method=NULL) +
        theme_classic() +
        theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20), # can't do margins for x and y axis titles at once
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        axis.title.y = element_text(margin=margin(r=5)),
        axis.title.x = element_text(margin=margin(t=5))) 
        # geom_ribbon(alpha=0.2) +
        # geom_line()
      
    }
    
    return(pl)
}


#####################################################################################
# create_grid <- function(dat, var){
intPredict <- function(model, var){

  # nm <- paste("intGrid_", var)
  # grid
  # assign(x = nm, value=marginaleffects::datagrid(newdata=ndGLM_scl, grid_type="balanced"))
  # return(nm)
  # model   <- with(dat,
  #                 glm(as.formula(paste0(resp, modList[[modnum]])),
  #                 family=binomial,
  #                 method=regMet
  #                 ))
  # intPred <- 
  # grid1 <- marginaleffects::datagrid(newdata=ndGLM_scl, grid_type="balanced")
  # grid1 <- get(nm)
  intPred <- marginaleffects::predictions(model, newdata=grid1)
}

# if(grepl("*", modList[top_num], fixed=TRUE)){
# var <- "obs_int"
# # intGridNew1 <- marginaleffects::datagrid(newdata=ndGLM_scl,


# if(grepl("*", modList[top_num], fixed=TRUE)){ # the defaults are for this model
#   # change the grid used based on which variables are in the interaction:
#   if(intVar2 == "obs_int"){
#     # gr = intGridObs 
#     grName = "intGridObs"
#     lab = "Final interval - centered (days)"
#   } else {
#     # gr = intGridAge
#     grName = "intGridAge"
#     lab = "Nest age - centered (days)"
#   }
# 
#   intPlot(impDat,                             # (interaction between species & nest_age)
#           # modnum=5)
#           modnum=top_num,
#           vars = c(intVar1, intVar2), # in this case, intVar1 should always be species
#           # grid1=gr,
#           grid1=get(grName),
#           x_lab=lab,
#           y_lab=yl)
#   # ,
#           # y_lab="Predicted probability of\n unknown field fate")                   
# }
# # g <- get(grName)
# # rm(g)
# 
# if(grepl("*", modList[sec_top], fixed=TRUE)){
#   # if("obs_int" %in% c(intVar1.2, intVar2.2)){
#   if(intVar2.2 == "obs_int"){
#     # gr = intGridObs
#     grName = "intGridObs"
#     lab = "Final interval - centered (days)"
#   } else {
#     # gr = intGridAge
#     grName = "intGridAge"
#     lab = "Nest age - centered (days)"
#   }
#   
#   intPlot(impDat,
#             # modnum=4, 
#             modnum=sec_top, 
#             # grid1=gr, 
#             grid1=get(grName), 
#             # vars=c("species", "obs_int"), 
#             vars=c(intVar1.2, intVar2.2), 
#             x_lab=lab,
#           y_lab = yl)
# }


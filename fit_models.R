
fit_models <- function(
    dataMod, resp, fam=binomial, met="glm.fit", modlist,iter=1000
    ){
  # library(stats)
  # library(brglm2)
  # library(bbmle)

  mods <- list()
  if(met=="brglm2"){
    cat("using brglmFit method")

    for(m in seq_along(modlist)){
      cat("\nm:", m)
      # cat("\nmodel:", paste0(resp,modlist))
      # cat("\nmodel:", paste0(resp,modlist[m]))
      cat("\nmodel:", paste(resp,modlist[m]))
      mods[[m]] <- stats::glm(
      as.formula(
        # paste0(
        paste(
          # resp, modlist
          # this is pasting resp in front of each value of modlist
          # then somehow the list of mod calls comes out OK on the other end
          # but it messes up the modnames() function
          # I guess it's because I didn't do seq_along, so m was not a number
          # and then the elements of the mods list were named after m?
          # but the formula was identical for each of them. the names were misleading.
          # I guess it just took the first legitimate formula in the vector
          # even though it was fed the entire vector each time...
          # see "nest_models/mod_calls_bad.png" for examples
          resp, modlist[m], sep=" "
      )  ),
        data=dataMod,
        family=fam,
        method=brglm2::brglm_fit,
        # control=brglmControl(maxit=1000)
        control=brglmControl(maxit=iter)

      )}
  }else if(met=="jeffreys"){
    cat("using Jeffreys's prior")
    for(m in seq_along(modlist)){
      cat("\nm:", m)
      cat("\nmodel:", paste(resp,modlist[m]))
      mods[[m]] <- stats::glm(
        as.formula(
          paste(
            resp, modlist[m], sep=" "
          )
        ),
        data=dataMod,
        family=fam,
        method=brglm2::brglmFit,
        type="MPL_Jeffreys",
        control=brglmControl(maxit=iter)
      )
    }
  } else if(met=="logistf"){
    cat("using logistf method")
    for(m in seq_along(modlist)){
      cat("\nm:", m)
      cat("\nmodel:", paste(resp,modlist[m]))
      mods[[m]] <- logistf::logistf(
      as.formula(
        # paste0(
        paste(
          resp, modlist[m], sep=" "
        )),
        data=dataMod
      )}

  } else{
    for(m in seq_along(modlist)){
      cat("\nm:", m)
      cat("\nmodel:", paste(resp,modlist[m]))
      mods[[m]] <- stats::glm(
      as.formula(
        # paste0(
        paste(
          resp, modlist[m], sep=" "
        )),
        data=dataMod,
        family=fam,
        method=met

      )}
  }
  return(mods)
}

  # } else if(met={
  #   for(m in seq_along(modlist)){
  #     mods[[m]] <- stats::glm(
  #     as.formula(
  #       # paste0(
  #       paste(
  #         resp, modlist[m], sep=" "
  #       )),
  #       data=dataMod,
  #       family=fam,
  #       method=met
  #
  #     )}

############## UNIVARIATE MODELS #################################################
univar_mod <- function(resp, modData, retFormat="summary" ){
  # browser()
  vars <- c("species", "nest_age", "cam_fate", "fdate", "obs_int")
  # vars <- c("species", "nest_age", "c.fate", "fdate", "obs_int")
  # vars <- c("species", "nest_age", "cfate", "fdate", "obs_int")
  umod_ <- list()
  for(v in seq_along(vars)){
    f <- paste0(resp," ~ ", vars[v])
    # cat(f,"\t")                             # print the models as they are called
    umod_[[v]] <- stats::glm(f, data=modData, family="binomial")
  }

  modsum <- list()

  for(i in 1:5){
    # modsum[[i]] <- summary(umod_[[i]])
    # modsum[[i]] <- tbl_summary(umod_[[i]]) # try to make it a gtsummary object
    modsum[[i]] <- gtsummary::tbl_regression(umod_[[i]]) # tbl_summary is for variables; try this instead
  }
  modtab <- modData %>%
    dplyr::select(c(resp, "species", "nest_age", "cam_fate", "fdate", "obs_int")) %>%
    # select(c(resp, "species", "nest_age", "c.fate", "fdate", "obs_int")) |>
    gtsummary::tbl_uvregression(
      method       = "glm",
      y            = resp,
      method.args  = list(family=binomial),
      exponentiate = TRUE
  )
  modtab2 <- gtsummary::as_gt(modtab)
  # if(retFormat=="summary") modinfo = modsum else modinfo = modtab
  # modsum is a list of summaries, doesn't work with tbl_merge
  if(retFormat=="summary") modinfo = modtab else modinfo = modtab2
  # modinfo = list(modsum, modtab)
  return(modinfo ) # not sure how to return one or the other - just return both
}


############## MOD NAMES #################################################
modnames <- function(
    # mods, null=TRUE, debug=FALSE
    mods, null=TRUE, debug=FALSE
    ){
  dp <- function(x) deparse(x, width.cutoff=150)
  if(is.list(mods)){
    modEQ <- paste(sapply(
      mods,
      function(x) stringr::str_extract(dp(x$formula), "(?<=~\\s).*")
      # function(x) stringr::str_extract(x$formula, "(?<=~\\s)(.*)  ")
    ), sep = ",")
  } else {
    modEQ <- stringr::str_extract(mods, "(?<=~\\s).*")
  }
  modEQ <- stringr::str_replace_all(modEQ, c(
    "species" = "SPECIES",
    "obs_int" = "FINAL_INTERVAL",
    "cam_fate" = "TRUE_FATE",
    "nest_age" = "NEST_AGE",
    "fdate" = "END_DATE",
    "1"     = "NULL"))
  # cat("mod equations before substitution:", modEQ)
  if(debug) message(paste("\nmod equations before substitution: ", modEQ))
  # print(con=stdout(), "mod equations before substitution:", modEQ)
  m <- stringr::str_match(modEQ, "(\\w+)\\s\\*\\s(\\w+)" )
  # don't want to replace wheree this is no match (just replaces with NA)
  # modEQ <- stringr::str_replace(modEQ[!is.na(m)], coll(m[,1]), coll((paste0(m[,2]," + ",m[,3]," + ",m[,2],":",m[,3]))))
  modEQ <- ifelse(!is.na(m[,1]),
                  str_replace(modEQ, coll(m[,1]), paste0(m[,2]," + ",m[,3]," + ", m[,2],":",m[,3])),
                  modEQ
  )
  # if(debug) message(paste("\nand again after str_match", m))
  if(debug) message(paste("\nand again after str_match/str_replace:", modEQ))

  return(modEQ)
}

dp <- function(x) deparse(x, width.cutoff=150)

############## MOD_DIFF #################################################

mod_diff <- function(mOne, mAlt, modtype){ # argument order matters - names should help

  difVList <- which(!(names(stats::coef(mOne))) %in% names(stats::coef(mAlt)))

  comp <- stats::anova(mOne, mAlt, test="LRT")
  cat("\nCompare to full model:\n\n")
  print(comp)

  cat("\n\nPercent difference in model coefficients (vs full model):\n\n")
  cat(names(stats::coef(mOne))[-difVList], "\n")
  cat((abs(stats::coef(mAlt)-stats::coef(mOne)[-difVList]) / stats::coef(mOne)[-difVList])) #
  cat("\n---------------------------------------------------------------")
}
#
############## ODDS RATIOS #################################################
odd_rat <- function(mods){
  oddrats <- paste(sapply(mods, function(x) exp(x$coefficients)), sep=' ')
  names(oddrats) <- paste0("OR_mod",seq_along(mods))
  return(oddrats)
}

############## AIC TABLE #################################################
aic_mod <- function(mods){
 modEQ = modnames(mods)
  # modEQ = sapply(mods, function(x) stringr::str_extract(dp(x$formula)))
  nestAIC <- tryCatch(
    {
      AICcmodavg::aictab(cand.set = mods, modnames = modEQ, sort = T) %>%
        dplyr::mutate(df = K-1) %>%
        cat("k=", K, "df=", df) %>%
        cat("; mod names=\n", paste(Modnames, sep="\n")) %>%
        dplyr::select(Modnames, ,AICc, Delta_AICc, df, AICcWt, LL)
        # rename()
    },
    error=function(e){

        cat("using bbmle; mod names=\n", paste(modEQ, sep="\n"))
        bbmle::AICctab(mods, weights=TRUE, logLik=TRUE, base=TRUE, nobs=121, mnames=modEQ) %>%
          as.data.frame() %>%
          dplyr::select(df, AICc, dAICc, weight, logLik, dLogLik)
    }
  )
   AICtab <- nestAIC %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.numeric), ~round(.x, digits=2))) %>%
    gt::gt(rownames_to_stub = TRUE) %>%
    gt::tab_header(title=paste("AIC scores")) %>%
    gt::tab_options(
      column_labels.font.weight = "bold"
    ) %>%
  return(AICtab)
}

#'
############## BIC TABLE #################################################
bic_mod <- function(mods){

  modEQ <- paste(sapply(mods,
                        function(x) stringr::str_extract(dp(x$call), "(?<=~\\s).+(?=,\\sf)")
                        ), sep = " ")

  nestBIC <- AICcmodavg::bictab(cand.set = mods, modnames = modEQ, sort = T) # include VIF?

  BICtab <- gt::gt(nestBIC) %>%
    gt::tab_header(title="BIC scores") %>%
    gt::tab_style(style =
                gt::cell_text(size="small"),
              locations = gt::cells_body())
  return(BICtab)
}

############## REGRESSION TABLE #################################################
tabRegr <- function(modNum, mods, resp,
                    vars=c('species', 'nest_age', 'obs_int', 'cam_fate', 'fdate'),
                    newNames=c('SPECIES', 'NEST_AGE', 'FINAL_INTERVAL', 'TRUE_FATE', 'END_DATE'),
                    suffix="",
                    debug=F){
  # browser()
  now = format(Sys.time(), "%m%d_%H%M%S") # seconds bc tables all generated in < 1 min
  nm       <- function(x) substitute(x)
  if (resp=="is_u"){
    r = "MARKED_UNKNOWN"
  } else if(resp=="HF_mis" | resp=="misclass"){
    r = "MISCLASSIFIED"
  }
  if(debug) cat("response is", r, "\n\n")
  modEQ <- paste(sapply(
    mods,
    function(x) stringr::str_extract(dp(x$call), "(?<=~\\s)(.*)(?=\")")
    ), sep = " ")
  modVar <- stringr::word(modEQ)# extracts words
  if(debug) cat("extracted mod var:", modVar, "\n\n")
  # vars <- c('species', 'nest_age', 'obs_int', 'cam_fate', 'fdate')
  # newName <- c('SPECIES', 'NEST_AGE', 'FINAL_INTERVAL', 'TRUE_FATE', 'END_DATE')
  varNames <- ifelse(vars %in% modVar, paste(vars,"~'", newNames,"'",sep=""), NA)
  varNames <- varNames[!is.na(varNames)]
  vn <- as.list(varNames) #maybe all of this is not necessary?
  vnames <- as.list(newNames)
  if(debug) print(vnames)
  names(vnames) <- vars
  if(debug) print(vnames)
  tabOR <- gtsummary::tbl_regression(mods[[modNum]],
                          label = vnames,
                          exponentiate=T)
  headr <- paste("Regression Summary for Model:",
                 modnames(list(x = mods[[modNum]]), null=FALSE),
                 sep="\n")


  now = format(Sys.time(), "%m%d_%H%M_")
  filename3  <- sprintf("reg_%s_%s.rtf", now, suffix)
  filename4  <- sprintf("reg_%s_%s.png", now, suffix)
  # tabOR %>% gt() %>% gtsave(filename=filename3, path="analysis/", vwidth=1200, vheight=800)
  tabOR %>%
    gtsummary::as_gt() %>%
    # gt::tab_header(title=paste(question)) %>%
    gt::tab_header(title=headr) %>%
    gt::gtsave(filename=filename3, path="analysis/", vwidth=1200, vheight=800)


  tabOR %>%
    gtsummary::as_gt() %>%
    # gt::tab_header(title=paste(question)) %>%
    gt::tab_header(title=headr) %>%
    gt::gtsave(filename=filename4, path="analysis/", vwidth=1200, vheight=800)

  # return(tabOR)
  return( filename4)
}

############## SAVE TABLES #################################################
save_tab <- function(tabName, dpi=(1800/6), dir, suffix="", rtf=FALSE){

  now = format(Sys.time(), "_%m%d_%H%M%S") # seconds bc tables all generated in < 1 min
  nm       <- function(x) substitute(x)
  # dpi      <- (1800/6)                          # img width (px) / desired img width (in)
  if(rtf){
    file_name <- paste0(nm(tab), suffix, now, ".rtf") # name of tab will always be tab within the function?
    gt::gtsave(tab, file_name, path=paste0(dir,"analysis/"))
  }

  file_name1 <- paste0(nm(tab), now, suffix, ".png")
  gt::gtsave(tab, file_name1, path=paste0(dir,"analysis/"))

  # file_name2 <- paste0(nm(tab), now, "_rounded.png")
  # tab2 <- tab %>% fmt_number(decimals=2)
  # gtsave(tab2, file_name2, path-"analysis/")
  return(file_name1)
}
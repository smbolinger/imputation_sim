
###  to perform the multiple imputations ###############################################

imputeDat <- function(impDat1, m = 10, return = "complete"){
  
  # impDat1 <- get(datNameImp)
  # impDat  <- mice::mice(impDat1, seed=613, m=m, printFlag=FALSE)
  impDat  <- mice::mice(impDat1, seed=613, method="rf", m=m, printFlag=FALSE)
  impDat2 <- mice::complete(impDat, "all") # generate a list of all the completed datasets
  
  print(plot(impDat))
  print(attributes(impDat))
  
  if(return=="complete") return(impDat2) else return(impDat)
}

#####################################################################################
# Stratified multiple imputation (to deal with interactions)
impStrat <- function(datNameImp, met, col_sel){
  # If you don't edit the predictor matrix, get "logged events" where there are linear dependencies
  # in this case, I think it's because is_u and HF_mis are correlated
  dat <- get(datNameImp)
  cat(">> predictor matrix for imputation of", datNameImp, ":\n")
  print(mice::make.predictorMatrix(dat))
  # col_sel <- c(prVars, resp) # columns to select, as strings
  datLT <- get(datNameImp) %>% filter(species=="LETE") %>% select(all_of(col_sel))
  impLT <- mice::mice(datLT, seed=613, m=30, method=met, printFlag=FALSE)
  
  datCO <- get(datNameImp) %>% filter(species=="CONI") %>% select(all_of(col_sel))
  impCO <- mice::mice(datCO, seed=613, m=30,method=met, printFlag=FALSE)
  
  # impInt <- mice::cbind(impLT, impCO)
  impInt <- mice::rbind(impLT, impCO)

  # imp_plot(impInt, plType = "density")
  # imp_plot(impInt, nest_age ~ fdate, plType="xy")
}
### other ways of dealing with interactions: 
###   transform, then impute; 
###   impute, then transform; 
###   passive imputation; 
###   smcfcs

#### transform, then impute ############################################################

# get(datNameImp) %>% select(species, cam_fate, fdate, obs_int, nest_age) %>% mice::md.pattern()
# van Buuren found that transform, then impute gives more biased results unless you know your data are MCAR
# 
# # interaction with obs_int
# imppDat1 <- get(datNameImp) %>% 
#   mutate(species = as.numeric(species=="LETE"),
#          sp_obsInt = species * obs_int) %>% 
#   # select(species, cam_fate, fdate, nest_age, obs_int) 
#   select(species, obs_int, cam_fate, fdate, nest_age, sp_obsInt) 
# 
# imppDat2 <- get(datNameImp) %>% 
#   mutate(species = as.numeric(species=="LETE"),
#          sp_age = species * nest_age) %>% 
#   # select(species, cam_fate, fdate, nest_age, obs_int) 
#   select(species, obs_int, cam_fate, fdate, nest_age, sp_age) 
# # impDat       <- imputeDat(impDat1 = imppDat, m=100, return="")
# # impDat2 <- impDat
# # impDat  <- mice::mice(imppDat, seed=613, m=100, printFlag=FALSE)
# # impDat  <- mice::mice(imppDat, seed=613, m=30, printFlag=FALSE)
# impDat1  <- mice::mice(imppDat1, seed=613, m=30, printFlag=FALSE)
# impDat2  <- mice::mice(imppDat2, seed=613, m=30, printFlag=FALSE)
# 
# impDat1$predictorMatrix
# impDat2$predictorMatrix
# May actually only require 25 imputations, but that will change all the numbers slightly (AIC
# rankings stay the same)
# impDat  <- mice::mice(imppDat, seed=613, m=25, printFlag=FALSE)

#####################################################################################
# fit the models using each of the m datasets
impFit <- function(impDat, modList){
  
    impMod <- vector(mode="list", length=length(modList))
    
    for (m in seq_along(modList)){
      impMod[[m]] = with(impDat, glm(as.formula(
        # paste0(resp, modList[[m]])), family=fam, method=brglm2::brglm_fit))
        
        # how come regMet works in the nestGLM package but not here??
        paste0(resp, modList[[m]])), family=fam, method = regMet, 
        # paste0(resp, modList[[m]])), family=fam, method =impMet,
        control=brglmControl(maxit=iter)
      ))
    }
    
    return(impMod)
    
}

#####################################################################################
# pool the results
poolImp <- function(impMod, modList){
  impModPool <- vector(mode="list", length=length(modList))
  
  for(m in seq_along(modList)){
    impModPool[[m]] = mice::pool(impMod[[m]])
    summary(impModPool[[m]])
    # impModPool[[m]] %>%
    #   gt::gt() # was this new? not working 12 aug 25
  }
  return(impModPool)
}


# > mice::is.mids
is.mids <-  function (x) 
{
  inherits(x, "mids")
}
# <bytecode: 0x0000029113e13780>
#   <environment: namespace:mice>

#####################################################################################

# > mice:::rm.whitespace
rm.whitespace <- function (string, side = "both") 
{
  side <- match.arg(side, c("left", "right", "both"))
  pattern <- switch(side, left = "^\\s+", right = "\\s+$", 
                    both = "^\\s+|\\s+$")
  sub(pattern, "", string)
}
# <bytecode: 0x00000291311b1000>
#   <environment: namespace:mice>

#####################################################################################

mice.theme <- function (transparent = TRUE, alpha.fill = 0.3) 
{
  filler <- function(transparent, alpha) {
    if (transparent) {
      return(c(grDevices::hcl(240, 100, 40, alpha), grDevices::hcl(0,  100, 40, alpha)))
    }
    return(c(grDevices::hcl(240, 100, 40), grDevices::hcl(0,  100, 40)))
  }
  if (missing(transparent)) 
    transparent <- supports.transparent()
  if (missing(alpha.fill)) 
    alpha.fill <- ifelse(transparent, 0.3, 0)
  list(superpose.symbol = list(col = mdc(1:2),
                               fill = filler(transparent,  alpha.fill),
                               pch = 1),
       fontsize = list(text=22),
       superpose.line = list(col = mdc(4:5),
       # superpose.line = list(col = alpha(mdc(4:5), 0.5),
                             lwd = 2
                             ),
       box.dot = list(col = mdc(1:2)),
       box.rectangle = list(col = mdc(4:5)), 
       box.symbol = list(col = mdc(1:2)),
       plot.symbol = list(col = mdc(1:2),  
                          fill = filler(transparent, alpha.fill),
                          pch = 1), 
       plot.line = list(col = mdc(4:5)),
       superpose.polygon = list(col = filler(transparent,  alpha.fill)),
       strip.background = list(col = "grey95"), 
       mice = list(flag = TRUE))
}
# <bytecode: 0x00000291311a8540>
#   <environment: namespace:mice>
#####################################################################################

# > mice::supports.transparent
supports.transparent <- function () 
{
  query <- grDevices::dev.capabilities("semiTransparency")$semiTransparency
  if (is.na(query)) {
    query <- FALSE
  }
  query
}
# <bytecode: 0x00000291311ac758>
  # <environment: namespace:mice>

#####################################################################################
mdc <- function (r = "observed", s = "symbol", transparent = TRUE,
                 cso = grDevices::hcl(240,  100, 40, 0.7),
                 csi = grDevices::hcl(0, 100, 40, 0.7),
                 csc = "gray50",
                 # clo = grDevices::hcl(240, 100, 40, 0.8),
                 clo = grDevices::hcl(240, 100, 40, 0.5),
                 # cli = grDevices::hcl(0,  100, 40, 0.8),
                 cli = grDevices::hcl(0,  100, 40, 0.5),
                 clc = "gray50") 
{
  if (missing(transparent)) {
    if (!supports.transparent()) {
      cso <- grDevices::hcl(240, 100, 40)
      csi <- grDevices::hcl(0, 100, 40)
      csc <- "black"
      clo <- grDevices::hcl(240, 100, 40)
      cli <- grDevices::hcl(0, 100, 40)
      clc <- "black"
    }
  }
  else if (!transparent) {
    cso <- grDevices::hcl(240, 100, 40)
    csi <- grDevices::hcl(0, 100, 40)
    csc <- "black"
    clo <- grDevices::hcl(240, 100, 40)
    cli <- grDevices::hcl(0, 100, 40)
    clc <- "black"
  }
  fallback <- grDevices::palette()[1]
  if (is.numeric(r)) {
    idx <- floor(r)
    idx[r < 1 | r > 6] <- 7
    myc <- c(cso, csi, csc, clo, cli, clc, fallback)[idx]
    return(myc)
  }
  rc <- pmatch(r, c("observed", "missing", "both"))
  sc <- pmatch(s, c("symbol", "line"))
  idx <- rc + (sc - 1) * 3
  idx[is.na(idx)] <- 7
  myc <- c(cso, csi, csc, clo, cli, clc, fallback)[idx]
  myc
}
# <bytecode: 0x000002913112cb68>
#   <environment: namespace:mice>

#####################################################################################

if(FALSE){
  x <- impInt
  na.groups = NULL
  groups = NULL
  panel = lattice::lattice.getOption("panel.densityplot")
  default.prepanel = lattice::lattice.getOption("prepanel.default.densityplot")
  allow.multiple = TRUE
  thicker = 2.5
  outer = TRUE
  drop.unused.levels = lattice::lattice.getOption("drop.unused.levels")
  subscripts = TRUE
  as.table = TRUE
  plot.points = FALSE
  subset = TRUE
  mayreplicate = TRUE
  theme = mice.theme()
  xl    = xlabels
  dots <- list()
  # plType="xy"
  data <- nest_age ~ fdate
  plType="density"
  # something in the list might be causing the xyplot not to work?
}

#####################################################################################

vlabels <- c("nest_age" = "NEST_AGE",
             "cam_fate" = "TRUE_FATE",
             "fdate"    = "END_DATE",
             "species"  = "SPECIES",
             "obs_int"  = "FINAL_INTERVAL")

# > mice:::densityplot.mids
# density_plot<- function (x, data, na.groups = NULL, groups = NULL, as.table = TRUE, 
# missing_plot<- function (x, data, plType, na.groups = NULL, groups = NULL, as.table = TRUE, 
imp_plot<- function (x, data, plType, na.groups = NULL, groups = NULL, as.table = TRUE, 
          plot.points = FALSE, theme = mice.theme(), mayreplicate = TRUE, 
          thicker = 2.5, allow.multiple = TRUE, outer = TRUE, drop.unused.levels = lattice::lattice.getOption("drop.unused.levels"), 
          panel = lattice::lattice.getOption("panel.densityplot"), 
          default.prepanel = lattice::lattice.getOption("prepanel.default.densityplot"), 
          ..., subscripts = TRUE, subset = TRUE, vl = vlabels) 
{
  call <- match.call()
  
  if (!is.mids(x)) 
    stop("Argument 'x' must be a 'mids' object")
  
  cd <- data.frame(complete(x, "long", include = TRUE))
  r <- as.data.frame(is.na(x$data))
  
  nagp <- eval(expr = substitute(na.groups), envir = r, enclos = parent.frame())
  if (is.expression(nagp)){ 
    nagp <- eval(expr = nagp, envir = r, enclos = parent.frame())
  }
  
  ngp <- eval(expr = substitute(groups), envir = cd, enclos = parent.frame())
  if (is.expression(ngp)) {
    ngp <- eval(expr = ngp, envir = cd, enclos = parent.frame())
    }
  groups <- ngp
  
  ss <- eval(expr = substitute(subset), envir = cd, enclos = parent.frame())
  if (is.expression(ss)) {
    ss <- eval(expr = ss, envir = cd, enclos = parent.frame())
    }
  subset <- ss
  
  dots <- list(...)
  
  #####################################################################################
  args <- list(panel = panel, default.prepanel = default.prepanel, 
               allow.multiple = allow.multiple, outer = outer, drop.unused.levels = drop.unused.levels, 
               subscripts = subscripts, as.table = as.table, plot.points = plot.points)
  #####################################################################################
  args_xy <- list(allow.multiple = allow.multiple, outer = outer,
               drop.unused.levels = drop.unused.levels, subscripts = subscripts,
               as.table=as.table)
  #####################################################################################
  
  vnames <- setdiff(names(cd), c(".id", ".imp"))
  allfactors <- unlist(lapply(cd[vnames], is.factor))
  
  if (missing(data)) {
    vnames <- vnames[!allfactors & x$nmis > 2 & x$nmis < 
                       nrow(x$data) - 1]
    formula <- as.formula(paste("~", paste(vnames, collapse = "+", 
                                           sep = ""), sep = ""))
  } else {
    formula <- data
  }
  
  form <- lattice::latticeParseFormula(model = formula, data = cd, 
                                       subset = subset, groups = groups, multiple = allow.multiple, 
                                       outer = outer, subscripts = TRUE, drop = drop.unused.levels)
  
  xnames <- unlist(lapply(strsplit(form$right.name, " \\+ "),  rm.whitespace))
  ynames <- unlist(lapply(strsplit(form$left.name, " \\+ "), rm.whitespace))
  
  nona <- is.null(call$na.groups)
  if (!is.null(call$groups) && nona) {
    
    gp <- call$groups
    
  } else {
    
    if (nona) {
      
      if (plType=="density"){
        for (i in seq_along(xnames)) {
          xvar <- xnames[i]
          select <- cd$.imp != 0 & !r[, xvar]
          cd[select, xvar] <- NA
        }
        gp <- rep.int(cd$.imp, length(xnames)) # replicate values
      } else if (plType=="xy"){
        na.df <- r[, ynames, drop = FALSE]
        gp <- unlist(lapply(na.df, rep.int, x$m + 1))
      }
      
    } else {
      
      if (plType=="density"){
        for (i in seq_along(xnames)) {
          xvar <- xnames[i]
          select <- cd$.imp != 0 & !nagp
          cd[select, xvar] <- NA
        }
        gp <- rep.int(cd$.imp, length(xnames))
      } else if (plType == "xy"){
        gp <- rep.int(nagp, length(ynames) * (x$m + 1))
      }
      
    }
  }
  
  if(plType=="density"){
    
    mustreplicate <- !(!is.null(call$groups) && nona) && mayreplicate
    if (mustreplicate) {
      # theme$superpose.line$col <- rep(theme$superpose.line$col[seq_len(2)], # sequence of length 2
      # theme$superpose.line$col[seq_len(2)] # sequence of length 2
                                      # c(1, x$m))
      # theme$superpose.line$col <- c("deepskyblue2", grDevices::colorRampPalette(c("coral","orangered4"))(x$m))
      # theme$superpose.line$col <- c("deepskyblue2", grDevices::colorRampPalette(c("lightsalmon","orangered4"))(x$m))
      # theme$superpose.line$col <- c("deepskyblue2", grDevices::colorRampPalette(c("burlywood1","orangered4"))(x$m))
      # theme$superpose.line$col <- c("deepskyblue2",alpha( grDevices::colorRampPalette(c("burlywood1","orangered4"))(x$m), 0.5))
      # theme$superpose.line$col <- c("deepskyblue2",alpha( grDevices::colorRampPalette(c("pink","orangered4"))(x$m), 0.4))
      # theme$superpose.line$col <- alpha(col, 0.5)
      theme$superpose.line$col <- c("deepskyblue2", alpha(rep(x="deeppink4", times=x$m), 0.3))
      theme$superpose.line$lwd <- rep(c(theme$superpose.line$lwd[1] * 
                                          thicker, theme$superpose.line$lwd[1]), c(1, x$m))
      theme$superpose.symbol$col <- rep(theme$superpose.symbol$col[seq_len(2)], 
                                        c(1, x$m))# repeat the first one 1x, and the 2nd one m times
      theme$superpose.symbol$pch <- c(NA, 49:(49 + x$m - 1))
    }
    
    if (is.null(call$xlab)) {
      args$xlab <- ""
      if (length(xnames) == 1) 
        # args$xlab <- xnames
        args$xlab <- vl[xnames]
    
    }
  } else if (plType == "xy"){
    args <- args_xy
    if (is.null(call$ylab)) {
    args$ylab <- ""
    if (length(ynames) == 1) {
      # args$ylab <- ynames
      args$ylab <- vl[ynames]
      }
    # added this to xy plot to change the x label:
    if (length(xnames) == 1) {
      # args$xlab <- xnames
      args$xlab <- vl[xnames]
    }
    }
    
  }
  
  if (is.null(call$scales)) {
    args$scales <- list()
    if (length(xnames) > 1) {
      args$scales <- list(x = list(relation = "free"), 
                          y = list(relation = "free"))
    }
    if (length(ynames) > 1) {
      args$scales <- list(x = list(relation = "free"), 
                          y = list(relation = "free"))
    }
  }
  
  args <- c(x = formula, data = list(cd), groups = list(gp), 
            args, dots, subset = call$subset)
  if(plType=="density"){
    tp <- do.call(lattice::densityplot, args)
    tp <- update(tp, par.settings = theme)
  }else if(plType=="xy"){
    # tp <- do.call(lattice::xyplot, args_xy)
    tp <- do.call(lattice::xyplot, args)
    tp <- update(tp, par.settings = theme)
  }else{
    cat("invalid plType")
    return(c("invalid pl Type:", plType))
  }
  return(tp)
}
# tp
# <bytecode: 0x0000029131ab47d0>
#   <environment: namespace:mice>



#####################################################################################
### Call the function
#####################################################################################

# density_plot(impInt)
# missing_plot(impInt, plType="density")
# imp_plot(impInt, plType = "density")
# missing_plot(impInt, nest_age ~ fdate, plType="xy") # need data argument for xy plot 
# missing_plot(impInt, nest_age ~ obs_int, plType="xy") # need data argument for xy plot 




library(data.table)
library(tidyverse)
library(brglm2)

data.file <- "dat_complete.csv"
homeDir   <- "$HOME/Projects/fate_glm"
modList   <- readLines("modList.txt") 
mods4sim  <- modList[c(1,8,16)]
names(mods4sim) <- c("m1", "m8", "m16")

#realDat   <- read.csv(data.file)
#realDat$species <- relevel(as.factor(realDat$species), ref="LETE")
#realDat$cam_fate <- relevel(as.factor(realDat$cam_fate), ref="H")

#cMat <- list()
#cMat[["CONI"]] <- realDat %>% filter(species=="CONI") %>% select(nest_age, obs_int, fdate) %>% cor()
#mList <- list()

cMat <- readRDS("cormat.rds")
mList <- readRDS("means.rds")
# this is a list of vectors or something:
betaList <- readRDS("betas.rds")

print(names(betaList))

# fprob <- c(0.1, 0.06, 0.13, 0.06, 0.53, 0.13)
# names(fprob) <- c("A", "H", "F", "D", "S", "Hu")
# fprob <- c( 0.53,0.1, 0.06, 0.13, 0.06, 0.13)
# names(fprob) <- c( "H","A", "F", "D", "S", "Hu")

sprob <- c(0.32, 0.68)
names(sprob) <- c("CONI", "LETE")

resp <- "is_u"
mod <- mods4sim[[1]]
modInd <- which(mods4sim == mod)
#print(modInd)
#grepl(resp, names(betas))
ma <- grepl(resp, names(betaList)) & grepl(modInd, names(betaList))
#print(betaList[ma])
betas <- unlist(betaList[ma])
#print(ma)
#bb <- betas[resp]

prList <- names(betas)
prList[1] <- 1
prList[length(prList)] <- str_replace(prList[length(prList)], pattern=":", replacement="*")
print(prList)

mkResp <- function(sDat, betas, form, debug=FALSE, xdebug=FALSE){
    dummySim <- model.matrix(form, data = sDat)
    #dummySim <- as.data.table(model.matrix(form, data = sDat))
    # if(xdebug) cat("\n>> and with dummy variables:\n") # ~*~*~*
    # print(class(dummySim))
    #if(debug) aprint(dummySim)
    #if(debug) print.data.table(dummySim)
    # if(debug) bprint(dummySim)
    #if(debug) print(dummySim, topn=5, trunc.cols=T, )
    # if(debug) cat("\n>>>> beta values, class:", class(betas), "\n\n") # ~*~*~*
    # if(debug) print(betas)
    # print(head(betas))
    # if(debug) cat("\n>>> multiplied by design matrix, class:", str(dummySim), class(dummySim), "\n")
    #if(debug) aprint(dummySim)
    eta <- dummySim %*% betas
    # if(debug) cat("\nequals eta:\n") # ~*~*~*
    # if(debug) bprint(eta)
    y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
    return(y)
}

mkSim <- function(resp_list, mod, s_size, cMat, mList, betas,fprob, sprob, prList, debug=FALSE, xdebug=FALSE){
    fates <- names(fprob)
    spp   <- names(sprob)
    if(debug) cat("\nfates:", fates, "& species:", spp) # ~*~*~*
    # if(debug){ # ~*~*~*
    #     cat("\ncoefficients, class:\n", class(betas))
    #     print(betas)
    # }
    simDat <- list()
    for (s in seq_along(spp)){
        sp <- spp[s]
        #simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=s_size*sprob[s], mu=mList[[sp]], Sigma=cMat[[sp]]))
        # cat("\n>>> species:", sp, "\t& means:", mList[[sp]], "\t& correlation matrix:\n") # ~*~*~*
        # print(cMat[[sp]])
        simDat[[sp]] <- as.data.table(MASS::mvrnorm(n=round(s_size*sprob[s]), mu=mList[[sp]], Sigma=cMat[[sp]]))
        # cat("\n>> adding camera fates\n") # ~*~*~*
        simDat[[sp]][,cam_fate := sample(fates, size=nrow(simDat[[sp]]), replace=T, prob=fprob)]
    }
    cat("\n>>> created sim data") # ~*~*~*
    names(simDat) <- c("CONI", "LETE")
    sDat <- data.table::rbindlist(simDat, idcol="species", use.names=T, fill=T)
    cat(" - merged sim data for the two species") # ~*~*~*
    sDat[,cam_fate := relevel(as.factor(sDat[,cam_fate]), ref="H")]
    sDat[,species  := relevel(as.factor(sDat[,species] ), ref="LETE")]
    cat(" - releved factors") # ~*~*~*
    form <- as.formula(mod)
    #dummySim <- model.matrix(mod, data = sDat)
    # if(xdebug) cat("\n>>>>> simulated data frame w/o dummy variables:\n") # ~*~*~*
    # if(xdebug) print(sDat)
    #sDat[,is_u := mkResp(sDat, as.matrix(betas[[1]]), form, debug=debug )]
    #sDat[,HF_mis := mkResp(sDat, as.matrix(betas[[2]]), form, debug=debug )]
    sDat[,is_u := mkResp(sDat, (betas[[1]]), form, debug=debug )]
    sDat[,HF_mis := mkResp(sDat, (betas[[2]]), form, debug=debug )]
    cat(' - added response variables\n') # ~*~*~*
    #for (r in seq_along(resp_list)){
    #    #resp <- mkResp(sDat, betas, form, debug=debug)
    #    sDat[,resp_list[r] := mkResp(sDat, betas[[r]], form, debug=debug)]
    #}

    #dummySim <- model.matrix(form, data = sDat)
    #if(debug) cat("\ndesign matrix, class:\n", str(dummySim))
    #if(debug) print(dummySim)
    #eta <- dummySim %*% betas
    #if(debug) cat("\nmultiplied by beta values:\n")
    #print(eta)
    #y <- rbinom(nrow(dummySim), size = 1, prob = binomial()$linkinv(eta)) # The outcome    
    #sDat[,is_u := y]
    #sDat[,is_u := resp]
    # if(xdebug){ # ~*~*~*
    #     fitReal2 <- summary(glm(is_u ~ nest_age * species + obs_int + cam_fate + fdate,
    #                        family=binomial,
    #                        data=sDat,
    #                        method=brglm2::brglmFit))
    #     cat("\n>>>> model fit for simulated data:\n")
    #     print(fitReal2)
    # }
    return(sDat)
}

#m <- mkSim(resp="is_u", mod=mods4sim[[1]], s_size=500, cMat=cMat, mList=mList, betas=betas, fprob=fprob, sprob=sprob, prList=prList, debug=TRUE)



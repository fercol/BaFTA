# ============================= CODE METADATA ================================ #
# PACKAGE: BaFTA
# AUTHOR: Fernando Colchero
# DATE CREATED: 2023-08-17
# DESCRIPTION: Functions to estimate average age-specific fertility from
#              alternative model with mixed effects for repeated individuals.
# COMMENTS: dataType: aggregated ("aggregated"), individual seasonal 
#           ("indivSimple"), individual continuous age with IBI 
#           ("indivExtended").
# ============================== START CODE ================================== #
# ======================================== #
# A) FUNCTIONS AVAILABLE TO THE USER: ==== 
# ======================================== #
# Main bafta function:
bafta <- function(object, ...) UseMethod("bafta")

bafta.default <- function(object, dataType = "aggregated",
                          model = "quadratic", minAge = NA, gestTime = NA, 
                          niter = 55000, burnin = 5001, thinning = 20, 
                          nsim, ncpus, UPDJUMP = TRUE, jumpSD = NULL, ...) {
  # Extract arguments:
  argList <- list(...)
  argNames <- names(argList)
  
  # BaFTAstart timer:
  BaFTAstart <- Sys.time()
  
  # Fertility models:
  models <- c("quadratic", "PeristeraKostaki", "ColcheroMuller", 
              "Hadwiger", "gamma", "beta", "skewNormal", 'gammaMixture',
              "HadwigerMixture", "skewSymmetric", "skewLogistic")
  nmods <- length(models)
  modList <- paste(paste("'", models, "'", sep = ""), collapse = ", ")
  
  # Data types:
  dTypes <- c("aggregated", "indivSimple", "indivExtended")
  dTypeList <- paste(paste("'", dTypes, "'", sep = ""), collapse = ", ")
  
  # Check model argument:
  if (!model %in% models) {
    stop(sprintf("Wrong model specification. Available models are: %s.\n", 
                 modList))
  } else {
    # Create fertility functions:
    FertFun <- function(beta, ...) UseMethod("FertFun")
    FertFun.matrix <- .DefineFertilityMatrix(modelFert = model)
    FertFun.numeric <- .DefineFertilityNumeric(modelFert = model)
  }
  
  # check dataType argument:
  if (!dataType %in% dTypes) {
    stop(sprintf("Wrong dataType specification. Available data types are: %s.\n", 
                 dTypeList))
  }
  
  # Logical for random effects in fertility:
  # if (grepl("indiv", dataType)) {
  #   RANDEFFU <- TRUE
  # } else {
  #   RANDEFFU <- FALSE
  # }
  RANDEFFU <- FALSE
  
  # Logical for random effects in IBI:
  if (dataType == "indivExtended") {
    RANDEFFV <- TRUE
  } else {
    RANDEFFV <- FALSE
  }
  
  # Create algorithm object:
  algObj <- .CreateAlgObj(model, dataType, minAge, gestTime, formula, niter, 
                          burnin, nsim, thinning, UPDJUMP, jumpSD)
  
  # Create data object:
  dataObj <- .CreateDataObj(object = object, algObj = algObj)
  
  # Create parameter object:
  parObj <- .BuildParObj(algObj = algObj, dataObj = dataObj)
  
  # Update parObj with user parameters:
  parObj <- .CreateUserPar(argList = argList, argNames = argNames, 
                           parObj = parObj, dataObj = dataObj)
  
  # Number of iterations kept for inference:
  keep <- seq(burnin, niter, thinning)
  nkeep <- length(keep)
  # keepind <- burnin:niter

  # Run jump sd:
  if (UPDJUMP) {
    cat("\nRunning sequence to find jump SDs... ")
    Start <- Sys.time()
    
    # Number of iterations for jump sds:
    niterJump <- 10000
    
    # Run jump SD sequence:
    outJump <- .RunMCMC(sim = 1, dataObj = dataObj, parObj = parObj, 
                        niter = niterJump, algObj = algObj, FertFun = FertFun,
                        FertFun.numeric = FertFun.numeric, 
                        FertFun.matrix = FertFun.matrix,
                        jumpSD = NULL, UPDJUMP = TRUE)
    
    # Extract jump SD vector:
    jumpSD <- outJump$jumps
    End <- Sys.time()
    cat("Done\n")
    compTime <- round(as.numeric(End-Start, units = units(End - Start)), 2)
    cat(sprintf("Total jump SDs computing time: %.2f %s.\n\n", compTime, 
                units(End - Start)))
  } else {
    if (is.null(jumpSD) | length(jumpSD) != parObj$pSamp) {
      stop(sprintf("Length of 'jumpSD' argument should be %s.\n", 
                   parObj$pSamp))
      
    }
  }
  
  # Variables to be loaded to the CPUS:
  # cpuVars <- c(".rtnorm", ".dtnorm", ".qtnorm", ".CalcLikeFert",
  #              ".CalcLikeFert.baftaAggr", ".CalcLikeFert.baftaIndSimp",
  #              ".CalcLikeFert.baftaIndExt", ".CalcPostTheta",
  #              ".CalcPostRandEffU", ".CalcPostRandEffV", ".CalcMHratio",
  #              ".SampleUSig", ".SampleVSig", "FertFun",
  #              "FertFun.numeric", "FertFun.matrix")

  
  cat("Multiple simulations started...\n\n") 
  
  # BaFTAstart parallel computing:
  sfInit(parallel = TRUE, cpus = ncpus)
  
  # Load BaFTA to CPUS:
  sfLibrary("BaFTA", character.only = TRUE,
            warn.conflicts = FALSE)
  # sfLibrary(BaFTA)
  # sfSource("pkg/R/bafta.R")
  
  # Load variables to CPUS:
  # sfExport(list = cpuVars)
  
  # Run MCMC in parallel:
  outMCMC <- sfClusterApplyLB(1:nsim, .RunMCMC, dataObj = dataObj, 
                              parObj = parObj, niter = niter, algObj = algObj,
                              FertFun = FertFun,
                              FertFun.numeric = FertFun.numeric, 
                              FertFun.matrix = FertFun.matrix,
                              jumpSD = jumpSD, UPDJUMP = FALSE)
  
  # Stop cluster:
  sfStop()
  
  cat("Simulations finished.\n")
  
  # Extract variables for coefficients and DIC:
  for (ic in 1:nsim) {
    if (ic == 1) {
      thetaMat <- outMCMC[[ic]]$theta[keep, ]
      likeMat <- outMCMC[[ic]]$likePost[keep, ]
      if (RANDEFFU) {
        # uSdVec <- outMCMC[[ic]]$uSd[keep]
        uMat <- outMCMC[[ic]]$u[keep, ]
      } else {
        uMat <- NA
      }
      if (RANDEFFV) {
        # vSdVec <- outMCMC[[ic]]$vSd[keep]
        vMat <- outMCMC[[ic]]$v[keep, ]
      } else {
        vMat <- NA
      }
      
    } else {
      thetaMat <- rbind(thetaMat, outMCMC[[ic]]$theta[keep, ])
      likeMat <- rbind(likeMat, outMCMC[[ic]]$likePost[keep, ])
      if (RANDEFFU) {
        # uSdVec <- c(uSdVec, outMCMC[[ic]]$uSd[keep])
        uMat <- rbind(uMat, outMCMC[[ic]]$u[keep, ])
      } 
      if (RANDEFFV) {
        # vSdVec <- c(vSdVec, outMCMC[[ic]]$vSd[keep])
        vMat <- rbind(vMat, outMCMC[[ic]]$v[keep, ])
      } 
    }
  }
  
  # Calculate convergence statistics:
  Conv <- .CalcPSRF(object = outMCMC, keep = keep, nsim = ncpus)
  
  # Effective sample size:
  Neff <- .CalcNeff(object = outMCMC, keep = keep, nsim = ncpus, Rhat = Conv)
  
  # Extract average parameters:
  idPars <- parObj$idSamp
  if (inherits(dataObj, "baftaIndExt")) {
    idPars <- c(idPars, max(idPars) + 1)
  }
  coeffs <- cbind(Mean = apply(thetaMat[, idPars],  2, mean), 
                  SD = apply(thetaMat[, idPars], 2, sd),
                  Lower = apply(thetaMat[, idPars], 2, quantile, 0.025),
                  Upper = apply(thetaMat[, idPars], 2, quantile, 0.975),
                  Neff = Neff[idPars], Rhat = Conv[idPars, "Rhat"])
  
  if (RANDEFFU) {
    # Include random effect standard error:
    # coeffs <- rbind(coeffs, uSd = c(Mean = mean(uSdVec), SD = sd(uSdVec),
    #                                 Lower = quantile(uSdVec, 0.025),
    #                                 Upper = quantile(uSdVec, 0.975),
    #                                 Neff = length(keep) * nsim,
    #                                 Rhat = 1))
    # Merge thetaMat with uSd:
    # thetaMat <- cbind(thetaMat, uSD = uSdVec)
    
    # Random effect quantile:
    uQuant <- cbind(Mean = apply(uMat, 1, mean, na.rm = TRUE), 
                    Lower = apply(uMat, 1, quantile, 0.025, na.rm = TRUE),
                    Upper = apply(uMat, 1, quantile, 0.975, na.rm = TRUE))
  } else {
    uSdVec <- NA
    uQuant <- NA
  }
  
  if (RANDEFFV) {
    # Include random effect standard error:
    # coeffs <- rbind(coeffs, vSd = c(Mean = mean(vSdVec), SD = sd(vSdVec),
    #                                 Lower = quantile(vSdVec, 0.025),
    #                                 Upper = quantile(vSdVec, 0.975),
    #                                 Neff = length(keep) * nsim,
    #                                 Rhat = 1))
    # Merge thetaMat with uSd:
    # thetaMat <- cbind(thetaMat, uSD = uSdVec)
    
    # Random effect quantile:
    vQuant <- cbind(Mean = apply(vMat, 1, mean, na.rm = TRUE), 
                    Lower = apply(vMat, 1, quantile, 0.025, na.rm = TRUE),
                    Upper = apply(vMat, 1, quantile, 0.975, na.rm = TRUE))
  } else {
    vSdVec <- NA
    vQuant <- NA
  }
  
  # Summary for unknown ages:
  if (dataObj$AgeUpdate) {
    # Index of values to be kept (after burn-in):
    keepAges <- algObj$burnin:algObj$niter
    
    # Age Estimation matrix:
    ageMat <- outMCMC[[1]]$age[keepAges, ]
    fullAgeMat <- matrix(dataObj$data$Age, nrow = nkeep, ncol = dataObj$n,
                         byrow = TRUE)
    fullAgeMat[, dataObj$idAgeUpd] <- t(outMCMC[[1]]$age[keep, ])
    for (isim in 2:algObj$nsim) {
      ageMat <- rbind(ageMat, outMCMC[[isim]]$age[keepAges, ])
      tempFull <- matrix(dataObj$data$Age, nrow = nkeep, ncol = dataObj$n,
                         byrow = TRUE)
      tempFull[, dataObj$idAgeUpd] <- t(outMCMC[[isim]]$age[keep, ])
      fullAgeMat <- rbind(fullAgeMat, tempFull)
    }
    
    # ID ages:
    agesID <- object$indID[dataObj$idFirstAge]
    
    # Summary of age estimation:
    ageQuant <- data.frame(indID = agesID, Mean = apply(ageMat, 2, mean),
                      SD = apply(ageMat, 2, sd),
                      Lower = apply(ageMat, 2, quantile, prob = 0.025, 
                                    names = FALSE),
                      Upper = apply(ageMat, 2, quantile, prob = 0.975, 
                                    names = FALSE),
                      userAge = object$Age[dataObj$idFirstAge])
    # Create full age matrix for predictive loss calculations:
    
  } else {
    ageQuant <- NA
    fullAgeMat <- NA
  }
  
  
  # Summary for unknown offspring:
  if (dataObj$offsUpdate) {
    # Index of values to be kept (after burn-in):
    keepAges <- algObj$burnin:algObj$niter
    
    # Age Estimation matrix:
    offsMat <- outMCMC[[1]]$offsOut[keepAges, ]
    fullOffsMat <- matrix(dataObj$data$nOffspring, nrow = nkeep, 
                          ncol = dataObj$n, byrow = TRUE)
    fullOffsMat[, dataObj$idOffsUpd] <- t(outMCMC[[1]]$offsOut[keep, ])
    
    for (isim in 2:algObj$nsim) {
      offsMat <- rbind(offsMat, outMCMC[[isim]]$offsOut[keepAges, ])
      tempFull <- matrix(dataObj$data$nOffspring, nrow = nkeep, 
                         ncol = dataObj$n, byrow = TRUE)
      tempFull[, dataObj$idOffsUpd] <- t(outMCMC[[isim]]$offsOut[keep, ])
      fullOffsMat <- rbind(fullOffsMat, tempFull)
      
    }
    
    # ID ages:
    offsID <- object$indID[dataObj$idOffsUpd]
    
    # Summary of age estimation:
    offsQuant <- data.frame(indID = offsID, Mean = apply(offsMat, 2, mean),
                       SD = apply(offsMat, 2, sd),
                       Lower = apply(offsMat, 2, quantile, prob = 0.025, 
                                     names = FALSE),
                       Upper = apply(offsMat, 2, quantile, prob = 0.975, 
                                     names = FALSE),
                       obsOffs = object$nOffspring[dataObj$idOffsUpd])
    
  } else {
    offsQuant <- NA
    fullOffsMat <- NA
  }
  
  # Calculate DIC:
  DIC <- .CalcDIC(likelihood = likeMat[, "Likelihood"], 
                  k = length(parObj$idSamp))
  
  # Y predicted:
  yPred <- .CalcYpred(dataObj = dataObj, thetaMat = thetaMat, uMat = uMat,
                      ageMat = fullAgeMat, offsMat = fullOffsMat, 
                      parObj = parObj, FertFun = FertFun, 
                      FertFun.numeric = FertFun.numeric)
  
  # Calculate predictive loss:
  PredLoss <- .CalcPredLoss(dataObj = dataObj, yPred = yPred)
  
  # Extract aggregated fertility and number of offspring:
  
  if (algObj$dataType == "aggregated") {
    aggrData <- object[which(object$Age >= dataObj$alpha), ]
    dx <- diff(aggrData$Age[1:2])
    aggrData$Age <- aggrData$Age + dx / 2
  } else {
    xag <- dataObj$alpha:ceiling(max(object$Age))
    nxag <- length(xag)
    indID <- unique(object$indID)
    indObs <- t(sapply(indID, function(iid) {
      idi <- which(object$indID == iid)
      agei <- range(floor(object$Age[idi]))
      return(agei)
    }))
    
    tempag <- t(sapply(1:(nxag - 1), function(ix) {
      idi <- which(indObs[, 2] >= xag[ix] & indObs[, 1] < xag[ix + 1])
      idb <- which(object$Age >= xag[ix] & object$Age < xag[ix + 1])
      nPariAll <- length(idi)
      nPariAvail <- length(idb)
      nOffs <- sum(object$nOffspring[idb])
      return(c(nParents = nPariAvail, nOffspring = nOffs, 
               nParentsAll = nPariAll))
    }))
    dx <- diff(xag[1:2])
    aggrData <- data.frame(Age = xag[-nxag] + dx / 2, tempag, 
                           Fertility = tempag[, 2] / tempag[, 1],
                           RealizedFert = tempag[, 2] / tempag[, 3])
  }
  
  # Calculate mean fertility and 95% credible intervals:
  xv <- seq(dataObj$alpha, ceiling(max(dataObj$data$Age)), 0.1) - 
    dataObj$alpha
  
  # Calculate estimated fertility from parameter posteriors:
  
  fertAll <- apply(thetaMat, 1, function(be) {
    fert <- FertFun(be, xv)
    return(fert)
  })
  
  # Calculate mean and quantiles:
  fertQuant <- cbind(Mean = apply(fertAll, 1, mean, na.rm = TRUE), 
                     Lower = apply(fertAll, 1, quantile, 0.025, na.rm = TRUE),
                     Upper = apply(fertAll, 1, quantile, 0.975, na.rm = TRUE))
  
  # Calculate predictive quantiles:
  if (algObj$dataType == "aggregated") {
    predQuant <- cbind(Mean = apply(yPred, 2, mean, na.rm = TRUE), 
                       Lower = apply(yPred, 2, quantile, 0.025, na.rm = TRUE),
                       Upper = apply(yPred, 2, quantile, 0.975, na.rm = TRUE))
    
  } else {
    xag <- dataObj$alpha:ceiling(max(object$Age))
    nxag <- length(xag)
    yPredAgg <- sapply(1:(nxag - 1), function(ix) {
      idb <- which(object$Age >= xag[ix] & object$Age < xag[ix + 1])
      if (length(idb) > 1) {
        nOffs <- apply(yPred[, idb], 1, sum)
      } else if (length(idb) == 1) {
        nOffs <- yPred[, idb]
      } else {
        nOffs <- rep(0, nrow(yPred))
      }
      return(nOffs)
    })
    predQuant <- cbind(Mean = apply(yPredAgg, 2, mean, na.rm = TRUE), 
                       Lower = apply(yPredAgg, 2, quantile, 0.025, 
                                     na.rm = TRUE),
                       Upper = apply(yPredAgg, 2, quantile, 0.975, 
                                     na.rm = TRUE))
  }
  
  # end timer:
  BaFTAend <- Sys.time()
  compTime <- sprintf("%s mins", signif(as.numeric(BaFTAend - BaFTAstart, 
                                                   units = "mins"), 2))
  timeDiff <- BaFTAend - BaFTAstart
  diffUnits <- units(timeDiff)
  if (diffUnits == "secs" & as.numeric(timeDiff) < 60) {
    tdiffUn <- round(as.numeric(BaFTAend - BaFTAstart, 
                                units = 'secs'), 2)
    compTime <- sprintf("%s secs", tdiffUn)
  } else if (diffUnits == "secs" & as.numeric(timeDiff) >= 60 | 
      diffUnits == "mins" & as.numeric(timeDiff) < 60) {
    tdiffUn <- round(as.numeric(BaFTAend - BaFTAstart, 
                                      units = 'mins'), 2)
    compTime <- sprintf("%s mins", tdiffUn)
  } else if (diffUnits == "mins" & as.numeric(timeDiff) >= 60 | 
             diffUnits == "hours") {
    tdiffUn <- round(as.numeric(BaFTAend - BaFTAstart, 
                                      units = 'hours'), 2)
    compTime <- sprintf("%s hours", tdiffUn)
  } 
  cat(sprintf("Total MCMC computing time: %s\n\n", compTime))
  
  # Settings:
  settings <- list(model = model, dataType = dataType, niter = niter, 
                   burnin = burnin, thinning = thinning, nsim = nsim, 
                   ncpus = ncpus, compTime = compTime)
  
  # store output:
  fullOut <- list(coefficients = coeffs, x = xv, fert = fertQuant, 
                  theta = thetaMat, uSd = uSdVec, 
                  likePost = likeMat, DIC = DIC, PredLoss = PredLoss, 
                  pred = predQuant, estAges = ageQuant, estOffs = offsQuant, 
                  aggrData = aggrData, runs = outMCMC, data = dataObj, 
                  settings = settings, keep = keep, params = parObj)
  
  class(fullOut) <- "bafta"
  return(fullOut)
}

# Print BaFTA:
print.bafta <- function(x, ...) {
  extraArgs <- list(...)
  if (length(extraArgs) > 0) {
    if (!is.element('digits', names(extraArgs))){
      digits <- 4
    } else {
      digits <- extraArgs$digits
    }
  } else {
    digits <- 4
  }
  # Call:
  cat("\nCall:\n")
  cat(paste("Model          \t: ", x$settings$model, "\n", sep = ""))
  cat(paste("Data type      \t: ", x$settings$dataType, "\n", sep = ""))
  cat(paste("Num. iterations\t: ", x$settings$niter, "\n", sep = ""))
  cat(paste("Burnin         \t: ", x$settings$burnin, "\n", sep = ""))
  cat(paste("Thinning       \t: ", x$settings$thinning, "\n", sep = ""))
  cat(paste("Number of sims.\t: ", x$settings$nsim, "\n", sep = ""))
  cat(paste("Computing time \t: ", x$settings$compTime, "\n", sep = ""))
  
  # Coefficients:
  cat("\nCoefficients:\n")
  print.default(x$coefficients, digits, ...)
  
  # Convergence:
  cat("\nConvergence:\n")
  if (all(!is.na(x$coefficients[, "Rhat"])) & 
      all(x$coefficients[, "Rhat"] < 1.05)) {
    convMessage <- "All parameter chains converged.\n"
  } else {
    convMessage <- "Some parameter chains did not converge.\n"
  }
  cat(convMessage)
  
  # DIC:
  cat("\nModel fit:\n")
  if (is.na(x$DIC[1])){
    cat("\nDIC not calculated.")
  } else {
    cat(sprintf("\nDIC = %s\n", round(x$DIC["DIC"], 2)))
  }
  
  # Predictive loss:
  cat("\nPredictive loss:\n")
  print.default(x$PredLoss, digits = 3)
  
}

# Summary BaFTA:
summary.bafta <- function(object, ...) {
  extraArgs <- list(...)
  if (length(extraArgs) > 0) {
    if (!is.element('digits', names(extraArgs))){
      digits <- 4
    } else {
      digits <- extraArgs$digits
    }
  } else {
    digits <- 4
  }
  # Call:
  cat("\nCall:\n")
  cat(paste("Model          \t: ", object$settings$model, "\n", sep = ""))
  cat(paste("Data type      \t: ", object$settings$dataType, "\n", sep = ""))
  cat(paste("Num. iterations\t: ", object$settings$niter, "\n", sep = ""))
  cat(paste("Burnin         \t: ", object$settings$burnin, "\n", sep = ""))
  cat(paste("Thinning       \t: ", object$settings$thinning, "\n", sep = ""))
  cat(paste("Number of sims.\t: ", object$settings$nsim, "\n", sep = ""))
  cat(paste("Computing time \t: ", object$settings$compTime, "\n", sep = ""))
  
  # Coefficients:
  cat("\nCoefficients:\n")
  print.default(object$coefficients, digits, ...)
  
  # Convergence:
  cat("\nConvergence:\n")
  if (all(!is.na(object$coefficients[, "Rhat"])) & 
      all(object$coefficients[, "Rhat"] < 1.05)) {
    convMessage <- "All parameter chains converged.\n"
  } else {
    convMessage <- "Some parameter chains did not converge.\n"
  }
  cat(convMessage)
  
  # DIC:
  cat("\nModel fit:\n")
  if (is.na(object$DIC[1])){
    cat("\nDIC not calculated.")
  } else {
    cat(sprintf("\nDIC = %s\n", round(object$DIC["DIC"], 2)))
  }
  
  # Predictive loss:
  cat("\nPredictive loss:\n")
  print.default(object$PredLoss, digits = 3)
  
}

# Plotting function:
plot.bafta <- function(x, type = "traces", ...) {
  argList <- list(...)
  plTypes <- c("traces", "density", "fertility", "predictive")
  if (!type %in% plTypes) {
    stop(sprintf("%s.\n%s %s.\n", "Wrong type of plot for object of class 'bafta'", 
                 "Available types are", 
                 paste(paste("'", plTypes, "'", sep = ""), collapse = ", ")))
  }
  if (type == "traces") {
    idkeep <- seq(1, x$settings$niter, x$settings$thinning)
    pSamp <- x$params$pSamp
    idSamp <- x$params$idSamp
    pName <- x$params$thetaName[idSamp]
    if (x$settings$dataType == "indivExtended") {
      pSamp <- pSamp + 1
      idSamp <- c(idSamp, max(idSamp) + 1)
      pName <- c(pName, "vSd")
    }
    nsim <- x$settings$nsim
    ylim <- sapply(1:pSamp, function(ipar) {
      range(sapply(1:nsim, function(isim) {
        range(x$runs[[isim]]$theta[, idSamp[ipar]], na.rm = TRUE)
      }), na.rm = TRUE)
    })
    pPars <- pSamp
    if (grepl("indiv", x$settings$dataType)) {
      # uyl <- range(sapply(1:nsim, function(isim) {
      #   range(x$runs[[isim]]$uSd)
      # }))
      # ylim <- cbind(ylim, uyl)
      # pPars <- pSamp + 1
      # pName <- c(pName, "uSd")
      # if (x$settings$dataType == "indivExtended") {
      #   vyl <- range(sapply(1:nsim, function(isim) {
      #     range(x$runs[[isim]]$vSd)
      #   }))
      #   ylim <- cbind(ylim, vyl)
      #   pPars <- pSamp + 1
      #   pName <- c(pName, "vSd")
      # }
    } else {
      pPars <- pSamp
    }
    if ("ylim" %in% names(argList)) {
      ylim <- argList$ylim
    }
    
    colnames(ylim) <- pName
    par(mfrow = c(ceiling(pPars / 2), 2), mar = c(4, 4, 3, 1))
    for (ipar in idSamp) {
      idp <- which(idSamp == ipar)
      plot(idkeep, x$runs[[1]]$theta[idkeep, ipar], type = 'l', 
           ylim = ylim[, idp], main = pName[idp], xlab = "Iteration", 
           ylab = "Parameter")
      for (ic in 2:nsim) {
        lines(idkeep, x$runs[[ic]]$theta[idkeep, ipar], col = ic)
      }
    }
    # if (grepl("indiv", x$settings$dataType)) {
      # plot(idkeep, x$runs[[1]]$uSd[idkeep], type = 'l',
      #      ylim = ylim[, ncol(ylim)], xlab = "Iteration", ylab = "Parameter",
      #      main = "uSd")
      # for (ic in 2:nsim) {
      #   lines(idkeep, x$runs[[ic]]$uSd[idkeep], col = ic)
      # }
      # if (x$settings$dataType == "indivExtended") {
      #   plot(idkeep, x$runs[[1]]$vSd[idkeep], type = 'l',
      #        ylim = ylim[, ncol(ylim)], xlab = "Iteration", ylab = "Parameter",
      #        main = "vSd")
      #   for (ic in 2:nsim) {
      #     lines(idkeep, x$runs[[ic]]$vSd[idkeep], col = ic)
      #   }
        
      # }
    # }
    
  } else if (type == "density") {
    idkeep <- seq(1, x$settings$niter, x$settings$thinning)
    pSamp <- x$params$pSamp
    idSamp <- x$params$idSamp
    pName <- x$params$thetaName[idSamp]
    if (x$settings$dataType == "indivExtended") {
      pSamp <- pSamp + 1
      idSamp <- c(idSamp, max(idSamp) + 1)
      pName <- c(pName, "vSd")
    }
    parMat <- x$theta[, idSamp]
    # if (grepl("indiv", x$settings$dataType)) {
    #   # parMat <- cbind(parMat, uSD = x$uSd)
    #   # idSamp <- c(idSamp, max(idSamp) + 1)
    #   # pSamp <- length(idSamp)
    #   # pName <- c(pName, "uSd")
    #   if (x$settings$dataType == "indivExtended") {
    #     parMat <- cbind(parMat, vSD = x$vSd)
    #     idSamp <- c(idSamp, max(idSamp) + 1)
    #     pSamp <- length(idSamp)
    #     pName <- c(pName, "vSd")
    #   }
    # }
    
    nsim <- x$settings$nsim
    pDens <- lapply(1:pSamp, function(ipar) {
      density(parMat[, ipar])
    })
    names(pDens) <- pName
    ylim <- sapply(1:pSamp, function(ipar) {
      c(0, max(pDens[[ipar]]$y))
    })
    if ("ylim" %in% names(argList)) {
      ylim <- argList$ylim
    }
    xlim <- sapply(1:pSamp, function(ipar) {
      range(pDens[[ipar]]$x)
    })
    if ("xlim" %in% names(argList)) {
      xlim <- argList$xlim
    }
    
    colnames(xlim) <- colnames(ylim) <- pName
    
    pCis <- lapply(pName, function(ipar) {
      which(pDens[[ipar]]$x >= x$coefficients[ipar, "Lower"] &
              pDens[[ipar]]$x <= x$coefficients[ipar, "Upper"])
    })
    par(mfrow = c(ceiling(pSamp / 2), 2), mar = c(4, 4, 3, 1))
    for (ipar in idSamp) {
      plot(xlim[, ipar], ylim[, ipar], col = NA, xlab = "Iteration", 
           ylab = "Posterior", main = pName[ipar])
      xx <- pDens[[ipar]]$x
      yy <- pDens[[ipar]]$y
      idcis <- pCis[[ipar]]
      ncis <- length(idcis)
      polygon(c(xx[idcis], rev(xx[idcis])), c(yy[idcis], rep(0, ncis)), 
              col = "orange", border = NA)
      lines(xx, yy, lwd = 2, col = 'dark red')
      
    }
  } else if (type == "fertility") {
    if ("showRealized" %in% names(argList)) {
      REALIZED <- TRUE
    } else {
      REALIZED <- FALSE
    }
    dat <- x$aggrData
    fertQuant <- x$fert
    idplf <- 1:min(apply(fertQuant, 2, function(ff) max(which(!is.na(ff)))))
    fertQuant <- fertQuant[idplf, ]
    xv <- x$x[idplf]
    alpha <- x$data$alpha
    ylim <- c(0, max(c(dat$Fertility, fertQuant[, "Upper"]), 
                     na.rm = TRUE))
    if ("ylim" %in% names(argList)) {
      ylim <- argList$ylim
    }
    xlim <- c(0, max(dat$Age) + alpha)
    # Offset for models with support x > 0:
    if (x$settings$model %in% c("gamma", "beta", "gammaMixture", "Hadwiger", 
                            "HadwigerMixture")) {
      xoffs <- 0.005
    } else {
      xoffs <- 0
    }
    
    
    par(mfrow = c(1, 1))
    plot(xlim, ylim, col = NA, xlab = "Age", ylab = "Fertility")
    polygon(alpha + c(xv, rev(xv)), 
            c(fertQuant[, "Lower"], rev(fertQuant[, "Upper"])),
            col = "orange", border = NA)
    lines(xv + alpha, fertQuant[, "Mean"], col = 'red', lwd = 2)
    lines(dat$Age + xoffs, dat$Fertility, type = 'b',
          lwd = 2)
    if ("RealizedFert" %in% colnames(dat) & REALIZED) {
      lines(dat$Age + xoffs, dat$RealizedFert, type = 'b',
            lwd = 2, col = 'grey80', lty = 2)
    }
  } else if (type == "predictive") {
    dat <- x$aggrData
    predQuant <- x$pred
    idna <- which(is.na(predQuant[, 1]))
    if (length(idna) > 0) {
      for (ina in idna) predQuant[ina, ] <- 0
    }
    alpha <- x$data$alpha
    ylim <- c(0, max(c(dat$nOffspring, predQuant[, "Upper"]), 
                     na.rm = TRUE))
    if ("ylim" %in% names(argList)) {
      ylim <- argList$ylim
    }
    xlim <- c(0, max(dat$Age))
    # Offset for models with support x > 0:
    if (x$settings$model %in% c("gamma", "beta", "gammaMixture", "Hadwiger", 
                                "HadwigerMixture")) {
      xoffs <- 0.005
    } else {
      xoffs <- 0
    }
    par(mfrow = c(1, 1))
    plot(xlim, ylim, col = NA, xlab = "Age", ylab = "Number of offspring")
    polygon(c(dat$Age, rev(dat$Age)), 
            c(predQuant[, "Lower"], rev(predQuant[, "Upper"])),
            col = "orange", border = NA)
    lines(dat$Age, predQuant[, "Mean"], col = 'red', lwd = 2)
    lines(dat$Age + xoffs, dat$nOffspring, col = 1, 
          lwd = 2, type = 'b')
  }
}

# Add min-max ages for unknown ages:
CalcMinMaxAge <- function(object, minAge, maxAge, unkAgeID = NULL) {
  sunID <- unique(object$indID)
  nind <- length(sunID)
  if (is.null(unkAgeID)) {
    idUnkAges <- 1:nind
  } else {
    if (!all(unkAgeID %in% sunID)) {
      stop("Some individual IDs in 'unkAgeID' not found in object.")
    }
    idUnkAges <- which(sunID %in% unkAgeID)
  }
  if (any(object$Age < minAge)) {
    stop("Some ages are earlier than minAge.")
  }
  if (any(object$Age > maxAge)) {
    stop("Some ages are earlier than minAge.")
  }
  MinAge <- MaxAge <- object$Age
  for (ii in idUnkAges) {
    idi <- which(object$indID == sunID[ii])
    ranAges <- range(object$Age[idi])
    deltaLow <- ranAges[1] - minAge
    deltaHigh <- maxAge - ranAges[2]
    MinAge[idi] <- object$Age[idi] - deltaLow
    MaxAge[idi] <- object$Age[idi] + deltaHigh
  }
  object <- data.frame(object, MinAge = MinAge, MaxAge = MaxAge)
  return(object)
}

# =========================== #
# B) INTERNAL FUNCTIONS: ==== 
# =========================== #
# ------------------------------------- #
# ---- Internal object management: ----
# ------------------------------------- #
# Algorithm object function:
.CreateAlgObj <- function(model, dataType, minAge, gestTime, formula, niter, 
                          burnin, nsim, thinning, UPDJUMP, jumpSD) {
  return(list(model = model, dataType = dataType, minAge = minAge, 
              gestTime = gestTime, formula = formula, niter = niter, 
              burnin = burnin, thinning = thinning, nsim = nsim, 
              minAge = minAge, UPDJUMP = UPDJUMP, jumpSD = jumpSD))
}

# Prepare data object:
.CreateDataObj <- function(object, algObj) {
  n <- nrow(object)
  # Aggretated data:
  if (algObj$dataType == "aggregated") {
    # Extract ages:
    alpha <- object$Age[which(object$Fertility > 0)[1]]
    idages <- which(object$Age >= alpha & object$nParents > 0)
    idpar <- which(object$nParents[which(object$Age >= alpha)] > 0)
    x <- object$Age[idages] - alpha
    y <- object$nOffspring
    dx <- diff(x[1:2])
    
    # Create data object:
    do <- list(data = object, idages = idages, alpha = alpha, x = x, y = y,
               xMax = max(x), n = n, dx = dx, idpar = idpar, AgeUpdate = FALSE,
               offsUpdate = FALSE)
    class(do) <- "baftaAggr"
  } else if (algObj$dataType == "indivSimple") {
    if (is.na(algObj$minAge)) {
      alpha <- min(object$Age)
    } else {
      alpha <- algObj$minAge
    }
    x <- object$Age - alpha
    dx <- min(diff(sort(unique(x))))

    if (any(x < 0)) {
      warning("Some ages occur before minAge.", 
              "These records have been excluded from the analysis.\n", 
              call. = FALSE)
      idincl <- which(x >= alpha)
      object <- object[idincl, ]
      x <- x[idincl]
    }
    y <- object$nOffspring
    rMat <- model.matrix(~ indID - 1, data = object)
    unID <- unique(object$indID)
    rMat <- rMat[, sprintf("indID%s", unID)]
    ni <- ncol(rMat)
    if (all(c("MinAge", "MaxAge") %in% colnames(object))) {
      idUpd <- which(object$Age != object$MinAge | 
                          object$Age != object$MaxAge)
      idAgeUpd <- which(unID %in% unique(object$indID[idUpd]))
      AgeUpdate <- TRUE
      nAgeUpd <- length(idAgeUpd)
      idFirstAge <- sapply(idAgeUpd, function(iid) {
        idi <- which(object$indID == unID[iid])
        return(idi[which(object$Age[idi] == min(object$Age[idi]))])
      })
    } else {
      idAgeUpd <- NA
      AgeUpdate <- FALSE
      nAgeUpd <- 0
      idFirstAge <- NA
    }
    if ("obsProp" %in% colnames(object)) {
      idOffsUpd <- which(object$obsProp < 1)
      offsUpdate <- TRUE
      nOffsUpd <- length(idOffsUpd)
      object$estOffspring <- object$nOffspring
      o <- object$nOffspring
    } else {
      idOffsUpd <- NA
      offsUpdate <- FALSE
      nOffsUpd <- 0
      o <- NA
    }
    do <- list(data = object, x = x, y = y, o = o, rMat = rMat, ni = ni, 
               xMax = max(x), n = n, dx = dx, alpha = alpha, 
               unID = unID, AgeUpdate = AgeUpdate, idAgeUpd = idAgeUpd,
               nAgeUpd = nAgeUpd, idFirstAge = idFirstAge, 
               offsUpdate = offsUpdate, idOffsUpd = idOffsUpd, 
               nOffsUpd = nOffsUpd)
    class(do) <- "baftaIndSimp"
  } else if (algObj$dataType == "indivExtended") {
    if (is.na(algObj$minAge)) {
      alpha <- floor(min(object$Age))
    } else {
      alpha <- algObj$minAge
    }
    if (is.na(algObj$gestTime)) {
      tau <- floor(min(object$IBI[which(object$First == 0)]))
    } else {
      tau <- algObj$gestTime
      if (any(object$IBI[which(object$First == 0)] < tau)) {
        warning("Specified gestation time is longer than the minimum IBI.\nGestation time was adjusted to min(IBI).")
        tau <- floor(min(object$IBI[which(object$First == 0)]))
      }
    }
    x <- object$Age - alpha
    # Adjust first age > 0 for certain models:
    if (algObj$model %in% c("gamma", "beta", "gammaMixture", "Hadwiger", 
                            "HadwigerMixture")) {
      x[which(x == 0)] <- 0.005
    }
    if (any(x < 0)) {
      warning("Some ages occur before minAge.", 
              "These records have been excluded from the analysis.\n", 
              call. = FALSE)
      idincl <- which(x >= alpha)
      object <- object[idincl, ]
      x <- x[idincl]
    }
    y <- object$nOffspring
    rMat <- model.matrix(~ indID - 1, data = object)
    unID <- unique(object$indID)
    rMat <- rMat[, sprintf("indID%s", unID)]
    ni <- ncol(rMat)
    z <- object$IBI
    nIBI <- c(rep(1, n) %*% (rMat * (1 - object$First)))
    idv <- which(nIBI > 0)
    idIBI <- rep(0, ni)
    idIBI[idv] <- 1
    niv <- sum(idIBI)
    z[which(object$First == 0)] <- z[which(object$First == 0)] - tau
    w <- object$Age - alpha
    if (all(c("MinAge", "MaxAge") %in% colnames(object))) {
      idUpd <- which(object$Age != object$MinAge | 
                       object$Age != object$MaxAge)
      idAgeUpd <- which(unID %in% unique(object$indID[idUpd]))
      AgeUpdate <- TRUE
      nAgeUpd <- length(idAgeUpd)
      idFirstAge <- sapply(idAgeUpd, function(iid) {
        idi <- which(object$indID == unID[iid])
        return(idi[which(object$Age[idi] == min(object$Age[idi]))])
      })
      
    } else {
      idAgeUpd <- NA
      AgeUpdate <- FALSE
      nAgeUpd <- 0
      idFirstAge <- NA
    }
    
    do <- list(data = object, x = x, y = y, rMat = rMat, ni = ni, 
               xMax = max(x), z = z, w = w, n = n, alpha = alpha, tau = tau,
               indv = idIBI, idv = idv, niv = niv, unID = unID, 
               AgeUpdate = AgeUpdate, idAgeUpd = idAgeUpd, nAgeUpd = nAgeUpd,
               idFirstAge = idFirstAge, offsUpdate = FALSE)
    class(do) <- "baftaIndExt"
  }
  return(do)
}

# Fertility parameters:
.SetDefaultBeta <- function(algObj, dataObj) {
  model <- algObj$model
  if (model == "quadratic") {
    nBe <- 3
    lowBe <- rep(0, 3)
    uppBe <- rep(Inf, 3)
    startBe <- c(1, 0.025, 10)
    priorMeanBe <- c(1, 0.025, 10)
    idSamp <- 1:3
  } else if (model == "PeristeraKostaki") {
    nBe <- 4
    lowBe <- c(0, 0, 0, 0)
    uppBe <- c(rep(Inf, 3), dataObj$xMax)
    startBe <- c(1, 5, 15, dataObj$xMax / 2)
    priorMeanBe <- c(1, 5, 15, dataObj$xMax / 2)
    idSamp <- 1:4
  } else if (model == "ColcheroMuller") {
    nBe <- 4
    lowBe <- c(0, 0, 0, -Inf)
    uppBe <- rep(Inf, 4)
    startBe <- c(1, 0.005, 0.0001, -1)
    priorMeanBe <- c(1, 0.005, 0.0001, -1)
    idSamp <- 1:4
  } else if (model == "Hadwiger") {
    nBe <- 3
    lowBe <- c(0, 0, 0)
    uppBe <- rep(Inf, 3)
    startBe <- c(15, 1, 20)
    priorMeanBe <- c(15, 1, 20)
    idSamp <- 1:3
  } else if (model == "gamma") {
    nBe <- 3
    lowBe <- c(0, 0, 0)
    uppBe <- rep(Inf, 3)
    startBe <- c(2, 2, 0.1)
    priorMeanBe <- c(2, 2, 0.1)
    idSamp <- 1:3
  } else if (model == "beta") {
    nBe <- 5
    lowBe <- rep(0, 5)
    uppBe <- rep(Inf, 5)
    xMax <- ceiling(dataObj$xMax)
    if (algObj$dataType %in% c("aggregated", "indivSimple")) xMax <- xMax + 1
    startBe <- c(10, 1.5, 2.5, 0, xMax)
    priorMeanBe <- c(10, 1.5, 2.5, 0, xMax)
    idSamp <- 1:3
  } else if (model == "skewNormal") {
    nBe <- 4
    lowBe <- rep(0, 4)
    uppBe <- rep(Inf, 4)
    startBe <- c(5, 20, 5, 5)
    priorMeanBe <- c(5, 20, 5, 5)
    idSamp <- 1:4
  } else if (model == "gammaMixture") {
    nBe <- 6
    lowBe <- rep(0, 6)
    uppBe <- c(Inf, 1, rep(Inf, 4))
    startBe <- c(20, 0.2, 1.5, 0.25, 4, 0.2)
    priorMeanBe <- c(20, 0.2, 1.5, 0.25, 4, 0.2)
    idSamp <- 1:6
  } else if (model == "HadwigerMixture") {
    nBe <- 6
    lowBe <- rep(0, 6)
    uppBe <- c(Inf, 1, rep(Inf, 4))
    startBe <- c(1, 1, 10, 1.25, 30)
    priorMeanBe <- c(1, 0.2, 1, 10, 1.25, 30)
    idSamp <- 1:6
  } else if (model == "skewSymmetric") {
    nBe <- 5
    lowBe <- c(rep(0, 4), -Inf)
    uppBe <- rep(Inf, 5)
    startBe <- c(5, 10, 20, 1, -0.5)
    priorMeanBe <- c(5, 10, 20, 1, -0.5)
    idSamp <- 1:5
  } else if (model == "skewLogistic") {
    nBe <- 5
    lowBe <- c(rep(0, 4), -Inf)
    uppBe <- rep(Inf, 5)
    startBe <- c(5, 10, 20, 1, -0.5)
    priorMeanBe <- c(5, 10, 20, 1, -0.5)
    idSamp <- 1:5
  }
  priorSdBe <- rep(5, nBe)
  nameBe <- sprintf("b%s", 1:nBe - 1)
  names(startBe) <- names(priorMeanBe) <- names(priorSdBe) <- nameBe
  defaultBeta  <- list(beta = startBe, priorMean = priorMeanBe, 
                       priorSD = priorSdBe, p = nBe, name = nameBe, low = lowBe, 
                       upp = uppBe, idSamp = idSamp, pSamp = length(idSamp))
  attr(defaultBeta, "model") = model
  return(defaultBeta)
}

# Define parameter object:
.BuildParObj <- function(algObj, dataObj) {
  # NA on random effects (to be replaced based on need):
  uSd <- NA
  uPrior1 <- NA
  uPrior2 <- NA
  vSd <- NA
  vPrior1 <- NA
  vPrior2 <- NA
  
  # Basic beta parameters:
  defBet <- .SetDefaultBeta(algObj = algObj, dataObj = dataObj)
  # ---------------------------------- #
  # Negative binomial instead of log-normal random effects
  # if (algObj$dataType %in% c("indivSimple", "indivExtended")) {
  #   uSd <- 0.5
  #   uPrior1 <- 0.01
  #   uPrior2 <- 0.1
  # } 
  theta <- c(defBet$beta, gamma = 0.1)
  thetaPriorMean <- c(defBet$priorMean, gamma = 1)
  thetaPriorSD <- c(defBet$priorSD, gamma = 5)
  thetaLower <- c(defBet$low, gamma = 0)
  thetaUpper <- c(defBet$upp, gamma = Inf)
  idSamp <- c(defBet$idSamp, defBet$p + 1)
  if (dataObj$AgeUpdate) {
    agePriorSD <- 2
  } else {
    agePriorSD <- NA
  }
  if (algObj$dataType == "indivSimple" & dataObj$offsUpdate) {
    offsPrior <- c(-5, 10)
  } else {
    offsPrior <- NA
  }
  if (algObj$dataType == "indivExtended") {
    etaStart <- 1 / mean(dataObj$z[which(dataObj$data$First == 0)])
    kappaStart <- 1 / mean(dataObj$w[which(dataObj$data$First == 1)])
    theta <- c(theta, eta = etaStart, kappa = kappaStart)
    thetaPriorMean <- c(thetaPriorMean, eta = 1, kappa = 1)
    thetaPriorSD <- c(thetaPriorSD, eta = 1, kappa = 1)
    thetaLower <- c(thetaLower, eta = 0, kappa = 0)
    thetaUpper <- c(thetaUpper, eta = Inf, kappa = Inf)
    idSamp <- c(defBet$idSamp, defBet$p + c(1:3))

    # Random effects for IBI:
    vSd <- 0.5
    vPrior1 <- 1
    vPrior2 <- 0.1
    
  } 
  if (dataObj$AgeUpdate) {
    idSdAge <- which(dataObj$data$indID %in% dataObj$unID[dataObj$idAgeUpd] & 
                       1:dataObj$n %in% dataObj$idFirstAge)
    agePriorSD <- mean(c(dataObj$data$MaxAge - 
                           dataObj$data$MinAge)[idSdAge] / 2) / 
      qnorm(0.999)
    # agePriorSD <- 2
  } else {
    agePriorSD <- NA
  }
  parList <- list(thetaStart = theta, thetaPriorMean = thetaPriorMean, 
                  thetaPriorSD = thetaPriorSD, thetaLower = thetaLower,
                  thetaUpper = thetaUpper, thetaName = names(theta),
                  p = length(theta), idSamp = idSamp, 
                  pSamp = length(idSamp), uSd = uSd, uPrior1 = uPrior1, 
                  uPrior2 = uPrior2, vSd = vSd, vPrior1 = vPrior1, 
                  vPrior2 = vPrior2, agePriorSD = agePriorSD, 
                  offsPrior = offsPrior)
  return(parList)
}

# Extract parameters specified by the user:
.CreateUserPar <- function(argList, argNames, parObj, dataObj) {
  # Verify if starting parameters are specified:
  if ("thetaStart" %in% argNames) {
    if (length(argList$thetaStart) != parObj$p) {
      stop(sprintf("Length of 'thetaStart' argument should be %s.\n", 
                   parObj$p), call. = FALSE)
    } else {
      parObj$thetaStart <- argList$thetaStart
      names(parObj$thetaStart) <- parObj$thetaName
    }
  }
  
  # Verify if mean priors are specified:
  if ("thetaPriorMean" %in% argNames) {
    if (length(argList$thetaPriorMean) != parObj$p) {
      stop(sprintf("Length of 'thetaPriorMean' argument should be %s.\n", 
                   parObj$p), call. = FALSE)
    } else {
      parObj$thetaPriorMean <- argList$thetaPriorMean
      names(parObj$thetaPriorMean) <- parObj$thetaName
    }
  }
  
  # Verify if SD priors are specified:
  if ("thetaPriorSD" %in% argNames) {
    if (length(argList$thetaPriorSD) != parObj$p) {
      stop(sprintf("Length of 'thetaPriorSD' argument should be %s.\n", 
                   parObj$p), call. = FALSE)
    } else {
      parObj$thetaPriorSD <- argList$thetaPriorSD
      names(parObj$thetaPriorSD) <- parObj$thetaName
    }
  }
  
  if ("agePriorSD" %in% argNames) {
    if (length(argList$agePriorSD) == 1 | 
        length(argList$agePriorSD) == nrow(dataObj$data)) {
      parObj$agePriorSD <- argList$agePriorSD
    } else {
      stop(sprintf("Argument 'agePriorSD' should be of length one or equal to the number of rows in the dataset (i.e., %s).\n", 
                   nrow(dataObj$data)), call. = FALSE)
    }
  }
  return(parObj)
}

# ------------------------ #
# ---- Distributions: ----
# ------------------------ #
# Truncated normal:
.rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  ru <- runif(n, Flow, Fup)
  rx <- .qtnorm(ru, mean, sd, lower = lower, upper = upper)
  return(rx)
}

.dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

.ptnorm <- function(q, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  p <- (pnorm(q, mean, sd) - pnorm(lower, mean, sd)) / 
    (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
  if (log) {
    p <- log(p)
  }
  return(p)
}

.qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  p2 <- (p) * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) + 
    pnorm(lower, mean, sd)
  q <- qnorm(p2, mean, sd)
  return(q)
}

# --------------------------- #
# ---- Fertility models: ----
# --------------------------- #
# For matrix of parameters:
.DefineFertilityMatrix <- function(modelFert = "quadratic") {
  if (modelFert == "quadratic") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * exp(-beta[, "b1"] * (x - beta[, "b2"])^2)
      return(fert)
    }
  } else if (modelFert == "PeristeraKostaki") {
    fertfun <- function(beta, x) {
      be1 <- x * 0 + beta[, "b1a"]
      be1[which(x > beta[, "b2"])] <- beta[, "b1b"]
      fert <- beta[, "b0"] * exp(-((x - beta[, "b2"]) / be1)^2)
      return(fert)
    }
  } else if (modelFert == "ColcheroMuller") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * exp(-beta[, "b1"] * (x - beta[, "b2"])^2 +
                                   beta[, "b3"] * 1/(x + 1))
      return(fert)
    }
  } else if (modelFert == "Hadwiger") {
    fertfun <- function(beta, x) {
      fert <- (beta[, "b0"] * beta[, "b1"])/beta[, "b2"] * 
        (beta[, "b2"]/x)^(3/2) * exp(-beta[, "b1"]^2 * 
                                       (beta[, "b2"] / x + 
                                          x / beta[, "b2"] - 2))
      return(fert)
    }
  } else if (modelFert == "gamma") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * dgamma(x, shape = beta[, "b1"], 
                                    rate = beta[, "b2"])
      return(fert)
    }
  } else if (modelFert == "beta") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * ((x - beta[, "b3"])^(beta[, "b1"] - 1) * 
                                (beta[, "b4"] - x)^(beta[, "b2"] - 1)) / 
        ((beta[, "b4"] - beta[, "b3"])^(beta[, "b1"] + beta[, "b2"] - 1) * 
           beta(beta[, "b1"], beta[, "b2"]))
      return(fert)
    }
  } else if (modelFert == "skewNormal") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * 2 * 1/beta[, "b1"] * 
        dnorm((x - beta[, "b2"]) / beta[, "b1"]) * 
        pnorm(beta[, "b3"] * ((x - beta[, "b2"]) / beta[, "b1"]))
      return(fert)
    }
  } else if (modelFert == "gammaMixture") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * (beta[, "b1"] * 
                                dgamma(x, shape = beta[, "b2"], 
                                       rate = beta[, "b3"]) +
                                (1 - beta[, "b1"]) * 
                                dgamma(x, shape = beta[, "b4"], 
                                       rate = beta[, "b5"]))
      return(fert)
    }
  } else if (modelFert == "HadwigerMixture") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * 
        (beta[, "b1"] * 
           (beta[, "b2"] / beta[, "b3"] * (beta[, "b3"]/x)^(3/2) * 
              exp(-beta[, "b2"]^2 * (beta[, "b3"] / x + 
                                       x / beta[, "b3"] - 2))) +
           (1 - beta[, "b1"]) * 
           (beta[, "b4"] / beta[, "b5"] * (beta[, "b5"]/x)^(3/2) * 
              exp(-beta[, "b4"]^2 * (beta[, "b5"] / x + 
                                       x / beta[, "b5"] - 2))))
      return(fert)
    }
  } else if (modelFert == "skewSymmetric") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * 2 * 1/beta[, "b1"] * 
        dnorm((x - beta[, "b2"]) / beta[, "b1"]) * 
        pnorm(beta[, "b3"] * ((x - beta[, "b2"]) / beta[, "b1"]) + 
                beta[, "b4"] * ((x - beta[, "b2"]) / beta[, "b1"])^3)
      return(fert)
    }
    
  } else if (modelFert == "skewLogistic") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * 2 * 1 / beta[, "b1"] * 
        ((exp(-(x - beta[, "b2"]) / beta[, "b1"])) / 
           ((1 + exp(-(x - beta[, "b2"]) / beta[, "b1"]))^2 *
              (1 + exp(-beta[, "b3"] * (x - beta[, "b2"]) / beta[, "b1"] - 
                         beta[, "b4"] * ((x - beta[, "b2"]) / 
                                           beta[, "b1"])^3))))
      return(fert)
    }
    
  }
  return(fertfun)
}

# For single vector of parameters:
.DefineFertilityNumeric <- function(modelFert = "quadratic") {
  if (modelFert == "quadratic") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2)
      return(fert)
    }
  } else if (modelFert == "PeristeraKostaki") {
    fertfun <- function(beta, x) {
      be1 <- x * 0 + beta["b1"]
      be1[which(x > beta["b3"])] <- beta["b2"]
      fert <- beta["b0"] * exp(-((x - beta["b3"]) / be1)^2)
      return(fert)
    }
  } else if (modelFert == "ColcheroMuller") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2 +
                                 beta["b3"] * 1/(x + 1))
      return(fert)
    }
  } else if (modelFert == "Hadwiger") {
    fertfun <- function(beta, x) {
      fert <- (beta["b0"] * beta["b1"])/beta["b2"] * 
        (beta["b2"]/x)^(3/2) * exp(-beta["b1"]^2 * (beta["b2"] / x + 
                                                      x / beta["b2"] - 2))
      return(fert)
    }
  } else if (modelFert == "gamma") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * dgamma(x, shape = beta["b1"], rate = beta["b2"])
      return(fert)
    }
  } else if (modelFert == "beta") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * ((x - beta["b3"])^(beta["b1"] - 1) * 
                              (beta["b4"] - x)^(beta["b2"] - 1)) / 
        ((beta["b4"] - beta["b3"])^(beta["b1"] + beta["b2"] - 1) * 
           beta(beta["b1"], beta["b2"]))
      return(fert)
    }
  } else if (modelFert == "skewNormal") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * 2 * 1/beta["b1"] * 
        dnorm((x - beta["b2"]) / beta["b1"]) * 
        pnorm(beta["b3"] * ((x - beta["b2"]) / beta["b1"]))
      return(fert)
    }
  } else if (modelFert == "gammaMixture") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * (beta["b1"] * 
                              dgamma(x, shape = beta["b2"], 
                                     rate = beta["b3"]) +
                              (1 - beta["b1"]) * dgamma(x, shape = beta["b4"],
                                                        rate = beta["b5"]))
      return(fert)
    }
  } else if (modelFert == "HadwigerMixture") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * 
        (beta["b1"] * (beta["b2"]/beta["b3"] * (beta["b3"]/x)^(3/2) * 
                         exp(-beta["b2"]^2 * 
                               (beta["b3"] / x + x / beta["b3"] - 2))) +
           (1 - beta["b1"]) * 
           (beta["b4"]/beta["b5"] * (beta["b5"]/x)^(3/2) * 
              exp(-beta["b4"]^2 * (beta["b5"] / x + x / beta["b5"] - 2))))
      return(fert)
    }
  } else if (modelFert == "skewSymmetric") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * 2 * 1/beta["b1"] * 
        dnorm((x - beta["b2"]) / beta["b1"]) * 
        pnorm(beta["b3"] * ((x - beta["b2"]) / beta["b1"]) + 
                beta["b4"] * ((x - beta["b2"]) / beta["b1"])^3)
      return(fert)
    }
  } else if (modelFert == "skewLogistic") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * 2 * 1 / beta["b1"] * 
        ((exp(-(x - beta["b2"]) / beta["b1"])) / 
           ((1 + exp(-(x - beta["b2"]) / beta["b1"]))^2 *
              (1 + exp(-beta["b3"] * (x - beta["b2"]) / beta["b1"] - 
                         beta["b4"] * ((x - beta["b2"]) / beta["b1"])^3))))
      return(fert)
    }
  }
  return(fertfun)
}

# ------------------------------------------ #
# ---- Likelihood, posterior, sampling: ----
# ------------------------------------------ #
# Likelihood function:
.CalcLikeFert <- function(dataObj, ...) UseMethod(".CalcLikeFert")

.CalcLikeFert.baftaAggr <- function(dataObj, pars, FertFun, 
                                    FertFun.numeric) {
  fert <- FertFun(beta = pars$theta, x = dataObj$x + dataObj$dx / 2)
  # ----------------------------------- #
  # Use negative binomial (2023-11-14):
  
  # lk <- dpois(dataObj$data$nOffspring[dataObj$idages], 
  #             lambda = dataObj$data$nParents[dataObj$idages] * fert, 
  #             log = TRUE)
  rho <- pars$theta["gamma"] / (fert + pars$theta["gamma"])
  lk <- dnbinom(x = dataObj$data$nOffspring[dataObj$idages],
                size = dataObj$data$nParents[dataObj$idages] *
                  pars$theta["gamma"], prob = rho, log = TRUE)
  # ---------------------------------- #
  return(lk)
}

.CalcLikeFert.baftaIndSimp <- function(dataObj, pars, FertFun, 
                                       FertFun.numeric) {
  fert <- FertFun(beta = pars$theta, x = dataObj$x + dataObj$dx / 2)
  rho <- pars$theta["gamma"] / (fert + pars$theta["gamma"])
  lk <- dnbinom(x = dataObj$y, 
                size = pars$theta["gamma"], prob = rho, log = TRUE)
  return(lk)
}

# Likelihood function:
.CalcLikeFert.baftaIndExt <- function(dataObj, pars, FertFun, 
                                      FertFun.numeric) {
  fert <- FertFun(beta = pars$theta, x = dataObj$x)
  rho <- pars$theta["gamma"] / (fert + pars$theta["gamma"])
  lk <- dnbinom(x = dataObj$y, size = pars$theta["gamma"],
                prob = rho, log = TRUE) +
    dexp(dataObj$z, rate = pars$theta["eta"] * 
           exp(c(dataObj$rMat %*% pars$v)), log = TRUE) * 
    (1 - dataObj$data$First) +
    dexp(dataObj$w, rate = pars$theta["kappa"], log = TRUE) * 
    dataObj$data$First
  return(lk)
}

# Posterior for Beta:
.CalcPostTheta <- function(pars, like, parObj) {
  post <- sum(like) + sum(.dtnorm(pars$theta, mean = parObj$thetaPriorMean, 
                                 sd = parObj$thetaPriorSD, 
                                 lower = parObj$thetaLower, 
                                 upper = parObj$thetaUpper, log = TRUE))
  return(post)
}

# Metropolis-Hastings ratio:
.CalcMHratio <- function(parsNow, parsNew, jumpSD, parObj, ip) {
  idj <- which(parObj$idSamp == ip)
  MHr <- .dtnorm(parsNow$theta[ip], mean = parsNew$theta[ip], sd = jumpSD[idj],
                lower = parObj$thetaLower[ip], upper = parObj$thetaUpper[ip], 
                log = TRUE) -
    .dtnorm(parsNew$theta[ip], mean = parsNow$theta[ip], sd = jumpSD[idj],
           lower = parObj$thetaLower[ip], upper = parObj$thetaUpper[ip], 
           log = TRUE)
  return(MHr)
}

# Posterior for random effects on fertility:
.CalcPostRandEffU <- function(dataObj, pars, like) {
  postu <- c(t(like * dataObj$rMat) %*% rep(1, dataObj$n)) + 
    dnorm(pars$u, mean = 0, sd = pars$uSd, log = TRUE)
  return(postu)
}

# Sampling fertility random effects variance:
.SampleUSig <- function(dataObj, pars, parObj) {
  u1 <- parObj$uPrior1 + dataObj$ni / 2
  u2 <- parObj$uPrior2 + 0.5 * sum(pars$u^2) 
  sig <- sqrt(1 / rgamma(1, u1, u2))
  return(sig)
}

# Posterior for random effects on IBI:
.CalcPostRandEffV <- function(dataObj, pars, like) {
  postv <- c(t(like * dataObj$rMat * (1 - dataObj$data$First)) 
             %*% rep(1, dataObj$n))[dataObj$idv] + 
    dnorm(pars$v[dataObj$idv], mean = 0, sd = pars$vSd, log = TRUE)
  return(postv)
}

# Sampling IBI random effects variance:
.SampleVSig <- function(dataObj, pars, parObj) {
  v1 <- parObj$vPrior1 + dataObj$niv / 2
  v2 <- parObj$vPrior2 + 0.5 * sum(pars$v[dataObj$idv]^2) 
  sig <- sqrt(1 / rgamma(1, v1, v2))
  return(sig)
}

# Sample unknown ages:
.SampleAges <- function(dataObj, ...) UseMethod(".SampleAges")

.SampleAges.baftaAggr <- function(dataObj, algObj) {
  return(dataObj)
}

.SampleAges.baftaIndSimp <- function(dataObj, algObj) {
  dataObjUpd <- dataObj
  dataNew <- dataObj$data
  for (iid in dataObj$idAgeUpd) {
    id <- which(dataNew$indID == dataObj$unID[iid])
    updRan <- c(dataNew$MinAge[id[1]], dataNew$MaxAge[id[1]]) - 
      dataNew$Age[id[1]]
    deltaAge <- round(.rtnorm(n = 1, mean = 0, sd = 1, lower = updRan[1],
                              upper = updRan[2]))
    dataNew$Age[id] <- dataObj$data$Age[id] + deltaAge
  }
  
  x <- dataNew$Age - dataObj$alpha
  dx <- min(diff(sort(unique(x))))
  # Adjust first age > 0 for certain models:
  if (algObj$model %in% c("gamma", "beta", "gammaMixture", "Hadwiger", 
                          "HadwigerMixture")) {
    x[which(x == 0)] <- 0.005
  }
  dataObjUpd$data <- dataNew
  dataObjUpd$x <- x
  return(dataObjUpd)
}

.SampleAges.baftaIndExt <- function(dataObj, algObj) {
  dataObjUpd <- dataObj
  dataNew <- dataObj$data
  for (iid in dataObj$idAgeUpd) {
    id <- which(dataNew$indID == dataObj$unID[iid])
    updRan <- c(dataNew$MinAge[id[1]], dataNew$MaxAge[id[1]]) - 
      dataNew$Age[id[1]]
    deltaAge <- .rtnorm(n = 1, mean = 0, sd = 1, lower = updRan[1],
                              upper = updRan[2])
    dataNew$Age[id] <- dataObj$data$Age[id] + deltaAge
  }
  
  x <- dataNew$Age - dataObj$alpha
  dx <- min(diff(sort(unique(x))))
  # Adjust first age > 0 for certain models:
  if (algObj$model %in% c("gamma", "beta", "gammaMixture", "Hadwiger", 
                          "HadwigerMixture")) {
    x[which(x == 0)] <- 0.005
  }
  dataObjUpd$data <- dataNew
  dataObjUpd$x <- x
  dataObjUpd$w <- dataNew$Age - dataObj$alpha
  
  return(dataObjUpd)
}

.CalcHastingsRatioAges <- function(dataObjNew, dataObjNow, algObj) {
  hrat <- rep(0, dataObj$nAgeUpd)
  ide <- 0
  for (iid in dataObj$idAgeUpd) {
    ide <- ide + 1
    id <- which(dataNew$indID == dataObj$unID[iid])
    deltaAge <- dataObjNew$x[id[1]] - dataObjNow$x[id[1]]
    updRan <- c(dataNew$MinAge[id[1]], dataNew$MaxAge[id[1]]) - 
      dataNew$Age[id[1]]
    hrat[ide] <- .dtnorm(deltaAge, mean = 0, sd = 1, lower = updRan[1],
                        upper = updRan[2]) / 
      .dtnorm(0, mean = deltaAge, sd = 1, lower = updRan[1],
              upper = updRan[2])
  }
  return(hrat)
}

# Prior of partially observed offspring:
.CalcPostOffs <- function(dataObj, parObj, like, FertFun, FertFun.numeric) {
  ages <- dataObj$x[dataObj$idOffsUpd]
  mu <- FertFun(parObj$thetaPriorMean, ages)
  po <- dataObj$data$obsProp[dataObj$idOffsUpd]
  o <- dataObj$o[dataObj$idOffsUpd]
  y <- dataObj$y[dataObj$idOffsUpd]
  priorOffs <- dnbinom(x = y - o, size = parObj$thetaPriorMean["gamma"],
                       mu = mu * (1 - po), log = TRUE)
  postOffs <- like[dataObj$idOffsUpd] + priorOffs
  return(postOffs)
}

# Sample unknown ages:
.SampleOffspring <- function(dataObj, ...) UseMethod(".SampleOffspring")

.SampleOffspring.baftaAggr <- function(dataObj, ...) {
  return(dataObj)
}

.SampleOffspring.baftaIndSimp <- function(dataObj, algObj) {
  dataObjUpd <- dataObj
  oObs <- dataObj$o
  yNew <- dataObj$y
  yNow <- dataObj$y
  yNew[dataObj$idOffsUpd] <- round(.rtnorm(n = dataObj$nOffsUpd,
                                           mean = yNow[dataObj$idOffsUpd],
                                           sd = 1, 
                                           lower = oObs[dataObj$idOffsUpd]))
  dataObjUpd$data$estOffspring <- yNew
  dataObjUpd$y <- yNew
  return(dataObjUpd)
}

.SampleOffspring.baftaIndExt <- function(dataObj, ...) {
  return(dataObj)
}

# Function to update jumps in MCMC:
.UpdateJumps <- function(updMat, jumps, iter, iterUpd = 50, 
                         updTarg = 0.25) {
  updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd  
  updRate[updRate == 0] <- 1e-2
  jumps <- jumps * updRate / updTarg
  return(jumps)
}


# =============== #
# ==== MCMC: ====
# =============== #
.RunMCMC <- function(sim, dataObj, parObj, niter, algObj, FertFun, 
                     FertFun.numeric, FertFun.matrix,
                     jumpSD = NULL, UPDJUMP = TRUE) {
  # Restart the random seed for each core:
  if (sim > 1) {
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
  }
  
  # Create sampling parameter object:
  parsNow <- list(theta = parObj$thetaStart)
  
  # Fertility random effects parameter:
  # if (grepl("indiv", algObj$dataType)) {
  #   RANDEFFU <- TRUE
  #   parsNow$u <- rep(0, dataObj$ni)
  #   parsNow$uSd <- 1
  # } else {
  #   RANDEFFU <- FALSE
  # }
  RANDEFFU <- FALSE
  
  # IBI random effects parameter:
  if (algObj$dataType == "indivExtended") {
    RANDEFFV <- TRUE
    parsNow$v <- rep(0, dataObj$ni)
    parsNow$vSd <- 1
  } else {
    RANDEFFV <- FALSE
  }
  
  # Start dataObj:
  dataObjNow <- dataObj
  
  # Calculate likelihood and posteriors:
  likeNow <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNow, 
                           FertFun = FertFun, 
                           FertFun.numeric = FertFun.numeric)
  postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                            parObj = parObj)
  if (is.na(postNow) | postNow == -Inf | postNow == Inf) {
    postNA <- TRUE
  } else {
    postNA <- FALSE
  }
  
  # Jitter parameters:
  while(postNA) {
    parsNow$theta[parObj$idSamp] <- 
      .rtnorm(n = 1, mean = parObj$thetaStart[parObj$idSamp], 
              sd = rep(0.1, length(parObj$idSamp)), 
              lower = parObj$thetaLower[parObj$idSamp],
              upper = parObj$thetaUpper[parObj$idSamp])
    
    # Calculate likelihood and posteriors:
    likeNow <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNow, 
                             FertFun = FertFun, 
                             FertFun.numeric = FertFun.numeric)
    postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                              parObj = parObj)
    if (is.na(postNow) | postNow == -Inf | postNow == Inf) {
      postNA <- TRUE
    } else {
      postNA <- FALSE
    }
  }
  
  # Calculate likelihood and posteriors:
  likeNow <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNow, 
                           FertFun = FertFun, 
                           FertFun.numeric = FertFun.numeric)
  postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                            parObj = parObj)
  
  # If ages need to be updated:
  if (dataObj$AgeUpdate) {
    likeIndNow <- c(t(likeNow) %*% dataObj$rMat)
    
    # Prior for proposed ages:
    priorNow <- .dtnorm(x = dataObjNow$data$Age, 
                        mean = dataObj$data$Age, sd = parObj$agePriorSD, 
                        lower = dataObj$data$MinAge, 
                        upper = dataObj$data$MaxAge, 
                        log = TRUE)
    priorNow[which(priorNow == Inf)] <- 0
    priorIndNow <- c(t(priorNow) %*% dataObj$rMat)
    
    # Posterior for proposed ages:
    postIndNow <- likeIndNow + priorIndNow
    
  }
  
  # If offspring number needs to be updated (partial observations):
  if (dataObj$offsUpdate) {
    postOffsNow <- .CalcPostOffs(dataObj = dataObjNow, parObj = parObj, 
                                 like = likeNow, FertFun = FertFun, 
                                 FertFun.numeric = FertFun.numeric)
  }
  
  
  if (RANDEFFU) {
    postUNow <- .CalcPostRandEffU(dataObj = dataObjNow, pars = parsNow,
                                  like = likeNow)
  }
  
  if (RANDEFFV) {
    postVNow <- .CalcPostRandEffV(dataObj = dataObjNow, pars = parsNow, 
                                  like = likeNow)
  }
  
  # Prepare jumps and outputs:
  if (UPDJUMP) {
    jumpSD <- rep(0.1, parObj$pSamp)
    jumpOut <- matrix(jumpSD, nrow = 1, ncol = parObj$pSamp,
                      dimnames = list(NULL, parObj$thetaName[parObj$idSamp]))
    updJumpIter <- seq(50, niter, 50)
    updateMat <- matrix(0, nrow = niter, ncol = parObj$pSamp,
                        dimnames = list(NULL, parObj$thetaName[parObj$idSamp]))
    parOut <- NA
    likePostOut <- NA
    uOut <- NA
    uSdOut <- NA
    vOut <- NA
    vSdOut <- NA
    ageOut <- NA
    offsOut <- NA
    
  } else {
    pAll <- parObj$p
    parNames <- parObj$thetaName
    if (inherits(dataObj, "baftaIndExt")) {
      pAll <- pAll + 1
      parNames <- c(parNames, "vSd")
    }
    parOut <- matrix(NA, nrow = niter, ncol = pAll,
                     dimnames = list(NULL, parNames))
    likePostOut <- matrix(NA, nrow = niter, ncol = 2,
                          dimnames = list(NULL, c("Likelihood", "Posterior")))
    parOut[1, parObj$thetaName] <- parsNow$theta
    if (inherits(dataObj, "baftaIndExt")) {
      parOut[1, "vSd"] <- parsNow$vSd
    }
    likePostOut[1, ] <- c(sum(likeNow), postNow)
    if (dataObj$AgeUpdate) {
      ageOut <- matrix(NA, nrow = niter, ncol = dataObj$nAgeUpd)
      ageOut[1, ] <- dataObjNow$data$Age[dataObj$idFirstAge]
    } else {
      ageOut <- NA
    }
    
    # If offspring number needs to be updated (partial observations):
    if (dataObj$offsUpdate) {
      offsOut <- matrix(NA, nrow = niter, ncol = dataObj$nOffsUpd)
      offsOut[1, ] <- dataObjNow$y[dataObj$idOffsUpd]
    } else {
      offsOut <- NA
    }
    
    if (RANDEFFU) {
      uOut <- matrix(NA, nrow = niter, ncol = dataObj$ni)
      uOut[1, ] <- parsNow$u
      uSdOut <- rep(NA, niter)
      uSdOut[1] <- parsNow$uSd
    } else {
      uOut <- NA
      uSdOut <- NA
    }
    if (RANDEFFV) {
      vOut <- matrix(NA, nrow = niter, ncol = dataObj$ni)
      vOut[1, ] <- parsNow$v
      # vSdOut <- rep(NA, niter)
      # vSdOut[1] <- parsNow$vSd
    } else {
      vOut <- NA
      # vSdOut <- NA
    }
  }
  
  # ----------- #
  # Start MCMC:
  # ----------- #
  for (iter in 2:niter) {
    for (ip in parObj$idSamp) {
      idj <- which(parObj$idSamp == ip)
      parsNew <- parsNow
      parsNew$theta[ip] <- .rtnorm(n = 1, mean = parsNow$theta[ip], 
                                   sd = jumpSD[idj],
                                   lower = parObj$thetaLower[ip], 
                                   upper = parObj$thetaUpper[ip])
      likeNew <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNew, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postNew <- .CalcPostTheta(pars = parsNew, like = likeNew, 
                                parObj = parObj)
      mhRatio <- .CalcMHratio(parsNow = parsNow, parsNew = parsNew, 
                              jumpSD = jumpSD, parObj = parObj, ip = ip)
      postRatio <- exp(postNew - postNow + mhRatio)
      
      if (!is.na(postRatio)) {
        if (postRatio > runif(n = 1)) {
          parsNow <- parsNew
          likeNow <- likeNew 
          postNow <- postNew
          if (UPDJUMP) updateMat[iter, idj] <- 1
        }        
      }
    }
    
    # Random effects fertility:
    if (RANDEFFU & iter / 2 == floor(iter / 2)) {
      postUNow <- .CalcPostRandEffU(dataObj = dataObjNow, pars = parsNow,
                                    like = likeNow)
      # Sample U values:
      parsNew <- parsNow
      parsNew$u <- rnorm(n = dataObj$ni, mean = parsNow$u, sd = 0.1)
      likeNew <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNew,
                               FertFun = FertFun,
                               FertFun.numeric = FertFun.numeric)
      postUNew <- .CalcPostRandEffU(dataObj = dataObjNow, pars = parsNew,
                                    like = likeNew)
      acceptRatio <- exp(postUNew - postUNow)
      ranU <- runif(dataObj$ni)
      idUpd <- which(acceptRatio > ranU)
      parsNow$u[idUpd] <- parsNew$u[idUpd]
      postUNow[idUpd] <- postUNew[idUpd]
      likeNow <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNow,
                               FertFun = FertFun,
                               FertFun.numeric = FertFun.numeric)
      postNow <- .CalcPostTheta(pars = parsNow, like = likeNow,
                                parObj = parObj)
      
      # Sample sigma:
      parsNow$uSd <- .SampleUSig(dataObj = dataObjNow, pars = parsNow,
                                 parObj = parObj)
      postUNow <- .CalcPostRandEffU(dataObj = dataObjNow, pars = parsNow,
                                    like = likeNow)
      
    }
    
    # Random effects IBI:
    if (RANDEFFV & iter / 2 == floor(iter / 2)) {
      postVNow <- .CalcPostRandEffV(dataObj = dataObjNow, pars = parsNow, 
                                    like = likeNow)
      # Sample U values:
      parsNew <- parsNow
      parsNew$v[dataObj$idv] <- rnorm(n = dataObj$niv, 
                                      mean = parsNow$v[dataObj$idv], sd = 0.2)
      likeNew <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNew, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postVNew <- .CalcPostRandEffV(dataObj = dataObjNow, pars = parsNew, 
                                    like = likeNew)
      acceptRatio <- exp(postVNew - postVNow)
      ranV <- runif(dataObj$niv)
      idUpd <- which(acceptRatio > ranV)
      parsNow$v[dataObj$idv[idUpd]] <- parsNew$v[dataObj$idv[idUpd]]
      postVNow[dataObj$idv[idUpd]] <- postVNew[dataObj$idv[idUpd]]
      likeNow <- .CalcLikeFert(dataObj = dataObjNow, pars = parsNow, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                                parObj = parObj)
      if (RANDEFFU) {
        postUNow <- .CalcPostRandEffU(dataObj = dataObjNow, pars = parsNow, 
                                      like = likeNow)
        
      }
      
      # Sample sigma:
      parsNow$vSd <- .SampleVSig(dataObj = dataObjNow, pars = parsNow, 
                                 parObj = parObj)
      postVNow <- .CalcPostRandEffV(dataObj = dataObjNow, pars = parsNow, 
                                    like = likeNow)
      
    } 
    
    # Sample unknown ages:
    if (dataObj$AgeUpdate) {
      # Update individual likelihood:
      likeIndNow <- c(t(likeNow) %*% dataObj$rMat)
      
      # Posterior for proposed ages:
      postIndNow <- likeIndNow + priorIndNow
      
      # Sample ages:
      dataObjNew <- .SampleAges(dataObj = dataObjNow, algObj = algObj)
      
      # Calculate likelihood and posteriors:
      likeNew <- .CalcLikeFert(dataObj = dataObjNew, pars = parsNow, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      
      # Update individual likelihood:
      likeIndNew <- c(t(likeNew) %*% dataObj$rMat)
      
      # Prior for proposed ages:
      priorNew <- .dtnorm(x = dataObjNew$data$Age, 
                          mean = dataObj$data$Age, sd = parObj$agePriorSD, 
                          lower = dataObj$data$MinAge, 
                          upper = dataObj$data$MaxAge, 
                          log = TRUE)
      priorNew[which(priorNew == Inf)] <- 0
      priorIndNew <- c(t(priorNew) %*% dataObj$rMat)
      
      # Posterior for proposed ages:
      postIndNew <- likeIndNew + priorIndNew
      postRatio <- exp(postIndNew - postIndNow)
      idUpd <- dataObj$idAgeUpd[which(postRatio[dataObj$idAgeUpd] > 
                                     runif(n = dataObj$nAgeUpd))]
      if (length(idUpd) > 0) {
        idUpdAll <- which(dataObj$data$indID %in% dataObj$unID[idUpd])
        postIndNow[idUpd] <- postIndNew[idUpd]
        priorIndNow[idUpd] <- priorIndNew[idUpd]
        likeIndNow[idUpd] <- likeIndNew[idUpd]
        priorNow[idUpdAll] <- priorNew[idUpdAll]
        likeNow[idUpdAll] <- likeNew[idUpdAll]
        dataObjNow$data$Age[idUpdAll] <- dataObjNew$data$Age[idUpdAll]
        dataObjNow$x[idUpdAll] <- dataObjNew$x[idUpdAll]
        if (inherits(dataObjNow, "baftaIndExt")) {
          dataObjNow$w[idUpdAll] <- dataObjNew$w[idUpdAll]
        }
        postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                                  parObj = parObj)
      }
    }
    
    # Sample unknown number of offspring:
    if (dataObj$offsUpdate) {
      # Update current prior and posterior:
      postOffsNow <- .CalcPostOffs(dataObj = dataObjNow, parObj = parObj, 
                                     like = likeNow, FertFun = FertFun, 
                                     FertFun.numeric = FertFun.numeric)

      # Sample new number of offspring:
      dataObjNew <- .SampleOffspring(dataObj = dataObjNow, algObj = algObj)
      
      # Calculate likelihood and posteriors:
      likeNew <- .CalcLikeFert(dataObj = dataObjNew, pars = parsNow, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      
      postOffsNew <- .CalcPostOffs(dataObj = dataObjNew, parObj = parObj, 
                                   like = likeNew, FertFun = FertFun, 
                                   FertFun.numeric = FertFun.numeric)
      mhRatio <- .dtnorm(x = dataObjNow$y[dataObj$idOffsUpd],
                           mean = dataObjNew$y[dataObj$idOffsUpd],
                           sd = 1, 
                           lower = dataObj$o[dataObj$idOffsUpd], log = TRUE) -
        .dtnorm(x = dataObjNew$y[dataObj$idOffsUpd],
                mean = dataObjNow$y[dataObj$idOffsUpd],
                sd = 1, 
                lower = dataObj$o[dataObj$idOffsUpd], log = TRUE)
      
      postRatio <- exp(postOffsNew - postOffsNow + mhRatio)
      idUpd <- which(postRatio > runif(n = dataObj$nOffsUpd))
      if (length(idUpd) > 0) {
        idUpdAll <- dataObj$idOffsUpd[idUpd]
        likeNow[idUpdAll] <- likeNew[idUpdAll]
        postOffsNow[idUpd] <- postOffsNew[idUpd]
        dataObjNow$y[idUpdAll] <- dataObjNew$y[idUpdAll]
        dataObjNow$data$estOffspring[idUpdAll] <- dataObjNew$y[idUpdAll]
        postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                                  parObj = parObj)
      }
    }

    # Update jumps:
    if (UPDJUMP) {
      if (iter %in% updJumpIter) {
        jumpSD <- .UpdateJumps(updMat = updateMat, jumps = jumpSD, 
                               iter = iter)
        jumpOut <- rbind(jumpOut, jumpSD)
      }
    } else {
      parOut[iter, parObj$thetaName] <- parsNow$theta
      if (inherits(dataObj, "baftaIndExt")) {
        parOut[iter, "vSd"] <- parsNow$vSd
      }
      likePostOut[iter, ] <- c(sum(likeNow), postNow)
      if (dataObj$AgeUpdate) {
        ageOut[iter, ] <- dataObjNow$data$Age[dataObj$idFirstAge]
      }
      if (dataObj$offsUpdate) {
        offsOut[iter, ] <- dataObjNow$y[dataObj$idOffsUpd]
      }
      if (RANDEFFU) {
        uOut[iter, ] <- parsNow$u
        uSdOut[iter] <- parsNow$uSd
      }
      if (RANDEFFV) {
        vOut[iter, ] <- parsNow$v
        # vSdOut[iter] <- parsNow$vSd
      }
    }
  }
  
  
  if (UPDJUMP) {
    nJumps <- nrow(jumpOut)
    idJumps <- floor(nJumps / 2):nJumps
    jumpSdFin <- apply(jumpOut[idJumps, ], 2, mean)
  } else {
    jumpSdFin <- jumpSD
    jumpOut <- NA
  }
  outList <- list(theta = parOut, likePost = likePostOut, u = uOut,
                  uSd = uSdOut, v = vOut, age = ageOut,
                  offsOut = offsOut, jumps = jumpSdFin, jumpMat = jumpOut) 
  return(outList)
}

# ======================= #
# ==== MCMC outputs: ====
# ======================= #
# Function to calculate convergence statistics 
# based on Gelman et al. (2014).
.CalcPSRF <- function(object, keep, nsim) {
  nthin <- length(keep)
  Means <- t(sapply(1:nsim, function(i) {
    apply(object[[i]]$theta[keep, ], 2, mean)
  }))
  Vars <- t(sapply(1:nsim, function(i) {
    apply(object[[i]]$theta, 2, var)
  }))
  meanall <- apply(Means, 2, mean)
  B <- nthin / (nsim - 1) * apply(t((t(Means) - meanall)^2), 2, sum)
  W <- 1 / nsim * apply(Vars, 2, sum)
  Varpl <- (nthin - 1) / nthin * W + 1 / nthin * B
  Rhat <- sqrt(Varpl / W)
  Rhat[which(Varpl == 0)] <- 1
  Rhat[which(Rhat < 1)] <- 1
  conv <- cbind(B, W, Varpl, Rhat)
  rownames(conv) <- colnames(Means)
  return(conv)
}

# Calculate effective sample size:
.CalcNeff <- function(object, keep, nsim, Rhat) {
  nthin <- length(keep)
  Varpl <- Rhat[, "Varpl"]
  Tlags <- min(c(200, nthin - 1))
  nTheta <- ncol(object[[1]]$theta)
  rhoHat <- matrix(NA, nrow = Tlags, ncol = nTheta)
  for (tt in 1:Tlags) {
    Vt <- 1 / (nsim * (nthin - tt)) * 
      apply(sapply(1:nsim, function(im) {
        bMat <- object[[im]]$theta[keep, ]
        apply(bMat, 2, function(bi) {
          sum((bi[1:(nthin - tt + 1)] - bi[tt:nthin])^2)
        })
      }), 1, sum)
    rhoHat[tt, ] <- 1 - Vt / Varpl
  }
  
  # Find T:
  Tthresh <- apply(rhoHat, 2, function(rhoi) {
    which(rhoi[-1] + rhoi[-Tlags] < 0)[1]
  })
  
  idna <- which(is.na(Tthresh))
  if (length(idna) > 0) Tthresh[idna] <- Tlags
  
  neff <- floor(sapply(1:nTheta, function(ip) {
    nsim * nthin / (1 + 2 * sum(rhoHat[1:Tthresh[ip], ip]))
  }))
  names(neff) <- rownames(Rhat)
  return(neff)
}

# Calculate DIC:
.CalcDIC <- function(likelihood, k) {
  L <- length(likelihood)
  Dm <- -2 * likelihood
  Dave <- mean(Dm)
  pD <- 1/2 * 1/(L-1) * sum((Dm - Dave)^2)
  DIC <- Dave + pD
  modSel <- c(Dave, pD, k, DIC)
  names(modSel) <- c("D.ave", "pD", "k", "DIC")
  return(modSel)
}

# Y predicted:
.CalcYpred <- function(dataObj, ...) UseMethod(".CalcYpred")
.CalcYpred.baftaAggr <- function(dataObj, thetaMat, uMat, ageMat, 
                                 offsMat, parObj, FertFun, 
                                 FertFun.numeric) {
  x <- dataObj$x
  
  # Offset:
  nx <- dataObj$data$nParents[dataObj$idages]
  
  # Ny:
  ny <- length(dataObj$idages)
  
  # Total adult ages:
  nytot <- length(which(dataObj$data$Age >= dataObj$alpha))
  
  # Calculate estimated fertility from parameter posteriors:
  yPred <- t(apply(thetaMat, 1, function(be) {
    fert <- FertFun(beta = be, x = x + 0.5)
    # ------------------------------ #
    # Negative binomial 2023-11-14:
    rho <- be["gamma"] / (fert + be["gamma"])
    yp <- rep(NA, nytot)
    yp[dataObj$idpar] <- rnbinom(n = ny, size = nx * be["gamma"], prob = rho)
    # ------------------------------ #
    # yp <- rpois(n = ny, lambda = nx * fert)
    return(yp)
  }))
  
  return(yPred)
}

.CalcYpred.baftaIndSimp <- function(dataObj, thetaMat, uMat, ageMat, 
                                    offsMat, parObj, FertFun, 
                                    FertFun.numeric) {
  # Calculate estimated fertility from parameter posteriors:
  nth <- nrow(thetaMat)
  if (dataObj$AgeUpdate) {
    yPred <- t(sapply(1:nth, function(ith) {
      be <- thetaMat[ith, ]
      x <- ageMat[ith, ] + 0.5
      fert <-  FertFun(beta = be, x = x)
      rho <- be["gamma"] / (fert + be["gamma"])
      yp <- rnbinom(n = dataObj$n, size = be["gamma"], prob = rho)
      return(yp)
    }))
  } else if (dataObj$offsUpdate) {
    ages <- dataObj$x[dataObj$idOffsUpd]
    mu <- FertFun(parObj$thetaPriorMean, ages)
    po <- dataObj$data$obsProp[dataObj$idOffsUpd]
    
    yPred <- t(sapply(1:nth, function(ith) {
      be <- thetaMat[ith, ]
      o <- offsMat[ith, ]
      fert <-  FertFun(beta = be, x = dataObj$x + 0.5)
      rho <- be["gamma"] / (fert + be["gamma"])
      yp <- rnbinom(n = dataObj$n, size = be["gamma"], prob = rho)
      op <- yp
      op[dataObj$idOffsUpd] <- yp[dataObj$idOffsUpd] - 
        rnbinom(n = dataObj$nOffsUpd, 
                size = parObj$thetaPriorMean["gamma"], 
                mu = mu * (1 - po))
      op[which(op < 0)] <- 0
      return(op)
    }))
    
  } else {
    yPred <- t(apply(thetaMat, 1, function(be) {
      fert <-  FertFun(beta = be, x = dataObj$x + 0.5)
      rho <- be["gamma"] / (fert + be["gamma"])
      yp <- rnbinom(n = dataObj$n, size = be["gamma"], prob = rho)
      return(yp)
    }))
  }
  return(yPred)
}

.CalcYpred.baftaIndExt <- function(dataObj, thetaMat, uMat, ageMat, 
                                   offsMat, parObj, FertFun, 
                                   FertFun.numeric) {
  # Calculate estimated fertility from parameter posteriors:
  nth <- nrow(thetaMat)
  if (dataObj$AgeUpdate) {
    yPred <- t(sapply(1:nth, function(ith) {
      be <- thetaMat[ith, ]
      x <- ageMat[ith, ]
      fert <-  FertFun(beta = be, x = x)
      rho <- be["gamma"] / (fert + be["gamma"])
      yp <- rnbinom(n = dataObj$n, size = be["gamma"], prob = rho)
      return(yp)
    }))
    
  } else {
    yPred <- t(apply(thetaMat, 1, function(be) {
      fert <-  FertFun(beta = be, x = dataObj$x)
      rho <- be["gamma"] / (fert + be["gamma"])
      yp <- rnbinom(n = dataObj$n, size = be["gamma"], prob = rho)
      return(yp)
    }))
    
  }
  return(yPred)
}

# Predictive loss:
.CalcPredLoss <- function(dataObj, ...) UseMethod(".CalcPredLoss")
.CalcPredLoss.baftaAggr <- function(dataObj, yPred) {
  y1 <- apply(yPred[, dataObj$idpar], 2, mean, na.rm = TRUE)
  y2 <- apply(yPred[, dataObj$idpar], 2, var, na.rm = TRUE)
  gm <- sum((dataObj$data$nOffspring[dataObj$idages] - y1)^2)
  pm <- sum(y2)
  dev <- gm + pm
  devMat <- matrix(c(gm, pm, dev), 1, 3,
                   dimnames = list("", c("Good. Fit", "Penalty", "Deviance")))
  return(devMat)
}

.CalcPredLoss.baftaIndSimp <- function(dataObj, yPred) {
  y1 <- apply(yPred, 2, mean, na.rm = TRUE)
  y2 <- apply(yPred, 2, var, na.rm = TRUE)
  gm <- sum((dataObj$data$nOffspring - y1)^2)
  pm <- sum(y2)
  dev <- gm + pm
  devMat <- matrix(c(gm, pm, dev), 1, 3,
                   dimnames = list("", c("Good. Fit", "Penalty", "Deviance")))
  return(devMat)
}

.CalcPredLoss.baftaIndExt <- function(dataObj, yPred) {
  y1 <- apply(yPred, 2, mean, na.rm = TRUE)
  y2 <- apply(yPred, 2, var, na.rm = TRUE)
  gm <- sum((dataObj$data$nOffspring - y1)^2)
  pm <- sum(y2)
  dev <- gm + pm
  devMat <- matrix(c(gm, pm, dev), 1, 3,
                   dimnames = list("", c("Good. Fit", "Penalty", "Deviance")))
  return(devMat)
}

# ============================= END OF CODE ================================== #
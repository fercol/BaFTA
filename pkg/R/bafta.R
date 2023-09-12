# ============================= CODE METADATA ================================ #
# PACKAGE: BaFTA
# AUTHOR: Fernando Colchero
# DATE:
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
  
  # Start timer:
  Start <- Sys.time()
  
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
  if (grepl("indiv", dataType)) {
    RANDEFFU <- TRUE
  } else {
    RANDEFFU <- FALSE
  }
  
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
  
  # Number of iterations kept for inference:
  keep <- seq(burnin, niter, thinning)
  
  # Verify if starting parameters are specified:
  if ("thetaStart" %in% argNames) {
    if (length(argList$thetaStart) != parObj$p) {
      stop(sprintf("Length of 'thetaStart' argument should be %s.\n", parObj$p))
    } else {
      parObj$thetaStart <- argList$thetaStart
      names(parObj$thetaStart) <- parObj$thetaName
    }
  }
  
  # Verify if mean priors are specified:
  if ("thetaPriorMean" %in% argNames) {
    if (length(argList$thetaPriorMean) != parObj$p) {
      stop(sprintf("Length of 'thetaPriorMean' argument should be %s.\n", parObj$p))
    } else {
      parObj$thetaPriorMean <- argList$thetaPriorMean
      names(parObj$thetaPriorMean) <- parObj$thetaName
    }
  }
  
  # Verify if SD priors are specified:
  if ("thetaPriorSD" %in% argNames) {
    if (length(argList$thetaPriorSD) != parObj$p) {
      stop(sprintf("Length of 'thetaPriorSD' argument should be %s.\n", 
                   parObj$p))
    } else {
      parObj$thetaPriorSD <- argList$thetaPriorSD
      names(parObj$thetaPriorSD) <- parObj$thetaName
    }
  }
  
  # Run jump sd:
  if (UPDJUMP) {
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

  # Start parallel computing:
  sfInit(parallel = TRUE, cpus = ncpus)
  
  # Load BaFTA to CPUS:
  sfLibrary("BaFTA", character.only = TRUE,
            warn.conflicts = FALSE)
  # sfLibrary(BaFTA)
  
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
  
  # Extract variables for coefficients and DIC:
  for (ic in 1:nsim) {
    if (ic == 1) {
      thetaMat <- outMCMC[[ic]]$theta[keep, ]
      likeMat <- outMCMC[[ic]]$likePost[keep, ]
      if (RANDEFFU) {
        uSdVec <- outMCMC[[ic]]$uSd[keep]
        uMat <- outMCMC[[ic]]$u[keep, ]
      } else {
        uMat <- NA
      }
      if (RANDEFFV) {
        vSdVec <- outMCMC[[ic]]$vSd[keep]
        vMat <- outMCMC[[ic]]$v[keep, ]
      } else {
        vMat <- NA
      }
      
    } else {
      thetaMat <- rbind(thetaMat, outMCMC[[ic]]$theta[keep, ])
      likeMat <- rbind(likeMat, outMCMC[[ic]]$likePost[keep, ])
      if (RANDEFFU) {
        uSdVec <- c(uSdVec, outMCMC[[ic]]$uSd[keep])
        uMat <- rbind(uMat, outMCMC[[ic]]$u[keep, ])
      } 
      if (RANDEFFV) {
        vSdVec <- c(vSdVec, outMCMC[[ic]]$vSd[keep])
        vMat <- rbind(vMat, outMCMC[[ic]]$v[keep, ])
      } 
    }
  }
  
  # Calculate convergence statistics:
  Conv <- .CalcPSRF(object = outMCMC, keep = keep, nsim = ncpus)
  
  # Extract average parameters:
  coeffs <- cbind(Mean = apply(thetaMat[, parObj$idSamp],  2, mean), 
                  SD = apply(thetaMat[, parObj$idSamp], 2, sd),
                  Lower = apply(thetaMat[, parObj$idSamp], 2, quantile, 0.025),
                  Upper = apply(thetaMat[, parObj$idSamp], 2, quantile, 0.975),
                  Rhat = Conv[parObj$idSamp, "Rhat"])
  
  if (RANDEFFU) {
    # Include random effect standard error:
    coeffs <- rbind(coeffs, uSd = c(Mean = mean(uSdVec), SD = sd(uSdVec),
                                    Lower = quantile(uSdVec, 0.025),
                                    Upper = quantile(uSdVec, 0.975),
                                    Rhat = 1))
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
    coeffs <- rbind(coeffs, vSd = c(Mean = mean(vSdVec), SD = sd(vSdVec),
                                    Lower = quantile(vSdVec, 0.025),
                                    Upper = quantile(vSdVec, 0.975),
                                    Rhat = 1))
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
  
  # Calculate DIC:
  DIC <- .CalcDIC(likelihood = likeMat[, "Likelihood"], 
                  k = length(parObj$idSamp))
  
  # Y predicted:
  yPred <- .CalcYpred(dataObj = dataObj, thetaMat = thetaMat, uMat = uMat,
                      FertFun = FertFun, FertFun.numeric = FertFun.numeric)
  
  # Calculate predictive loss:
  PredLoss <- .CalcPredLoss(dataObj = dataObj, yPred = yPred)
  
  # Extract aggregated fertility and number of offspring:
  
  if (algObj$dataType == "aggregated") {
    aggrData <- object[which(object$Age >= dataObj$alpha), ]
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
    aggrData <- data.frame(Age = xag[-nxag], tempag, 
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
  End <- Sys.time()
  compTime <- sprintf("%s mins", signif(as.numeric(End - Start, 
                                                   units = "mins"), 2))
  
  # Settings:
  settings <- list(model = model, dataType = dataType, niter = niter, 
                   burnin = burnin, thinning = thinning, nsim = nsim, 
                   ncpus = ncpus, compTime = compTime)
  
  # store output:
  fullOut <- list(coefficients = coeffs, x = xv, fert = fertQuant, 
                  theta = thetaMat, uSd = uSdVec, vSd = vSdVec, 
                  likePost = likeMat, DIC = DIC, PredLoss = PredLoss, 
                  pred = predQuant, aggrData = aggrData,
                  runs = outMCMC, data = dataObj, settings = settings,
                  keep = keep, params = parObj)
  
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
    nsim <- x$settings$nsim
    ylim <- sapply(1:pSamp, function(ipar) {
      range(sapply(1:nsim, function(isim) {
        range(x$runs[[isim]]$theta[, ipar], na.rm = TRUE)
      }), na.rm = TRUE)
    })
    if (grepl("indiv", x$settings$dataType)) {
      uyl <- range(sapply(1:nsim, function(isim) {
        range(x$runs[[isim]]$uSd)
      }))
      ylim <- cbind(ylim, uyl)
      pPars <- pSamp + 1
      pName <- c(pName, "uSd")
      if (x$settings$dataType == "indivExtended") {
        vyl <- range(sapply(1:nsim, function(isim) {
          range(x$runs[[isim]]$vSd)
        }))
        ylim <- cbind(ylim, vyl)
        pPars <- pPars + 1
        pName <- c(pName, "vSd")
      }
    } else {
      pPars <- pSamp
    }
    if ("ylim" %in% names(argList)) {
      ylim <- argList$ylim
    }
    
    colnames(ylim) <- pName
    par(mfrow = c(ceiling(pPars / 2), 2), mar = c(4, 4, 3, 1))
    for (ipar in idSamp) {
      plot(idkeep, x$runs[[1]]$theta[idkeep, ipar], type = 'l', 
           ylim = ylim[, ipar], main = pName[ipar], xlab = "Iteration", 
           ylab = "Parameter")
      for (ic in 2:nsim) {
        lines(idkeep, x$runs[[ic]]$theta[idkeep, ipar], col = ic)
      }
    }
    if (grepl("indiv", x$settings$dataType)) {
      plot(idkeep, x$runs[[1]]$uSd[idkeep], type = 'l',
           ylim = ylim[, ncol(ylim)], xlab = "Iteration", ylab = "Parameter",
           main = "uSd")
      for (ic in 2:nsim) {
        lines(idkeep, x$runs[[ic]]$uSd[idkeep], col = ic)
      }
      if (x$settings$dataType == "indivExtended") {
        plot(idkeep, x$runs[[1]]$vSd[idkeep], type = 'l',
             ylim = ylim[, ncol(ylim)], xlab = "Iteration", ylab = "Parameter",
             main = "vSd")
        for (ic in 2:nsim) {
          lines(idkeep, x$runs[[ic]]$vSd[idkeep], col = ic)
        }
        
      }
    }
    
  } else if (type == "density") {
    pSamp <- x$params$pSamp
    idSamp <- x$params$idSamp
    pName <- x$params$thetaName[idSamp]
    parMat <- x$theta[, idSamp]
    if (grepl("indiv", x$settings$dataType)) {
      parMat <- cbind(parMat, uSD = x$uSd)
      idSamp <- c(idSamp, max(idSamp) + 1)
      pSamp <- length(idSamp)
      pName <- c(pName, "uSd")
      if (x$settings$dataType == "indivExtended") {
        parMat <- cbind(parMat, vSD = x$vSd)
        idSamp <- c(idSamp, max(idSamp) + 1)
        pSamp <- length(idSamp)
        pName <- c(pName, "vSd")
      }
    }
    
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
    idages <- which(object$Fertility > 0)
    alpha <- object$Age[which(object$Fertility > 0)[1]]
    idages <- which(object$Age >= alpha)
    x <- object$Age[idages] - alpha
    
    # Adjust first age > 0 for certain models:
    if (algObj$model %in% c("gamma", "beta", "gammaMixture", "Hadwiger", 
                     "HadwigerMixture")) {
      x[which(x == 0)] <- 0.005
    }

    # Create data object:
    do <- list(data = object, idages = idages, alpha = alpha, x = x,
               xMax = max(x), n = n)
    class(do) <- "baftaAggr"
  } else if (algObj$dataType == "indivSimple") {
    if (is.na(algObj$minAge)) {
      alpha <- min(object$Age)
    } else {
      alpha <- algObj$minAge
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
    ni <- ncol(rMat)
    do <- list(data = object, x = x, y = y, rMat = rMat, ni = ni, 
               xMax = max(x), n = n, alpha = alpha)
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
    ni <- ncol(rMat)
    z <- object$IBI
    nIBI <- c(rep(1, n) %*% (rMat * (1 - object$First)))
    idv <- which(nIBI > 0)
    idIBI <- rep(0, ni)
    idIBI[idv] <- 1
    niv <- sum(idIBI)
    z[which(object$First == 0)] <- z[which(object$First == 0)] - tau
    w <- object$Age - alpha
    do <- list(data = object, x = x, y = y, rMat = rMat, ni = ni, 
               xMax = max(x), z = z, w = w, n = n, alpha = alpha, tau = tau,
               indv = idIBI, idv = idv, niv = niv)
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
    startBe <- c(10, 1.5, 2.5, 0, dataObj$xMax)
    priorMeanBe <- c(10, 1.5, 2.5, 0, dataObj$xMax)
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
  priorSdBe <- rep(1, nBe)
  nameBe <- sprintf("b%s", 1:nBe - 1)
  names(startBe) <- nameBe
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
  if (algObj$dataType %in% c("indivSimple", "indivExtended")) {
    uSd <- 0.5
    uPrior1 <- 0.01
    uPrior2 <- 0.1
  } 
  if (algObj$dataType == "indivExtended") {
    etaStart <- 1 / mean(dataObj$z[which(dataObj$data$First == 0)])
    gammaStart <- 1 / mean(dataObj$w[which(dataObj$data$First == 1)])
    theta <- c(defBet$beta, eta = etaStart, gamma = gammaStart)
    thetaPriorMean <- c(defBet$priorMean, eta = 1, gamma = 1)
    thetaPriorSD <- c(defBet$priorSD, eta = 1, gamma = 1)
    thetaLower <- c(defBet$low, 0, 0)
    thetaUpper <- c(defBet$upp, Inf, Inf)
    idSamp <- c(defBet$idSamp, defBet$p + c(1:2))

    # Random effects for IBI:
    vSd <- 0.5
    vPrior1 <- 0.01
    vPrior2 <- 0.1
    
  } else {
    theta <- defBet$beta
    thetaPriorMean <- defBet$priorMean
    thetaPriorSD <- defBet$priorSD
    thetaLower <- defBet$low
    thetaUpper <- defBet$upp
    idSamp <- defBet$idSamp
  }
  parList <- list(thetaStart = theta, thetaPriorMean = thetaPriorMean, 
                  thetaPriorSD = thetaPriorSD, thetaLower = thetaLower,
                  thetaUpper = thetaUpper, thetaName = names(theta),
                  p = length(theta), idSamp = idSamp, 
                  pSamp = length(idSamp), uSd = uSd, uPrior1 = uPrior1, 
                  uPrior2 = uPrior2, vSd = vSd, vPrior1 = vPrior1, 
                  vPrior2 = vPrior2)
  return(parList)
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
  fert <- FertFun(beta = pars$theta, x = dataObj$x)
  lk <- dpois(dataObj$data$nOffspring[dataObj$idages], 
              lambda = dataObj$data$nParents[dataObj$idages] * fert, 
              log = TRUE)
  return(lk)
}

.CalcLikeFert.baftaIndSimp <- function(dataObj, pars, FertFun, 
                                       FertFun.numeric) {
  fert <- FertFun(beta = pars$theta, x = dataObj$x) * 
    exp(c(dataObj$rMat %*% pars$u))
  lk <- dpois(dataObj$data$nOffspring, lambda = fert, log = TRUE)
  return(lk)
}

# Likelihood function:
.CalcLikeFert.baftaIndExt <- function(dataObj, pars, FertFun, 
                                      FertFun.numeric) {
  fert <- FertFun(beta = pars$theta, x = dataObj$x) * 
    exp(c(dataObj$rMat %*% pars$u))
  lk <- dpois(dataObj$data$nOffspring, lambda = fert, log = TRUE) +
    dexp(dataObj$z, rate = pars$theta["eta"] * 
           exp(c(dataObj$rMat %*% pars$v)), log = TRUE) * 
    (1 - dataObj$data$First) +
    dexp(dataObj$w, rate = pars$theta["gamma"], log = TRUE) * 
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
  MHr <- .dtnorm(parsNow$theta[ip], mean = parsNew$theta[ip], sd = jumpSD[ip],
                lower = parObj$thetaLower[ip], upper = parObj$thetaUpper[ip], 
                log = TRUE) -
    .dtnorm(parsNew$theta[ip], mean = parsNow$theta[ip], sd = jumpSD[ip],
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
  if (grepl("indiv", algObj$dataType)) {
    RANDEFFU <- TRUE
    parsNow$u <- rep(0, dataObj$ni)
    parsNow$uSd <- 1
  } else {
    RANDEFFU <- FALSE
  }
  
  # IBI random effects parameter:
  if (algObj$dataType == "indivExtended") {
    RANDEFFV <- TRUE
    parsNow$v <- rep(0, dataObj$ni)
    parsNow$vSd <- 1
  } else {
    RANDEFFV <- FALSE
  }
  
  # Jitter parameters:
  if (sim > 1) {
    parsNow$theta[parObj$idSamp] <- 
      .rtnorm(n = 1, mean = parObj$thetaStart[parObj$idSamp], 
              sd = rep(0.1, length(parObj$idSamp)), 
              lower = parObj$thetaLower[parObj$idSamp],
              upper = parObj$thetaUpper[parObj$idSamp])
  }
  
  # Calculate likelihood and posteriors:
  likeNow <- .CalcLikeFert(dataObj = dataObj, pars = parsNow, FertFun = FertFun, 
                           FertFun.numeric = FertFun.numeric)
  postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                            parObj = parObj)
  if (RANDEFFU) {
    postUNow <- .CalcPostRandEffU(dataObj = dataObj, pars = parsNow, 
                                 like = likeNow)
  }
  
  if (RANDEFFV) {
    postVNow <- .CalcPostRandEffV(dataObj = dataObj, pars = parsNow, 
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
    
  } else {
    parOut <- matrix(NA, nrow = niter, ncol = parObj$p,
                     dimnames = list(NULL, parObj$thetaName))
    likePostOut <- matrix(NA, nrow = niter, ncol = 2,
                          dimnames = list(NULL, c("Likelihood", "Posterior")))
    parOut[1, ] <- parsNow$theta
    likePostOut[1, ] <- c(sum(likeNow), postNow)
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
      vSdOut <- rep(NA, niter)
      vSdOut[1] <- parsNow$vSd
    } else {
      vOut <- NA
      vSdOut <- NA
    }
    
    
  }
  
  for (iter in 2:niter) {
    for (ip in parObj$idSamp) {
      idj <- which(parObj$idSamp == ip)
      parsNew <- parsNow
      parsNew$theta[ip] <- .rtnorm(n = 1, mean = parsNow$theta[ip], 
                                   sd = jumpSD[idj],
                                   lower = parObj$thetaLower[ip], 
                                   upper = parObj$thetaUpper[ip])
      likeNew <- .CalcLikeFert(dataObj = dataObj, pars = parsNew, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postNew <- .CalcPostTheta(pars = parsNew, like = likeNew, 
                                parObj = parObj)
      mhRatio <- .CalcMHratio(parsNow = parsNow, parsNew = parsNew, 
                              jumpSD = jumpSD, parObj = parObj, ip = ip)
      postRatio <- exp(postNew - postNow + mhRatio)
      # ========== #
      # Debugging SEC 2b:
      if (DEBUG & iter == 2) {
        cat("Sec 2b\n", file = logFile, append = TRUE)
      }
      # ========== #
      
      if (!is.na(postRatio)) {
        if (postRatio > runif(n = 1)) {
          parsNow <- parsNew
          likeNow <- likeNew 
          postNow <- postNew
          if (UPDJUMP) updateMat[iter, idj] <- 1
        }        
      }
    }
    
    # Random effects:
    if (RANDEFFU & iter / 2 == floor(iter / 2)) {
      postUNow <- .CalcPostRandEffU(dataObj = dataObj, pars = parsNow, 
                                    like = likeNow)
      # Sample U values:
      parsNew <- parsNow
      parsNew$u <- rnorm(n = dataObj$ni, mean = parsNow$u, sd = 0.1)
      likeNew <- .CalcLikeFert(dataObj = dataObj, pars = parsNew, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postUNew <- .CalcPostRandEffU(dataObj = dataObj, pars = parsNew, 
                                    like = likeNew)
      acceptRatio <- exp(postUNew - postUNow)
      ranU <- runif(dataObj$ni)
      idUpd <- which(acceptRatio > ranU)
      parsNow$u[idUpd] <- parsNew$u[idUpd]
      postUNow[idUpd] <- postUNew[idUpd]
      likeNow <- .CalcLikeFert(dataObj = dataObj, pars = parsNow, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                                parObj = parObj)
      
      # Sample sigma:
      parsNow$uSd <- .SampleUSig(dataObj = dataObj, pars = parsNow, 
                                 parObj = parObj)
      postUNow <- .CalcPostRandEffU(dataObj = dataObj, pars = parsNow, 
                                    like = likeNow)
      
    } 
    
    # Random effects:
    if (RANDEFFV & iter / 2 == floor(iter / 2)) {
      postVNow <- .CalcPostRandEffV(dataObj = dataObj, pars = parsNow, 
                                    like = likeNow)
      # Sample U values:
      parsNew <- parsNow
      parsNew$v[dataObj$idv] <- rnorm(n = dataObj$niv, 
                                      mean = parsNow$v[dataObj$idv], sd = 0.1)
      likeNew <- .CalcLikeFert(dataObj = dataObj, pars = parsNew, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postVNew <- .CalcPostRandEffV(dataObj = dataObj, pars = parsNew, 
                                    like = likeNew)
      acceptRatio <- exp(postVNew - postVNow)
      ranV <- runif(dataObj$niv)
      idUpd <- which(acceptRatio > ranV)
      parsNow$v[dataObj$idv[idUpd]] <- parsNew$v[dataObj$idv[idUpd]]
      postVNow[dataObj$idv[idUpd]] <- postVNew[dataObj$idv[idUpd]]
      likeNow <- .CalcLikeFert(dataObj = dataObj, pars = parsNow, 
                               FertFun = FertFun, 
                               FertFun.numeric = FertFun.numeric)
      postNow <- .CalcPostTheta(pars = parsNow, like = likeNow, 
                                parObj = parObj)
      postUNow <- .CalcPostRandEffU(dataObj = dataObj, pars = parsNow, 
                                    like = likeNow)
      
      # Sample sigma:
      parsNow$vSd <- .SampleVSig(dataObj = dataObj, pars = parsNow, 
                                 parObj = parObj)
      postVNow <- .CalcPostRandEffV(dataObj = dataObj, pars = parsNow, 
                                    like = likeNow)
      
    } 
    
    # Update jumps:
    if (UPDJUMP) {
      if (iter %in% updJumpIter) {
        jumpSD <- .UpdateJumps(updMat = updateMat, jumps = jumpSD, 
                               iter = iter)
        jumpOut <- rbind(jumpOut, jumpSD)
      }
    } else {
      parOut[iter, ] <- parsNow$theta
      likePostOut[iter, ] <- c(sum(likeNow), postNow)
      if (RANDEFFU) {
        uOut[iter, ] <- parsNow$u
        uSdOut[iter] <- parsNow$uSd
      }
      if (RANDEFFV) {
        vOut[iter, ] <- parsNow$v
        vSdOut[iter] <- parsNow$vSd
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
                  uSd = uSdOut, v = vOut, vSd = vSdOut, jumps = jumpSdFin, 
                  jumpMat = jumpOut)
  
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
.CalcYpred.baftaAggr <- function(dataObj, thetaMat, uMat, FertFun, 
                                 FertFun.numeric) {
  x <- dataObj$x
  
  # Offset:
  offs <- dataObj$data$nParents[dataObj$idages]
  
  # Ny:
  ny <- length(dataObj$idages)
  
  # Calculate estimated fertility from parameter posteriors:
  yPred <- t(apply(thetaMat, 1, function(be) {
    fert <- FertFun(beta = be, x = x)
    yp <- rpois(n = ny, lambda = offs * fert)
    return(yp)
  }))
  
  return(yPred)
}

.CalcYpred.baftaIndSimp <- function(dataObj, thetaMat, uMat, FertFun, 
                                    FertFun.numeric) {
  # Calculate estimated fertility from parameter posteriors:
  nth <- nrow(thetaMat)
  yPred <- t(sapply(1:nth, function(ith) {
    the <- thetaMat[ith, ]
    fert <-  FertFun(beta = the, x = dataObj$x) * 
      exp(c(dataObj$rMat %*% uMat[ith, ]))
    yp <- rpois(n = dataObj$n, lambda = fert)
    return(yp)
  }))
  
  return(yPred)
}

.CalcYpred.baftaIndExt <- function(dataObj, thetaMat, uMat, FertFun, 
                                   FertFun.numeric) {
  # Calculate estimated fertility from parameter posteriors:
  nth <- nrow(thetaMat)
  yPred <- t(sapply(1:nth, function(ith) {
    the <- thetaMat[ith, ]
    fert <-  FertFun(beta = the, x = dataObj$x) * 
      exp(c(dataObj$rMat %*% uMat[ith, ]))
    yp <- rpois(n = dataObj$n, lambda = fert)
    return(yp)
  }))
  
  return(yPred)
}

# Predictive loss:
.CalcPredLoss <- function(dataObj, ...) UseMethod(".CalcPredLoss")
.CalcPredLoss.baftaAggr <- function(dataObj, yPred) {
  y1 <- apply(yPred, 2, mean, na.rm = TRUE)
  y2 <- apply(yPred, 2, var, na.rm = TRUE)
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

# ============= #
# ==== END ====
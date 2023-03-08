# ============================= CODE METADATA ================================ #
# PACKAGE: BaFTA
# AUTHOR: Fernando Colchero
# DATE:
# DESCRIPTION: Functions to estimate average age-specific fecundity from
#              alternative model with mixed effects for repeated individuals.
# COMMENTS: dataType: aggregated ("aggr"), individual yearly ("indYear"), 
#                     individual continuous ("indCont").
# ============================== START CODE ================================== #
# ======================================== #
# A) FUNCTIONS AVAILABLE TO THE USER: ==== 
# ======================================== #

# A.1) Data check function:
# ------------------------- #


# A.2) main bafta function:
# ------------------------- #
bafta <- function(object, ...) UseMethod("bafta")

bafta.default <- function(object, dataType = "indYear", formula = NULL, 
                          model = "quadratic", minAge = 0, niter = 11000, 
                          burnin = 1001, thinning = 20, nsim = 1, 
                          parallel = FALSE, ncpus = 2, 
                          updateJumps = TRUE, ...) {
  
  # ------------------- #
  # General data setup:
  # ------------------- #
  # Algorithm information:
  algObj <- .CreateAlgObj(model, dataType, formula, niter, burnin,
                          thinning, updateJumps, nsim, minAge)
  
  # Dataset object:
  dataObj <- .CreateDataObj(object, algObj)
  
  # Covariate object:
  covObj <- .CreateCovObj(object, dataObj, algObj)
  
  
}

# A.3) plotting bafta outputs:
# ---------------------------- #
plot.bafta <- function(x, plot.type = "traces", minAge = 0) {
  op <- par(no.readonly = TRUE)
  if (plot.type == "traces") {
    np <- x$settings['np']
    randEffs <- x$settings["randEffs"]
    ns <- ifelse(randEffs == 1, 1, 0)
    parmatfull <- x$parmatfull
    parmat <- x$parmat
    parnames <- expression(a[0], a[1], a[2], a[3], a[4], a[5])[1:np]
    parnames <- c(parnames, expression(sigma))
    niter <- x$settings["niter"]
    # Visualize results:
    par(mfrow = c(np + ns, 2)) 
    for (pp in 1:(np + ns)) {
      xlim <- c(1, niter)
      ylim <- range(parmatfull[, pp])
      plot(xlim, ylim, col = NA, main = parnames[pp], xlab = "", ylab = "")
      for (al in 1:nsim) {
        idtr <- 1:niter + niter * (al - 1)
        lines(parmatfull[idtr, pp], col = al)
      }
      abline(h = x$coeff[pp, 1], col = 2)
      
      plot(density(parmat[, pp]), main = parnames[pp], xlab = "", ylab = "")
      abline(v = x$coeff[pp, 1], col = 2)
    }
  } else {
    # Extract fecundity model:
    fecmod <- x$settings["model"]
    fec <- fecMods[[sprintf("mod%s", fecmod)]]$fun
    np <- fecMods[[sprintf("mod%s", fecmod)]]$np
    yave <- x$avefec
    xd <- 1:length(yave) - 1 + minAge
    fecCI <- x$fecICs
    par(mfrow = c(1, 1))
    ylim <- c(0, max(yave, fecCI))
    par(mfrow = c(1, 1))
    plot(range(xd), ylim, col = NA, xlab = "Age", ylab = "Ave. fec.")
    polygon(x = c(xd, rev(xd)), c(fecCI[, 2], rev(fecCI[, 3])), 
            col = adjustcolor("dark green", alpha.f = 0.25), border = NA)
    lines(xd, yave, col = 'dark red', lwd = 4)
    lines(xd, fecCI[, 1], col = 'dark green', lwd = 2)
    
  }
}

# A.4) Printing bafta outputs:
# ---------------------------- #
print.bafta <- function(x, ...) {
  cat("Settings:\n")
  print(x$settings)
  
  cat("Coefficients:\n")
  print(x$coefficients, ...)
  
  cat("\nDIC:\n")
  print(x$modSel)
}

# =========================== #
# B) INTERNAL FUNCTIONS: ==== 
# =========================== #
# B.1) Manage user inputs:
# ------------------------ #
# Algorithm object function:
.CreateAlgObj <- function(model, dataType, formula, niter, burnin,
                          thinning, updateJumps, nsim, minAge) {
  return(list(model = model, dataType = dataType, formula = formula, 
              niter = niter, burnin = burnin, thinning = thinning, 
              updJump = updateJumps, nsim = nsim, minAge = minAge))
}

# Prepare data object:
.CreateDataObj <- function(object, algObj) {
  n <- nrow(object)
  
}


# Prepare data object:
.CreateDataObj <- function(object, algObj) {
  classDataObj <- c("bastacmr", "ageUpd")
  # Data Object for CMR data type:
  if (algObj$dataType == "CMR") {
    dataObj <- list()
    # Extract study year sequence and length:
    dataObj$study <- algObj$start:algObj$end
    dataObj$studyLen <- length(dataObj$study)
    
    # Number of observations:
    dataObj$n <- nrow(object)
    
    # Recapture matrix:
    Y <- as.matrix(object[, 1:dataObj$studyLen + 3])
    dataObj$Y <- Y
    dataObj$Y[Y > 1] <- 1
    colnames(dataObj$Y) <- dataObj$study
    
    # NOTE: addition of censTime for studies where individuals are censored
    #       before the end of the study (e.g. two studies of different 
    #       duration).
    # Find possible times of censoring before the study end:
    censTime <- rep(algObj$end, dataObj$n)
    for (tt in 1:dataObj$studyLen) {
      idCens <- which(Y[, tt] > 1)
      censTime[idCens] <- dataObj$study[tt]
    }
    dataObj$censTime <- censTime
    
    # Birth - death matrix:
    bd <- as.matrix(object[, 2:3])
    dataObj$bi <- bd[, 1]
    dataObj$di <- bd[, 2]
    bi0 <- which(dataObj$bi == 0 | is.na(dataObj$bi))
    if (length(bi0) > 0) {
      dataObj$idNoB <- bi0
      dataObj$updB <- TRUE
      dataObj$nUpdB <- length(bi0)
    } else {
      dataObj$updB <- FALSE
      dataObj$nUpdB <- 0
    }
    di0 <- which(dataObj$di == 0 | is.na(dataObj$di))
    if (length(di0) > 0) {
      dataObj$idNoD <- di0
      dataObj$updD <- TRUE
      dataObj$nUpdD <- length(di0)
    } else {
      dataObj$updD <- FALSE
      dataObj$nUpdD <- 0
    }
    
    if (!dataObj$updB & !dataObj$updD) {
      classDataObj[2] <- "noAgeUpd"
      dataObj$updA <- FALSE
    } else {
      dataObj$idNoA <- sort(unique(c(dataObj$idNoB, dataObj$idNoD)))
      dataObj$nUpdA <- length(dataObj$idNoA)
      dataObj$updA <- TRUE
      
      # NOTE: addition of min-max birth and death (2022-05-17):
      if ("Min.Birth" %in% colnames(object)) {
        dataObj$minBirth <- object[, "Min.Birth"]
        idMinB <- which(!is.na(dataObj$minBirth[dataObj$idNoB]))
        if (length(idMinB) > 0) {
          dataObj$idMinB <- dataObj$idNoB[idMinB]
          dataObj$updMinB <- TRUE
        } else {
          dataObj$idMinB <- NA
          dataObj$updMinB <- FALSE
        }
      } else {
        dataObj$minBirth <- rep(NA, dataObj$n)
        dataObj$idMinB <- NA
        dataObj$updMinB <- FALSE
      }
      if ("Max.Birth" %in% colnames(object)) {
        dataObj$maxBirth <- object[, "Max.Birth"]
        idMaxB <- which(!is.na(dataObj$maxBirth[dataObj$idNoB]))
        if (length(idMaxB) > 0) {
          dataObj$idMaxB <- dataObj$idNoB[idMaxB]
          dataObj$updMaxB <- TRUE
        } else {
          dataObj$idMaxB <- NA
          dataObj$updMaxB <- FALSE
        }
      } else {
        dataObj$maxBirth <- rep(NA, dataObj$n)
        dataObj$idMaxB <- NA
        dataObj$updMaxB <- FALSE
      }
      if ("Min.Death" %in% colnames(object)) {
        dataObj$minDeath <- object[, "Min.Death"]
        idMinD <- which(!is.na(dataObj$maxDeath[dataObj$idNoD]))
        if (length(idMinD) > 0) {
          dataObj$idMinD <- dataObj$idNoD[idMinD]
          dataObj$updMinD <- TRUE
        } else {
          dataObj$idMinD <- NA
          dataObj$updMinD <- FALSE
        }
      } else {
        dataObj$minDeath <- rep(NA, dataObj$n)
        dataObj$idMinD <- NA
        dataObj$updMinD <- FALSE
      }
      if ("Max.Death" %in% colnames(object)) {
        dataObj$maxDeath <- object[, "Max.Death"]
        idMaxD <- which(!is.na(dataObj$maxDeath[dataObj$idNoD]))
        if (length(idMaxD) > 0) {
          dataObj$idMaxD <- dataObj$idNoD[idMaxD]
          dataObj$updMaxD <- TRUE
        } else {
          dataObj$idMaxD <- NA
          dataObj$updMaxD <- FALSE
        }
        
      } else {
        dataObj$maxDeath <- rep(NA, dataObj$n)
        dataObj$idMaxD <- NA
        dataObj$updMaxD <- FALSE
      }
      
      # 4.1.2 Calculate first and last time observed 
      #       and total number of times observed:
      ytemp <- t(t(dataObj$Y) * dataObj$study)
      dataObj$lastObs <- c(apply(ytemp, 1, max))
      ytemp[ytemp == 0] <- 10000
      dataObj$firstObs <- c(apply(ytemp, 1, min))
      dataObj$firstObs[dataObj$firstObs == 10000] <- 0
      dataObj$oi <- dataObj$Y %*% rep(1, dataObj$studyLen)
      
      # 4.1.3 Define study duration:
      dataObj$Tm <- matrix(dataObj$study, dataObj$n, 
                           dataObj$studyLen, byrow = TRUE)
      fii <- dataObj$firstObs
      id1 <- which(dataObj$bi > 0 & dataObj$bi >= algObj$start)
      fii[id1] <- dataObj$bi[id1] + 1
      fii[dataObj$bi > 0 & dataObj$bi < algObj$start]  <- algObj$start
      lii <- dataObj$lastObs
      id2 <- which(dataObj$di > 0 & dataObj$di <= dataObj$censTime)
      lii[id2] <- dataObj$di[id2] - 1
      idCens <- which(dataObj$di > 0 & dataObj$di > dataObj$censTime)
      lii[idCens] <- dataObj$censTime[idCens]
      dataObj$obsMat <- .BuildAliveMatrix(fii, lii, dataObj)
      dataObj$obsMat[lii == 0 | fii == 0, ] <- 0
      classDataObj[2] <- "ageUpd"
    }
    dataObj$Dx <- 1 
    classDataObj[1] <- "bastacmr"
  }
  
  # Data Object for Census data type:
  else {
    n <- nrow(object)
    # Calculate Julian times:
    bi <- round(as.numeric(as.Date(object$Birth.Date, 
                                   format = "%Y-%m-%d")) / 
                  365.25, 2) + 1970
    bil <- round(as.numeric(as.Date(object$Min.Birth.Date, 
                                    format = "%Y-%m-%d")) /
                   365.25, 2) + 1970
    biu <- round(as.numeric(as.Date(object$Max.Birth.Date, 
                                    format = "%Y-%m-%d")) /
                   365.25, 2) + 1970
    firstObs <- round(as.numeric(as.Date(object$Entry.Date, 
                                         format = "%Y-%m-%d")) /
                        365.25, 2) + 1970
    lastObs <- round(as.numeric(as.Date(object$Depart.Date, 
                                        format = "%Y-%m-%d")) /
                       365.25, 2) + 1970
    
    # Entry and departure types:
    entryType <- as.character(object$Entry.Type)
    departType <- as.character(object$Depart.Type)
    
    # Censored individuals:
    idCens <- which(departType %in% c("O", "C"))
    nCens <- length(idCens)
    
    # Birth times to be estimated:
    idNoBirth <- which(bi != bil | bi != biu)
    nNoBirth <- length(idNoBirth)
    if (nNoBirth == 0) {
      updB <- FALSE
      classDataObj[2] <- "noAgeUpd"
    } else {
      updB <- TRUE
      classDataObj[2] <- "ageUpd"
    }
    
    # Create data object:
    dataObj <- list(bi = bi, bil = bil, biu = biu, firstObs = firstObs, 
                    lastObs = lastObs, idCens = idCens,
                    idNoB = idNoBirth, nUpdB = nNoBirth, updB = updB, 
                    updD = FALSE, idNoD = NA, nUpdD = 0,
                    idNoA = idNoBirth, nUpdA = nNoBirth, updA = updB, n = n, 
                    nCens = nCens, studyLen = NA)
    classDataObj[1] <- "bastacensus"
  }
  class(dataObj) <- classDataObj
  
  # Output:
  return(dataObj)
}

# Create initial age object:
.CreateAgeObj <- function(dataObj, algObj) {
  ageObj <- list()
  bi <- dataObj$bi
  if (dataObj$updB & inherits(dataObj, "bastacmr")) {
    bi[dataObj$idNoB] <- dataObj$firstObs[dataObj$idNoB] - 
      sample(6:1, size = dataObj$nUpdB, replace = TRUE)
  }
  if (inherits(dataObj, "bastacmr")) {
    di <- dataObj$di
    if (dataObj$updD) {
      di[dataObj$idNoD] <- apply(cbind(bi[dataObj$idNoD], 
                                       dataObj$lastObs[dataObj$idNoD]),
                                 1, max) + 
        sample(6:1, size = dataObj$nUpdD, replace = TRUE)
    }
    age <- di - bi
    ageTr <- algObj$start - bi
    ageTr[ageTr < 0] <- 0
  } else {
    di <- dataObj$lastObs
    age <- dataObj$lastObs - bi
    ageTr <- dataObj$firstObs - bi
    ageTr[ageTr < 0] <- 0
  }
  
  # indicator for uncensored:
  indUncens <- rep(1, dataObj$n)
  if (inherits(dataObj, "bastacensus")) {
    indUncens <- rep(1, dataObj$n)
    indUncens[dataObj$idCens] <- 0
  }
  
  # Recalculate ages based on minAge:
  # ---------------------------------
  if (algObj$minAge > 0) {
    # Ages and indicator after min age:
    ageAftMa <- age - algObj$minAge
    ageAftMa[ageAftMa < 0] <- 0
    
    # Ages and indicator before min age:
    ageBefMa <- age
    indBefMa <- rep(0, dataObj$n)
    indBefMa[ageBefMa < algObj$minAge] <- 1
    ageBefMa[age >= algObj$minAge] <- algObj$minAge
    
    # Ages at truncation and indicator after min age:
    ageTrAftMa <- ageTr - algObj$minAge
    ageTrAftMa[ageTrAftMa < 0] <- 0
    
    # Ages at truncation and indicator before min age:
    ageTrBefMa <- ageTr
    ageTrBefMa[ageTr >= algObj$minAge] <- algObj$minAge
  } else {
    ageBefMa <- ageTrBefMa <- indBefMa <- rep(0, dataObj$n)
    ageAftMa <- age
    ageTrAftMa <- ageTr
  }
  
  
  # Create alive matrix for bastacmr:
  # ---------------------------------
  if (inherits(dataObj, "bastacmr")) {
    firstObs <- c(apply(cbind(algObj$start, bi + 1), 1, max))
    lastObs <- c(apply(cbind(algObj$end, dataObj$censTime, di), 1, min))
    alive <- .BuildAliveMatrix(firstObs, lastObs, dataObj)
  } else {
    alive <- NA
  }
  
  # Fill-in matrices for age object:
  # --------------------------------
  ageObj$ages <- data.frame(birth = bi, death = di, age = age, ageTr = ageTr,
                            ageAft = ageAftMa, ageBef = ageBefMa, 
                            truAft = ageTrAftMa, truBef = ageTrBefMa)
  ageObj$inds <- data.frame(ageBef = indBefMa, uncens = indUncens)
  ageObj$alive <- alive
  
  # Assign class to age object:
  # ---------------------------
  minAgeClass <- ifelse(algObj$minAge > 0, "minAge", "noMinAge")
  if (inherits(dataObj, "bastacmr")) {
    class(ageObj) <- c("agecmr", minAgeClass)
  } else {
    class(ageObj) <- c("agecensus", minAgeClass)
  }
  return(ageObj)
}


# B.2) Create and manage covariates:
# ---------------------------------- #


# B.3) Fecundity parameters:
# -------------------------- #
# Fertility:
.SetDefaultBeta <- function(beta, modelFert = "quadratic") {
  if (is.null(beta)) {
    stop("Missing 'beta' parameter vector or matrix.\n", call. = FALSE)
  } 
  if (modelFert == "quadratic") {
    nBe <- 3
    startBet <- c(-2, 0.01, 0.001)
    priorMean <- c(-2, 0.01, 0.001)
    priorSd <- c(2, 0.5, 0.5)
    lowBet <- c(0, 0, 0)
    jumpBet <- c(1, 0.2, 0.2)
    jitter <- c(1, 0.2, 0.2)
  } else if (modelFert == "PeristeraKostaki") {
    nBe <- 4
    startBet <- c(1, 20, 10, 20)
    priorMean <- c(1, 10, 10, 10)
    priorSd <- c(2, 5, 5, 5)
    lowBet <- c(0, 0, 0, 0)
    jumpBet <- c(1, 1, 1, 1)
    jitter <- c(1, 5, 5, 5)
  } else if (modelFert == "ColcheroMuller") {
    nBe <- 4
    startBet <- c(-2, 0.01, 0.001, -4)
    priorMean <-c(-2, 0.01, 0.001, -4)
    priorSd <- c(2, 0.5, 0.5, 2)
    lowBet <- c(0, 0, 0, -Inf)
    jumpBet <- c(1, 0.2, 0.2, 1)
    jitter <- c(1, 0.2, 0.2, 1)
  } else if (modelFert == "Hadwiger") {
    nBe <- 3
    lowBe <- c(0, 0, 0)
    startBet <- c(1, 2, 10)
    priorMean <-c(1, 2, 10)
    priorSd <- c(2, 2, 5)
    lowBet <- c(0, 0, 0)
    jumpBet <- c(0.2, 0.2, 1)
    jitter <- c(0.2, 0.2, 1)
  } else if (modelFert == "gamma") {
    nBe <- 3
    startBet <- c(2, 2, 0.1)
    priorMean <-c(1, 1, 0.5)
    priorSd <- c(2, 2, 1)
    lowBet <- c(0, 0, 0)
    jumpBet <- c(0.2, 0.2, 0.1)
    jitter <- c(0.5, 0.5, 0.2)
  }
  nameBet <- sprintf("b%s", 1:nBe - 1)
  if (is.matrix(beta)) {
    nbeUser <- ncol(beta)
    clBe <- "matrix"
    stBe <- "columns"
  } else if (is.numeric(beta)) {
    nbeUser <- length(beta)
    clBe <- "vector"
    stBe <- "elements"
  } else {
    stop("The beta parameters should either be of class matrix", 
         " or a numeric vector.\n", call. = FALSE)
  }
  if (nbeUser != nBe) {
    stop(sprintf("The beta %s should have %s %s.\n", clBe, nBe, stBe),
         call. = FALSE)
  } else {
    if (is.matrix(beta)) {
      colnames(beta) <- nameBe
    } else {
      names(beta) <- nameBe
    }
  }
  # check if beta parameters conform to their support:
  if (is.matrix(beta)) {
    BETLOW <- all(sapply(1:nBe, function(bi) {
      bl <- all(beta[, bi] >= lowBe[bi])
    }))
  } else {
    BETLOW <- all(beta >= lowBe)
  }
  if(!BETLOW) {
    stop(sprintf("Some beta parameters are below their lower bound.\n %s.\n",
                 paste(sprintf("min(%s) = %s", nameBe, lowBe), 
                       collapse = ", ")),
         call. = FALSE)
    
  }
  defaultBeta  <- list(length = nTh, start = startTh, jump = jumpTh, 
                       priorMean = priorMean, priorSd = priorSd, name = nameTh, 
                       low = lowTh, jitter = jitter)
  attr(defaultBeta, "model") = algObj$model
  return(defaultBeta)
}

# B.4) Fecundity models:
# ---------------------- #
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
      be1 <- x * 0 + beta["b1a"]
      be1[which(x > beta["b2"])] <- beta["b1b"]
      fert <- beta["b0"] * exp(-((x - beta["b2"]) / be1)^2)
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
  } 
  return(fertfun)
}

# Calculate Fertility:
.CalcFert <- function(beta, x, modelFert = "quadratic", checkBeta = TRUE) {
  
  # Verify fertility model:
  .VerifyFertMod(modelFert = modelFert)
  
  # Extract beta attributes:
  if (checkBeta) {
    betaAttr <- .SetBeta(beta, modelFert = modelFert)
    beta <- betaAttr$beta
  }
  
  # a) Fertility method:
  .CalcFert <- function(beta, ...) UseMethod(".CalcFert")
  .CalcFert.matrix <- .DefineFertilityMatrix(modelFert = modelFert)
  .CalcFert.numeric <- .DefineFertilityNumeric(modelFert = modelFert)
  
  # b) Calculate Fertility:
  fertfun <- .CalcFert(beta, x)
  return(fertfun)
}


# Binomial or Poisson likelihood:
like <- function(ga, sig, u, xi, yi, zi, np, binom = TRUE, fec, 
                 randEffs = FALSE, forPars = TRUE) {
  fi <- fec(xi, ga) * c(exp(zi %*% u))
  if (binom) {
    lka <- dbinom(yi, 1, fi, log = TRUE)
  } else {
    lka <- dpois(yi, fi, log = TRUE)
  }
  if (randEffs) {
    lkr <- dnorm(u, mean = 0, sd = sig, log = TRUE)
  } else {
    lkr <- 0
  }
  if (forPars) {
    lk <- sum(lka) + sum(lkr)
  } else {
    lk <- c(t(lka * zi) %*% rep(1, n)) + lkr
  }
  return(list(full = lk, a = sum(lka)))
}

# Calculate age at maximum fecundity:
CalcAgeMaxFec <- function(b, fecmod = 1, xstart = 10, minAge = 0) {
  # Extract fecundity model and parameters:
  feclist <- fecMods[[sprintf("mod%s", fecmod)]]
  fec <- feclist$fun
  np <- feclist$np
  bfull <- b
  if (fecmod == 1) {
    b <- bfull[-1]
    an <- rep(0, 10)
    x0 <- xstart
    dd <- 1
    ii <- 0
    while (dd > 0.0001 & ii <= 25) {
      ii <- ii + 1
      x <- x0
      x1 <- x - (x^3 - x^2 * ((b[1] - 4 * b[2]) / (2 * b[2])) - 
                   x * ((b[1] - b[2]) / b[2]) - (b[1] - b[3]) / (2 * b[2])) / 
        (3 * x^2 - x * 2 * ((b[1] - 4 * b[2]) / (2 * b[2])) - 
           ((b[1] - b[2]) / b[2]))
      dd <- abs(x1 - x0)
      x0 <- x1
    }
    if (dd < 0.0001) {
      xm <- x1
      fecmax <- fec(xm, bfull[1:np])
      xm <- xm + minAge
    } else {
      xm <- NA
      fecmax <- NA
    }
  } else if (fecmod == 0) {
    b <- bfull[-1]
    xm <- b[1] / (2 * b[2]) 
    fecmax <- fec(xm, bfull[1:np])
    xm <- xm + minAge
    dd <- NA
    ii <- NA
  } else if (fecmod == 3) {
    dd <- 1
    dx <- 0.1
    xv <- seq(0, b[3] * 5, dx)
    ii <- 0
    while (dd > 0.000001 & ii <= 10) {
      ii <- ii + 1
      ddv <- abs(-2 * b[1] * b[2] * (xv - b[3]) * exp(-b[2] * (xv - b[3])^2) * 
                  (1 + exp(-b[4] * (xv - b[5]))) + b[1] * b[4] * 
                  exp(-b[2] * (xv - b[3])^2) *
                  exp(-b[4] * (xv - b[5])))
      idxmax <- which(ddv == min(ddv))
      x1 <- xv[idxmax]
      dd <- ddv[idxmax]
      dx <- dx / ii
      xv <- x1 + seq(- 3 * dx, 3 * dx, dx)
    }
    if (dd < 0.000001) {
      xm <- x1
      fecmax <- fec(xm, bfull[1:np])
      xm <- xm + minAge
    } else {
      xm <- NA
      fecmax <- NA
    }
  } else {
    xm <- NA
    fecmax <- NA
    dd <- NA
    ii <- NA
  }
  return(c(xm = xm, fmax = fecmax, dd = dd, niter = ii))
}

# Truncated normal:
.rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  ru <- runif(n, Flow, Fup)
  rx <- qnorm(ru, mean, sd)
  return(rx)
}

.dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

# Function to update jumps in MCMC:
.UpdateJumps <- function(jumps, updMat, iter, iterUpd = 100, updTarg = 0.25) {
  updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd  
  updRate[updRate == 0] <- 1e-2
  jumps <- jumps * updRate / updTarg
  return(jumps)
}

# MCMC function:
.RunMCMC <- function(sim, fecmod = 1, niter = 11000, burnin = 1001, 
                    updJump = TRUE, jumps = NULL, randEffs = FALSE, 
                    jumpInt = 50) {
  # Restart the random seed for each core:
  if (sim > 1) {
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
  }
  
  # Extract fecundity model and parameters:
  feclist <- fecMods[[sprintf("mod%s", fecmod)]]
  fec <- feclist$fun
  np <- feclist$np
  
  # Fecundity parameter names:
  parnames <- sprintf("a%s", 1:np)
  
  # Priors:
  parmean <- feclist$mean
  parsd <- rep(20, np)
  
  # Lower bound of parameter support:
  parlow <- feclist$low
  
  # Initial values:
  parini <- feclist$ini
  
  # Random effects parameters:
  if (randEffs) {
    signow <- 0.1
    unow <- rnorm(ni, mean = 0, sd = signow)
    s1 <- 0.01
    s2 <- 0.1
  } else {
    signow <- 0
    unow <- 0
    zi <- 0
  }
  
  # Jump sds:
  if (updJump) {
    parjump <- rep(0.1, np)
    # Jump matrix setup:
    parjumpMat <- parjump
    updMat <- matrix(0, nrow = burnin, ncol = np)
  } else {
    parjump <- jumps
    parjumpMat <- NA
  }
  
  # ga, u, xi, yi, zi, binom = TRUE, fec
  # Start parameters and calculate conditional posterior:
  parnow <- rtnorm(1, mean = parini, sd = abs(parini) * 0.2, lower = parlow)
  likenow <- like(parnow, signow, unow, xi, yi, zi, np, binom = binom, fec, 
                  randEffs = randEffs)
  postnow <- likenow$a + sum(dtnorm(parnow, parmean, parsd, lower = parlow, 
                                  log = TRUE)) 
  while (postnow == -Inf | is.na(postnow)) {
    parnow <- rtnorm(1, mean = parini, sd = abs(parini) * 0.2, lower = parlow)
    likenow <- like(parnow, signow, unow, xi, yi, zi, np, binom = binom, fec, 
                    randEffs = randEffs)
    postnow <- likenow$a + 
      sum(dtnorm(parnow, parmean, parsd, lower = parlow, log = TRUE))
  }
  
  # Results matrices:
  parmat <- matrix(0, niter, np, dimnames = list(NULL, parnames))
  parmat[1, ] <- parnow
  postmat <- matrix(0, niter, 3, dimnames = list(NULL, c("Post", "Likefull",
                                                         "likea")))
  postmat[1, ] <- c(postnow, likenow$full, likenow$a)
  if (randEffs) {
    postUnow <- like(parnow, signow, unow, xi, yi, zi, np, binom = binom, fec, 
                     forPars = FALSE)$full
    umat <- matrix(0, niter, ni)
    umat[1, ] <- unow
    sigvec <- rep(0, niter); sigvec[1] <- signow
  } else {
    umat <- NA
    sigvec <- NA
  }
  
  # ==============
  # 4) MCMC run
  # ============
  for (iter in 2:niter) {
    for (pp in 1:(np)) {
      parnew <- parnow
      parnew[pp] <- rtnorm(n = 1, mean = parnow[pp], sd = parjump[pp], 
                           lower = parlow[pp])
      likenew <- like(parnew, signow, unow, xi, yi, zi, np, binom = binom, fec, 
                      randEffs = randEffs) 
      postnew <- likenew$a + sum(dtnorm(parnew, parmean, parsd, 
                                           lower = parlow, log = TRUE)) 
      if (!is.na(postnew)) {
        hastRatio <- dtnorm(parnow[pp], mean = parnew[pp], sd = parjump[pp], 
                            lower = parlow[pp]) / 
          dtnorm(parnew[pp], mean = parnow[pp], sd = parjump[pp], 
                 lower = parlow[pp])
        acceptRatio <- exp(postnew - postnow) * hastRatio
        if (acceptRatio > runif(1)) {
          parnow <- parnew
          likenow <- likenew
          postnow <- postnew
          if (updJump) {
            if (iter < burnin) updMat[iter, pp] <- 1
          }
        }
      }
    }
    
    # Random effects:
    if (randEffs) {
      # Sample U values:
      postUnow <- like(parnow, signow, unow, xi, yi, zi, np, binom = binom, 
                       fec, forPars = FALSE, randEffs = randEffs)$full
      unew <- rnorm(ni, mean = unow, sd = 0.5)
      postUnew <- like(parnow, signow, unew, xi, yi, zi, np, binom = binom, 
                       fec, forPars = FALSE, randEffs = randEffs)$full 
      acceptRatio <- exp(postUnew - postUnow)
      ranU <- runif(ni)
      idUpd <- which(acceptRatio > ranU)
      unow[idUpd] <- unew[idUpd]
      postUnow[idUpd] <- postUnew[idUpd]
      likenow <- like(parnow, signow, unow, xi, yi, zi, np, binom = binom, fec, 
                      randEffs = randEffs)
      postnow <- likenow$a + 
        sum(dtnorm(parnow, parmean, parsd, lower = parlow, log = TRUE))
      umat[iter, ] <- unow
      
      # Sample sigma:
      u1 <- s1 + ni / 2
      u2 <- s2 + 0.5 * sum(unow^2) 
      signow <- 1 / rgamma(1, u1, u2)
      sigvec[iter] <- signow
    } 
    
    # Fill up output matrices:
    parmat[iter, ] <- parnow
    postmat[iter, ] <- c(postnow, likenow$full, likenow$a)
    
    # Update jumps:
    if (updJump) {
      if(iter %in% c(1:((burnin - 1) / jumpInt) * jumpInt)) {
        parjump <- .UpdateJumps(parjump, updMat, iter, jumpInt, 0.25)
        parjumpMat <- rbind(parjumpMat, parjump)
        if (iter == burnin) {
          idjmean <- round(nrow(parjumpMat) / 2):nrow(parjumpMat)
          parjump <- apply(parjumpMat[idjmean, ], 2, mean)
        }
      }
    }
  }
  output <- list(params = parmat, postLike = postmat, jumps = parjump, 
                 jumpmat = parjumpMat, umat = umat, sigma = sigvec, 
                 settings = c(fecmod = fecmod, randEffs = randEffs))
  class(output) <- "baftamcmc"
  return(output)
}


# Function to extract data from parallel MCMC:
ExtractMCMC <- function(out) {
  
  # output list:
  outlist <- list()
  fecmod <- out[[1]]$settings["fecmod"]
  randEffs <- out[[1]]$settings["randEffs"]
  
  # Find values to keep for summary statistics:
  keep <- ceiling(seq(burnin, niter, thinning))
  
  # Extract fecundity model:
  fec <- fecMods[[sprintf("mod%s", fecmod)]]$fun
  np <- fecMods[[sprintf("mod%s", fecmod)]]$np
  
  # Extract MCMC parameter values:
  parmatfull <- out[[1]]$params
  if (randEffs == 1) parmatfull <- cbind(parmatfull, sigma = out[[1]]$sigma)
  parmat <- parmatfull[keep, ]
  likePost <- out[[1]]$postLike[keep, ]
  for (i in 2:nsim) {
    pmfulli <- out[[i]]$params
    if (randEffs == 1) pmfulli <- cbind(pmfulli, sigma = out[[i]]$sigma)
    parmatfull <- rbind(parmatfull, pmfulli)
    parmat <- rbind(parmat, pmfulli[keep, ])
    likePost <- rbind(likePost, out[[i]]$postLike[keep, ])
  }
  
  # Extract results:
  # Coefficients and their SE and CIs:
  coeffs <- cbind(Mean = apply(parmat, 2, mean), SE = apply(parmat, 2, sd),
                  "2.5%" = apply(parmat, 2, quantile, 0.025),
                  "97.5%" = apply(parmat, 2, quantile, 0.975))
  
  # Fecundity CIs:
  fecMat <- t(apply(parmat, 1, function(aa) fec(xd, aa)))
  fecCI <- cbind(Mean = apply(fecMat, 2, mean), 
                 "2.5%" = apply(fecMat, 2, quantile, 0.025),
                 "97.5%" = apply(fecMat, 2, quantile, 0.975))
  
  if (nsim > 1) {
    nthin <- length(keep)
    idSims <- rep(1:nsim, each = nthin)
    Means <- apply(parmat, 2, function(x) 
      tapply(x, idSims, mean))
    Vars <- apply(parmat, 2, function(x) 
      tapply(x, idSims, var))
    meanall <- apply(Means, 2, mean)
    B <- nthin / (nsim - 1) * apply(t((t(Means) - meanall)^2), 2, sum)
    W <- 1 / nsim * apply(Vars, 2, sum)
    Varpl <- (nthin - 1) / nthin * W + 1 / nthin * B
    Rhat <- sqrt(Varpl / W)
    Rhat[Varpl==0] <- 1
    conv <- cbind(B, W, Varpl, Rhat)
    rownames(conv) <- colnames(parmat)
    coeffs <- cbind(coeffs, conv[, 'Rhat'])
    colnames(coeffs) <- c(colnames(coeffs)[-ncol(coeffs)], "PotScaleReduc")
    idnconv <- which(conv[, 'Rhat'] > 1.1)
    if (length(idnconv) == 0) {
      # DIC:
      Dave <- mean(- 2 * likePost[, 'likea'])
      pD <- 1/2 * var(-2 * likePost[, 'likea'])
      DIC <- pD + Dave
      Dmode <- Dave - 2 * pD
      k <- np
      modSel <- c(Dave, Dmode, pD, k, DIC)
      names(modSel) <- c("D.ave", "D.mode", "pD", "k", "DIC")
      convmessage <- "All parameters converged properly.\n"
    } else {
      modSel <- NA
      convmessage <- "Convergence not reached for some parameters.\n"
    }
  } else {
    modSel <- NA
    convmessage <- "Convergence not calculated due to\ninsuficcient number of simulations.\n"
  }
  
  # Age of maximum fecundity:
  xstart <- coeffs["b2", 1] / (2 * coeffs["b3", 1])
  xstart <- coeffs["b2", 1] / (2 * coeffs["b3", 1])
  xmvec <- t(apply(parmat, 1, CalcAgeMaxFec, fecmod = fecmod, xstart = xstart,
                   minAge = minAge))
  idnotkeep <- lapply(1:4, function(ii) {
    pv <- parmat[, ii]
    pci <- coeffs[ii, c("2.5%", "97.5%")]
    idout <- which(pv < pci[1] | pv > pci[2])
    return(idout)
  })
  
  idnot <- unique(unlist(idnotkeep))
  xmq <- t(apply(xmvec[-idnot, c(1, 2)], 2, function(xx) {
    c(Mean = mean(xx, na.rm = TRUE), 
      SE = sd(xx, na.rm = TRUE), 
      quantile(xx, c(0.025, 0.975), na.rm = TRUE))}))
  rownames(xmq) <- c("Age", "MaxFec")
  
  # Outputs:
  outlist$coefficients <- coeffs
  outlist$fecICs <- fecCI
  outlist$modSel <- modSel
  outlist$convergence <- conv
  outlist$settings <- c(niter = niter, burnin = burnin,
                        thinning = thinning, nsim = nsim, fecmod,
                        np = np, randEffs)
  outlist$avefec <- sapply(xd, function(xx) {
    idx <- which(xii == xx)
    ymean <- mean(yii[idx])
    return(ymean)
  })
  outlist$maxFec <- xmq
  outlist$jumpsd <- out[[1]]$jumps
  names(outlist$jumpsd) <- colnames(parmat)[colnames(parmat) != "sigma"]
  outlist$keep <- keep
  outlist$parmat <- parmat
  outlist$parmatfull <- parmatfull
  outlist$likepost <- likePost
  if (randEffs == 1) {
    umat <- out[[1]]$umat[keep, ]
    for (sim in 2:nsim) {
      umat <- rbind(umat, out[[sim]]$umat[keep, ])
    }
    usumar <- cbind(Mean = apply(umat, 2, mean), 
                    "2.5%" = apply(umat, 2, quantile, 0.025),
                    "97.5%" = apply(umat, 2, quantile, 0.975))
    outlist$randeff <- usumar
  } else {
    outlist$randeff <- NA
  }
  class(outlist) <- "baftaOut"
  return(outlist)
}


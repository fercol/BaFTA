# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2023-02-21
# DESCRIPTION: Code to simulate different levels of age-specific fertility 
#              datasets to test with package.
# ================================ CODE START ================================ #
# ======================== #
# ==== GENERAL SETUP: ====
# ======================== #
# Load libraries:
library(paramDemo)

# ==================== #
# ==== FUNCTIONS: ====
# ==================== #
# Function to construct Leslie matrix:
FillMatr <- function(p, f, n) {
  idcol <- 1:n
  idrow <- c(2:n, n)
  idpa <- (idcol - 1) * n + idrow
  Aa <- matrix(0, n, n)
  Aa[1, ] <- f
  Aa[idpa] <- p
  return(Aa)
}

# Likelihood function:
like <- function(params) {
  params <- abs(params)
  bet <- params[1:3]
  et <- params[4]
  ga <- params[5]
  gxl <- CalcFert(beta = bet, x = xi - alpha, 
                  modelFert = "gamma")
  lk <- dpois(yit, lambda = gxl, log = TRUE) + 
    dexp(zit, rate = et, log = TRUE) * (1 - first) + 
    dexp(wit, rate = ga, log = TRUE) * first
  return(-sum(lk))
}

# ============================== #
# ==== PARAMETRIC SETTINGS: ====
# ============================== #
# Initial population size:
Nini <- 100

# Number of years to simulate:
Nyears <- 44
Tmax <- Nyears

# Siler mortality parameters:
theta <- c(a0 = -1, a1 = 2, c = 0.0001, b0 = -6, b1 = 0.15)

# Gamma fertility parameters:
beta <- c(b0 = 4, b1 = 2, b2 = 0.15)

# Eta inter-birth interval parameter:
eta <- 2

# Gamma parameter for age at first birth:
gamma <- 2

# Age at maturity:
alpha <- 4

# Gestation and weaning period:
gest <- 0.5
wean <- 0.5
tau <- gest + wean

# Vector of ages:
dx <- 0.001
x <- seq(0, 100, dx)
z <- seq(0, 5, dx)
w <- seq(0, 5, dx)

# Standard deviation of individual effects:
uSd <- 0.25

# ============================ #
# ==== VISUAL INSPECTION: ====
# ============================ #
# Calculate mortality and survival:
Sx <- CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")
idsx <- which(Sx >= 0.001)
xv <- x[idsx]
Sx <- Sx[idsx]
mux <- CalcMort(theta = theta, x = xv, model = "GO", shape = "bathtub")

# Calculate age-specific fertility:
gx <- CalcFert(beta = beta, x = xv[which(xv >= alpha)] - alpha, 
               modelFert = "gamma")

# Calculate interbirth interval (after gestation and weaning):
fz <- dexp(z, rate = eta)

# Calculate interbirth interval (after gestation and weaning):
fw <- dexp(w, rate = gamma)

# Age at minimum mortality:
xmin <- xv[which(mux == min(mux))]

# Plot vital rates:
par(mfrow = c(3, 2), mar = c(4, 4, 1, 1))
plot(xv, mux, type = 'l', col = 'dark red', lwd = 4, xlab = "", 
     ylab = "Mortality")
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)
plot(xv, Sx, type = 'l', col = 'dark red', lwd = 4, xlab = "", 
     ylab = "Survival")
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)
plot(xv[which(xv >= alpha)], gx, type = 'l', col = 'dark red', lwd = 4, 
     xlab = "Age", ylab = "Fertility", xlim = c(0, max(xv)), 
     ylim = c(0, max(gx)*exp(uSd)))
for (ii in c(-1, 1)) {
  lines(xv[which(xv >= alpha)], gx * exp(uSd * ii), col = 'dark red')
}
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)

plot(w + alpha, fw, type = 'l', col = 'dark red', lwd = 4, 
     xlim = c(0, max(w) + alpha), xlab = "Age at first birth", 
     ylab = "PDF age first birth")
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)

plot(z + tau, fz, type = 'l', col = 'dark red', lwd = 4, 
     xlim = c(0, max(z) + tau), xlab = "Interbirth interval", ylab = "PDF IBI")
abline(v = gest, col = 'red')
abline(v = tau, col = 'red')

# =================================== #
# ==== SIMULATE POPULATION DATA: ====
# =================================== #
# ---------------------- #
# Simulate initial pop.:
# ---------------------- #
# Transition matrix:
Sx <- CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")
gx <- x * 0
gx[which(x >= alpha)] <- CalcFert(beta = beta, 
                                  x = x[which(x >= alpha)] - alpha, 
                                  modelFert = "gamma")
xd <- 0:ceiling(max(xv))
idx <- which(x %in% xd)
idx2 <- which(x %in% c(xd + 1))
px <- Sx[idx2] / Sx[idx]
A <- FillMatr(p = px, f = gx[idx], n = length(px))

# Eigen analysis:
eiA <- eigen(x = A)

# Asymptotic population growth rate:
lambda <- Re(eiA$values[1])

# Right-eigen vector:
v <- Re(eiA$vectors[, 1])
v <- abs(v) / sum(abs(v))

# Plot asymptotic age-structure vs survival:
par(mfrow = c(1, 1), mar = c(4, 4, 1, 4))
plot(x[idx], v, type = 'l', xlab = "Age", ylab = "Assymptotic age-structure",
     lwd = 4, axes = FALSE)
lines(xd, Sx[idx] * v[1], col = 2, lwd = 4)
axis(side = 4, at = seq(0, 1, 0.2) * v[1], labels = seq(0, 1, 0.2), col = 2,
     col.axis = 2, las = 2)
Axis(range(v), side = 2, las = 2)
mtext("Survival", side = 4, line = 2, col = 2)
Axis(range(xd), side = 1)

# Cumulative age structure:
cumw <- cumsum(v) / sum(v)
plot(xd, cumw)

# Simulate ages at death:
ri <- runif(n = Nini, min = min(cumw), max = max(cumw))
ageIni <- xd[findInterval(ri, cumw)]

# Plot histogram of X:
par(mfrow = c(1, 1))
hist(ageIni, freq = FALSE)

# ---------------------------- #
# Simulate single life course:
# ---------------------------- #
# ages for which Sx in [0, 1]:
idsx <- 1:which(Sx == 0)[1]
xpop <- x[idsx]
Sxpop <- Sx[idsx]
Fxpop <- 1 - Sx
idgx <- which(x >= alpha & which(Sx > 0)[1] + 1)
xg <- x[idgx]
gxpop <- gx[idgx]

# Death date:
Xini <- sapply(1:Nini, function(i) {
  idi <- which(xpop >= ageIni[i])
  ri <- runif(n = 1, Fxpop[idi[1]], 1)
  xii <- xpop[idi[findInterval(ri, Fxpop[idi])]]
  return(xii)
})

# Birth date:
Bini <- 0 - ageIni
Dini <- Bini + Xini

# Longevity table:
longdat <- data.frame(ID = 1:Nini, birthDate = Bini, entryDate = 0, 
                      departDate = Dini, departType = rep("D", Nini), 
                      ageEntry = ageIni, ageDeath = Xini, 
                      parent = rep(NA, Nini))

# Start individual counter:
i <- 0

# Logical to continue with loop:
iCONT <- TRUE
Start <- Sys.time()
while (iCONT) {
  # Update ind. counter:
  i <- i + 1
  
  if (nrow(longdat) < i) {
    iCONT <- FALSE
  } else {
    # Extract dates and ages for ind. i:
    Ai <- longdat$ageEntry[i]
    if (Ai > alpha) {
      ti <- longdat$entryDate[i]
      xit <- Ai
    } else {
      ti <- longdat$entryDate[i] + alpha - Ai
      xit <- alpha
    }
    Bi <- longdat$birthDate[i]
    Di <- longdat$departDate[i]
    Xi <- longdat$ageDeath[i]
    ibi <- NA
    
    if (ti <= Tmax) {
      # Start reproductive events counter:
      irep <- 0
      
      while (xit < Xi & ti < Tmax) {
        # Age at time ti:
        if (xit == alpha) {
          wi <- rexp(n = 1, rate = gamma)
          xit <- xit + wi
          ti <- ti + wi
          r1i <- 1
        } else {
          ibi <- tau + rexp(n = 1, rate = eta)
          xit <- xit + ibi
          ti <- ti + ibi
          r1i <- 0
        }
        
        if (xit < Xi) {
          # Birth event:
          idgxi <- findInterval(xit, xg)
          yit <- rpois(n = 1, lambda = gxpop[idgxi])
          reprtemp <- data.frame(ID = i, age = xit, offsp = yit, IBI = ibi, 
                                 time = ti, stage = "A", event = "R", 
                                 first = r1i)
          irep <- irep + 1
          if (i == 1 & irep == 1) {
            reprdat <- reprtemp
          } else {
            reprdat <- rbind(reprdat, reprtemp)
          }
          
          # Offspring longevity:
          if (yit > 0) {
            for (j in 1:yit) {
              Bj <- ti
              Ej <- ti
              rj <- runif(n = 1)
              Xj <- xpop[findInterval(rj, Fxpop)]
              Dj <- Bj + Xj
              if (Dj > Tmax) {
                Dj <- Tmax
                depj <- "C"
              } else {
                depj <- "D"
              }
              longtemp <- data.frame(ID = max(longdat$ID) + 1, 
                                     birthDate = Bj, entryDate = Ej, 
                                     departDate = Dj, departType = depj, 
                                     ageEntry = 0, ageDeath = Xj, parent = i)
              longdat <- rbind(longdat, longtemp)
            }
          }
        }
      }
      # Stop loop if all remaining individuals were born after Tmax - alpha:
      iCONT <- !all(longdat$birthDate[i:nrow(longdat)] > Tmax - alpha) 
    }
  }
}
End <- Sys.time()

print(End - Start)

# ========================================= #
# ==== EXTRACT POPULATION INFORMATION: ====
# ========================================= #
# Year vector:
years <- 0:(Nyears - alpha)
nt <- length(years)

# Population size:
Nt <- sapply(years, function(it) {
  length(which(longdat$entryDate < it + 1 & longdat$departDate >= it))
})

# Sizes by age:
agev <- 0:max(ceiling(longdat$ageDeath))
nage <- length(agev)
Mt <- matrix(0, nage, nt, dimnames = list(as.character(agev), NULL))
for (it in 1:nt) {
  idt <- which(longdat$entryDate <= years[it] + 1 & 
                 longdat$departDate > years[it])
  aget <- floor(years[it] - longdat$birthDate[idt])
  idb <- which(longdat$birthDate[idt] > years[it])
  aget[idb] <- floor(longdat$birthDate[idt[idb]] - years[it])
  taget <- table(aget)
  Mt[names(taget), it] <- c(taget)
}

# Realized age-specific fertility:
ageAds <- agev[which(agev >= alpha)]
nagead <- length(ageAds)

parents <- allParents <- offspring <- ageAds * 0
for (ia in 1:nagead) {
  ida <- which(reprdat$age >= ageAds[ia] & reprdat$age <= ageAds[ia] + 1 &
                 reprdat$event == "R")
  offspring[ia] <- sum(reprdat$offsp[ida])
  idp <- which(longdat$ageEntry <= ageAds[ia] + 1 &
                 (longdat$departDate - longdat$birthDate) >= ageAds[ia])
  allParents[ia] <- length(idp)
  parents[ia] <- length(unique(reprdat$ID[ida]))
}

# ======================= #
# ==== PLOT RESULTS: ====
# ======================= #
# ---------------------------- #
# ---- Population trends: ----
# ---------------------------- #
par(mfrow = c(2, 1))
plot(years, Nt, type = 'l', col = 'dark red', xlab = "Time", ylab = "N")

# Population Growth:
plot(years[-1], Nt[-1] / Nt[-nt], type = 'l', col = 'dark red', xlab = "Time",
     ylab = "Population growth")
abline(h = lambda, col = 'orange', lwd = 4, lty = 2)


# ------------------------------------ #
# ----- Changes in Age structure: ----
# ------------------------------------ #
Sxt <- Sx[which(x %in% agev)]
xt <- x[which(x %in% agev)]
pSxt <- Sxt / sum(Sxt)
par(mfrow = c(5, 1))
for (it in c(1, 10, 20, 30, 40)) {
  plot(agev, Mt[, it] / Nt[it], type = 'h', col = 'dark red', lwd = 4, 
       xlab = "", ylab = "")
  lines(xt, pSxt , col = 'orange', lwd = 2)
  lines(x[idx], v, col = 'grey50', lwd = 2)
}

# ------------------------------ #
# ---- Realized fertility: -----
# ------------------------------ #
dage <- diff(ageAds[1:2])
mar <- c(4, 4, 1, 1)
layout(mat = cbind(c(1, 2), c(1, 3)), widths = c(1, 1), heights = c(1, 0.5))

# Observed vs actual fertility:
par(mar = mar)
plot(ageAds + dage / 2, offspring / parents, type = 'b', pch = 19, 
     col = 'dark red', lty = 2, xlab = "Age", 
     ylab = "Average number of offspring")
lines(ageAds + dage / 2, offspring / allParents, type = 'b', 
      col = 'dark green', pch = 15, lty = 2)
lines(x, gx, col = 'orange', lwd = 2)
legend("topright", legend = c("Theoretical", "Only avail. parents", "Realized"),
       pch = c(NA, 19, 15), col = c("orange", "dark red", "dark green"),
       lwd = c(2, 1, 1), lty = c(1, 2, 2), bty = 'n')

# IBI
hist(reprdat$IBI[which(reprdat$first == 0)] - tau, freq = FALSE,  
     xlab = "IBI", ylab = "Density", main = "IBI")
lines(z, fz, col = 'dark red', lwd = 2)

# First birth
hist(reprdat$age[which(reprdat$first == 1)] - alpha, freq = FALSE,  
     xlab = "Time to first birth", ylab = "Density", main = "First birth")
lines(w, fw, col = 'dark red', lwd = 2)


# ======================================= #
# ==== MAXIMUM LIKELIHOOD INFERENCE: ====
# ======================================= #
# Subset data:
nr <- nrow(reprdat)
indID <- unique(reprdat$ID)
nind <- length(indID)
ndat <- 200
sampid <- sample(x = 1:nind, size = ndat, replace = FALSE)
idsub <- which(reprdat$ID %in% indID[sampid])
idlon <- which(longdat$ID %in% indID[sampid])

parentSub <- allParentSub <- offspringSub <- ageAds * 0
for (ia in 1:nagead) {
  ida <- which(reprdat$age[idsub] >= ageAds[ia] & 
                 reprdat$age[idsub] <= ageAds[ia] + 1)
  offspringSub[ia] <- sum(reprdat$offsp[idsub][ida])
  idp <- which(longdat$ageEntry[idlon] <= ageAds[ia] + 1 &
                 c(longdat$departDate - longdat$birthDate)[idlon] >= ageAds[ia])
  allParentSub[ia] <- length(idp)
  parentSub[ia] <- length(unique(reprdat$ID[idsub][ida]))
}

# Reproductive output:
xi <- reprdat$age[idsub]
yit <- reprdat$offsp[idsub]

# IBI:
zit <- reprdat$IBI[idsub] - tau
zit[which(is.na(zit))] <- 0

# First birth time:
wit <- reprdat$age[idsub] - alpha
first <- reprdat$first[idsub]

# Original parameters:
params <- c(beta, eta, gamma)

# Starting parameters:
startpars <- c(b0 = 2, b1 = 1, b2 = 0.2, eta = 1.5, gamma = 1)

# Compare likelihood:
like(params)
like(startpars)

# Fit likelihood function:
out <- optim(startpars, like)
params
out$par

dage <- diff(ageAds[1:2])
mar <- c(4, 4, 1, 1)

# Plot results:
layout(mat = cbind(c(1, 2), c(1, 3)), widths = c(1, 1), heights = c(1, 0.5))

# Observed vs actual fertility:
par(mar = mar)
plot(ageAds + dage / 2, offspring / parents, type = 'b', pch = 19, 
     col = adjustcolor('dark red', alpha.f = 0.35), lty = 2, xlab = "Age", 
     ylab = "Average number of offspring")
lines(ageAds + dage / 2, offspring / allParents, type = 'b', 
      col = adjustcolor('dark green', alpha.f = 0.35), pch = 15, lty = 2)

lines(ageAds + dage / 2, offspringSub / parentSub, type = 'b', pch = 19, 
      col = 'dark red', lty = 2)
lines(ageAds + dage / 2, offspringSub / allParentSub, type = 'b', 
      col = 'dark green', pch = 15, lty = 2)
lines(x, gx, col = 'orange', lwd = 2)
lines(x + alpha, CalcFert(beta = out$par[1:3], x = x, modelFert = "gamma"), 
      col = "purple", lwd = 2)
legend("topright", legend = c("Theoretical", "Only avail. parents", "Realized",
                              "Estimated"),
       pch = c(NA, 19, 15, NA), col = c("orange", "dark red", "dark green",
                                                "purple"),
       lwd = c(2, 1, 1, 2), lty = c(1, 2, 2, 1), bty = 'n')

# IBI
hist(reprdat$IBI[which(reprdat$first == 0)] - tau, freq = FALSE,  
     xlab = "IBI", ylab = "Density", main = "IBI")
lines(z, fz, col = 'dark red', lwd = 2)
lines(z, dexp(x = z, rate = out$par[4]), col = 'purple', lwd = 2)

# First birth
hist(reprdat$age[which(reprdat$first == 1)] - alpha, freq = FALSE,  
     xlab = "Time to first birth", ylab = "Density", main = "First birth")
lines(w, fw, col = 'dark red', lwd = 2)
lines(w, dexp(x = w, rate = out$par[5]), col = 'purple', lwd = 2)

print(cbind(Real = params, Est = signif(out$par, 3)))

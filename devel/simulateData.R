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

# Beta parameter for age at first birth:
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
z <- seq(0, 10, dx)

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

# Age at minimum mortality:
xmin <- xv[which(mux == min(mux))]

# Plot vital rates:
par(mfrow = c(4, 1), mar = c(4, 4, 1, 1))
plot(xv, mux, type = 'l', col = 'dark red', lwd = 4, xlab = "", 
     ylab = "Mortality")
abline(v = xmin, col = 'orange')
plot(xv, Sx, type = 'l', col = 'dark red', lwd = 4, xlab = "", 
     ylab = "Survival")
abline(v = xmin, col = 'orange')
plot(xv[which(xv >= alpha)], gx, type = 'l', col = 'dark red', lwd = 4, 
     xlab = "Age", ylab = "Fertility", xlim = c(0, max(xv)), 
     ylim = c(0, max(gx)*exp(uSd)))
for (ii in c(-1, 1)) {
  lines(xv[which(xv >= alpha)], gx * exp(uSd * ii), col = 'dark red')
}
abline(v = xmin, col = 'orange')

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
w <- Re(eiA$vectors[, 1])
w <- abs(w) / sum(abs(w))

# Plot assymptotic age-structure vs survival:
par(mfrow = c(1, 1), mar = c(4, 4, 1, 4))
plot(x[idx], w, type = 'l', xlab = "Age", ylab = "Assymptotic age-structure",
     lwd = 4, axes = FALSE)
lines(xd, Sx[idx] * w[1], col = 2, lwd = 4)
axis(side = 4, at = seq(0, 1, 0.2) * w[1], labels = seq(0, 1, 0.2), col = 2,
     col.axis = 2, las = 2)
Axis(range(w), side = 2, las = 2)
mtext("Survival", side = 4, line = 2, col = 2)
Axis(range(xd), side = 1)

# Cumulative age structure:
cumw <- cumsum(w) / sum(w)
plot(xd, cumw)

# Simulate ages at death:
ri <- runif(n = Nini, min = min(cumw), max = max(cumw))
ageIni <- xd[findInterval(ri, cumw)]

# Plot histogram of X:
par(mfrow = c(1, 1))
hist(ageIni, freq = FALSE)

# Cumulative age structure:
cumw <- cumsum(w) / sum(w)
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
while (iCONT) {
  # Update ind. counter:
  i <- i + 1
  
  if (nrow(longdat) < i) {
    CONT <- FALSE
  } else {
    # Extract dates and ages for ind. i:
    Ai <- longdat$ageEntry[i]
    if (Ai > alpha) {
      ti <- longdat$entryDate[i]
      xt <- Ai
    } else {
      ti <- longdat$entryDate[i] + alpha - Ai
      xt <- alpha
    }
    Bi <- longdat$birthDate[i]
    Di <- longdat$departDate[i]
    Xi <- longdat$ageDeath[i]
    ibi <- NA
    
    # Start reproductive events counter:
    irep <- 0
    
    while (xt < Xi) {
      # Age at time ti:
      if (xt == alpha) {
        wi <- rexp(n = 1, rate = gamma)
        xt <- xt + wi
        ti <- ti + wi
        r1i <- 1
      } else {
        ibi <- tau + rexp(n = 1, rate = eta)
        xt <- xt + ibi
        ti <- ti + ibi
        r1i <- 0
      }
      
      # Birth event:
      idgxi <- findInterval(xt, xg)
      yit <- rpois(n = 1, lambda = gxpop[idgxi])
      reprtemp <- data.frame(ID = i, age = xt, offsp = yit, IBI = ibi, 
                             time = ti, stage = "A", event = "R", first = r1i)
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
    # Stop loop if all remaining individuals were born after Tmax - alpha:
    iCONT <- !all(longdat$birthDate[i:nrow(longdat)] > Tmax - alpha) 
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

realFert <- parents <- offspring <- ageAds * 0
for (ia in 1:nagead) {
  ida <- which(reprdat$age >= ageAds[ia] & reprdat$age <= ageAds[ia] + 1 &
                 reprdat$event == "R")
  offspring[ia] <- sum(reprdat$offsp[ida])
  idp <- which(longdat$ageEntry <= ageAds[ia] & 
                 (longdat$departDate - longdat$birthDate) >= ageAds[ia])
  parents[ia] <- length(idp)
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
  lines(x[idx], w, col = 'grey50', lwd = 2)
}

# ------------------------------ #
# ---- Realized fertility: -----
# ------------------------------ #
par(mfrow = c(1, 1))
plot(ageAds, offspring / parents, pch = 19, col = 'dark red', xlab = "Age",
     ylab = "Average number of offspring")
lines(x, gx, col = 'orange', lwd = 2)


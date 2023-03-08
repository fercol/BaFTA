# ============================= CODE METADATA ================================ #
# FILE NAME: testFecMods03.R
# DATE CREATED: June 02, 2021
# DESCRIPTION: Develop example of tractable age-specific fecundity models that
#              can be fitted with standard maximum likelihood methods.
# AUTHOR: Fernando Colchero
# NOTES: - Found alternative model that combines a quadratic exponential
#          with a logistic models.
#        - 2021-10-02: Incorporate interbirth intervals to estimate the proba-
#        bility that an individual of age breeds given its last breeding.
# ============================== START CODE ================================== #
# ===================== #
# 0) GENERAL SETUP ====
# ===================== #
library(snowfall)
library(RColorBrewer)
setwd("~/FERNANDO/PROJECTS/1.ACTIVE/FecundEstimation/")

# ============
# 1) FUNCTIONS
# ============
# Fecundity models:
fecmods <- list()

# model 1:
fecmods$m1 <- list()
fecmods$m1$mod <- function(x, a) {
  exp(a[1] + a[2] * x - a[3] * x^2)
}
fecmods$m1$inp <- c(-1, 0.00001, 0.0001)
fecmods$m1$col <- "#E41A1C"
fecmods$m1$lty <- 1

# model 2:
fecmods$m2 <- list()
fecmods$m2$mod <- function(x, a) {
  b1 <- a[1] / a[2]
  b2 <- a[5] / a[6]
  mx <- b1 * (1 - exp(-a[2] * x)) + a[3] * x + a[4] * x^2 + 
    b2 * (exp(a[6] * x) - 1)
  return(mx)
}
fecmods$m2$inp <- c(0.16, 1.75, 0.09, - 0.007, 0.0000005, 0.8)
fecmods$m2$col <- "#377EB8"
fecmods$m2$lty <- 2

# model 3:
fecmods$m3 <- list()
fecmods$m3$mod <- function(x, a) {
  a[1] * (- exp(-a[2] * x) + exp(-a[3] * (x - a[4])^2) + exp(a[5] * x))
}
fecmods$m3$inp <- c(0.4, 0.25, 0.002, 10, 0.00002)
fecmods$m3$col <- "#4DAF4A"
fecmods$m3$lty <- 4

# model 4:
fecmods$m4 <- list()
fecmods$m4$mod <- function(x, a) {
  a[1] * exp(-a[2] * (x - a[3])^2) / (1 + exp(-(a[4] * (x - a[3]))))
}
fecmods$m4$inp <- c(1, 0.001, 20, 0.1)
fecmods$m4$col <- "#984EA3"
fecmods$m4$lty <- 5

# model 5:
fecmods$m5 <- list()
fecmods$m5$mod <- function(x, a) {
  exp(a[1] + a[2] * x - a[3] * x^2 + a[4] * 1/(x + 1))
}
fecmods$m5$inp <- c(-1, 0.00001, 0.0001, -1)
fecmods$m5$col <- "#A65628"
fecmods$m5$lty <- 6

# Number of models:
nmods <- length(fecmods)

# Least squares function:
lsf <- function(a, ffun, dat) {
  f <- ffun(x, a)
  lsq <- sum((rfec - f)^2)
  return(lsq)
}

# Optimization function:
myopt <- function(a, lsf, ffun) {
  out1 <- optim(a, lsf, ffun = ffun)
  iter <- 0
  while(out1$convergence != 0 | iter < 100) {
    iter <- iter + 1
    out1 <- optim(out1$par, lsf, ffun = ffun)
  }
  return(out1)
}


# ========================= #
# 2) READ IN DATASETS: ====
# ========================= #
lion <- read.table("data/lionfec.txt", header = TRUE, sep = "\t")
baboon <- read.table("data/baboonfec.txt", header = TRUE, sep = "\t")
bison <- read.csv("data/datarepro.csv", header = TRUE)

# List objects:
dat <- list()

# Data list:
dat$ba <- list(x = baboon$Age[baboon$Age >= 5] - 5, minx = 5, 
               fec = baboon$Fecundity[baboon$Age >= 5], sp = "Baboon")
dat$li <- list(x = lion$Age[lion$Age >= 2] - 2, minx = 2, 
               fec = lion$Fecundity[lion$Age >= 2], sp = "Lion")
dat$bi <- list(x = bison$Age - bison$Age[1], minx = 2, 
               fec = bison$R, sp = "Bison")
spname <- names(dat)
spfull <- c(ba = "Baboon", li = "Lion", bi = "Bison")

# =====================
# 3) FIT MODELS TO DATA
# =====================
mod <- list()

for (sp in spname) {
  x <- dat[[sp]]$x
  rfec <- dat[[sp]]$fec
  mod[[sp]] <- list()
  for (mm in 1:nmods) {
    out <- myopt(fecmods[[mm]]$inp, lsf, fecmods[[mm]]$mod)
    mod[[sp]][[sprintf("m%s", mm)]] <- out
    cat(sprintf("%s - model %s: %s\n", sp, mm, out$convergence))
  }
}

# Print results to the screen:
for (sp in spname) {
  cat("-------")
  cat(sprintf("%s:\n", sp))
  for (i in 1:nmods) {
    cat(sprintf("RMSE%s = %s\n", i, round(sqrt(mod[[sp]][[i]]$value / 
                                            length(dat[[sp]]$x)), 3)))
  }
}

# Plotting settings:
maxAge <- sapply(spname, function(sp) max(dat[[sp]]$x + dat[[sp]]$minx))
maxAge <- (ceiling(max(maxAge) / 5) * 5)
layout(matrix(c(rep(1, 3), 0, 2:5), 4, 2), widths = c(0.25, 1), 
       heights = c(rep(1, 3), 0.25))
par(mar = c(2, 0, 1, 0))
plot(c(0, 1), c(0, 1), col = NA, axes = FALSE, xlab = "", ylab = "")
text(0.5, 0.5, "Age-specific fecundity", srt = 90, cex = 2)

par(mar = c(2, 2, 1, 0.5))
ii <- 0
for (sp in spname) {
  ii <- ii + 1
  plot(dat[[sp]]$x + dat[[sp]]$minx, dat[[sp]]$fec, pch = 19, cex = 1.5,
       xlim = c(0, maxAge), ylim = c(0, 1.2), xlab = "",
       ylab = "", axes = FALSE)
  text(maxAge * 0.05, 1.2, sprintf("%s) %s", letters[ii], spfull[sp]), font = 2,
       cex = 1.5)
  Axis(c(0, maxAge), side = 1, labels = NA, tcl = 0.5, pos = 0)
  Axis(c(0, 1), side = 2, tcl = 0.5, las = 2, pos = 0)
  for (mm in 1:nmods) {
    lines(dat[[sp]]$x + dat[[sp]]$minx, 
          fecmods[[mm]]$mod(dat[[sp]]$x, mod[[sp]][[mm]]$par), 
          col = fecmods[[mm]]$col, lty = fecmods[[mm]]$lty, lwd = 3)
    
  }
}
legend('topright', names(fecmods), 
       col = sapply(1:5, function(mm) fecmods[[mm]]$col),
       lty = sapply(1:5, function(mm) fecmods[[mm]]$lty), lwd = 3,
       bty = 'n', seg.len = 8)
par(mar = c(0, 2, 0, 0.5))
plot(c(0, maxAge), c(0, 1), col = NA, axes = FALSE, xlab = "", ylab = "")
text(maxAge * 0.5, 0.75, "Age (years)", cex = 2)
Axis(c(0, maxAge), side = 3, lwd = NA, line = 0)



# Fancy plot:
PDF <- FALSE
cols <- c('dark red', 'dark blue', 'dark green')
cold <- 'grey80'
ltys <- c(1, 2, 4)
modpl <- c(1, 2, 3)
leglabs <- c("Data", "Null Model", 
             sprintf("Alt. mod. %s", 1:(length(modpl) - 1)))
ylim <- c(0, max(sapply(1:3, function(i) max(dat[[i]]$fec)))) * c(1, 1.1)
xlim <- c(0, maxAge)
if (PDF) {
  pdf("results/fecMods02.pdf", width = 10, height = 4.5)
}
layout(matrix(c(1, 0, 2:((length(dat) * 2) + 1)), 2, length(dat) + 1), 
       widths = c(0.25, rep(1, length(dat))), heights = c(1, 0.35))
# Y-axis:
par(mar = c(0.5, 0, 2, 0))
plot(c(0, 1), ylim, col = NA, xlab = "", ylab = "", axes = FALSE)
Axis(ylim, side = 2, pos = 1.2, las = 2, lwd = NA)
text(0.25, mean(ylim), "Average fecundity", srt = 90, cex = 2)
for (sp in 1:length(dat)) {
  # Fecundity plot:
  par(mar = c(0.5, 0.5, 2, 2))
  #xlim <- c(0, max(dat[[i]]$x + dat[[i]]$minx) + 2)
  plot(xlim, ylim, axes = FALSE, col = NA,
       xlab = "", ylab = "")
  Axis(ylim, side = 2, labels = NA, pos = 0)
  Axis(xlim, side = 1, labels = NA, pos = 0)
  mtext(sprintf("%s) %s", letters[sp], dat[[sp]]$sp), side = 3, line = 0,
        at = xlim[1] + diff(xlim) / 10, adj = 0)
  # points(dat[[sp]]$x + dat[[sp]]$minx, dat[[sp]]$fec, col = cold, pch = 19,
  #        cex = 1)
  # lines(dat[[sp]]$x + dat[[sp]]$minx, dat[[sp]]$fec, type ='h', 
  #       col = cold)
  for (mm in modpl) {
    lines(dat[[sp]]$x + dat[[sp]]$minx, 
          fecmods[[mm]]$mod(dat[[sp]]$x, mod[[sp]][[mm]]$par), 
          col = fecmods[[mm]]$col, lty = fecmods[[mm]]$lty, lwd = 3)
  }
  points(dat[[sp]]$x + dat[[sp]]$minx, dat[[sp]]$fec, col = cold, pch = 19,
         cex = 1)
  lines(dat[[sp]]$x + dat[[sp]]$minx, dat[[sp]]$fec, type ='h', 
        col = cold)
  par(mar = c(0, 0.5, 0, 2))
  plot(xlim, c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
  Axis(xlim, side = 1, lwd = NA, pos = 1.2)
  if (sp == 2) text(mean(xlim), 0.5, "Age (years)", cex = 2)
  if (sp == 3) {
    legend('right', leglabs, lwd = c(NA, rep(4, length(modpl))), 
           pch = c(19, rep(NA, length(modpl))), pt.cex = c(2, NA, NA),
           col = c(cold, sapply(modpl, function(mm) fecmods[[mm]]$col)), 
           lty = c(NA, sapply(modpl, function(mm) fecmods[[mm]]$lty)), 
           bty = 'n', seg.len = 10, cex = 1.25)
  }
}
if (PDF) dev.off()
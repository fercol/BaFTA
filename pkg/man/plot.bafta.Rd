\name{plot.bafta}
\alias{plot.bafta}
\alias{print.bafta}
\alias{summary.bafta}

\title{
 Summarizing and plotting Bayesian Fertility Trajectory Analysis (BaFTA) model outputs
}

\description{
 These functions are all generic methods for class \code{bafta}.
}

\usage{
  \method{plot}{bafta}(x, type = "traces", \dots)
  \method{summary}{bafta}(object, \dots)
  \method{print}{bafta}(x, \dots)
}

\arguments{
  \item{x }{An object of class \code{bafta}.}
  \item{object }{An object of class \code{bafta}.}
  \item{type }{Characther string for the type of BaFTA plot to be produced. Options are \dQuote{\code{traces}} (the default), \dQuote{\code{density}}, \dQuote{\code{fertility}}, and \dQuote{\code{predictive}} (see details).}
  \item{\dots }{Additional arguments passed to functions \code{print}, \code{summary} and \code{plot} (see details).}
}

\details{
For objects of class \code{bafta}, the \code{print} and \code{summary} functions display on the console four summary elements describing the model and its results, namely \code{Call}, \code{Coefficients}, \code{Convergence} and \code{Model fit}. The section \code{call} describes the fertility model used (e.g. quadratic, PeristeraKostaki, Gamma, etc.), the data type, and general settings of the Bayesian implementation (e.g., number of iterations, burn-in sequence, number of simulations, etc.). The \code{Model fit} section includes the \code{DIC} and \code{Predictive loss} values. 

The function \code{plot} can produce four different types of plots specified with argument \code{type}, namely \dQuote{\code{traces}} for the MCMC parameter traces which faciliates visual inspection of their convergence, \dQuote{\code{density}} for the parameter posterior densities calculated from the remaining MCMC iterations after burn-in and thinning, \dQuote{\code{fertility}} for the estimated age-specific fertility (i.e., average number of offspring produced per individual of a given age) with their 95\% credible intervals under the actual fertility (i.e., as from the data), and \dQuote{\code{predictive}} for the predicted total number of offspring used to calculate predictive loss as proposed by Gelfand and Gosh (1998), including the 95\% predictive intervals and the actual total number of offspring as in the data.
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{bafta}}}

\examples{
# Load data:
data("lionFert", package = "BaFTA")

# Run quadratic model:
out <- bafta(object = lionFert, model = "quadratic", 
                dataType = "aggregated", nsim = 2, ncpus = 2, 
                niter = 100, burnin = 11, UPDJUMP = FALSE, 
                jumpSD = rep(0.1, 4))

# Print results to the console:
summary(out)

# Plot traces:
plot(out)

# Plot posterior densities:
plot(out, type = "density")

# Plot estimated fertility:
plot(out, type = "fertility")

# Plot predicted number of offspring:
plot(out, type = "predictive")

}

\keyword{methods}


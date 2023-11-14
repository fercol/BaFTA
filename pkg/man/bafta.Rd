\name{bafta}
\alias{bafta}
\alias{bafta.default}

\title{
  Parametric Bayesian estimation of age-specific fertility for aggregated and individual level data
}

\description{
  This function performs multiple Markov Chain Monte Carlo (MCMC) simulations for the Bayesian estimation of age-specific fertility for aggregated and individual level data. Fertility parameters are estimated, allowing the user to test a range of fertility patterns, following the approach in Colchero (\emph{in prep.}).
}

\usage{
bafta(object, \dots)

\method{bafta}{default} (object, dataType = "aggregated", model = "quadratic",  
        minAge = NA, gestTime = NA, niter = 55000, burnin = 5001, 
        thinning = 20, nsim, ncpus, UPDJUMP = TRUE, 
        jumpSD = NULL,  ...)}

\arguments{
	\item{object }{Data frame with either aggregated or individual level fertility data (see details)}

	\item{dataType }{Character string specifying the type of data. Options are \dQuote{\code{aggregated}} (default), \dQuote{\code{indivSimple}}, and \dQuote{\code{indivExtended}}.}

	\item{model }{Age-specific fertility model. Options are \dQuote{\code{quadratic}} (default), \dQuote{\code{PeristeraKostaki}}, \dQuote{\code{ColcheroMuller}}, \dQuote{\code{Hadwiger}}, \dQuote{\code{gamma}}, \dQuote{\code{beta}}, \dQuote{\code{skewNormal}}, \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} (see details)}

	\item{minAge }{Numeric value for the earliest age at maturity (default \code{NA}). If the value is \code{NA}, then the earliest age in the dataset is taken.}

	\item{gestTime }{Numeric value for the minimum time of gestation and, if applicable, weaning (default \code{NA}), only necessary for \code{dataType} = \dQuote{\code{indivExtendend}}. If the value is \code{NA}, then the earliest time between reproductive bouts is taken.}

  \item{niter }{The total number of MCMC steps.}
  
  \item{burnin }{The number of iterations for the burn in (see \code{details}).}
  
  \item{thinning }{The number of skipped MCMC steps to minimize serial autocorrelation (see \code{details}).}

  \item{nsim }{An integer for the number of MCMC simulations to be run.}

	\item{ncpus }{An integer for the number of CPUS to be used for parallel computing.}

	\item{UPDJUMP }{Logical specifying if the jump standard deviations should be estimated.}

	\item{jumpSD }{Numerical vector with the jump standard deviations per parameter. Default is \code{NULL} when \code{UPDJUMP = TRUE}.}

	\item{ ... }{Additional arguments for prior means, prior standard deviations, etc. (see details).}

}

\details{
\code{1) DATA TYPES:}

The input \code{object} required by BaFTA needs to be constructed differently whether the data are aggregated, individual level with seasonal reproduction, or individual level with variable IBI. 

\code{1.a) Aggregated data:}

If \code{dataType} = \dQuote{\code{aggregated}}, then the data frame requires the following columns: \dQuote{\code{Age}} with discrete ages; \dQuote{\code{nParents}} with the number of available adults per age; \dQuote{\code{nOffspring}} with the number of offspring produced in that age interval; \dQuote{\code{Fertility}} with the age-specific fertility, calculated as \code{Fertility = nOffspring / nParents}. 

\code{1.b) Individual data with discrete age intervals:}

If \code{dataType} = \dQuote{\code{indivSimple}}, then the data are individual level data but with discrete age intervals, commonly used when reproduction is seasonal, or when only the realized fertility is necessary, for instance, when modeling population dynamics. In this case, the data needed consist of a data frame with individual rows per reproductive event, and with the following columns: \dQuote{\code{indID}} with character string of individual IDs of the parents; \dQuote{\code{Age}} with the discrete ages of the parents at the time of birth; \dQuote{\code{nOffspring}} with the number of offspring produced during that reproductive event by the given parent.

\code{1.c) Individual data with continuous ages and variable IBIs:}

If \code{dataType} = \dQuote{\code{indivExtended}}, then the data are individual level data but with continuous ages, commonly used when reproduction is not seasonal or when interbirth intervals (IBIs) are variable. In this case, the data needed consist of a data frame with individual rows per reproductive event, and with the following columns: \dQuote{\code{indID}} with character string of individual IDs of the parents; \dQuote{\code{Age}} with the discrete ages of the parents at the time of birth; \dQuote{\code{nOffspring}} with the number of offspring produced during that reproductive event by the given parent, including failed events (i.e., \code{nOffspring = 0}); \dQuote{\code{IBI}} with the time since last reproduction in years.

2) AGE-SPECIFIC FERTILITY MODELS:

The age-specific fertility models correspond to the expected number of offspring produced by adults of a given age. Therefore, for a random variable \eqn{Y_{x}} with realizations \eqn{y_x} for the number of offspring produced by adults of age \eqn{x}, we have that \eqn{E(Y_x) = g(x | \theta)}, where \eqn{g: \mathbb{R}_{\ge 0} \rightarrow \mathbb{R}_{\ge 0}} is a non-negative smooth fertility function and \eqn{\theta} is a vector of parameters to be estimated. The functional forms of function \eqn{g} fall in two main categories, namely polynomial and distributional models.

2.a) Polynomial:

Of the models available in BaFTA, the \dQuote{\code{quadratic}}, \dQuote{\code{PeristeraKostaki}} (Peristera and Kostaki 2007), and \dQuote{\code{ColcheroMuller}} (Colchero  \emph{et al.} 2021, Muller \emph{et al.} 2020) fall within the polynomial category. Both, the \dQuote{\code{PeristeraKostaki}} and \dQuote{\code{ColcheroMuller}} are non-symmetric around the age at maximum fertility.

2.b) Distributional:

The distributional models are of the form

\deqn{g(x | \theta) = R f(x | \theta_1),}

where \eqn{f: \mathbb{R}_{\ge 0} \rightarrow \mathbb{R}_{\ge 0}} is a probability density function, \eqn{R > 0} is a parameter for the total fertility rate, and \eqn{\theta_1 \in \theta} is a vector of parameters. The
\dQuote{\code{Hadwiger}} (Hadwiger 1940), \dQuote{\code{gamma}} (Hoem \emph{et al.} 1981), \dQuote{\code{beta}}  (Hoem \emph{et al.} 1981), \dQuote{\code{skewNormal}} (Mazzuco and Scarpa 2011, 2015), \dQuote{\code{gammaMixture}} (Hoem \emph{et al.} 1981), \dQuote{\code{HadwigerMixture}} (Chandola \emph{et al.} 1991), \dQuote{\code{skewSymmetric}} (Mazzuco and Scarpa 2011, 2015), and \dQuote{\code{skewLogistic}} (Asili \emph{et al.} 2014) all fall in this category. Notably, the \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} are appropriate when fertility might be bimodal (Hoem \emph{et al.} 1981, Chandola \emph{et al.} 1999, Mazzuco and Scarpa 2011, 2015, Asili \emph{et al.} 2014). 

%3) 

}

\value{
\item{coefficients }{A matrix with estimated coefficients (i.e. mean values per parameter on the thinned sequences after burnin), which includes standard errors, upper and lower 95\% credible intervals, and the potential scale reduction factor (Gelman \emph{et al.} 2013) for convergence (convergence is achieved for values below 1.05).}

\item{x }{Vector of ages to plot estimated fertility and predicted number of offspring.}

\item{fert }{Data frame with the estimated age-specific fertility including the upper and lower 95\% credible intervals.} 

\item{theta }{Matrix of the converged estimated parameters.} 

\item{uSd }{Vector of converged estimated standard deviation of random effects for fertility (only produced for \code{dataType}  \dQuote{\code{indivSimple}} and \dQuote{\code{indivExtended}}).} 

\item{vSd }{Vector of converged estimated standard deviation of random effects for IBI (only produced for \code{dataType} = \dQuote{\code{indivExtended}}).} 

\item{likePost }{Matrix of likelihood and posterior for the converged sequence.} 

\item{DIC }{Deviance information criterion (DIC) as a measure of model fit (Spiegelhalter \emph{et al.} 2002, Celeux \emph{et al.} 2006). Small differences between values should only be used as a reference (see comments in Spiegelhalter \emph{et al.} 2002). } 

\item{PredLoss }{Predictive loss (PL) as measure of model fit as proposed by Gelfand \emph{et al.} (1998). } 

\item{pred }{Data frame per age with the predicted number of offspring produced by adults of the given age including the mean and upper and lower 95\% predictive intervals (used for plotting the predictive estimates).} 

\item{aggrData }{Data frame with the aggregated data, aggregated by age for plotting. For \code{dataType} = \dQuote{\code{aggregated}} it corresponds to the data in \code{object}.} 

\item{runs }{List with the output for the individual MCMC runs.} 

\item{data }{List with the information in \code{object} and other variables relevant for modelling.} 

\item{settings }{List object indicating the model, the data type, the number of iterations for each MCMC, the burn in sequence, the thinning interval, the number of simulations, the number of CPUS, and the total computing time in minutes.} 

\item{keep }{Numerical vector of indices for the kept sequence for each MCMC chain.} 

\item{params }{List object with general information on the parameters such as initial values, prior means, prior sd, etc.}
}

\references{
Asili, S., Rezaei, S. & Najjar, L. (2014) Using Skew-Logistic Probability Density Function as a Model for Age-Specific Fertility Rate Pattern. \emph{BioMed Research International}, 2014, 790294.

Azzalini, A. (1985) A Class of Distributions Which Includes the Normal Ones. \emph{Scandinavian Journal of Statistics}, 2, 171-178.

Chandola, T., Coleman, D.A. & Hiorns, R.W. (1999) Recent European fertility patterns: Fitting curves to `distorted' distributions. \emph{Population Studies}, 53, 317-329.

Celeux, G., Forbes, F., Robert, C.P. & Titterington, D.M. (2006) Deviance information criteria for missing data models. \emph{Bayesian Analysis}, 1, 651-673.

Colchero, F. (\emph{In prep.}) Inference on age-specific fertility in ecology and evolution. Learning from other disciplines and improving the state of the art. 

Colchero, F., Eckardt, W. & Stoinski, T. (2021) Evidence of demographic buffering in an endangered great ape: Social buffering on immature survival and the role of refined sex-age classes on population growth rate. \emph{Journal of Animal Ecology}, 90, 1701-1713.

Gelfand, A. & Ghosh, S. (1998) Model Choice: A Minimum Posterior Predictive Loss Approach. \emph{Biometrika}, 85, 1-11.

Hadwiger, H. (1940) Eine analytische Reproduktionssunktion fur biologische Gesamtheiten. \emph{Scandinavian Actuarial Journal}, 1940, 101-113.

Hoem, J.M., Madien, D., Nielsen, J.L., Ohlsen, E.M., Hansen, H.O. & Rennermalm, B. (1981) Experiments in modelling recent Danish fertility curves. \emph{Demography}, 18, 231-244.

Mazzuco, S. & Scarpa, B. (2011) Fitting age-specific fertility rates by a skew- symmetric probability density function. (Working paper 10), University of Padua.

Mazzuco, S. & Scarpa, B. (2015) Fitting age-specific fertility rates by a flexible generalized skew normal probability density function. \emph{Journal of the Royal Statistical Society: Series A}, 178, 187-203.

Muller, M. N., Blurton Jones, N. G, Colchero, F., Thompson, M. E., Enigk, D. K. (2020) Sexual dimorphism in chimpanzee (\emph{Pan troglodytes schweinfurthii}) and human age-specific fertility. \emph{Journal of human evolution}, 144, 102795.

Peristera, P. & Kostaki, A. (2007) Modeling fertility in modern populations. \emph{Demographic Research}, 16, 141-194.

Spiegelhalter, D., Best, N., Carlin, B. & Linde, A.V.D. (2002) Bayesian measures of model complexity and fit. \emph{Journal Of The Royal Statistical Society Series B}, 64, 583-639.

}
\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{summary.bafta}}, \code{\link{print.bafta}}, \code{\link{plot.bafta}} to visualise summary outputs for objects of class \dQuote{\code{bafta}}. 
}

\examples{
# Load data:
data("lionFert", package = "BaFTA")

# Run quadratic model:
out <- bafta(object = lionFert, model = "quadratic", 
                dataType = "aggregated", nsim = 4, ncpus = 4)

# Print results to the console:
out

# Plot traces:
plot(out)

# Plot posterior densities:
plot(out, type = "density")

# Plot estimated fertility:
plot(out, type = "fertility")

# Plot predicted number of offspring:
plot(out, type = "predictive")

}

\keyword{htest,models}


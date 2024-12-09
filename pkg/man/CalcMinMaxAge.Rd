\name{CalcMinMaxAge}
\alias{CalcMinMaxAge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating minimum and maximum possible ages for individual level data.
}
\description{
This function is inteded to help the user put bounds to ages in case these are uncertain or unknown for some individuals in the dataset.
}
\usage{
CalcMinMaxAge(object, minAge, maxAge, unkAgeID = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object }{Data frame with individual level fertility data (see details)}
  \item{minAge }{Numerical value for the minimum age, commonly the minimum age at sexual maturity}
  \item{maxAge }{Numerical value for the maximum age, commonly the maximum recorded lifespan of the species}
  \item{unkAgeID }{IDs of individuals with uncertain ages. If \code{NULL}, then all idividuals are assumed to have unknonw ages.}
}
\details{
\code{1) DATA TYPES:}

The input \code{object} required by BaFTA needs to be constructed differently whether the data are aggregated, individual level with seasonal reproduction, or individual level with variable IBI. The function \code{CalcMinMaxAge} allows users to put bounds to the ages of individuals for which these might be uncertain. Here we specify to types of datasets that the function works on.

\code{1.a) Individual data with discrete age intervals:}

If \code{dataType} = \dQuote{\code{indivSimple}}, then the data are individual level data but with discrete age intervals, commonly used when reproduction is seasonal, or when only the realized fertility is necessary, for instance, when modeling population dynamics. In this case, the data needed consist of a data frame with individual rows per reproductive event, and with the following columns: \dQuote{\code{indID}} with character string of individual IDs of the parents; \dQuote{\code{Age}} with the discrete ages of the parents at the time of birth; \dQuote{\code{nOffspring}} with the number of offspring produced during that reproductive event by the given parent.

\code{1.b) Individual data with continuous ages and variable IBIs:}

If \code{dataType} = \dQuote{\code{indivExtended}}, then the data are individual level data but with continuous ages, commonly used when reproduction is not seasonal or when interbirth intervals (IBIs) are variable. In this case, the data needed consist of a data frame with individual rows per reproductive event, and with the following columns: \dQuote{\code{indID}} with character string of individual IDs of the parents; \dQuote{\code{Age}} with the discrete ages of the parents at the time of birth; \dQuote{\code{nOffspring}} with the number of offspring produced during that reproductive event by the given parent, including failed events (i.e., \code{nOffspring = 0}); \dQuote{\code{IBI}} with the time since last reproduction in years.

\code{2) CALCULATION OF MIN-MAX AGES}

The function calculates the range of ages for 
}
\value{
The function returns the original data frame with two additional numerical columns, \dQuote{\code{MinAge}} and \dQuote{\code{MaxAge}}.
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bafta}} for inference on age-specific fertility, \code{\link{summary.bafta}}, \code{\link{print.bafta}}, \code{\link{plot.bafta}} to visualise summary outputs for objects of class \dQuote{\code{bafta}}. 
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

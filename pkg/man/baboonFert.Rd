\name{baboonFert}
\alias{baboonFert}
\docType{data}
\title{
Aggregated age-specific fertility of female  olive babbons
}
\description{
Age specific fertility of olive baboons as reported by Packer \emph{et al.} (1998).}
\usage{data("baboonFert")}
\format{
  A data frame with 27 observations on the following 4 variables.
  \describe{
    \item{\code{Age}}{a numeric vector of discrete ages}
    \item{\code{nParents}}{a numeric vector of the number of available parents of the given age}
    \item{\code{nOffspring}}{a numerical vector of the number of offspring produced by parents of the given age}
    \item{\code{Fertility}}{a numerical vector of the average age-specific fertility}
  }
}
\details{
The fertility data was digitized from Fig. 1c in Packer \emph{et al.} (1998), while the number of females (\dQuote{\code{nParents}}) was obtained from Fig. 1a. The number of offspring (\dQuote{\code{nOffspring}}) was calculated as \code{floor(nParents * Fertility)}.

For inference on age-specific fertility on this data using package BaFTA, see Colchero (\emph{in prep.}).
}
\references{
Colchero, F. (\emph{In prep.}) Inference on age-specific fertility in ecology and evolution. Learning from other disciplines and improving the state of the art. 

Packer, C., Tatar, M. & Collins, A. (1998) Reproductive cessation in female mammals. \emph{Nature}, 392, 807-811.
}
\examples{
# Load data:
data("baboonFert", package = "BaFTA")

# First six rows:
head(baboonFert)
}
\keyword{datasets}

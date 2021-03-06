\name{l1pmreg}
\alias{l1pmreg}
\alias{quantpmreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Piecewise monotone nonparameric quantile regression}
\description{
Applies the generalized taut string method to quantile regression.   
}
\usage{
quantpmreg(y, beta = 0.5, squeezing.factor = 0.5, verbose = FALSE, localsqueezing = TRUE, DYADIC = TRUE, thr.const = 2, extrema.nr = -1, bandwidth = -1,SETTOMEAN = FALSE, method = 1, ...)
l1pmreg(y, beta=0.5, squeezing.factor = 0.5, verbose = FALSE, localsqueezing = TRUE, DYADIC = TRUE, thr.const = 2, extrema.nr = -1, bandwidth = -1,SETTOMEAN = FALSE, method = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{observed values (ordered by value of independent variable)}
  \item{beta}{quantile. The default is 0.5 which corresponds to 
the robust taut string.}
  \item{squeezing.factor}{The amount of decrement applied to the bandwidthes}
  \item{verbose}{logical, if T progress (for each iteration) is illustrated grahically}
  \item{localsqueezing}{logical, if T (default) the bandwidth is changed
locally.}
  \item{DYADIC}{logical, if T (default) the multiresolution criterion
is only verified on intervals with dyadic endpoints.}
  \item{thr.const}{smoothing parameter for the multiresoultion criterion
(should be approximately 2)}
  \item{extrema.nr}{if set to a positive integer an approximation with
the specified number of local extreme values is calculated}
  \item{bandwidth}{if set to a positive value the specified bandwidth is
used instead of the multiresolution criterion.}
  \item{SETTOMEAN}{logical, if T (default) the value of the taut string
approximation at local extreme values is set to the mean or median of the observations
on the interval where the extremum is taken.}
  \item{method}{The method used which can be 1 (quantile regression),
2 (usual taut string), 3 (logistic regression) and 4 (Poisson regression) }
  \item{\dots}{Passed to the plot command if verbose=T}
}
\value{
A list with components
  \item{y}{The approximation of the given data}
  \item{lambda}{The final values of lambda}
  \item{nmax}{Number of local extreme values}
}
\references{D\"umbgen, L. and Kovac, A. (2003) Extensions of smoothing via taut strings}     
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{pmreg}},\code{\link{frun}}}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ robust }% __ONLY ONE__ keyword per line


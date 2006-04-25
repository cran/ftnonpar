\name{genpmreg}
\alias{genpmreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Piecewise monotone regression with generalised taut strings}
\description{
Applies the generalised taut string method to one-dimensional data.
}
\usage{
genpmreg(y, beta = 0.5, squeezing.factor = 0.5, verbose = FALSE, localsqueezing = TRUE, DYADIC = TRUE, thr.const = 2, extrema.nr = -1, bandwidth = -1, SETTOMEAN = FALSE, method = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{observed values (ordered by value of independent variable)}
  \item{beta}{If method=1 specifies the quantile}
  \item{squeezing.factor}{The amount of decrement applied to the bandwidthes}
  \item{verbose}{logical, if T progress (for each iteration) is illustrated grahically}
  \item{localsqueezing}{logical, if T (default) the bandwidth is changed
locally.}
  \item{DYADIC}{If T checks the multiresolution criterion only on
dyadic intervals, otherwise all intervals are checked}
  \item{thr.const}{smoothing parameter for the multiresolution criterion
(should be approximately 2.5)}
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
  \item{lambda}{Values for lambda used}
  \item{nmax}{Number of local extreme values}
}
\references{D\"umbgen, L. and Kovac, A. (2006) Extensions of smoothing via taut strings}
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{pmreg}}}
\examples{
data(djdata)
par(mfrow=c(2,2))
plot(djblocks,col="grey")
lines(genpmreg(djblocks,verbose=FALSE,method=2,thr.const=2.5)$y,col="red")
title("Usual taut string method")
ind <- sample(1:length(djblocks),300)
djblocks[ind] <- djblocks[ind]+rnorm(length(ind),0,100)
plot(djblocks,col="grey")
lines(genpmreg(djblocks,verbose=FALSE,method=2)$y,col="red")
title("Usual taut string method with outliers")
plot(djblocks,col="grey")
lines(genpmreg(djblocks,verbose=FALSE,method=1)$y,col="red")
title("Robust taut string method with outliers")
plot(djblocks,col="grey",ylim=c(-10,10))
lines(genpmreg(djblocks,verbose=FALSE,method=1)$y,col="red")
title("Again robust taut string method, different scaling of y-axis")
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ robust }% __ONLY ONE__ keyword per line

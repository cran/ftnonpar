\name{pmreg}
\alias{pmreg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Piecewise monotone regression with taut strings}
\description{
Applies the taut string method to one-dimensional data.
}
\usage{
pmreg(y, thr.const=2.3, verbose=FALSE, extrema.nr=-1, bandwidth=-1, 
sigma=-1, localsqueezing=TRUE, squeezing.factor=0.5, tolerance=1e-08,
extrema.mean=TRUE,DYADIC=TRUE,dyad.factor=1.1,POSTISO=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{observed values (ordered by value of independent variable)}
  \item{thr.const}{smoothing parameter for the multiresolution criterion
(should be approximately 2.3)}
  \item{verbose}{logical, if T progress (for each iteration) is illustrated grahically}
  \item{extrema.nr}{if set to a positive integer an approximation with
the specified number of local extreme values is calculated}
  \item{bandwidth}{if set to a positive value the specified bandwidth is 
used instead of the multiresolution criterion.}
  \item{sigma}{if set to a positive value sigma the standard deviation is set to sigma and not estimated from the data}
  \item{localsqueezing}{logical, if TRUE (default) the bandwidth is changed
locally.}
  \item{squeezing.factor}{The amount of decrement applied to the bandwidthes}
  \item{tolerance}{Accuracy used for the determination of the bandwidth
when extrema.nr is greater than 0.}
  \item{extrema.mean}{logical, if TRUE (default) the value of the taut string
approximation at local extreme values is set to the mean of the observations
on the interval where the extremum is taken.}
  \item{DYADIC}{If TRUE the multiresolution constraints are only checked on dyadic intervals.}
  \item{dyad.factor}{If the multiresolution constraints are checked on dyadic intervals, dyad.factor determines the ratio between the lengths of two subsequent level (default is 1.1).}
  \item{POSTISO}{If TRUE (default) any bias caused by local squeezing is
removed by applying isotonic and isotonic regression between each two
local extreme values.}
}
\value{
A list with components
  \item{y}{The approximation of the given data}
  \item{sigma}{Standard deviation used}
  \item{widthes}{Bandwidth used}
  \item{nmax}{Number of local extreme values}
  \item{knotsind}{Indices of knot points}
  \item{knotsy}{y-koordinates of knots of the taut string}
}
\references{Davies, P. L. and Kovac, A. (2001) Local Extremes, Runs, Strings and Multiresolution (with discussion) Annals of Statistics. 29. p1-65}
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{mintvmon}},\code{\link{l1pmreg}},\code{\link{pmden}},\code{\link{pmspec}}}

\examples{
data(djdata)
pmreg(djdoppler,verbose=TRUE)$n
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line


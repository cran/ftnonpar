\name{pmlogreg}
\alias{pmlogreg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Piecewise monotone logistic regression with taut strings}
\description{
Applies the taut string method to binary data.
}
\usage{
pmlogreg(y, thr.const=2.5, verbose=FALSE, extrema.nr=-1, bandwidth=-1, 
localsqueezing=TRUE, squeezing.factor=0.5, tolerance=0.001,extrema.mean=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{observed values (ordered by value of independent variable)}
  \item{thr.const}{smoothing parameter for the multiresolution criterion
(should be approximately 2.5)}
  \item{verbose}{logical, if T progress (for each iteration) is illustrated grahically}
  \item{extrema.nr}{if set to a positive integer an approximation with
the specified number of local extreme values is calculated}
  \item{bandwidth}{if set to a positive value the specified bandwidth is 
used instead of the multiresolution criterion.}
  \item{localsqueezing}{logical, if T (default) the bandwidth is changed
locally.}
  \item{squeezing.factor}{The amount of decrement applied to the bandwidthes}
  \item{tolerance}{Accuracy used for the determination of the bandwidth
when extrema.nr is greater than 0.}
  \item{extrema.mean}{logical, if T (default) the value of the taut string
approximation at local extreme values is set to the mean of the observations
on the interval where the extremum is taken.}
}
\value{
A list with components
  \item{y}{The approximation of the given data}
  \item{widthes}{Bandwidth used}
  \item{nmax}{Number of local extreme values}
  \item{knotsind}{Indices of knot points}
  \item{knotsy}{y-koordinates of knots of the taut string}
}
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{l1pmreg}},\code{\link{pmden}},\code{\link{pmspec}}}

\examples{
aaa<-rbinom(1024,1,0.5+0.5*sin(seq(0,10*pi,len=1024)))
pmlogreg(aaa,verbose=TRUE)$n
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line


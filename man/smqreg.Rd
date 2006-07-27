\name{smqreg}
\alias{smqreg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Smooth piecewise monotone regression with taut strings}
\description{
Applies the smooth taut string method to one-dimensional data.
}
\usage{
smqreg(y, thr.const=2.5, verbose=FALSE, bandwidth=-1, 
sigma=-1, localsqueezing=TRUE, squeezing.factor=0.5, DYADIC=TRUE,
firstlambda=100,smqeps=1/length(y),fsign=double(0),gensign=TRUE,
tolerance = 1e-12,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{observed values (ordered by value of independent variable)}
  \item{thr.const}{smoothing parameter for the multiresolution criterion
(should be approximately 2.5)}
  \item{verbose}{logical, if T progress (for each iteration) is illustrated grahically}
  \item{bandwidth}{if set to a positive value the specified bandwidth is 
used instead of the multiresolution criterion.}
  \item{sigma}{if set to a positive value sigma the standard deviation is set to sigma and not estimated from the data}
  \item{localsqueezing}{logical, if T (default) the bandwidth is changed
locally.}
  \item{squeezing.factor}{The amount of decrement applied to the bandwidthes}
  \item{DYADIC}{If TRUE the multiresolution constraints are only checked on dyadic intervals.}
  \item{firstlambda}{Initial value of lambda's for local or global squeezing.}
  \item{smqeps}{Distance between the (equally-spaced) time points.}
  \item{fsign}{Monotonicity constraints, vector of size n-1 of -1,0 and 1's. If fsign[i]==1, then fhat[i+1]>= fhat[i]. If fsign[i]==-1, then fhat[i+1]<=f[i].
Otherwise no constraint at this position.}
  \item{gensign}{If TRUE the taut string method is used to automatically 
produce suitable monotonicity constraints.}
  \item{tolerance}{Precision for the nested intervals for solving the 
minimisation problem.}
  \item{\dots}{Passed to the plot command if verbose=T.}
}
\value{
A list with components
  \item{y}{The approximation of the given data}
  \item{nmax}{Number of local extreme values}
  \item{sigma}{Standard deviation used}
}
\references{Kovac, A. (2006) Smooth functions and local extreme values. Technical Report}
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{pmreg},\link{mintvmon}},\code{\link{l1pmreg}},\code{\link{pmden}},\code{\link{pmspec}}}

\examples{
data(djdata)
par(mfrow=c(2,2))
plot(djblocks,col="grey")
lines(smqreg(djblocks)$y,col="red")
plot(djbumps,col="grey")
lines(smqreg(djbumps)$y,col="red")
plot(djheavisine,col="grey")
lines(smqreg(djheavisine)$y,col="red")
plot(djdoppler,col="grey")
lines(smqreg(djdoppler)$y,col="red")
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line


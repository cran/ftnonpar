\name{smdenreg}
\alias{smdenreg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Piecewise monotone density estimation with smooth taut strings}
\description{
Applies the smooth taut string method to one-dimensional data.
}
\usage{
smdenreg(x, verbose = FALSE, bandwidth=-1, maxkuipnr=19,asympbounds=FALSE, 
squeezing.factor=0.9, firstlambda=10,smqeps=1/length(x),fsign=double(0),
gensign=TRUE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{observed values}
  \item{verbose}{logical, if T progress (for each iteration) is illustrated grahically}
  \item{bandwidth}{if set to a positive value the specified bandwidth is
used instead of the automatic criterion based on generalized Kuiper metrics.}
  \item{maxkuipnr}{The order of the highest generalized Kuiper metric used for the
automatic choice of the bandwidth}
  \item{asympbounds}{If set to T asymptotic bounds derived from a Brownian Bridge are used
for the Kuiper criterion. Otherwise simulated bounds for various sample sizes are interpolated
for the size of the data x}
  \item{squeezing.factor}{The amount of decrement applied to the bandwidthes}
  \item{firstlambda}{Initial value of lambda's for global squeezing.}
  \item{smqeps}{Distance between the (equally-spaced) time points.}
  \item{fsign}{Monotonicity constraints, vector of size n-1 of -1,0 and 1's. If fsign[i]==1, then fhat[i+1]>= fhat[i]. If fsign[i]==-1, then fhat[i+1]<=f[i].
Otherwise no constraint at this position.}
  \item{gensign}{If TRUE the taut string method is used to automatically
produce suitable monotonicity constraints.}
  \item{\dots}{Passed to the plot command if verbose=T.}
}
\value{
  \item{x}{The sorted data}
  \item{y}{values of the density approximation between the observations}
  \item{nmax}{Number of local extreme values}
  \item{trans}{taut string at the observations, should look like uniform
noise}
}
\references{Kovac, A. (2006) Smooth functions and local extreme values. Technical Report}
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{pmreg}},\code{\link{l1pmreg}},\code{\link{pmspec}}}

\examples{
y <- rclaw(500)
hist(y,col="lightgrey",40,freq=FALSE)
lines(smdenreg(y),col="red")
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line


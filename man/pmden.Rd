\name{pmden}
\alias{pmden}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Piecewise monotone density estimation with taut strings}
\description{
Applies the taut string method to one-dimensional data.
}
\usage{
pmden(x, DISCR=FALSE,verbose = FALSE, bandwidth=-1, extrema.nr = -1, accuracy = mad(x)/1000,  extrema.mean = TRUE,maxkuipnr=19,asympbounds=FALSE, tolerance = 0.001)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{observed values}
  \item{DISCR}{logical, if T a discrete density is fitted}
  \item{verbose}{logical, if T progress (for each iteration) is illustrated grahically}
  \item{bandwidth}{if set to a positive value the specified bandwidth is
used instead of the automatic criterion based on generalized Kuiper metrics.}
  \item{extrema.nr}{if set to a positive integer an approximation with
the specified number of local extreme values is calculated}
  \item{accuracy}{Precision of the data. Handling of identical observations
depends on this parameter.}
  \item{extrema.mean}{logical, if T the value at the local extrema is changed to the mean frequency of observations on that interval}
  \item{maxkuipnr}{The order of the highest generalized Kuiper metric used for the
automatic choice of the bandwidth}
  \item{asympbounds}{If set to T asymptotic bounds derived from a Brwonian Bridge are used
for the Kuiper criterion. Otherwise simulated bounds for various sample sizes are interpolated
for the size of the data x}
  \item{tolerance}{Accuracy used for the determination of the bandwidth
when extrema.nr is greater than 0.}
}
\value{
  \item{y}{values of the density approximation between the observations}
  \item{widthes}{bandwidth used for the taut string approximation}
  \item{nmax}{number of local extreme values}
  \item{ind}{indices of knots points of the taut string}
  \item{trans}{taut string at the observations, should look like uniform
noise}
}
\references{Davies, P. L. and Kovac, A. (2003) Densities, Spectral Densities and Modality}
\author{Arne Kovac \email{Arne.Kovac@uni-essen.de}}
\seealso{\code{\link{pmreg}},\code{\link{l1pmreg}},\code{\link{pmspec}}}

\examples{
aaa <- rclaw(500)
pmden(aaa,verb=TRUE)$n
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line


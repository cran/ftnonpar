\name{kkuip}
\alias{kkuip}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Kuiper Metrics}
\description{
Calculates the generalized k-th Kuiper metric
}
\usage{
kkuip(x, k=1)
}
\arguments{
  \item{x}{Data vector}
  \item{k}{Number of intervals to be used for the generalized Kuiper
metric. The usual Kuiper metric is obtained by the default value 1.}
}
\value{
  \item{metric}{Value of the metric}
  \item{a     }{The left borders of the intervals where the maximum is obtained.}
  \item{b     }{The right borders of the intervals where the maximum is obtained.}
}
\references{Davies, P. L. and Kovac, A. (2001) Densities, Spectral Densities and Modality}
\author{Arne Kovac \email{Arne.Kovac@uni-essen.de}}
\seealso{\code{\link{pmden}}}

\examples{
aaa <- rclaw(500)
kkuip(aaa,9)
}
\keyword{univar  }% at least one, from doc/KEYWORDS

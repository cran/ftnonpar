\name{rtennormal}
\alias{rtennormal}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mixture of ten normal distributions}
\description{
Generates a sample from a mixture of ten normal distributions.}
\usage{
rtennormal(n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{Sample size}
}
\author{Arne Kovac \email{Arne.Kovac@uni-essen.de}}
\seealso{\code{\link{rclaw}}}

\examples{
aaa <- rtennormal(500)
kkuip(aaa,9)
}
\keyword{ distribution  }% at least one, from doc/KEYWORDS

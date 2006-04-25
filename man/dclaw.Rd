\name{dclaw}
\alias{dclaw}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{The Claw Distribution}
\description{
Generates a sample from the claw distribution.}
}
\usage{
dclaw(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Vector of points where the claw density is evaluated.}
}
\references{Marron, J. S. and Wand, M. P. (1992) Exact mean integrated squared error}
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{rclaw}}}

\examples{
plot(dclaw(seq(-3,3,len=1000)),ty="l")
}
\keyword{ distribution  }% at least one, from doc/KEYWORDS

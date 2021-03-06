\name{rclaw}
\alias{rclaw}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{The Claw Distribution}
\description{
Generates a sample from the claw distribution.}
\usage{
rclaw(n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{Sample size}
}
\references{Marron, J. S. and Wand, M. P. (1992) Exact mean integrated squared error}
\author{Arne Kovac \email{A.Kovac@bristol.ac.uk}}
\seealso{\code{\link{dclaw}}}

\examples{
aaa <- rclaw(500)
kkuip(aaa,9)
}
\keyword{ distribution  }% at least one, from doc/KEYWORDS

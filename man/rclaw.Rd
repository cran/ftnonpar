\name{rclaw}
\alias{rclaw}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{The Claw Distribution}
\description{
Generates a sample from the claw distribution.}
}
\usage{
rclaw(n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{Sample size}
}
J.S. Marron and M.P. Wand.
\newblock Exact mean integrated squared error.
\newblock {\em Annals of Statistics}, 20:712--736, 1992.
\references{Marron, J. S. and Wand, M. P. (1992) Exact mean integrated squared error}
\author{Arne Kovac \email{Arne.Kovac@uni-essen.de}}
\seealso{\code{\link{dclaw}}}

\examples{
aaa <- rclaw(500)
kkuip(aaa,9)
}
\keyword{ distribution  }% at least one, from doc/KEYWORDS

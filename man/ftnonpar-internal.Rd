\name{ftnonpar-internal}
\alias{genmrcheck}
\alias{genstring}
\alias{kuipdiffbounds}
\alias{kuipdiffbounds.x}
\alias{multiwdwr}
\alias{tautstring}
\title{Internal ftnonpar functions}
\description{Internal ftnonpar functions}
\usage{
genmrcheck(res,thresh=-1,sigma=1,DYADIC=FALSE,beta=0.5,method=1,thr.const=2)
genstring(y,lambda,beta=0.5,method=1)
kuipdiffbounds
kuipdiffbounds.x
multiwdwr(y, thresh,firstwidth=1)
tautstring(ttt, fdist, y.low, y.up, y1 = 0.5 * (y.low[1] + y.up[1]),
    yn = 0.5 * (y.low[length(x)] + y.up[length(x)]),extrmean=TRUE)
}
\details{These are not to be called by the user}
\keyword{internal}

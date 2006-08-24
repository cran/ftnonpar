\name{ftnonpar-internal}
\alias{genmrcheck}
\alias{genstring}
\alias{settomedian}
\alias{settomean}
\alias{kuipdiffbounds}
\alias{kuipdiffbounds.x}
\alias{multiwdwr}
\alias{gendensgn}
\alias{gensign}
\alias{nondymwdr}
\alias{smqden}
\alias{smqnew}
\alias{findmod}
\alias{tautstring}
\alias{isoanti}
\title{Internal ftnonpar functions}
\description{Internal ftnonpar functions}
\usage{
genmrcheck(y,yhat,thresh=-1,sigma=1,DYADIC=FALSE,beta=0.5,method=1,thr.const=2,schwelle = 1e+09)
genstring(y,lambda,beta=0.5,method=1)
kuipdiffbounds
kuipdiffbounds.x
multiwdwr(y, thresh,firstwidth=1,dyad.factor=1.1)
tautstring(ttt, fdist, y.low, y.up, y1 = 0.5 * (y.low[1] + y.up[1]),
    yn = 0.5 * (y.low[length(x)] + y.up[length(x)]),extrmean=TRUE)
settomedian(f, y)
settomean(f, y)
gensign(y,...)
gendensgn(x,...)
smqden(x,lambda,eps=1,fsign=double(0))
smqnew(y,lambda,eps=1,fsign,tolerance = 1e-12)
nondymwdr(y, thresh, firstwidth = 1)
findmod(y,tol=1e-05)
isoanti(x,fhat)
}
\details{These are not to be called by the user}
\keyword{internal}

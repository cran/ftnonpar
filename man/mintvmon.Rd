\name{mintvmon}
\alias{mintvmon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimization of total variation}
\description{
Finds a function vector which minimizes the total variation of 
the function or a derivative under multiresolution constraints
and monotonicity and convexity constraints.
}
\usage{
mintvmon(y, sigma = -1, DYADIC = T, thresh = -1, method = 2, MONCONST = T, CONVCONST = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{observed values (ordered by value of independent variable).}
  \item{sigma}{if set to a positive value the standard deviation is set to sigma and not estimated from the data}
  \item{DYADIC}{logical, if T (default) the multiresolution constraints are only verifeid on intervals with dyadic endpoints}
  \item{thresh}{if set to a positive value other thresholds for the multiresolution criterion than the default sqrt(2*log(n))*sigma can be used.}
  \item{method}{Number of derivative the total variation of which is minimzed. Possible values are 0,1,2. Higher values lead to numerical inconsistencies.}
  \item{MONCONST}{logical, if T (default) additional monotonicty constraints are gathered from minimzing the total variation of f. Makes only sense, if method is 1 or 2.} 
  \item{CONVCONST}{logical, if T (default) additional convexity constraints are gathered from minimzing the total variation of f'. Makes only sense, if method is 2.}

}
\value{
A list with components
  \item{y        }{The approximation of the given data}           
  \item{derivsign}{Vector of 1 and -1, monotonicty constraints used if
MONCONST was true}
  \item{secsign  }{Vector of 1 and -1, convexity constraints used if
CONVCONST was true}
  \item{jact     }{Left endpoints of active multiresolution constraints
for the final approximation}
  \item{kact     }{Right endpoints of active multiresolution constraints
for the final approximation}
  \item{signact  }{Vector of 1 and -1, gives for each active multiresolution
constraints, if the residuals on that interval attain upper or lower bound}
  \item{pl       }{Left endpoint of piecewise constant intervals of the 
derivative of f being minmized}
  \item{pr       }{Right  endpoint of piecewise constant intervals of the
derivative of f being minmized}
}
\references{Kovac, A. (2003) Minimizing Total Variation under Multiresolution Conditions}
\author{Arne Kovac}

\seealso{\code{\link{pmreg}}}
\examples{
data(djdata)
djdoppler.tv0 <- mintvmon(djdoppler,method=0)
djdoppler.tv1 <- mintvmon(djdoppler,method=1)
djdoppler.tv2 <- mintvmon(djdoppler)
par(mfrow=c(2,2))
plot(djdoppler,col="lightgrey")
plot(djdoppler,col="lightgrey")
lines(djdoppler.tv0$y,col="blue")
plot(djdoppler,col="lightgrey")
lines(djdoppler.tv1$y,col="green")
plot(djdoppler,col="lightgrey")
lines(djdoppler.tv2$y,col="red")
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line


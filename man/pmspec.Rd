\name{pmspec}
\alias{pmspec}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Piecewise monotone spectral density approximation with taut strings}
\description{
Applies the taut string method to spectral densities.
}
\usage{
pmspec(x, pks=0, alpha=0.9, sqzf=0.9, mult=F, lcl=F, ln=F, fig = F, pow=10^(-2))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{data}
  \item{pks}{ if pks=0 then the number of peaks is determined
    automatically. If pks = k >0 then a density with n paeks is returned \code{pks}}
  \item{alpha}{ If pks=0 then on Gaussian test beds the number of peaks of the
    true density is at least the returned value with an asymptotic
    probability of at least alpha. If pks > 0 then the parameter is not
    operational \code{alpha} } 
  \item{sqzf}{ Squeeze factor for the taut string \code{sqzf}}
  \item{mult}{ If mult = T then rescaled empirical density is compared
    with the exponential distribution only on intervals forming a
    multiresolution scheme. If T = F then all intervals are used. If the
    sample size of the data is less than 512 then T = F is default \code{mult}}
  \item{lcl}{If lcl = T then string is calculated using local
    squeezing. If lcl = F then global squeezing is used.\code{lcl}}
  \item{fig}{If fig = T then the data are automticall plotted \code{fig}}
  \item{ln}{If fig = T and ln = T then the densities are plotted on a
    log scale. If ln = F then no transformation is performed \code{ln}}  
  \item{pow}{ Peaks with a power less than pow * total power will be ignored}
}
\details{}
\value{
  \item{edf   }{Emprical spectral density}
  \item{df    }{String density}
  \item{pks   }{Number of peaks}
  \item{ll    }{Lower boundary for string}
  \item{uu    }{Upper boundary for string}
  \item{str   }{String}
}
\references{Davies, P.L. snd Kovac, A. (2001) Local extremes, runs,
  strings and multiresolution, Annals of Statistics 29 (1) 1-65.
  Davies, P. L. and Kovac, A. (2002) Densities, Spectral Densities and Modality}
\author{Laurie Davies \email{Laurie.Davies@uni-essen.de}}
\seealso{\code{\link{pmreg}},\code{\link{l1pmreg}},\code{\link{pmden}}}

\examples{
}

\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ ts     }% __ONLY ONE__ keyword per line


\name{djdata}
\docType{data}
\alias{djdoppler}
\alias{djheavisine}
\alias{djblocks} 
\alias{djbumps}   
\title{Donoho-Johnstone test signals}
\description{
Four samples, each of size 2048 from Donoho and Johnstone's famous
signals commonly used as test beds for smoothing methods.
}
\usage{data(djdata)}
\format{ Four one-dimensional data sets, each of size 2048 }
\source{ Donoho, D. L. et al (1995) Wavelet Shrinkage: Asymptopia?
\emph{Journal of the Royal Statistical Society}
  series B, \bold{57}, 301--337.
}
\examples{
data(djdata)
par(mfrow=c(2,2))
plot(djdoppler)
plot(djheavisine)
plot(djblocks)  
plot(djbumps)
}
\keyword{datasets}

\name{balloon}
\docType{data}
\alias{balloon}   
\title{Data from a weather balloon}
\description{
The data consist of 4984 observations taken from a balloon about 30 kilometres
above the surface of the earth. The outliers are caused by the fact that the
balloon slowly rotates, causing the ropes from which the measuring instrument
is suspended to cut off the direct radiation from the sun. 
}
\usage{data(djdata)}
\format{ A one-dimensional data set of size 4984 }
\source{ Davies, L. and Gather, U. (1993), The Identification of Multiple
   Outliers. 
\emph{JASA}
  \bold{88}, 782--801.
}
\examples{
data(balloon)
plot(balloon)
}
\keyword{datasets}

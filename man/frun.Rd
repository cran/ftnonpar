\name{frun}
\alias{frun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Runs and local extremes}
\description{
Calculation of bounds for functions such that the residuals
satisfy a run criterion}
\usage{
frun(y, ..., alpha = 0.5, r = 0, mr = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the data}
  \item{\dots}{an optional argument which specifies
the approximate positions of the local extreme values. 
These should be consisten with the run length otherwise 
the result will be nonsensical. Should you wish to use
this option then you should first run the macro
without it. The item xb of the output list gives 
the acceptable limits of the local extreme values. You
can then specify the positions within these limits.}
  \item{alpha}{Quantile determining the acceptable run length.}
  \item{r}{Acceptable run length: Overrides alpha if not 0}
  \item{mr}{mr=0 minimizes the run length consistent with the
number of local extreme values found for the specified 
run length. mr=1 disables the option.}
}
\value{
  \item{l1    }{lower bounds}
  \item{u1    }{upper bounds}
  \item{l2    }{lower bound with specified extremes: the default choices for
the positions of the local extreme values are the mid-points
of the intervals specified by xb above.}
  \item{u2    }{upper bound with specified extremes}
  \item{f     }{function between l2 and u2 satisfying run condition}
  \item{xb    }{bounds for location of extremes: the position of the ith
extreme value lies between xb[2*i-1] and xb[2*i]}
  \item{nx    }{number of extremes}
  \item{r     }{run length: may differ from specified if mr=1}
}
\references{Davies, P. L. (1995) Data features. Statistica Neerlandica 49,185-245.

Davies, P. L. and Kovac, A. (2001) Local Extremes, Runs, Strings and Multiresolution (with discussion) Annals of Statistics. 29. p1-65}

\author{Laurie Davies \email{Laurie.Davies@uni-essen.de}}
\note{IN GENERAL THE MEAN OF THE BOUNDS l2 AND u2 (l2+u2)/2
GIVES A BETTER REGRESSION FUNCTION THAN f. HOWEVER THIS
FUNCTION IS INFINITE AT THE TWO ENDPOINTS AND AT LOCAL EXTREME
VALUES. IN THESE INTERVALS IT CAN BE REPLACED BY ANY VALUES
WHICH DO NOT ALTER THE NUMBER OF LOCAL EXTREME VALUES. THE
MEDIAN OF THE y-VALUES IN THESE INTERVALS IS A REASONABLE 
CHOICE. }

\seealso{\code{\link{mintvmon}},\code{\link{pmreg}},\code{\link{l1pmreg}}}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line


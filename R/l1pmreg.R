"quantpmreg" <-
function(y,beta=0.5,squeezing.factor=0.5,verbose=FALSE,localsqueezing=TRUE,DYADIC=TRUE,thr.const=2,extrema.nr = -1,bandwidth= -1)
{
n<-length(y)
sigma <- 1

if (bandwidth < 0)
        firstlambda <- 2^(floor(log(length(y), base = 2)) - 1)
    else firstlambda <- bandwidth

lambda <- rep(firstlambda,n-1)
currprecision <- firstlambda

while(1<2)
  {
  tmp <- genstring(y,lambda,beta,1)
  y.string <- tmp$y
  y.mr <- genmrcheck(y-y.string,sigma=sigma,beta=beta,method=1,DYADIC=DYADIC,thr.const=thr.const)
  if(verbose)
    {
    par(mfrow=c(1,2))
    plot(y,col="lightgrey")
    lines(y.string,col="red")
    plot(lambda,ylim=c(-max(lambda),max(lambda)),ty="l")
    lines(-lambda)
    lines(cumsum(sign(y-y.string)),col="red")
    print(c("lambda=",min(lambda)))
    print("Press Enter")
    readline()
    }
  if (bandwidth > 0)
     break
  if (extrema.nr > 0) {
    if (tmp$kext > extrema.nr)
       lambda <- lambda + currprecision
            if (currprecision < 0.5) {
                if (tmp$kext <= extrema.nr)
                  break
            }
            else {
                currprecision <- currprecision/2
                lambda <- lambda - currprecision
            }
        }
        else {

  if(sum(y.mr)<0.1)
    break
  if(!localsqueezing)
    lambda <- lambda*squeezing.factor
  else
    lambda[y.mr[-1] | y.mr[-n]] <- lambda[y.mr[-1] | y.mr[-n]] * squeezing.factor
}

  }
list(y=y.string,lambda=lambda,nmax=tmp$kext)
}
"genstring" <-
function(y,lambda,beta=0.5,method=1)
{
n <- length(y)
if(length(lambda)==1)
  lambda <- c(rep(lambda,n-1),0)
else
if(length(lambda)==n-1)
  lambda <- c(lambda,0)
if(method == 1)
  {
  tmp <- sort(y)
  eps <- min(c(0.00001*(max(y)-min(y)),min(tmp[-1] - tmp[-n])/4))
  if(eps < 1e-36)
    {
    if(mad(y[-1]-y[-n])>0)
      y <- y + rnorm(n,0,0.001*mad(y[-1]-y[-n]))
    else
      y <- y + rnorm(n,0,1e-12)
    tmp <- sort(y)
    eps <- min(tmp[-1] - tmp[-n])/4
    }
  }
else
  eps <- 0

print(c("eps is",eps))
tmp <- .C("genstring",y=as.double(y),as.integer(length(y)),as.double(lambda),
as.double(beta),as.integer(method),as.double(eps),kext=as.integer(0),PACKAGE="ftnonpar")

list(y=tmp$y,kext=tmp$kext)
}
"genmrcheck" <-
function(res,thresh=-1,sigma=1,DYADIC=FALSE,beta=0.5,method=1,thr.const=2)
{
if(method == 1)
  {
  res[res> 1e-04] <- beta
  res[res< -1e-04] <- -1+beta
  sigma <- sqrt(beta*(1-beta))

  currq <- pnorm(sqrt(thr.const*log(length(res))))
  tmp <- ceiling(qbinom(1-currq,1:length(res),beta) )
  upperthresh <- -tmp+beta*(1:length(res))
  tmp <- floor(qbinom(currq,1:length(res),beta)) 
  lowerthresh <- -tmp+beta*(1:length(res))
  }
else
  {
  upperthresh <- sqrt(thr.const*log(length(res)))*sigma*sqrt(1:length(res))
  lowerthresh <- -sqrt(thr.const*log(length(res)))*sigma*sqrt(1:length(res))
  }

  .C("genmrcheck",res=as.double(c(0,cumsum(res))),as.integer(length(res)),as.double(lowerthresh),as.double(upperthresh),as.integer(DYADIC),PACKAGE="ftnonpar")$res[1:length(res)]
}


"l1pmreg" <-
function(y,squeezing.factor=0.5,verbose=FALSE,localsqueezing=TRUE,DYADIC=TRUE,thr.const=2,extrema.nr = -1,bandwidth= -1)
{
quantpmreg(y,0.5,squeezing.factor,verbose,localsqueezing,DYADIC,thr.const,extrema.nr,bandwidth)
}

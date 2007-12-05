#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <R_ext/Utils.h>
     
void R_CheckUserInterrupt(void),fastrank(),sortx();

int dcomp(double *a, double *b)
  {
  if(*a>*b) return 1;
  else if(*a==*b) return 0;
  else return -1;
  }

void dqsort(double *x, int n)
  {
  qsort(x,n,sizeof(double),(void *)(*dcomp));
  }

double genM(int left, int right, double x, int sign, double *y, double beta, int method, double eps,int *orign, int *logn, double *sortedx, int *zwpo)
{
int i,l,n,copyl;
double sum,*tmpy,out,f,fprime,outnew;

n = right - left + 1;

if(method == 1) /* quantile regression */
  {
  if(sign == 1)
    l = (int)(0.1 + floor(1.5+0.5*(x + 2.0*(n - 1)*beta)));
  else
    l = (int)(0.1 + ceil(0.5+0.5*(x + 2.0*(n - 1)*beta)));
  
  if(l<1)
    out = -1e+38;
  else if(l>n)
    out = 1e+38;
  else
    {
    copyl=l;
    fastrank(sortedx,orign,logn,&copyl,&left,&right,zwpo,&out);
    out = out + (x + 2.0*(n-1)*beta-2*(l-1))*eps;

/*

    tmpy = malloc(n*sizeof(double));
    for(i=0;i<n;i++)
      tmpy[i] = y[left+i-1];
    dqsort(tmpy,n);

    out = tmpy[l-1] + (x + 2.0*(n-1)*beta-2*(l-1))*eps;
if(fabs(out-outnew)>1e-05)
printf("%f %f %d left:%d right:%d l:%d\n",out,outnew,sign,left,right,l);
    free(tmpy);
*/
    }
  }
else
if(method == 2) /* usual taut string */
  {
  sum=x;
  for(i=0;i<n;i++)
    sum += y[left+i-1];
  out = sum/n;
  }
/*
else 
if(method == 3)*/ /* logistic regression */
/*
  {
  sum = 0;
  for(i=0;i<n;i++)
    sum+=y[left+i-1];

  if(sum<0.001)
    {
    f=0;
    if(x<n)
      out = 0;
    else
      out = 1 - n/x;
    }
  else
  if(sum>n-0.001)
    {
    f=0;
    if(x > -n)
      out = 1;
    else
      out=-n/x;
    }
  else
    {
    out = 0.5;
    f = 100;
    while(fabs(f)>1e-08)
      {
      f = -sum/out + (n-sum)/(1-out)-x;
      fprime = sum/(out*out) + (n-sum)/((1-out)*(1-out));
      if(out - f/fprime<1e-16)
        out = 0.5*out;
      else
      if(out - f/fprime>1-1e-16)
        out = 0.5*(1+out);
      else
        out = out-f/fprime;
      }
    }
  }
*/
else /* Poisson regression */
  {
  sum=x;
  for(i=0;i<n;i++)
    sum += y[left+i-1];
  out = sum/n;
  }
 
return(out);
}

void genstring(double *y, int *n, double *lambda, double *beta, int *method, double *eps, int *kext)
  {
  int *Lowerleft, *Lowerright, *Upperleft, *Upperright, logn=ceil(log((double)(*n))/log(2.0))+1.0,*zwpo;
  double *lhat, *uhat, *L, *U,genM(),*sortedx;
  int a,b,c,K,i,j,lastsign=0;

  Lowerleft = malloc(*n*sizeof(int));
  Lowerright = malloc(*n*sizeof(int));
  Upperleft = malloc(*n*sizeof(int));
  Upperright = malloc(*n*sizeof(int));
  lhat = malloc(*n*sizeof(double));
  uhat = malloc(*n*sizeof(double));
  U = malloc(*n*sizeof(double));
  L = malloc(*n*sizeof(double));

  if(*method == 1)
    {
    sortedx = malloc(*n * logn * sizeof(double));
    zwpo = malloc(logn*sizeof(int));
    sortx(y,sortedx,n,&logn);
    zwpo[0]=1;
    for(i=1;i<logn;i++)
      zwpo[i]= zwpo[i-1]*2;
    }

  /* Part 1 (Initialization) */
  Lowerleft[0] = 1;
  Lowerright[0] = 1;
  Upperleft[0] = 1;
  Upperright[0] = 1;
  lhat[0] = genM(1,1,-lambda[0],-1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
  uhat[0] = genM(1,1,lambda[0],1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
  L[0] = -lambda[0];
  U[0] = lambda[0];
  a = 1;
  b = 1;
  c = 0;
  *kext=0;

  /* Part 2 (Induction) */
  for(K=2;K<=*n;K++)
    {
    /* Step 1 (K), Update of (Lower, lhat, L) */
    Lowerleft[a] = K;
    Lowerright[a] = K;
    lhat[a] = genM(K,K,-lambda[K-1]-L[a-1],-1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
    L[a] = -lambda[K-1];
    a++;

    while((a > c+1) && (lhat[a-1] >= lhat[a-2]))
      {
      Lowerright[a-2] = Lowerright[a-1];
      if(a>2)
        {
        lhat[a-2] = genM(Lowerleft[a-2],Lowerright[a-2],-lambda[K-1]-L[a-3],-1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
        L[a-2]=-lambda[K-1];
        }
      else
        {
        lhat[0] = genM(Lowerleft[a-2],Lowerright[a-2],-lambda[K-1],-1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
        L[0]=-lambda[K-1];
        }
      a--;
      }

    /* Step 2 (K), Update of (Upper, uhat, U) */
    Upperleft[b] = K;
    Upperright[b] = K;
    uhat[b] = genM(K,K,lambda[K-1]-U[b-1],1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
    U[b] = lambda[K-1];
    b++;

    while((b > c+1) && (uhat[b-1] <= uhat[b-2]))
      {
      Upperright[b-2] = Upperright[b-1];
      if(b>2)
        {
        uhat[b-2] = genM(Upperleft[b-2],Upperright[b-2],lambda[K-1]-U[b-3],1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
        U[b-2]=lambda[K-1];
        }
      else
        {
        uhat[0] = genM(Upperleft[b-2],Upperright[b-2],lambda[K-1],1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
        U[0]=lambda[K-1];
        }
      b--;
      }

    /* Step 3 (K), Comparing lhat[c+1] and uhat[c+1] */
    if((a == c+1) && (lhat[a-1] > uhat[c]))
      while(lhat[a-1] > uhat[c])
        {
        Lowerleft[c] = Upperleft[c];
        Lowerleft[c+1] = Upperright[c]+1;
        Lowerright[c+1] = Lowerright[a-1];
        Lowerright[c] = Upperright[c];
        lhat[c]=uhat[c];
        lhat[c+1]=genM(Lowerleft[a],Lowerright[a],-lambda[K-1]-U[c],-1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
        L[c]=U[c];
        L[c+1]=-lambda[K-1];
        a++;
        c++;
        if(lastsign==-1) (*kext)++;
        lastsign=1;
        }
    else
    if((b == c+1) && (uhat[b-1] < lhat[c]))
      while(uhat[b-1] < lhat[c])
        {
        Upperleft[c] = Lowerleft[c];
        Upperleft[c+1] = Lowerright[c]+1;
        Upperright[c+1] = Upperright[b-1];
        Upperright[c] = Lowerright[c];
        uhat[c]=lhat[c];
        uhat[c+1]=genM(Upperleft[b],Upperright[b],lambda[K-1]-L[c],1,y,*beta,*method,*eps,n,&logn,sortedx,zwpo);
        U[c]=L[c];
        U[c+1]=lambda[K-1];
        b++;
        c++;
        if(lastsign==1) (*kext)++;
        lastsign=-1;
        }
    }

  /* Part 3 (Definition of fhat) */

  if(c > 0)
    for(j=1;j<=c;j++)
      for(i=Lowerleft[j-1];i<=Lowerright[j-1];i++)
        y[i-1] = lhat[j-1];
  if(c == 0)
    for(i=1;i<=*n;i++)
      y[i-1] = (lhat[c]+uhat[c])/2;
  else
  if(L[c-1] == lambda[Lowerright[c-1]])
    for(i=Lowerleft[c];i<=*n;i++)
      y[i-1] = lhat[c];
  else
    for(i=Lowerleft[c];i<=*n;i++)
      y[i-1] = uhat[c];

  if(*method == 1)
    {
    free(sortedx);
    free(zwpo);
    }
  free(Lowerleft);
  free(Lowerright);
  free(Upperleft);
  free(Upperright);
  free(lhat);
  free(uhat);
  free(U);
  free(L);
  }


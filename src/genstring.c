#include "math.h"
#include "stdlib.h"
#include "stdio.h"

int dcomp(const double *a, const double *b)
  {
  if(*a>*b) return 1;
  else if(*a==*b) return 0;
  else return -1;
  }

void dqsort(double *x, int n)
  {
  qsort(x,n,sizeof(double),(void *)(*dcomp));
  }

double genM(int left, int right, double x, int sign, double *y, double beta, int method, double eps)
{
int i,l,n;
double sum,*tmpy,out;

n = right - left + 1;

if(method == 1)
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
    tmpy = malloc(n*sizeof(double));
    for(i=0;i<n;i++)
      tmpy[i] = y[left+i-1];
    dqsort(tmpy,n);

    out = tmpy[l-1] + (x + 2.0*(n-1)*beta-2*(l-1))*eps;
    free(tmpy);
    }
  }
else
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
  int *Lowerleft, *Lowerright, *Upperleft, *Upperright;
  double *lhat, *uhat, *L, *U,genM();
  int a,b,c,K,i,j,lastsign=0;

  Lowerleft = malloc(*n*sizeof(int));
  Lowerright = malloc(*n*sizeof(int));
  Upperleft = malloc(*n*sizeof(int));
  Upperright = malloc(*n*sizeof(int));
  lhat = malloc(*n*sizeof(double));
  uhat = malloc(*n*sizeof(double));
  U = malloc(*n*sizeof(double));
  L = malloc(*n*sizeof(double));

  /* Part 1 (Initialization) */
  Lowerleft[0] = 1;
  Lowerright[0] = 1;
  Upperleft[0] = 1;
  Upperright[0] = 1;
  lhat[0] = genM(1,1,-lambda[0],-1,y,*beta,*method,*eps);
  uhat[0] = genM(1,1,lambda[0],1,y,*beta,*method,*eps);
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
    lhat[a] = genM(K,K,-lambda[K-1]-L[a-1],-1,y,*beta,*method,*eps);
    L[a] = -lambda[K-1];
    a++;

    while((a > c+1) && (lhat[a-1] >= lhat[a-2]))
      {
      Lowerright[a-2] = Lowerright[a-1];
      if(a>2)
        {
        lhat[a-2] = genM(Lowerleft[a-2],Lowerright[a-2],-lambda[K-1]-L[a-3],-1,y,*beta,*method,*eps);
        L[a-2]=-lambda[K-1];
        }
      else
        {
        lhat[0] = genM(Lowerleft[a-2],Lowerright[a-2],-lambda[K-1],-1,y,*beta,*method,*eps);
        L[0]=-lambda[K-1];
        }
      a--;
      }

    /* Step 2 (K), Update of (Upper, uhat, U) */
    Upperleft[b] = K;
    Upperright[b] = K;
    uhat[b] = genM(K,K,lambda[K-1]-U[b-1],1,y,*beta,*method,*eps);
    U[b] = lambda[K-1];
    b++;

    while((b > c+1) && (uhat[b-1] <= uhat[b-2]))
      {
      Upperright[b-2] = Upperright[b-1];
      if(b>2)
        {
        uhat[b-2] = genM(Upperleft[b-2],Upperright[b-2],lambda[K-1]-U[b-3],1,y,*beta,*method,*eps);
        U[b-2]=lambda[K-1];
        }
      else
        {
        uhat[0] = genM(Upperleft[b-2],Upperright[b-2],lambda[K-1],1,y,*beta,*method,*eps);
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
        lhat[c+1]=genM(Lowerleft[a],Lowerright[a],-lambda[K-1]-U[c],-1,y,*beta,*method,*eps);
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
        uhat[c+1]=genM(Upperleft[b],Upperright[b],lambda[K-1]-L[c],1,y,*beta,*method,*eps);
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

  free(Lowerleft);
  free(Lowerright);
  free(Upperleft);
  free(Upperright);
  free(lhat);
  free(uhat);
  free(U);
  free(L);
  }

void genmrcheck(double *intres, int *n, double *lowerthresh, double *upperthresh, int *DYADIC)
  {
  int i,j,k,length,*ergebnis,currright;

  ergebnis = malloc(*n*sizeof(int));

  for(i=0;i<*n;i++)
    ergebnis[i]=-1;

  length=1;
  while(length <= *n)
    {
    j=0;
    while(j < *n)
      {
      k=j+length-1;
      if(k>=*n)
        k = *n-1;

      if((intres[k+1]-intres[j] > 1e-06+upperthresh[length-1])||(intres[k+1]-intres[j] +1e-06 < lowerthresh[length-1]))
        {
        ergebnis[j]=k;
        }

      if(*DYADIC)
        j += length;
      else
        j++;
      }
    if(*DYADIC)
      length *= 2;
    else
      length++;
    }

  currright = -1;
  for(i=0;i<*n;i++)
    {
    if(ergebnis[i]>currright)
      currright = ergebnis[i];
    if(i <= currright)
      intres[i]=1;
    else
      intres[i]=0;
    }

  free(ergebnis);
  }


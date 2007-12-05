#include "math.h"
#include "stdlib.h"
#include "stdio.h"

/* *n is length of lower, upper, etc not string */

void tautstring(double *fdist, double *t, double *lower, double *upper, double *y1, double *yn, int *n, double *string, int *knotsind, double *knotst, double *knotsy, int *nknots,int *nmax,int *EXTRMEAN)
{
double newmaxderiv,newminderiv,maxderiv=1e+38, minderiv=-maxderiv;
int actind=2, maxind=1, minind=1, i,j,lastbound=0;
int *knotssign,lastsign;


*nmax=0;
knotssign=calloc(*n,sizeof(int));
lastbound=0;

knotssign[0]=lastbound;
knotsind[0]=1;
knotsy[0]=*y1; 
knotst[0]=t[0]; 
*nknots=1;

while(actind<=*n)
  {
  if(actind<*n)
    {
    newmaxderiv = (upper[actind-1]-knotsy[*nknots-1])/(t[actind-1]-(knotst[*nknots-1]));
    newminderiv = (lower[actind-1]-knotsy[*nknots-1])/(t[actind-1]-(knotst[*nknots-1]));
    }
  else
    {
/*
    *yn=0.5*(upper[*n-1]+lower[*n-1]);
*/
    newmaxderiv = (*yn-knotsy[*nknots-1])/(t[actind-1]-(knotst[*nknots-1]));
    newminderiv = (*yn-knotsy[*nknots-1])/(t[actind-1]-(knotst[*nknots-1]));
    }
  if(newminderiv > maxderiv)
    {
      if(lastbound==-1) (*nmax)++;
      knotssign[*nknots]=1;
      knotsind[*nknots] = maxind;
      knotsy[*nknots] = upper[maxind-1];
      knotst[*nknots]=t[maxind-1]; 
      (*nknots)++;
    actind = maxind;
    lastbound=1;
    maxderiv = 1e+38;
    minderiv = -1e+38;
    }
  else
  if(newmaxderiv < minderiv)
    {
      if(lastbound==1) (*nmax)++;
      knotsind[*nknots] = minind;
      knotssign[*nknots]=-1;
      knotsy[*nknots] = lower[minind-1];
      knotst[*nknots]=t[minind-1]; 
      (*nknots)++;
    maxderiv = 1e+38;
    minderiv = -1e+38;
    lastbound=-1;
    actind = minind;
    }
  else
    {
    if(newmaxderiv < maxderiv)
      {
      maxderiv = newmaxderiv;
      maxind = actind;
      }
    if(newminderiv > minderiv)
      {
      minderiv = newminderiv;
      minind = actind;
      }
    if(actind==*n)
      if(lastbound!=0)
	{
        knotsind[*nknots] = actind;
        knotsy[*nknots] = *yn;
        knotst[*nknots]=t[actind-1]; 
        /*knotssign[*nknots]=lastsign; */
        knotssign[*nknots]=0; 
        (*nknots)++;
	}
      else
	{
        lastbound=-1;
        knotsind[0]=1;
        knotssign[0]=-1;
        knotsy[0]=*y1; 
        knotst[0]=t[0]; 
        actind=1;
        minind=1;
        maxind=1;
	}
    }
  actind++;
  }
 
  for(i=0;i<*nknots-1;i++)
      for(j=knotsind[i];j<knotsind[i+1];j++)
        if(knotssign[i]==knotssign[i+1])
          string[j-1]=(knotsy[i+1]-knotsy[i])/(knotst[i+1]-knotst[i]);
        else
          if(*EXTRMEAN)
            string[j-1]=(fdist[knotsind[i+1]-1]-fdist[knotsind[i]-1])/(knotst[i+1]-knotst[i]);
          else
            string[j-1]=(knotsy[i+1]-knotsy[i])/(knotst[i+1]-knotst[i]);
free(knotssign);
}

void multiwdwr(double *y, int *n, double *thresh, int *firstwidth, double *dyadfactor)
{
int j, leftind, rightind;
double *ysum, actwidth;

ysum=malloc((*n+1)*sizeof(double));
ysum[0]=0;
for(j=1;j<=*n;j++)
  ysum[j]=ysum[j-1]+y[j-1];
for(j=0;j<*n;j++)
  y[j]=0.0;

for(actwidth=*firstwidth;actwidth<=*n;actwidth*=(*dyadfactor))
  for(leftind=0,rightind=actwidth;leftind<*n;leftind=rightind,rightind+=actwidth)
    {
    if(rightind>*n) rightind= *n;
    if(fabs((ysum[rightind]-ysum[leftind])/sqrt((double)(rightind-leftind)))> *thresh)
      for(j=leftind;j<rightind;j++)
        y[j]=1.0;
    }

free(ysum);
}


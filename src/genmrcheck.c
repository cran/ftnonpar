#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
 
void qbinsum(double *p, double *probs, double *pout, int *k)
  {
  int i,j;

  pout[0] = 1-probs[0];
  pout[1] = probs[0];

  if(*k>=2)
  for(i=1;i<*k;i++)
    {
    pout[i+1] = pout[i] * probs[i];
    for(j=i;j>0;j--)
      pout[j] = pout[j-1] * probs[i] + pout[j] * (1-probs[i]);
    pout[0] = pout[0] * (1-probs[i]);
    }

  for(i=1;i<=*k;i++)
    pout[i] = pout[i-1] + pout[i];  

  if(*p <= pout[0])
    *p = 0.0;
  else
  if(*p >= 1.0)
    *p = *k;
  else
    for(i=1;i<=*k;i++)
      if(*p <= pout[i])
        {
        *p = i;
        break;
        }
  }

void genmrcheck(double *y, double *yhat, int *n, int *method, int *DYADIC, double *siglevel, int *schwelle,double *beta, double *sigma)
  {
  int i,j,k,length,*ergebnis,currright,currlen;
  double *ysum,*yhatsum,qlow,qup,currc,*wspace,mu,deltalow,deltaup,curry,currp,*uppersum,*lowersum,upperp,lowerp,normthresh;

  normthresh = qnorm(*siglevel,(double)0,(double)1,(int)1,(int)0)*(*sigma);

  if(*method == 1)
    {
    uppersum     = malloc((*n+1)*sizeof(double));
    lowersum     = malloc((*n+1)*sizeof(double));

    uppersum[0]=0;
    lowersum[0]=0;

    for(i=1;i<=*n;i++)
      {
      if(yhat[i-1]-y[i-1]>-1e-04)
        uppersum[i]=uppersum[i-1]+1.0;
      else
        uppersum[i]=uppersum[i-1];
      if(yhat[i-1]-y[i-1]>1e-04)
        lowersum[i]=lowersum[i-1]+1.0;
      else
        lowersum[i]=lowersum[i-1];
      }
    }
  else
    {
    ysum     = malloc((*n+1)*sizeof(double));
    yhatsum     = malloc((*n+1)*sizeof(double));
    wspace     = malloc((*n+1)*sizeof(double));

    ysum[0]=0;
    yhatsum[0]=0;

    for(i=1;i<=*n;i++)
      {
      ysum[i]=ysum[i-1]+y[i-1];
      yhatsum[i]=yhatsum[i-1]+yhat[i-1];
      }
    }

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

      if(*method == 1)
        {
        upperp = pbinom(uppersum[k+1] - uppersum[j],(double)(k-j+1.0),*beta, (int)1, (int)0);
        lowerp = pbinom(lowersum[k+1] - lowersum[j]-1e-02,(double)(k-j+1.0),*beta, (int)1, (int)0);

        if((lowerp > *siglevel)||(upperp < 1.0 - *siglevel))
/*
{
printf("%d %d %f(%f) %f(%f)\n",j,k,uppersum[k+1] - uppersum[j],upperp,lowersum[k+1] - lowersum[j]-1e-02,lowerp);
*/
          ergebnis[j]=k;
/*
}
*/
        }
      else
      if(*method == 2)
        {
        curry = ysum[k+1] - ysum[j]-yhatsum[k+1]+yhatsum[j];

        if(fabs(curry) > sqrt(k-j+1.0)*normthresh)
          ergebnis[j]=k;
        }
      else
      if(*method == 3)
        {
        currlen = k+1-j;
 
        if(currlen <= *schwelle)
          {
          qup  = *siglevel;
          qlow = 1-*siglevel;
          qbinsum(&qlow,&(yhat[j]),wspace,&currlen);
          qbinsum(&qup ,&(yhat[j]),wspace,&currlen);
          }
        else
          {
          mu = yhatsum[k+1]-yhatsum[j];

          deltalow = sqrt(-2.0*log(1 - *siglevel)/mu);
          deltaup  = sqrt(-3.0*log(1 - *siglevel)/mu);
          if(deltalow>1)
            qlow = 0;
          else
            qlow = (1.0 - deltalow) * mu;

          if(deltaup >1)
            qup  = 1e+38;
          else
            qup  = (1.0 + deltalow) * mu;
          }

        if((ysum[k+1]-ysum[j]>qup+1e-08)||(ysum[k+1]-ysum[j]<qlow-1e-08))
          ergebnis[j]=k;
        }
      else
      /* if(*method == 4) */
        {
        upperp = ppois(ysum[k+1] - ysum[j],yhatsum[k+1] - yhatsum[j], (int)1, (int)0);
        lowerp = ppois(ysum[k+1] - ysum[j]-0.2,yhatsum[k+1] - yhatsum[j], (int)1, (int)0);

        if((lowerp > *siglevel)||(upperp < 1.0 - *siglevel))
{
/*
printf("%d %d y:%f yhat:%f\n",j,k,ysum[k+1] - ysum[j],yhatsum[k+1] - yhatsum[j]);
*/
          ergebnis[j]=k;
}
/*
        currp = ppois(ysum[k+1] - ysum[j],yhatsum[k+1] - yhatsum[j], (int)1, (int)0);
        

        if((currp > *siglevel)||(currp < 1.0 - *siglevel))
{
printf("%d %d y:%f yhat:%f\n",j,k,ysum[k+1] - ysum[j],yhatsum[k+1] - yhatsum[j]);
          ergebnis[j]=k;
}
*/

        }

      if(*DYADIC)
        j += length;
      else
        {
        j++;
        if(j+length-1>=*n)
          j=*n;
        }
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
      y[i]=1;
    else
      y[i]=0;
    }

  free(ergebnis);
  if(*method == 1)
    {
    free(uppersum);
    free(lowersum);
    }
  else
    {
    free(ysum);
    free(yhatsum);
    free(wspace);
    }
  }



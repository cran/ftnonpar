#include "math.h"
#include "stdlib.h"

void easymax(),difficultmax(),kkuip();

void kkuip(double *x, int *n, int *kmax, double *norm, int *a, int *b)
  {
  double min,max,abmax,newmax;
  int i,k,mini,maxi,maxa,maxb,newa,newb,einordnen ;

  min=1e+38;
  max=-1e+38;
  for(i=0;i<*n;i++)
    {
    if(x[i]>max)
      {
      max=x[i];
      maxi=i;
      }
    if(x[i]<min)
      {
      min=x[i];
      mini=i;
      }
    }
  if(mini<maxi)
    {
    a[0]=mini;
    b[0]=maxi;
    }
  else
    {
    b[0]=mini;
    a[0]=maxi;
    }
  norm[0]=fabs(x[b[0]]-x[a[0]]);
  for(k=1;k<*kmax;k++)
    {
    abmax=-1e+38;
    einordnen =-1;
    if(a[0]>0)
      {
      easymax(x,n,0,a[0],&maxa,&maxb,&abmax);
      einordnen =0;
      }
    for(i=0;i<k;i++)
      {
      difficultmax(x,n,a[i],b[i],&newa,&newb,&newmax);
      if(newmax>abmax)
	{
	maxa=newa;
	maxb=newb;
	abmax=newmax;
	einordnen =i;
	}
      if(i<k-1)
	{
        easymax(x,n,b[i],a[i+1],&newa,&newb,&newmax);
        if(newmax>abmax)
	  {
	  maxa=newa;
	  maxb=newb;
	  abmax=newmax;
	  einordnen =i+1;
	  }
	}
      else
	{
        easymax(x,n,b[i],*n-1,&newa,&newb,&newmax);
        if(newmax>abmax)
	  {
	  maxa=newa;
	  maxb=newb;
	  abmax=newmax;
	  einordnen =i+1;
	  }
	}
     }
   if(einordnen <0)
     exit(0L);
   else
   if(einordnen ==k)
     {
     a[k]=maxa;
     b[k]=maxb;
     norm[k]=norm[k-1]+fabs(x[maxb]-x[maxa]);
     }
   else
     {
     for(i=k;i>=einordnen+1;i--)
       {
       a[i]=a[i-1];
       b[i]=b[i-1];
       }
     if(maxa<a[einordnen])
       {
       a[einordnen]=maxa;
       b[einordnen]=maxb;
       norm[k]=norm[k-1]+fabs(x[maxb]-x[maxa]);
       }
     else
       {
       b[einordnen]=maxa;
       a[einordnen+1]=maxb;
       norm[k]=norm[k-1]-fabs(x[a[einordnen]]-x[b[einordnen+1]])+fabs(x[a[einordnen]]-x[b[einordnen]])+fabs(x[a[einordnen+1]]-x[b[einordnen+1]]); 
       }
     }
    }
  }

void easymax(double *x, int *n,int left,int right,int *erga,int *ergb,double *ergmax)
  {
    double min,max;
  int i,mini,maxi;

  min=1e+38;
  max=-1e+38;
  for(i=left;i<=right;i++)
    {
    if(x[i]>max)
      {
      max=x[i];
      maxi=i;
      }
    if(x[i]<min)
      {
      min=x[i];
      mini=i;
      }
    }
  if(mini<maxi)
    {
    *erga=mini;
    *ergb=maxi;
    }
  else
    {
    *ergb=mini;
    *erga=maxi;
    }
  *ergmax=max-min;
  }

void difficultmax(double *x, int *n,int left,int right,int *erga,int *ergb,double *ergmax)
  {
    double min,max,*mins,*maxs;
  int i,maxi,*minis,*maxis;

  mins=malloc(*n*sizeof(double));
  maxs=malloc(*n*sizeof(double));
  minis=malloc(*n*sizeof(int));
  maxis=malloc(*n*sizeof(int));

  min=1e+38;
  max=-1e+38;

  if(x[left]<x[right])
    {
      maxs[left]=x[left];
      maxis[left]=left;
      for(i=left+1;i<=right;i++)
        if(x[i]>maxs[i-1])
          {
          maxs[i]=x[i];
          maxis[i]=i;
          }
        else
          {
          maxs[i]=maxs[i-1];
          maxis[i]=maxis[i-1];
          }
      mins[right]=x[right];
      minis[right]=right;
      for(i=right-1;i>=left;i--)
        if(x[i]<mins[i+1])
          {
          mins[i]=x[i];
          minis[i]=i;
          }
        else
          {
          mins[i]=mins[i+1];
          minis[i]=minis[i+1];
          }
      max=-1e+38;
      for(i=left;i<=right;i++)
        if(maxs[i]-mins[i]>max)
          {
	    max=maxs[i]-mins[i];
            maxi=i;
          }
      *erga=maxis[maxi];
      *ergb=minis[maxi];
      *ergmax=fabs(x[*erga]-x[left])+fabs(x[right]-x[*ergb])-fabs(x[right]-x[left]);
    }
  else
    {
      mins[left]=x[left];
      minis[left]=left;
      for(i=left+1;i<=right;i++)
        if(x[i]<mins[i-1])
          {
          mins[i]=x[i];
          minis[i]=i;
          }
        else
          {
          mins[i]=mins[i-1];
          minis[i]=minis[i-1];
          }
      maxs[right]=x[right];
      maxis[right]=right;
      for(i=right-1;i>=left;i--)
        if(x[i]>maxs[i+1])
          {
          maxs[i]=x[i];
          maxis[i]=i;
          }
        else
          {
          maxs[i]=maxs[i+1];
          maxis[i]=maxis[i+1];
          }
      max=-1e+38;
      for(i=left;i<=right;i++)
        if(maxs[i]-mins[i]>max)
          {
	    max=maxs[i]-mins[i];
            maxi=i;
          }
      *erga=minis[maxi];
      *ergb=maxis[maxi];
      *ergmax=fabs(x[*erga]-x[left])+fabs(x[right]-x[*ergb])-fabs(x[right]-x[left]);
    }
  free(minis);
  free(maxis);
  free(mins);
  free(maxs);
}

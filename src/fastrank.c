#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void sortx(double *x, double *sortedx, int *n, int *logn)
  {
  int zwpo,i,j,k1,k2,k,oldind1,oldind2;

  memcpy(sortedx,x,*n*sizeof(double));
  i=1;
  zwpo=1;
  while(i<*logn)
    {
    zwpo=zwpo*2;
    i++;

    for(j=0;j<*n;j+=zwpo)
      {
      k1=0;
      k2=zwpo/2;
      oldind1 = (i-2)*(*n) + j + k1;
      oldind2 = (i-2)*(*n) + j + k2;

      for(k=0;(k<zwpo)&&(j+k<*n);k++)
        if((k1==zwpo/2)||((k2<zwpo)&&(j+k2<*n)&&(sortedx[oldind2]<sortedx[oldind1])))
          {
          sortedx[(i-1)*(*n)+j+k]=sortedx[oldind2];
          k2++;
          oldind2 = (i-2)*(*n) + j + k2;
          }         
        else
          {
          sortedx[(i-1)*(*n)+j+k]=sortedx[oldind1];
          k1++;
          oldind1 = (i-2)*(*n) + j + k1;
          }         
      }
    }
  }

void fastrank(double *sortedx, int *n, int *logn, int *actxrank, int *left, int *right, int *zwpo, double *result)
  {
  double minq,maxq,res,**z;
  int i,j,currzwpo,countz,actn,STOPIT=-1,*zlen,*xranks;
  int maxind,minind,sumxranks;

  z = malloc(*n * sizeof(double *));
  zlen = malloc(*n * sizeof(int));
  xranks = malloc(*n * sizeof(int));

  if(*actxrank<1) 
    {
    *result=-1e+38;
    return;
    }
  if(*actxrank>*right-*left+1) 
    {
    *result=1e+38;
    return;
    }


  i=*left-1;
  j=*right-1;
  countz=0;

  /* Partition the given interval into subintervals with
     dyadic end points */
  currzwpo=0;
  while(i+zwpo[currzwpo]<j+1)
    {
    if(i%zwpo[currzwpo+1]>0)
      {
      z[countz]=&(sortedx[*n*currzwpo+i]);
      zlen[countz]=zwpo[currzwpo];

/*
printf("countz:%d, zlen:%d, currzwpo: %d, i:%d, j:%d\n", countz,zlen[countz],currzwpo,i,j);
if(zlen[countz] == 2)
  printf("%f %f\n",z[countz][0],z[countz][1]);
if(zlen[countz] == 2)
  printf("%f %f\n",sortedx[*n*currzwpo+i],sortedx[*n*currzwpo+i+1]);
*/
      i+=zwpo[currzwpo];
      countz++;
      }     
    currzwpo++;
    }
  while(i<=j)
    {
    if(i+zwpo[currzwpo]<=j+1)
      {
      z[countz]=&(sortedx[*n*currzwpo+i]);
      zlen[countz]=zwpo[currzwpo];
/*
printf("countz:%d, zlen:%d, currzwpo: %d, i:%d, j:%d\n", countz,zlen[countz],currzwpo,i,j);
if(zlen[countz] == 2)
  printf("%f %f\n",z[countz][0],z[countz][1]);
*/
      i+=zwpo[currzwpo];
      countz++;
      }     
    currzwpo--;
    }

  actn=*right-*left+1;




  while(STOPIT<0)
  {

/*
puts("HALLO");
printf("actn: %d, countz: %d, actxrank: %d\n",actn,countz,*actxrank);
for(i=0;i<countz;i++)
  {
  for(j=0;j<zlen[i];j++)
    printf("%f ",z[i][j]);
  puts("");
  }
*/

  /* If there is only one field, compute the rank immediately */
  if(countz==1)
    {
    res=z[0][*actxrank-1];
    STOPIT=1;
    }
  else
  if(*actxrank==1)
    {
    minq=1e+38;
    for(i=0;i<countz;i++)
      if(z[i][0]<minq)
        minq=z[i][0];
    res=minq;
    STOPIT=1;
    }
  else
  if(*actxrank==actn)
    {
    maxq=-1e+38;
    for(i=0;i<countz;i++)
      if(z[i][zlen[i]-1]>maxq)
        maxq=z[i][zlen[i]-1];
    res=maxq;
    STOPIT=1;
    }
  else
  /* If actxrank is smaller than or equal to countz, then everything is done here.
     In each step actxrank decreases by 1, while countz can only dropy by 2, if
     actxrank < actxrank, and will remain constant or drop by 1 otherwise. */
  if(*actxrank<=countz)
    {
    minq=1e+38;
    maxq=-1e+38;
    for(i=0;i<countz;i++)
      {
      if(z[i][0]<minq)
        {
        minq=z[i][0];
        minind=i;
        }
      if(z[i][0]>maxq)
        {
        maxq=z[i][0];
        maxind=i;
        }
      }
    if(*actxrank<countz)
      {
      actn-=zlen[maxind];
      if(maxind<countz-1)
        {
        z[maxind]=z[countz-1];
        zlen[maxind]=zlen[countz-1];
        }
      if(minind==countz-1) minind=maxind;
      countz--;
      }
    if(zlen[minind]==1)
      {
      if(minind<countz-1)
        {
        z[minind]=z[countz-1];
        zlen[minind]=zlen[countz-1];
        }
      countz--;
      }
    else
      {
      zlen[minind]--;
      z[minind]=(&(z[minind][1]));
      }  
    (*actxrank)--;
    actn--;
    }
  else
  /* Do the same if the rank is too large */
  if(*actxrank>actn-countz)
    {
    minq=1e+38;
    maxq=-1e+38;
    for(i=0;i<countz;i++)
      {
      if(z[i][zlen[i]-1]<minq)
        {
        minq=z[i][zlen[i]-1];
        minind=i;
        }
      if(z[i][zlen[i]-1]>maxq)
        {
        maxq=z[i][zlen[i]-1];
        maxind=i;
        }
      }
    if(*actxrank>actn-countz+1)
      {
      *actxrank-=zlen[minind];
      actn-=zlen[minind];
      if(minind<countz-1)
        {
        z[minind]=z[countz-1];
        zlen[minind]=zlen[countz-1];
        }
      if(maxind==countz-1) maxind=minind;
      countz--;
      }
    actn--;
    if(zlen[maxind]==1)
      {
      if(maxind<countz-1)
        {
        z[maxind]=z[countz-1];
        zlen[maxind]=zlen[countz-1];
        }
      countz--;
      }
    else
      zlen[maxind]-=1;
    }
  /* And now the general case */
  else
    {
    sumxranks=0;
    for(i=0;i<countz;i++)
      {
      xranks[i]=1+zlen[i]*(*actxrank-countz)/actn;
      sumxranks+=xranks[i];
      }
    i=0;
    while(sumxranks<*actxrank)
      {
      if(xranks[i]<zlen[i])
        {
        xranks[i]++;
        sumxranks++;
        }
      i++;
      if(i==countz) i=0;  
      }
    minq=1e+38;
    maxq=-1e+38;
    for(i=0;i<countz;i++)
      {
      if(z[i][xranks[i]-1]<minq)
        {
        minq=z[i][xranks[i]-1];
        minind=i;
        }
      if(z[i][xranks[i]-1]>maxq)
        {
        maxq=z[i][xranks[i]-1];
        maxind=i;
        }
      }
    actn-=(zlen[maxind]-xranks[maxind]);
    zlen[maxind]=xranks[maxind];
    *actxrank-=xranks[minind];
    actn-=xranks[minind];
    if(zlen[minind]==xranks[minind])
      {
      if(minind<countz-1)
        {
        z[minind]=z[countz-1];
        zlen[minind]=zlen[countz-1];
        }
      countz--;
      }
    else
      {
      zlen[minind]-=xranks[minind];
      z[minind]=&(z[minind][xranks[minind]]);
      }  
    }
  }
  *result=res;
  free(xranks);
  free(zlen);
  free(z);
  }


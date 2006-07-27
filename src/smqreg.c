#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <R_ext/Utils.h>

void R_CheckUserInterrupt(void);

void smqnew(double *y, double *f, int *n, double *lambda, double *eps, int *fsign,double *tol)
  {
  int i,l,k,j;
  double z,currsum;
  double ls[*n],us[*n];

  for(l=0;l<*n-1;l++)
    {
    ls[l]=y[l]-100;
    us[l]=y[l]+100;
    }

  for(l=0;l<*n-1;l++)
    {
    if((l>0)&&(fsign[l-1]==1)&&(ls[l]<f[l-1]))
      ls[l] = f[l-1];
    if((l>0)&&(fsign[l-1]==0)&&(us[l]>f[l-1]))
      us[l] = f[l-1];

if(ls[l]>us[l]) puts("ERROR");

    while(us[l]-ls[l]>*tol)
      {
      f[l] = 0.5*(us[l]+ls[l]);
      k = l+1;
      currsum = 0;
      for(j=0;j<k;j++)
        currsum+=(f[j]-y[j]);
      while(k<*n)
        {
        if((fsign[k-1]==1)&&(currsum<0))
          f[k]=f[k-1];
        else if((fsign[k-1]==0)&&(currsum>0))    
          f[k]=f[k-1];
        else
          { 
          z = currsum*2/ lambda[k-1];

          if(z>=1)
            break;        
          else
          if(z<=-1)
            break;       
          else
            {
            if(z>0)
              f[k] = f[k-1] + sqrt(eps[k-1]*z*z/(1-z*z));
            else
              f[k] = f[k-1] - sqrt(eps[k-1]*z*z/(1-z*z));
            }
          }
        currsum+=(f[k]-y[k]);
        k++;
        }
      if(currsum>0)
        {
        for(i=l;i<k;i++)
          if(f[i]<us[i])
            us[i] = f[i];
        }
      else
        {
        for(i=l;i<k;i++)
          if(f[i]>ls[i])
            ls[i] = f[i];
        }
      }
    }
  } 


void smqden(double *x, double *xeval, double *f, int *n, double *lambda, double *eps, int *fsign)
  {
  int i,l,k,j;
  double z,currsum;
  double ls[*n],us[*n];

  for(l=0;l<*n-2;l++)
    {
    ls[l]=0;                 
    us[l]=100;          
    }

  for(l=0;l<*n-2;l++)
    {
    if((l>0)&&(fsign[l-1]==1)&&(ls[l]<f[l-1]))
      ls[l] = f[l-1];
    if((l>0)&&(fsign[l-1]==0)&&(us[l]>f[l-1]))
      us[l] = f[l-1];

if(ls[l]>us[l]) puts("ERROR");

    while(us[l]-ls[l]>1e-12)
      {
      R_CheckUserInterrupt();
      f[l] = 0.5*(us[l]+ls[l]);
      k = l+1;
      currsum = 0;
      for(j=0;j<k;j++)
        currsum+=( f[j]*(x[j+1]-x[j])  - 1.0/(*n-1.0) )  ; /*  before (f[j]-y[j]); */
      while(k<*n-1)
        {
        if((fsign[k-1]==1)&&(currsum<0))
          f[k]=f[k-1];
        else if((fsign[k-1]==0)&&(currsum>0))    
          f[k]=f[k-1];
        else
          { 
          z = currsum*2/ lambda[k-1];

          if(z>=1)
            break;
          else
          if(z<=-1)
            break;
          else
            {
            if(z>0)
              f[k] = f[k-1] + sqrt((x[k]-x[k-1])*eps[k-1]*z*z/(1-z*z));
            else
              f[k] = f[k-1] - sqrt((x[k]-x[k-1])*eps[k-1]*z*z/(1-z*z));
            }
          }
        currsum+=( f[k]*(x[k+1]-x[k])  -      1.0/(*n-1.0) )  ; /*  before (f[k]-y[k]);*/
        k++;    
        }

      if(currsum>0)
        {
        for(i=l;i<k;i++)
          if(f[i]<us[i])
            us[i] = f[i];
        }
      else
        {
        for(i=l;i<k;i++)
          if(f[i]>ls[i])
            ls[i] = f[i];
        }



      }
    }
  } 

void nondymwdwr(double *y, int *n, double *thresh, int *firstwidth)
{
int j, actwidth, leftind, rightind;
double *ysum;
                                                                                
ysum=malloc((*n+1)*sizeof(double));
ysum[0]=0;
for(j=1;j<=*n;j++)
  ysum[j]=ysum[j-1]+y[j-1];
for(j=0;j<*n;j++)
  y[j]=0.0;
                                                                                
for(actwidth=1;actwidth<=*n;actwidth++)
  for(leftind=0,rightind=actwidth;leftind<*n;leftind=rightind,rightind+=actwidth)
    {
    if(rightind>*n) rightind= *n;
    if(fabs((ysum[rightind]-ysum[leftind])/sqrt((double)(rightind-leftind)))> *thresh)
      for(j=leftind;j<rightind;j++)
        y[j]=1.0;
    }
                                                                                
free(ysum);
}
                                                                                


#include "math.h"
#include "stdlib.h"
#include "stdio.h"

void yams(double *y, int *n, double *thresh, int *method)
  {
  double *fsum, *ysum, *f,newmaxf,newminf,maxf,minf,tmpmaxf,tmpminf;
  int i,j,k,m,lastsign,maxind,minind,lastind;
  int *indsum;

  ysum = malloc((*n+1)*sizeof(double));
  fsum = malloc((*n+1)*sizeof(double));
  f = malloc(*n*sizeof(double));
  indsum = malloc((*n+1)*sizeof(int));

  ysum[0] = 0.0;
  fsum[0] = 0.0;
  indsum[0] = 0.0;

  for(i=0;i<*n;i++)
    {
    f[i] = 0.0;
    fsum[i+1] = 0.0;
    indsum[i+1] = 1;
    ysum[i+1] = ysum[i] + y[i];
    }
  for(m=0;m <= *method;m++)
    for(i=1;i <= *n;i++)
      indsum[i] = indsum[i-1]+indsum[i];

  lastsign = 0;
  maxind = 0;
  minind = 0;
  lastind = 0;
 
  maxf = y[0] + *thresh;
  minf = y[0] - *thresh;

  k = 1;

  while(k < *n)
    {
    newmaxf = 1e+38;
    newminf = -1e+38;
    for(j=0;j<=k;j++)
      {
      tmpmaxf=(ysum[k+1]-ysum[j]-fsum[k+1]+fsum[j]+sqrt((double)k-j+1)* *thresh)/(indsum[k+1]-indsum[j]);
      tmpminf=(ysum[k+1]-ysum[j]-fsum[k+1]+fsum[j]-sqrt((double)k-j+1)* *thresh)/(indsum[k+1]-indsum[j]);

      if(tmpmaxf < newmaxf)
        newmaxf = tmpmaxf;

      if(tmpminf > newminf)
        newminf = tmpminf;
      }
    if(newminf > maxf)
      {
      lastsign = 1;
      for(i=lastind;i<=maxind;i++)
        f[i] = maxf;
      if(*method == 0)
        for(i=lastind;i< *n;i++)
          fsum[i+1] = fsum[i] + f[i];
      else
        {
        fsum[0]=0.0;
        for(i=0;i < *n;i++)
          fsum[i+1] = f[i];
        for(m=0;m <= *method;m++)
          for(i=1;i <= *n;i++)
            fsum[i] = fsum[i-1]+fsum[i];
        }
      indsum[0]=0.0;
      for(i=0;i < *n;i++)
        if(i< maxind+1)
          indsum[i+1] = 0;
        else
          indsum[i+1] = 1;
      for(m=0;m <= *method;m++)
        for(i=1;i <= *n;i++)
          indsum[i] = indsum[i-1]+indsum[i];
      k = maxind;
      maxf = 1e+38;
      minf = -1e+38;
      lastind = maxind + 1;
      }
    else
    if(newmaxf < minf)
      {
      lastsign = -1;
      for(i=lastind;i<=minind;i++)
        f[i] = minf;
      if(*method == 0)
        for(i=lastind;i< *n;i++)
          fsum[i+1] = fsum[i] + f[i];
      else
        {
        fsum[0]=0.0;
        for(i=0;i < *n;i++)
          fsum[i+1] = f[i];
        for(m=0;m <= *method;m++)
          for(i=1;i <= *n;i++)
            fsum[i] = fsum[i-1]+fsum[i];
        }
      indsum[0]=0.0;
      for(i=0;i < *n;i++)
        if(i< minind+1)
          indsum[i+1] = 0;
        else
          indsum[i+1] = 1;
      for(m=0;m <= *method;m++)
        for(i=1;i <= *n;i++)
          indsum[i] = indsum[i-1]+indsum[i];
      k = minind;
      maxf = 1e+38;
      minf = -1e+38;
      lastind = minind + 1;
      }
    else
      {
      if(newmaxf < maxf)
        {
        maxf = newmaxf;
        maxind = k;
        }
      if(newminf > minf)
        {
        minf = newminf;
        minind = k;
        }
      }
    k++;
    }

  if(lastsign == 1)
    for(i=lastind;i<*n;i++)
      f[i]=minf;
  else
    for(i=lastind;i<*n;i++)
      f[i]=maxf;

  for(i=0;i<*n;i++)
    y[i] = f[i];
  
  free(indsum);
  free(ysum);
  free(fsum);
  free(f);
  }

void mrcheck(double *intres, int *n, double *thresh, int *jact, int *kact, int *signact, int *nact)
  {
  int j,k;

  *nact=0;
  for(k=0;k<*n;k++)
    for(j=0;j<=k;j++)
      if(fabs(fabs(intres[k+1]-intres[j])-sqrt(k-j+1.0)*(*thresh))<1e-08)
        {
        if(*nact == *n)
          {
          *nact = -1;
          return;
          }
        jact[*nact] = j;
        kact[*nact] = k;
        if(intres[k+1]-intres[j] > 0)
          signact[*nact] = -1;
        else
          signact[*nact] = 1;
        (*nact)++;
        }
  }

void calcepsnewmr(double *intdirec, double *intres, int *n, double *thresh,
int *jact, int *kact, int *signact, int nact,
double *eps, int *newmrj, int *newmrk,int *DYADIC)
  {
  int i,j,k,length;
  double direcsum,curreps;

  *eps = 1e+38;
  length=1;
  while(length <= *n)
    {
    j=0;
    while(j < *n)
      {
      k=j+length-1;
      if(k>=*n)
        k = *n-1;
      direcsum = intdirec[k+1]-intdirec[j];
      if(direcsum > 1e-13)
        curreps = (intres[k+1]-intres[j]+sqrt(k-j+1.0)*(*thresh))/direcsum;
      else
      if(direcsum < -1e-13)
        curreps = (intres[k+1]-intres[j]-sqrt(k-j+1.0)*(*thresh))/direcsum;
      else
        curreps = 1e+38;
      if((curreps > 1e-13)&&(curreps < *eps))
        {
        for(i=0;i<nact;i++)
          if((signact[i] < 2)&&(jact[i] == j)&&(kact[i]==k))
            i=nact+10;
        if(i < nact + 2)
          {
          *eps = curreps;
          *newmrj = j;
          *newmrk = k;


          }
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
  }


void ludecomp (double *A, int *myindx, int *n, int maxband)
  {
  int i, imax, j, k,tmpind,*indx,l;
  double big, dum;
  
  indx=malloc(*n*sizeof(int));
  for(i=0;i<*n;i++)
    myindx[i]=i;
  for(j = 0; j < *n; j++)
    {
    big = 0.0;
    for(i = j; (i < *n)&&(i<=j+maxband); i++)
      if(fabs(A[i+j*(*n)]) >= big)
                {
                big = fabs(A[i+j*(*n)]);
                imax = i;
                }
    if(big == 0.0)
      return; 
    if(j != imax)
        for(k = 0; k < *n; k++)
            {
            dum = A[imax+k*(*n)];
            A[imax+k*(*n)] = A[j+k*(*n)];
            A[j+k*(*n)] = dum;
            }
    indx[j]=imax;    
    for(i = j + 1; (i < *n)&&(i<=j+maxband); i++)
      {
      A[i+j*(*n)] = A[i+j*(*n)] / A[j+j*(*n)];
      for(l = j + 1; l < *n; l++)
        A[i+l*(*n)] -= A[i+j*(*n)] * A[j+l*(*n)];
        
      }
    }
    
  for(i=*n-1;i>=0;i--)
    {
    tmpind=myindx[i];
    myindx[i]=myindx[indx[i]];
    myindx[indx[i]]=tmpind;
    }
  free(indx);
  }

void linsolve (double *A, double *b, int *myindx, int *n, int TRANSPOSE)
  {
  int i, j;
  double sum;
  double *bnew;

  bnew = malloc((*n)*sizeof(double));

  if(!TRANSPOSE)
    {
    for(i = 0; i < *n; i++)
      bnew[myindx[i]] = b[i];
    for(i = 0; i < *n; i++){
        sum = bnew[i];
        for(j = 0; j <= i - 1; j++)
          sum -= A[i+j*(*n)] * b[j];
        b[i] = sum;
      }
    for(i = *n - 1; i >= 0; i--){
        sum = b[i];
        for(j = i + 1; j < *n; j++){
            sum -= A[i+j*(*n)] * b[j];
        }
        b[i] = sum / A[i+i*(*n)];
      }
    }
  else
    {
    for(i = 0; i < *n; i++)
      {
        bnew[i] = b[i];
        for(j = 0; j <= i - 1; j++)
          bnew[i] -= A[j+i*(*n)] * bnew[j];
        bnew[i] = bnew[i] / A[i+i*(*n)];
      }
    for(i = *n - 1; i >= 0; i--)
        for(j = i + 1; j < *n; j++)
          bnew[i] -= A[j+i*(*n)] * bnew[j];
      
    for(i = 0; i < *n; i++)
      b[i] = bnew[myindx[i]];
    }
  free(bnew);
  }



void mintvmon(int *n, double *f, int *derivsign, int *secsign, int *method, int *DYADIC, double *thresh, double *y, int *MONCONS, int *CONVCONS, int *jact, int *kact, int *signact,int *outnact,int *piecesleft, int *piecesright)
  {
  int i,j,currI,currsign,m,DROPMR,colind,currpos,poseqjump,null=0,eins=1;
  int newmrj,newmrk,lpl,k,jumpsign,currind,maxband,Bsize,NULLJUMP,NULLNEWMON,NULLNEWCONV;
  int nact,*myindx,newmrpos,*actmonconsts,nmonconsts,posnewmonconst,*newpiecesleft,*newpiecesright;
  int *ISACTIVE,*ISACTIVECONV,newsign,zaehler,*actconvconsts,nconvconsts,posnewconvconst;
  double epseqjump,epsnewmonconst,minmrderiv,minjumpderiv,*direc,*realf,*feps,*realeps,*fderiv,*epsderiv,epsnewconvconst;
  double tmpeps,epsnewmr,currtv,*intres,*inty,*Blu,*v,*pi,*knotsdirec,*jumpderiv,*zufall;
  long int *B,*Bnew,*vint;
  char dummy[80];

  ISACTIVE   = malloc(*n * sizeof(int));
  ISACTIVECONV = malloc(*n * sizeof(int));
  myindx     = malloc(*n * sizeof(int));
  newpiecesleft     = malloc(*n * sizeof(int));
  newpiecesright    = malloc(*n * sizeof(int));
  actmonconsts= malloc(*n * sizeof(int));
  actconvconsts= malloc(*n * sizeof(int));
  realf      = malloc(*n * sizeof(double));
  fderiv     = malloc(*n * sizeof(double));
  epsderiv   = malloc(*n * sizeof(double));
  realeps    = malloc(*n * sizeof(double));
  vint       = malloc((*n+1) * sizeof(long int));
  intres     = malloc((*n+1) * sizeof(double));
  inty       = malloc((*n+1) * sizeof(double));
  pi         = malloc(*n * sizeof(double));
  feps       = malloc(*n * sizeof(double));
  knotsdirec = malloc(*n * sizeof(double));
  zufall     = malloc(*n * sizeof(double));
  v          = malloc((*n+1) * sizeof(double));
  direc      = malloc((*n+1) * sizeof(double));
  jumpderiv  = malloc((*n+1) * sizeof(double));

  currsign = 1;
  for(i=0;i<*n;i++)
    {
    zufall[i]=currsign*sin((double)i);
    currsign = currsign*1;
    }

  nmonconsts = 0;
  nconvconsts = 0;

  for(i=0;i<*n;i++)
    f[i] = y[i];
  for(i=0;i<*n;i++)
    ISACTIVECONV[i] = 0;
  for(i=0;i<*n;i++)
    ISACTIVE[i] = 0;



  if(!(*MONCONS)&&(!(*CONVCONS)))
    yams(f, n, thresh, method);
  else
  if((*MONCONS)&&(!(*CONVCONS)))
    {
    mintvmon(n, f, derivsign, secsign, &null, DYADIC, thresh, y, &null, &null, jact, kact, signact, outnact,piecesleft,piecesright);

    lpl=1;
    piecesleft[0] = 1;
    for(i=1;i<*n;i++)
      if(fabs(f[i]-f[i-1])>1e-13)
        {
        piecesleft[lpl] = i+1;
        piecesright[lpl-1] = i;
        lpl++;
        }
    piecesright[lpl-1]= (*n);

    /* Calculate monotonicity constraints */
 
    if(lpl == 1)
      {
      }
    else
      {
      if(f[piecesleft[1]-1] - f[piecesleft[0]-1] > 0)
        currsign = 1;
      else
        currsign = -1;

      k=0;
      i=0;
      while(i<lpl-1)
        {
        while((i<lpl-1)&&(currsign * (f[piecesleft[i+1]-1] - f[piecesleft[i]-1]) > 0))
          i++;
        if(i == lpl-1)
          while(k<*n)
            {
            derivsign[k] = currsign;
            k++;
            }
        else
          while(k<0.5*(piecesleft[i]+piecesright[i])-1)
            {
            derivsign[k] = currsign;
            k++;
            }
        currsign = -currsign;
        }
      }

    for(i=0;i<lpl;i++)
      for(m=0;m<*method;m++)
        if(piecesright[i]-piecesleft[i]>m)
          {
          ISACTIVE[piecesleft[i]-1+m] = 1;
          actmonconsts[nmonconsts] = piecesleft[i]-1+m;
          nmonconsts++;
          }


    for(m=1;m<=*method;m++)
      for(i=*n-1;i>0;i--)
        f[i] = f[i]-f[i-1];
    }
  else
  if(!(*MONCONS)&&(*CONVCONS))
    { 
    mintvmon(n, f, derivsign, secsign, &eins, DYADIC, thresh, y, &null, &null, jact, kact, signact, outnact,piecesleft,piecesright);
    for(i=*n-1;i>0;i--)
      f[i] = f[i]-f[i-1];

    lpl=1;
    piecesleft[0] = 1;
    for(i=1;i<*n;i++)
      if(fabs(f[i]-f[i-1])>1e-13)
        {
        piecesleft[lpl] = i+1;
        piecesright[lpl-1] = i;
        lpl++;
        }
    piecesright[lpl-1]= (*n);

    /* Calculate convexity constraints */
 
    if(lpl == 1)
      {
      }
    else
      {
      if(f[piecesleft[1]-1] - f[piecesleft[0]-1] > 0)
        currsign = 1;
      else
        currsign = -1;

      k=0;
      i=0;
      while(i<lpl-1)
        {
        while((i<lpl-1)&&(currsign * (f[piecesleft[i+1]-1] - f[piecesleft[i]-1]) > 0))
          i++;
        if(i == lpl-1)
          while(k<*n)
            {
            secsign[k] = currsign;
            k++;
            }
        else
          while(k<0.5*(piecesleft[i]+piecesright[i])-1)
            {
            secsign[k] = currsign;
            k++;
            }
        currsign = -currsign;
        }
      }

    for(i=0;i<lpl;i++)
      for(m=0;m<*method-1;m++)
        if(piecesright[i]-piecesleft[i]>m)
          {
          ISACTIVECONV[piecesleft[i]-1+m] = 1;
          actconvconsts[nconvconsts] = piecesleft[i]-1+m;
          nconvconsts++;
          }

    for(m=1;m<=*method-1;m++)
      for(i=*n-1;i>0;i--)
        f[i] = f[i]-f[i-1];
    }
  else
    { 
    mintvmon(n, f, derivsign, secsign, &eins, DYADIC, thresh, y, &eins, &null, jact, kact, signact, outnact,piecesleft,piecesright);
    for(i=*n-1;i>0;i--)
      f[i] = f[i]-f[i-1];

    /* Calculate convexity constraints */
 
    /*if(lpl == 1)*/
    if(*outnact == 1)
      {
      }
    else
      {
      k=0;
      while(fabs(f[piecesleft[k+1]-1] - f[piecesleft[k]-1])<1e-13) k++;

      if(f[piecesleft[k]-1] - f[piecesleft[k]-1] > 1e-13)
        currsign = 1;
      else
        currsign = -1;

      lpl = *outnact;
      k=0;
      i=0;
      while(i<lpl-1)
        {
        while((i<lpl-1)&&(currsign * (f[piecesleft[i+1]-1] - f[piecesleft[i]-1]) >= -1e-13))
          i++;
        if(i == lpl-1)
          while(k<*n)
            {
            secsign[k] = currsign;
            k++;
            }
        else
          while(k<0.5*(piecesleft[i]+piecesright[i])-1)
            {
            secsign[k] = currsign;
            k++;
            }
        currsign = -currsign;
        }
      }


    lpl=0;
    for(i=0;i<*outnact;i++)
      {
      newpiecesleft[lpl]=piecesleft[i];
      newpiecesright[lpl]=piecesright[i];
      lpl++;

      for(m=0;m<*method-1;m++)
        if(piecesright[i]-piecesleft[i]>m)
          {
          ISACTIVECONV[piecesleft[i]-1+m] = 1;
          actconvconsts[nconvconsts] = piecesleft[i]-1+m;
          nconvconsts++;
          newpiecesright[lpl]=newpiecesright[lpl-1];
          newpiecesleft[lpl]=newpiecesleft[lpl-1]+1;
          newpiecesright[lpl-1]=newpiecesleft[lpl-1];
          lpl++;
          }
      }
    for(i=0;i<lpl;i++)
      {
      piecesleft[i]=newpiecesleft[i]; 
      piecesright[i]=newpiecesright[i]; 
      }
    for(m=1;m<=*method-1;m++)
      for(i=*n-1;i>0;i--)
        f[i] = f[i]-f[i-1];
    }

  if((!(*MONCONS))||(!(*CONVCONS)))
    {
    lpl=1;
    piecesleft[0] = 1;
    for(i=1;i<*n;i++)
      if(fabs(f[i]-f[i-1])>1e-13)
        {
        piecesleft[lpl] = i+1;
        piecesright[lpl-1] = i;
        lpl++;
        }
    piecesright[lpl-1]= (*n);
    }

  intres[0]=0;
  inty[0]=0;
  for(i=1;i<=*n;i++)
    {
    intres[i]=f[i-1];
    inty[i]=inty[i-1]+y[i-1];
    }
  for(m=1;m<=*method;m++)
    for(i=1;i<=*n;i++)
      intres[i] = intres[i-1]+intres[i];
  for(i=1;i<=*n;i++)
    intres[i] = y[i-1]-intres[i];
  for(i=1;i<=*n;i++)
    intres[i] = intres[i-1]+intres[i];


  if(!(*MONCONS)&&(!(*CONVCONS)))
    mrcheck(intres, n, thresh, jact, kact, signact, &nact);
  else
    nact = *outnact;



nact += nmonconsts;
nact += nconvconsts;

if(nact != lpl)
  {
  return;
  }

/* Combine multiresolution and monotonicity constraints */

k = nmonconsts;
for(j = nact-nconvconsts-1 ; j >= 0 ; j--)
  {
  if((k == 0 ) || ((j >= k) && ((kact[j-k] > actmonconsts[k-1]+1)||((kact[j-k] == actmonconsts[k-1]+1)&&(jact[j-k] > actmonconsts[k-1])))))
    {
    jact[j] = jact[j-k];
    kact[j] = kact[j-k];
    signact[j] = signact[j-k];
    }
  else
    {
    jact[j] = actmonconsts[k-1];
    kact[j] = actmonconsts[k-1]+1;
    signact[j] = 2;
    k--;
    }
  }

/* Combine with convexity constraints */

k = nconvconsts;
for(j = nact-1 ; j >= 0 ; j--)
  {
  if((k == 0 ) || ((j >= k) && ((kact[j-k] > actconvconsts[k-1]+1)||((kact[j-k] == actconvconsts[k-1]+1)&&(jact[j-k] > actconvconsts[k-1])))))
    {
    jact[j] = jact[j-k];
    kact[j] = kact[j-k];
    signact[j] = signact[j-k];
    }
  else
    {
    jact[j] = actconvconsts[k-1]-1;
    kact[j] = actconvconsts[k-1]+1;
    signact[j] = 3;
    k--;
    }
  }

/* Calculate B */

Bsize=nact+10;
Blu     = malloc((Bsize)*(Bsize) * sizeof(double));
B     = malloc((Bsize)*(Bsize) * sizeof(long int));

for(j=0;j<nact;j++)
  {
  for(i=0;i<*n;i++)
    if((i<jact[j])||(i>kact[j]))
      vint[i] = 0;
    else
      if(signact[j] < 2)
        vint[i] = signact[j];
      else
      if(signact[j] == 2)
        {
        if(i == jact[j])
          vint[i] = derivsign[jact[j]];
        else
          vint[i] = -derivsign[jact[j]];
        }
      else
        {
        if(i == jact[j])
          vint[i] = -secsign[jact[j]+1];
        else
        if(i == jact[j]+1)
          vint[i] = 2*secsign[jact[j]+1];
        else
        if(i == jact[j]+2)
          vint[i] = -secsign[jact[j]+1];
        }
    for(m=1;m<=*method+1;m++)
      for(i=*n-2;i>=0;i--)
        vint[i] = vint[i]+vint[i+1];
  for(i=0;i<nact;i++)
    B[i+j*nact] = vint[piecesleft[i]-1];

  }
zaehler = 0;


while(1<2)
  {
  zaehler++;
  /* Calculate LU decomposition of B */
  maxband=0;
  for(j=0;j<nact;j++)
    for(i=0;i<nact;i++)
      {
      if(B[i+j*nact]!=0)
        if(i-j>maxband)
          maxband=i-j;
      Blu[i+j*nact] = B[i+j*nact];
      }

  ludecomp(Blu,myindx,&nact,maxband);

  /* calculate f, realf and intres */
if(zaehler!=1)
  {
  for(i=0;i<nact;i++)
    if(signact[i] < 2)
      v[i] = signact[i]*(inty[kact[i]+1]-inty[jact[i]])+sqrt(kact[i]-jact[i]+1.0)*(*thresh);
    else
      v[i] = 0;
  linsolve(Blu,v,myindx,&nact,1);
  for(i=1;i<nact;i++)
    v[i] = v[i]+v[i-1];

  for(i=0;i<nact;i++)
    for(j=piecesleft[i];j<=piecesright[i];j++)
      f[j-1] = v[i]; 
  }
  for(i=0;i<nact;i++)
    if(signact[i] < 2)
      v[i] = 0;
    else
    if(signact[i] == 2)
      v[i] = -derivsign[jact[i]]*zufall[jact[i]];
    else
      v[i] = -secsign[jact[i]+1]*zufall[jact[i]+1];
  linsolve(Blu,v,myindx,&nact,1);

  for(i=1;i<nact;i++)
    v[i] = v[i]+v[i-1];

  for(i=0;i<nact;i++)
    for(j=piecesleft[i];j<=piecesright[i];j++)
      feps[j-1] = v[i]; 


  pi[0] = 0;

  for(i=1;i<nact;i++)
    if(f[piecesleft[i]-1]>f[piecesleft[i-1]-1]+1e-08)
      pi[i] = 1;
    else
    if(f[piecesleft[i]-1]<f[piecesleft[i-1]-1]-1e-08)
      pi[i] = -1;
    else
    if(feps[piecesleft[i]-1]>feps[piecesleft[i-1]-1])
      pi[i] = 1;
    else
      pi[i] = -1;

  for(i=0;i<*n;i++)
    {
    realf[i] = f[i];
    realeps[i] = feps[i];
    fderiv[i] = f[i];
    epsderiv[i] = feps[i];
    }
  for(m=1;m<=*method;m++)
    for(i=1;i<*n;i++)
      {
      realf[i] = realf[i]+realf[i-1];
      realeps[i] = realeps[i]+realeps[i-1];
      }
  for(m=1;m<=*method-1;m++)
    for(i=1;i<*n;i++)
      {
      fderiv[i] = fderiv[i]+fderiv[i-1];
      epsderiv[i] = epsderiv[i]+epsderiv[i-1];
      }
  for(i=1;i<=*n;i++)
    intres[i] = intres[i-1]+y[i-1]-realf[i-1];

  /* Calc total variation of f */
  currtv=0;
  for(i=1;i<*n;i++)
    currtv += fabs(f[i]-f[i-1]);
/*
  printf("zaehler: %d -- TV is %f -- nact is %d -- maxband is %d\n",zaehler,currtv,nact,maxband);





if(*MONCONS && *CONVCONS)
  {
  for(i=0;i<*n-1;i++)
    if(derivsign[i]*fderiv[i+1]<-1e-07)
      {
      printf("MON %d hurt %3.15f\n",i,fderiv[i+1]);
      gets(dummy);
      }

  for(i=0;i<*n-1;i++)
    if(secsign[i]*f[i+1]<-1e-07)
      {

      printf("COINV %d hurt %3.15f\n",i,f[i+1]);
      gets(dummy);
      }
  }
*/  



  linsolve(Blu,pi,myindx,&nact,0);

  /* minmrderiv is the minimum derivative for leaving out an active multiresolution
     or monotonicity constraint, minjumpderiv is the minimum derivative for
     introducing a jump */

  minmrderiv = 1e+38;
  minjumpderiv = 1e+38;

  for(i=0;i<nact;i++)
    if(-pi[i] < minmrderiv)
      {
      minmrderiv = -pi[i];
      currI = i;
      }


  if(minmrderiv < -1e-07)
    {
    for(i=0;i<nact;i++)
      knotsdirec[i] = 0;
    knotsdirec[currI] = -1;
    
    DROPMR = 1;
    }
  else
    {
    for(i=0;i<*n;i++)
      jumpderiv[i] = 0.0;

    for(j=0;j<nact;j++)
      if(signact[j] < 2)
        {
        if(jact[j]-1>=0)
          jumpderiv[jact[j]-1]-=signact[j]*pi[j];
        jumpderiv[kact[j]]+=signact[j]*pi[j];
        }

    for(m=1;m<=2;m++) 
      for(i=*n-2;i>=0;i--)
        jumpderiv[i] = jumpderiv[i]+jumpderiv[i+1];

    for(j=0;j<nact;j++)
      if(signact[j] == 2)
        jumpderiv[kact[j]] -= derivsign[jact[j]] * pi[j];

    if(*method!=0)
      for(i=*n-2;i>=0;i--)
        jumpderiv[i] = jumpderiv[i]+jumpderiv[i+1];

    for(j=0;j<nact;j++)
      if(signact[j] == 3)
        jumpderiv[jact[j]+2] -= secsign[jact[j]+1] * pi[j];

    for(m=1;m<=*method-1;m++)
      for(i=*n-2;i>=0;i--)
        jumpderiv[i] = jumpderiv[i]+jumpderiv[i+1];
  
    /* Look at jumps on non-constant pieces first */
    for(i=0;i<nact;i++)
      for(j=piecesleft[i];j<piecesright[i];j++)
        if(fabs(realf[j]-realf[j-1])>1e-08)
          {
          if(jumpderiv[j]>0.0)
            {
            jumpderiv[j] *= -1.0;
            jumpsign = -1;
            }
          else
            jumpsign = 1;

          if(j > 0)
            jumpderiv[j] += 1.0;
          if(jumpderiv[j] < minjumpderiv)
            {
            minjumpderiv = jumpderiv[j];
            colind = i + 1;
            currpos = j;
            currsign = jumpsign;
            }
          }

    /* If still no direction found, look at jumps on constant pieces */
    if(minjumpderiv > -1e-07)
    for(i=0;i<nact;i++)
      for(j=piecesleft[i];j<piecesright[i];j++)
        if(fabs(realf[j]-realf[j-1])<=1e-08)
          {
          if(jumpderiv[j]>0.0)
            {
            jumpderiv[j] *= -1.0;
            jumpsign = -1;
            }
          else
            jumpsign = 1;

          if(j > 0)
            jumpderiv[j] += 1.0;
          if(jumpderiv[j] < minjumpderiv)
            {
            minjumpderiv = jumpderiv[j];
            colind = i + 1;
            currpos = j;
            currsign = jumpsign;
            }
          }

    if(minjumpderiv < -1e-07)
      {
      for(i=0;i<=*n;i++)
        v[i] = 0;
      for(i=piecesleft[colind-1] ; i<=currpos; i++)
        v[i] = 1;
      for(m=1;m<=*method;m++)
        for(i=1;i<=*n;i++)
          v[i] = v[i-1]+v[i];
      for(i=0;i<nact;i++)
        if(signact[i] == 2)
          knotsdirec[i] = -currsign * derivsign[jact[i]]*(v[jact[i]+1]-v[kact[i]+1]);
        else
        if(signact[i] == 3)
          knotsdirec[i] = currsign * secsign[jact[i]+1]*(v[jact[i]+1]-2*v[jact[i]+2]+v[jact[i]+3]);

      for(i=1;i<=*n;i++)
        v[i] = v[i-1]+v[i];

      for(i=0;i<nact;i++)
        if(signact[i] < 2)
          knotsdirec[i] = -currsign * signact[i]*(v[kact[i]+1]-v[jact[i]]);

      DROPMR = 0;
      }
    else /* no downward direction, we are finished */
      break;
    }

  /* Calculate direction */
  
  linsolve(Blu,knotsdirec,myindx,&nact,1);
  for(i=1;i<nact;i++)
    knotsdirec[i] = knotsdirec[i]+knotsdirec[i-1];
  direc[0] = 0;
  for(i=0;i<nact;i++)
    for(j=piecesleft[i];j<=piecesright[i];j++)
      if(!DROPMR && (i == colind -1) && (j<=currpos))
        direc[j] = knotsdirec[i]+currsign;
      else 
        direc[j] = knotsdirec[i]; 

  /* Calculate next tie in f */

  epseqjump = 1e+38;
  NULLJUMP = 0;

  for(i=0; i< nact-1; i++)
    {
    tmpeps = -(f[piecesleft[i+1]-1] - f[piecesright[i]-1])/(direc[piecesleft[i+1]] - direc[piecesright[i]]);
    if((tmpeps>0)&&(fabs(f[piecesleft[i+1]-1] - f[piecesright[i]-1])>1e-08)&&(!NULLJUMP)&&(tmpeps < epseqjump))
      {
      epseqjump = tmpeps;
      poseqjump = i+1;
      }
    else
    if((fabs(f[piecesleft[i+1]-1] - f[piecesright[i]-1])<=1e-08)&&(fabs(direc[piecesleft[i+1]] - direc[piecesright[i]])> 1e-05))
      {
      tmpeps = -(feps[piecesleft[i+1]-1] - feps[piecesright[i]-1])/(direc[piecesleft[i+1]] - direc[piecesright[i]]);
      if((tmpeps > 0)&&((tmpeps < epseqjump)||(!NULLJUMP)))
        {
        NULLJUMP = 1;
        epseqjump = tmpeps;
        poseqjump = i+1;
        }
      }
    }


  for(m=1;m<=*method-1;m++)
    for(i=1;i<=*n;i++)
      direc[i] = direc[i-1]+direc[i];
  
  epsnewconvconst = 1e+38;
  NULLNEWCONV = 0;

  if(*CONVCONS)
    {
    /* Calculate next tie in fderiv (active convexity condition) */

    for(i=1;i< *n-1;i++)
      if(!ISACTIVECONV[i])
        {
        tmpeps = -(fderiv[i+1] - fderiv[i])/(direc[i+2] - direc[i+1]);
        if((tmpeps>0)&&(fabs(fderiv[i+1]-fderiv[i])>1e-08)&&(!NULLNEWCONV)&&(tmpeps<epsnewconvconst))
          {
          epsnewconvconst = tmpeps;
          posnewconvconst = i;
          }
        else if((fabs(fderiv[i+1]-fderiv[i])<=1e-08)&&((direc[i+2]-direc[i+1]) * secsign[i]<-1e-05))
          { 
          tmpeps = -(epsderiv[i+1] - epsderiv[i]-zufall[i])/(direc[i+2] - direc[i+1]);
          if(((tmpeps < epsnewconvconst)||(!NULLNEWCONV)))
            {
            NULLNEWCONV = 1;
            epsnewconvconst = tmpeps;
            posnewconvconst = i;
            }
          }
        }
    }

if(*method!=0)
  for(i=1;i<=*n;i++)
    direc[i] = direc[i-1]+direc[i];

  epsnewmonconst = 1e+38;
  NULLNEWMON = 0;

  if(*MONCONS)
    {
    /* Calculate next tie in realf (active monotonicity condition) */

    for(i=0;i< *n-1;i++)
      if(!ISACTIVE[i])
        {
        tmpeps = -(realf[i+1] - realf[i])/(direc[i+2] - direc[i+1]);
        if((tmpeps>0)&&(fabs(realf[i+1]-realf[i])>1e-08)&&(!NULLNEWMON)&&(tmpeps<epsnewmonconst))
          {
          epsnewmonconst = tmpeps;
          posnewmonconst = i;
          }
        else if((fabs(realf[i+1]-realf[i])<=1e-08)&&((direc[i+2]-direc[i+1]) * derivsign[i]<-1e-05))
          { 
          tmpeps = -(realeps[i+1] - realeps[i]-zufall[i])/(direc[i+2] - direc[i+1]);
          if(((tmpeps < epsnewmonconst)||(!NULLNEWMON)))
            {
            NULLNEWMON = 1;
            epsnewmonconst = tmpeps;
            posnewmonconst = i;
            }
          }
        }
    }


  for(i=1;i<=*n;i++)
    direc[i] = direc[i-1]+direc[i];

  /* Calculate next active multiresolution condition */

  calcepsnewmr(direc, intres, n, thresh, jact, kact, signact, nact, &epsnewmr, &newmrj, &newmrk, DYADIC);

if(NULLNEWCONV | NULLJUMP | NULLNEWMON)
  {
  epsnewmr = 1e+38;
  if(!NULLNEWCONV)
    epsnewconvconst = 1e+38;
  if(!NULLNEWMON)
    epsnewmonconst = 1e+38;
  if(!NULLJUMP)
    epseqjump = 1e+38;
  }

  if((epseqjump < epsnewmr)&&(epseqjump < epsnewmonconst)&&(epseqjump < epsnewconvconst)) /* leave out some jump */
    {
    for(i=poseqjump-1;i<nact-1;i++)
      piecesright[i] =piecesright[i+1];
    for(i=poseqjump;i<nact-1;i++)
      piecesleft[i] =piecesleft[i+1];
    lpl = nact - 1;

    for(j=0;j<nact;j++)
      {
      for(i=0;i<poseqjump;i++)
        B[i+j*lpl] = B[i+j*nact]; 
      for(i=poseqjump;i<lpl;i++)
        B[i+j*lpl] = B[i+1+j*nact] ;
      }

    }
  else
    lpl = nact;

  if(!DROPMR) /* add some jump */
    {
    for(currind=0;piecesright[currind]<=currpos;currind++) 
      ;
    for(j=lpl;j>=currind+1;j--)
      piecesright[j]=piecesright[j-1];
    for(j=lpl;j>currind+1;j--)
      piecesleft[j]=piecesleft[j-1];
    piecesright[currind] = currpos;
    piecesleft[currind+1] = currpos+1;

    for(j=nact;j>=0;j--)
      {
      for(i=lpl;i>currind;i--)
        B[i+j*(lpl+1)] = B[i-1+j*lpl]; 
      for(i=currind;i>=0;i--)
        B[i+j*(lpl+1)] = B[i+j*lpl]; 
      }
    for(i=0;i<=*n;i++)
      vint[i] = 0;
    vint[currpos+1] = 1;
    for(m=1;m<=*method+1;m++)
      for(i=1;i<=*n;i++)
        vint[i] = vint[i]+vint[i-1];

    for(j=0;j<nact;j++)
      if(signact[j] == 2)
        B[currind+1+j*(lpl+1)] = derivsign[jact[j]]*(vint[jact[j]+1] - vint[kact[j]+1]);     
      else
      if(signact[j] == 3)
        B[currind+1+j*(lpl+1)] = -secsign[jact[j]+1]*(vint[jact[j]+1] - 2*vint[jact[j]+2] + vint[jact[j]+3]);     
 
    for(i=1;i<=*n;i++)
      vint[i] = vint[i]+vint[i-1];
    for(j=0;j<nact;j++)
      if(signact[j] < 2)
        B[currind+1+j*(lpl+1)] = signact[j]*(vint[kact[j]+1]-vint[jact[j]]);


    lpl++;
    }

  if(DROPMR) /* Drop some active condition */
    {
    if(signact[currI] == 2)
      {
      ISACTIVE[jact[currI]]=0;
      nmonconsts--;
      }
    else
    if(signact[currI] == 3)
      {
      ISACTIVECONV[jact[currI]+1]=0;
      nconvconsts--;
      }
    for(j=currI;j<nact-1;j++)
      {
      jact[j]=jact[j+1];
      kact[j]=kact[j+1];
      signact[j]=signact[j+1];
      }
    nact--;

    for(j=currI;j<nact;j++)
      for(i=0;i<lpl;i++)
        B[i+j*lpl] = B[i+(j+1)*lpl] ;

    }
  
  if((epseqjump >= epsnewmr) || (epseqjump >= epsnewmonconst) || (epseqjump >= epsnewconvconst)) /* new active condition */
    {
    for(i=0;i<=*n;i++)
      vint[i] = 0;
    if((epsnewmonconst > epsnewmr)&&(epsnewconvconst > epsnewmr))
      {
      if(direc[newmrk+1]-direc[newmrj] > 0)
        newsign = 1;
      else
        newsign = -1;
      for(i=newmrj+1 ; i<=newmrk+1; i++)
        vint[i] = newsign; 
      }
    else
    if(epsnewconvconst > epsnewmonconst)
      {
      ISACTIVE[posnewmonconst] = 1;
      newmrj = posnewmonconst;
      newmrk = posnewmonconst+1;
      newsign = 2;
      vint[posnewmonconst+1] = derivsign[newmrj];
      vint[posnewmonconst+2] = -derivsign[newmrj];
      nmonconsts++;
      }
    else
      {
      ISACTIVECONV[posnewconvconst] = 1;
      newmrj = posnewconvconst-1;
      newmrk = posnewconvconst+1;
      newsign = 3;
      if(posnewconvconst>0)
        vint[posnewconvconst] = -secsign[newmrj+1];
      vint[posnewconvconst+1] = 2*secsign[newmrj+1];
      vint[posnewconvconst+2] = -secsign[newmrj+1];
      nconvconsts++;
      }
    for(newmrpos=0;(newmrpos<nact)&&(kact[newmrpos]<=newmrk)&&((kact[newmrpos]+1<newmrk)||(jact[newmrpos]<=newmrj));newmrpos++) 
      ;
    for(j=nact;j>newmrpos;j--)
      {
      jact[j]=jact[j-1];
      kact[j]=kact[j-1];
      signact[j]=signact[j-1];
      }
    jact[newmrpos] = newmrj;
    kact[newmrpos] = newmrk;
    signact[newmrpos] = newsign;
    for(m=1;m<=*method+1;m++)
      for(i=*n-1;i>=1;i--)
        vint[i] = vint[i]+vint[i+1];
    for(j=nact;j>newmrpos;j--)
      for(i=(lpl-1);i>=0;i--)
        B[i+j*lpl] = B[i+(j-1)*lpl]; 
    for(i=0;i<lpl;i++)
      B[i+newmrpos*lpl] = vint[piecesleft[i]];
    nact++;
    }

  /* if B uses all the allocated memory, increase memory */

  if(nact == Bsize)
    {
    Bsize += 10;
    if(Bsize > *n)
      Bsize = (*n);

    Bnew  = malloc((Bsize)*(Bsize) * sizeof(long int));
    for(i=0;i<nact;i++)
      for(j=0;j<nact;j++)
        Bnew[i+j*nact] = B[i+j*nact];
    free(Blu);
    free(B);
    Blu = malloc((Bsize)*(Bsize) * sizeof(double));
    B = Bnew;
    }
  }


  for(i=0;i<*n;i++)
    f[i] = realf[i];

  *outnact = nact;

  free(pi);
  free(v);
  free(vint);
  free(knotsdirec);
  free(direc);
  free(myindx);
  free(jumpderiv);
  free(intres);
  free(Blu);
  free(B);
  free(inty);
  free(realf);
  free(feps);
  free(actmonconsts);
  free(actconvconsts);
  free(ISACTIVE);
  free(ISACTIVECONV);
  free(realeps);
  free(zufall);
  free(fderiv);
  free(epsderiv);
  free(newpiecesleft);
  free(newpiecesright);
  }


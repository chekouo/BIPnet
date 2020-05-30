#include <stdbool.h>
#include <gsl/gsl_sf_psi.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
 #include <gsl/gsl_statistics.h>
 #include <gsl/gsl_math.h>
#include "header.h"

 double mean(int n,double * x){
	    int i;
	       double me=0;
	          for (i=0;i<n;i++)
			       me+=x[i];
		     return me/n;
		      }
_Bool compar(int n,_Bool *u,_Bool *v){
	   int i;
	      for (i=0;i<n;i++){
		           if (u[i]!=v[i])
				          return 0;
			      }
	         return 1;
		  }

 _Bool **  UniqueModel(int nbrsample, int p, _Bool ** rhosample,int * modelidx,int *countd1){

	 int i;int j=0;
	    int countd=0;
	       modelidx[countd]=0;
	          for (i=0;i<nbrsample;i++){
			       for (j=0;j<countd;j++){
				            _Bool resultX =compar(p,rhosample[i], rhosample[modelidx[j]]);
					           if (resultX==1)
							                break;
						        }
			            if (j==countd){
					           modelidx[countd]=i;
						          countd++;
							       }
				         }
		     *countd1=countd;


		        _Bool ** UniqModel=malloc(countd*sizeof(_Bool*));
			   for (i=0;i<countd;i++){
				        UniqModel[i]=malloc(p*sizeof(_Bool));
					     for (j=0;j<p;j++){
						            UniqModel[i][j]=rhosample[modelidx[i]][j];
							         }
					        }
			      printf("\n\n\nNumber of models is= %d",countd);
			         printf("\n");
				    return UniqModel;
 }


void save1d(char *filename,int n,double *data)
         {
         FILE *qfile;
         qfile=fopen(filename,"wb");
          if(qfile==NULL){
                printf("\nUnable to open file.\n");
                      exit(1);
                        } else {
        int i;
        for (i = 0; i < n; i++) {
        fprintf(qfile, "%f ", data[i]);
                }
                fclose(qfile);
                 }
}


 void save2d(char *filename,int n,int p,double ** data)
         {
         FILE *qfile;
         qfile=fopen(filename,"wb");
          if(qfile==NULL){
                printf("\nUnable to open file.\n");
                      exit(1);
                        } else {
        int i, j;
        for (i = 0; i < n; i++) {
              for (j = 0; j < p; j++) {
        fprintf(qfile, "%lf ", data[i][j]);
        }
        fprintf(qfile,"\n");
                }
                fclose(qfile);
                 }
}

int multivariate_gaussian (const gsl_rng * r,
                               const gsl_vector * mu,
                               const gsl_matrix * L,
                               gsl_vector * result)
{
/*
 * L     matrix resulting from the Cholesky decomposition of
 *  *      the inverse of   variance-covariance matrix Sigma^-1 = L L^T (dimension d x d)
*/
  const size_t M = L->size1;
  const size_t N = L->size2;

  if (M != N)
    {
      GSL_ERROR("requires square matrix", GSL_ENOTSQR);
    }
  else if (mu->size != M)
    {
      GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
    }
  else if (result->size != M)
    {
      GSL_ERROR("incompatible dimension of result vector", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < M; ++i)
        gsl_vector_set(result, i, gsl_ran_ugaussian(r));

gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, L, result);     
 //gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, L, result);
      gsl_vector_add(result, mu);

      return GSL_SUCCESS;
    }
}





int logmultigaussianT(const gsl_vector * x, const gsl_vector * y,
                                       const gsl_matrix * L,
                                       double * result,double *quadF,
                                       gsl_vector * work){
  const size_t M = L->size1;
  const size_t N = x->size;
   size_t i;
      double quadForm=0;        /* (x)' Sigma^{-1} (x) */
      double logSqrtDetSigma=0; /* log [ sqrt(|Sigma|) ] */

      /* compute: work = x - mu*/ 
      for (i = 0; i < M; ++i)
        {
          double xi = gsl_vector_get(y, i);
          //double mui = gsl_vector_get(mu, i);
          gsl_vector_set(work, i, xi);
        }


      /* compute: work = L^{-1} * (x - mu) */
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L, work);

      /* compute: quadForm = (x - mu)' Sigma^{-1} (x - mu) */
      gsl_blas_ddot(work, work, &quadForm);
      double x2=0;
      gsl_blas_ddot(x, x, &x2);
quadForm=x2-quadForm;
      /* compute: log [ sqrt(|Sigma|) ] = sum_i log L_{ii} */
      for (i = 0; i < M; ++i)
        {
          double Lii = gsl_matrix_get(L, i, i);
          logSqrtDetSigma += log(Lii);
        }

*quadF=quadForm;
      *result = -0.5*quadForm - logSqrtDetSigma - 0.5*N*log(2.0*M_PI);
 return GSL_SUCCESS;

}




int logmultigaussian(const gsl_vector * x,
                                       const gsl_vector * mu, 
                                       const gsl_matrix * L,
                                       double * result,double *quadF,
                                       gsl_vector * work){
  const size_t M = L->size1;
  const size_t N = L->size2;
  if (M != N)
    {
      GSL_ERROR("requires square matrix", GSL_ENOTSQR);
    }
  else if (mu->size != M)
    {
      GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
    }
  else if (x->size != M)
    {
      GSL_ERROR("incompatible dimension of quantile vector", GSL_EBADLEN);
    }
  else if (work->size != M)
    {
      GSL_ERROR("incompatible dimension of work vector", GSL_EBADLEN);
    }
  else
{
size_t i;
      double quadForm=0;        /* (x - mu)' Sigma^{-1} (x - mu) */
      double logSqrtDetSigma=0; /* log [ sqrt(|Sigma|) ] */

      /* compute: work = x - mu */
      for (i = 0; i < M; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double mui = gsl_vector_get(mu, i);
          gsl_vector_set(work, i, xi - mui);
        }

      /* compute: work = L^{-1} * (x - mu) */
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L, work);

      /* compute: quadForm = (x - mu)' Sigma^{-1} (x - mu) */
      gsl_blas_ddot(work, work, &quadForm);

      /* compute: log [ sqrt(|Sigma|) ] = sum_i log L_{ii} */
      for (i = 0; i < M; ++i)
        {
          double Lii = gsl_matrix_get(L, i, i);
          logSqrtDetSigma += log(Lii);
        }
*quadF=quadForm;
      *result = -0.5*quadForm - logSqrtDetSigma - 0.5*M*log(2.0*M_PI);
 return GSL_SUCCESS;
}
}

double inverseGaussian(gsl_rng * r, double mu, double lambda) {
double v=gsl_ran_gaussian (r, 1);  // sample from a normal distribution with a mean of 0 and 1 standard deviation
        double y = v*v;
        double x = mu + (mu*mu*y)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*y + mu*mu*y*y);
double test=gsl_ran_flat (r, 0, 1);    // sample from a uniform distribution between 0 and 1
        if (test <= (mu)/(mu + x))
               return x;
        else
               return (mu*mu)/x;
}
void readBoolVector(char *filename, int nRows,  _Bool * data )
{

   FILE *fp = fopen (filename, "r");
   if (fp==NULL)
   {
      printf("We can't open the file (%s).\n", filename);
      exit(1);
   }
   else
   { int iR;
      for (iR = 0; iR < nRows; ++iR )
      {
int x;
          fscanf(fp, "%d" , &x );
_Bool bb=(x!=0);
data[iR]=bb;
                      }

                              fclose(fp);
}
}


void sort(int n,double *x,int *idx)
{
int i,j;
double a;
int id;
for (i = 0; i < n; i++)
idx[i]=i;
for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            if (x[i] <= x[j])
            {
                a =  x[i];
             id=idx[i];
                idx[i]=idx[j];
                x[i] = x[j];
                idx[j]=id;
                x[j] = a;
            }
        }
    }

}



double auc(int n, double * esti,_Bool class[n]){
double fpr[n+2],tpr[n+2];
double auc1=0;
int P=0;//P=positive instances
int i,j;
double esti1[n];
for (i=0;i<n;i++){
esti1[i]=esti[i];
if (class[i]==1)
P+=1;
}
int idx[n];
sort(n,esti1,idx);
fpr[n+1]=1;tpr[n+1]=1;
fpr[0]=0;tpr[0]=0;
for (i=n;i>=1;--i){
double af=0;double at=0;
for (j=0;j<n;j++){
if (esti[j]>esti1[i-1]){
if (class[j]==0){
af+=1;
}
else {
at+=1;
} } }
tpr[i]=at/P;
fpr[i]=af/(n-P);
auc1+=(fpr[i+1]-fpr[i])*(tpr[i+1]+tpr[i]);
}
auc1+=(fpr[1]-fpr[0])*(tpr[1]+tpr[0]);
auc1=0.5*(auc1);
return auc1;
}




void readBoolArray(char *filename, int nRows, int nCols, _Bool ** data )
{

   FILE *fp = fopen (filename, "r");
   if (fp==NULL)
   {
      printf("We can't open the file (%s).\n", filename);
      exit(1);
   }
   else
   { int iR,iC;
      for (iR = 0; iR < nRows; ++iR )
      {
         for (iC = 0; iC < nCols; ++iC )
         {
int x;
          fscanf(fp, "%d" , &x );
_Bool bb=(x!=0);
data[iR][iC]=bb;
                     }
                      }

                              fclose(fp);
}
}



void readDoubleArray(char *filename, int nRows, int nCols, double ** data )
{

   FILE *fp = fopen (filename, "r");
   if (fp==NULL)
   {
      printf("We can't open the file (%s).\n", filename);
      exit(1);
   }
   else
   { int iR,iC;
      for (iR = 0; iR < nRows; ++iR )
      {
         for (iC = 0; iC < nCols; ++iC )
         {
            fscanf(fp, "%lf" , &data[iR][iC] );
                     }
                           }
                              fclose(fp);
                             }
}
double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        double **m;

        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++)
   {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}

void nrerror(char error_text[])
{
        printf("Utils run-time error...\n");
        printf("%s\n",error_text);
        printf("...now exiting to system...\n");
        exit(1);
}
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));

        free((char*) (m+nrl));
}
void findc(int n,_Bool R[n],int a,int * IDX, int *nx)
{
int ii_data[n];
int idx = 0;
int  ii = 0;
_Bool  exitg2 = 0;
_Bool guard2=0;
  while ((exitg2 == 0) && (ii < n)) {
    guard2 = 0;
    if (R[ii]!= a) {
      idx++;
      ii_data[idx - 1] = ii;
      if (idx >= n) {
        exitg2 = 1;
      } else {
        guard2 = 1;
      }
    } else {
      guard2 = 1;
    }

    if (guard2 == 1) {
      ii++;
    }
  }

int loop_ub=idx;
 for (idx = 0; idx < loop_ub; idx++) {
    IDX[idx] = ii_data[idx];
  }
*nx=loop_ub;

}

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        _Bool **m;

        m=(_Bool **) malloc( (nrh-nrl+1)*sizeof(_Bool*));
        if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++)
   {
                m[i]=(_Bool *) malloc((nch-ncl+1)*sizeof(_Bool));
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch)
{
        int i;

        for(i=nrh;i>=nrl;i--) free((_Bool*) (m[i]+ncl));

        free((_Bool*) (m+nrl));
}


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
void SampleIntercept(gsl_rng * rr,int n, int r, double * intercept, double* sigma2, double sigma20,double ** U, double ** A, double **y){
int i,l;
double meany=0;
for (i = 0; i < n; i++){
double ua=0;
for (l=0;l<r;l++){
ua+=U[i][l]*A[l][0];
}
meany+=y[i][0]-ua;
}
meany=meany/n;
double invsig2=n/sigma2[0]+1/sigma20;
*intercept=(n/(invsig2*sigma2[0]))*meany+sqrt(1/invsig2)*gsl_ran_ugaussian(rr);

}
void logPostGam(double *logpo,int * IndVar,int Np,int r, int n,int* P,double *** Tau, double ** U,double *** X,double **s2,_Bool ** rho,_Bool *** Gam,double** qv,double* q){
double logpost=0;
int m,j,l;
for (m=0;m<Np;m++){
	for (j=0;j<P[m];j++){
double logp;double quad;
	logGausQuadForm(j,r, n,P[m], Tau[m], U,X[m], s2[m][j],&quad,Gam[m][j],&logp);
	logpost+=logp;
	for (l=0;l<r;l++){
	if (IndVar[m]!=2)
	logpost+=Gam[m][j][l]*log(qv[m][l])+(1-Gam[m][j][l])*log(1-qv[m][l]);
	}
	}
	//if (m<Np-1){
	if (IndVar[m]==0){
		for (l=0;l<r;l++){
logpost+=rho[m][l]*log(q[m])+(1-rho[m][l])*log(1-q[m]);
	}
	}
}
*logpo=logpost;
}

void logGausQuadForm(int j,int r, int n,int p,double ** Tau, double ** U,double ** X,double s2,double* quadForm,_Bool * Gam,double * loggauss){
int i,s,s1;
int NZ1[r];
int nx1=0;
findc(r,Gam,0,NZ1, &nx1);

double result=0;double quadF=0;
if (nx1>0){
gsl_vector *work1 =gsl_vector_alloc (nx1);
double * Sigma1=malloc(nx1*nx1*sizeof(double));
for (s=0;s<nx1;s++){
for (s1=0;s1<=s;s1++){
double a=0;
 for (i=0;i<n;i++){
 //if (Gam[NZ1[s]]*Gam[NZ1[s1]]==1)
  a+=U[i][NZ1[s]]*U[i][NZ1[s1]];
}   
  Sigma1[s*nx1+s1]=Sigma1[s1*nx1+s]=a;

}
Sigma1[s*nx1+s]+=(1.0/Tau[NZ1[s]][j]);
//printf("q1=%lf\n",Tau[NZ1[s]][j]);
}
gsl_matrix_view m1  = gsl_matrix_view_array (Sigma1, nx1,nx1);
//printf("Error on Gam\n");

gsl_linalg_cholesky_decomp (&m1.matrix);
gsl_vector *xi =gsl_vector_alloc (n);
for (i = 0; i < n; i++)
    {
      gsl_vector_set (xi, i, X[i][j]/sqrt(s2));
    }

gsl_vector *y =gsl_vector_alloc (nx1);
double sumT=0;
for (s=0;s<nx1;s++){
double a=0;                            
//if (Gam[NZ1[s]]==1){
sumT+=log(Tau[NZ1[s]][j]);
for (i=0;i<n;i++){
a+=U[i][NZ1[s]]*X[i][j]/sqrt(s2);
}
//}
gsl_vector_set (y, s, a);
}
logmultigaussianT(xi, y,  &m1.matrix,&result,&quadF, work1);
result=result-0.5*sumT-(n/2.0)*log(s2);
gsl_vector_free (y);gsl_vector_free (work1);gsl_vector_free (xi);
free(Sigma1); 
} else {
for (i = 0; i < n; i++){
quadF+=pow(X[i][j],2);
}
quadF=quadF/s2;
result=-(n/2.0)*log(s2)- 0.5*n*log(2.0*M_PI)-0.5*quadF;
}

/*gsl_ran_multivariate_gaussian_log_pdf (xj, mu, &m.matrix, result, work);*/
*quadForm=quadF;
*loggauss=result;

}

void sampleGam(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X,double *s2,double* quadForm,_Bool ** Gam,double * loggauss,double q){
int j,s,l;
for (j=0;j<p;j++){
for (l=0;l<r;l++) Gam[j][l]=rho[l]*Gam[j][l];
}

double phi=0.5;
int NZ1[r];
int nx1=0;
findc(r,rho,0,NZ1, &nx1);
if (nx1>0){
for (j=0;j<p;j++){
//for (l=0;l<r;l++) Gam[j][l]=rho[l]*Gam[j][l];
_Bool Gamnew2[r];
for (l=0;l<r;l++) Gamnew2[l]=0;
double logpostold=0;double logpostnew=0;
_Bool Gamold1[nx1];
for (s=0;s<nx1;s++){
Gamold1[s]=Gam[j][NZ1[s]];
}
_Bool Gamnew1[nx1];
proposal(nx1,Gamold1,Gamnew1,phi, rr);
for (s=0;s<nx1;s++){
Gamnew2[NZ1[s]]=Gamnew1[s];
logpostold+=Gamold1[s]*log(q)+(1-Gamold1[s])*log(1-q);
logpostnew+=Gamnew1[s]*log(q)+(1-Gamnew1[s])*log(1-q);
}

/*Logpost new*/
double loggaussnew=0;double quadForm1=0;
logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quadForm1,Gamnew2,&loggaussnew);
//void logGausQuadForm(int j,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X,double s2,double* quadForm,_Bool * Gam,double * loggauss)
logpostnew+=loggaussnew;
logpostold+=loggauss[j];
double u=gsl_ran_flat (rr, 0, 1);
if ((log(u)<logpostnew-logpostold)){
for (l=0;l<r;l++){
Gam[j][l]=Gamnew2[l];
}
loggauss[j]=loggaussnew;
quadForm[j]=quadForm1;
}
}
}
}

double  logpost(int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,_Bool ** Gam,double * loggauss){
double logpostrho=0;
int l,j;
for (j=0;j<p;j++){
logGausQuadForm(j,r, n,p,Tau, U,X,s2[j],&quadForm[j],Gam[j],&loggauss[j]);
logpostrho+= loggauss[j];
}
for (l=0;l<r;l++){
logpostrho+=rho[l]*log(q)+(1-rho[l])*log(1-q);
}
return logpostrho;
}


void rho1(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,_Bool** Gam,double *loggauss){
int l,j;
double logpostold=0;
double logprior=0;
for (l=0;l<r;l++){
logprior+=rho[l]*log(q)+(1-rho[l])*log(1-q);
}
for (j=0;j<p;j++){
logpostold+=loggauss[j];
}
logpostold+=logprior;
double * quadForm1=malloc(p*sizeof(double));
double * loggauss1=malloc(p*sizeof(double));
/*double logpostold=logpost(r, n,p, rho,Tau,U,X, q,s2,quadForm1,Gam);*/
double phi=0.5;
_Bool rhonew[r];
proposal(r,rho,rhonew,phi, rr);
double logpostnew=logpost(r, n,p, rhonew,Tau,U,X, q,s2,quadForm1,Gam,loggauss1);
double u=gsl_ran_flat (rr, 0, 1);
if (log(u)<logpostnew-logpostold)
{
for (l=0;l<r;l++){
rho[l]=rhonew[l];
}
for (j=0;j<p;j++){
quadForm[j]=quadForm1[j];
loggauss[j]=loggauss1[j];
}
} 
}

void SamplerhoGamma(gsl_rng * rr,int r, int n,int IndVar,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double* q1,double q2,double* s2,double* quadForm,_Bool** Gam,double *loggauss){
int l,j;
_Bool *rhonew=malloc(r*sizeof(_Bool));
_Bool Gamnew[p][r];
for (l=0;l<r;l++){
for (j=0;j<p;j++){Gamnew[j][l]=Gam[j][l];}
rhonew[l]=rho[l];
}

if (IndVar==1){
for (l=0;l<r;l++){
double quad1=0;double quad2=0;
double loggaussold=0; double loggaussnew=0;
Gam[0][l]=0;
logGausQuadForm(0,r, n,p, Tau,  U,X,s2[0],&quad1,Gam[0],&loggaussold);
Gamnew[0][l]=1;
logGausQuadForm(0,r, n,p,Tau,  U,X,s2[0],&quad2,Gamnew[0],&loggaussnew);
double rat=loggaussnew+log(q2)-loggaussold-log(1-q2);
double uni=gsl_ran_flat (rr, 0, 1);
//printf("RAT=%.2lf",rat);
if (log(uni/(1-uni))<rat){
Gamnew[0][l]=Gam[0][l]=rho[l]=1;
loggauss[0]=loggaussnew;
quadForm[0]=quad2;
} else {
Gamnew[0][l]=Gam[0][l]=rho[l]=0;
loggauss[0]=loggaussold;
quadForm[0]=quad1;
}
}

} else {

for (l=0;l<r;l++){
double logq=0;
if (rho[l]==0){
rhonew[l]=1;
double logpostnew=0;
double logpostold=0;
double * quadForm1=malloc(p*sizeof(double));
double * loggausnew1=malloc(p*sizeof(double));
double quad1=0; double quad2=0;
double loggaussold=0; double loggaussnew=0;
for (j=0;j<p;j++){
double logqj=0;double logmqj=0;
logpostold+=loggauss[j];
Gam[j][l]=0;
logGausQuadForm(j,r, n,p, Tau,  U,X,s2[j],&quad1,Gam[j],&loggaussold);
Gamnew[j][l]=1;
logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad2,Gamnew[j],&loggaussnew);
double rat=-loggaussold+loggaussnew+log(q1[l])-log(1-q1[l]);
//printf("%lf ",rat);
double x1=loggaussnew+log(q1[l]);
double x2=loggaussold+log(1-q1[l]);
double maxq=MAX(x1,x2);
double den=maxq+log(exp(x1-maxq)+exp(x2-maxq));
logqj=x1-den;
logmqj=x2-den;

double uni=gsl_ran_flat (rr, 0, 1);
if ((log(uni/(1-uni))<rat)){
Gam[j][l]=1;Gamnew[j][l]=1;
logpostnew+=loggaussnew+log(q1[l]);
loggausnew1[j]=loggaussnew;
logq+=logqj;
quadForm1[j]=quad2;
} else {
Gam[j][l]=0;Gamnew[j][l]=0;
logq+=logmqj;
logpostnew+=loggaussold+log(1-q1[l]);
quadForm1[j]=quad1;
loggausnew1[j]=loggaussold;
}
}
logpostnew+=log(q2);
logpostold+=log(1-q2);
double un=gsl_ran_flat (rr, 0, 1);
//printf("%lf ",logpostnew-logpostold-logq);
double rat1=logpostnew-logpostold-logq;
if (log(un)<rat1){
rho[l]=1;
for (j=0;j<p;j++){
quadForm[j]=quadForm1[j];
loggauss[j]=loggausnew1[j];
}
} else {
rho[l]=0;
for (j=0;j<p;j++)
Gam[j][l]=Gamnew[j][l]=0;
}
rhonew[l]=rho[l];
free(quadForm1);free(loggausnew1);
} else {
rhonew[l]=0;
double logpostnew=0;
double logpostold=0;
double * quadForm1=malloc(p*sizeof(double));
double * loggausnew1=malloc(p*sizeof(double));
double quad1;
double loggaussnew=0; double loggaussn=0;
double logpq=0;
for (j=0;j<p;j++){
logpostold+=loggauss[j];
Gamnew[j][l]=1;
logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad1,Gamnew[j],&loggaussn);
Gamnew[j][l]=0;
logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad1,Gamnew[j],&loggaussnew);
if (p!=1){
double x1=loggaussn+log(q1[l]);
double x2=loggaussnew+log(1-q1[l]);
double maxq=MAX(x1,x2);
double den=maxq+log(exp(x1-maxq)+exp(x2-maxq));
double logqj=x1-den;
double logmqj=x2-den;
logq+=Gam[j][l]*log(q1[l])+(1-Gam[j][l])*log(1-q1[l]);
logpq+=Gam[j][l]*logqj+(1-Gam[j][l])*logmqj;
}
logpostnew+=loggaussnew;
loggausnew1[j]=loggaussnew;
quadForm1[j]=quad1;
Gamnew[j][l]=Gam[j][l];
}
logpostnew+=log(1-q2)+logpq;
logpostold+=log(q2)+logq;
double un=gsl_ran_flat (rr, 0, 1);
double rat1=logpostnew-logpostold;
if (log(un)<rat1){
rho[l]=0;
for (j=0;j<p;j++){
Gamnew[j][l]=Gam[j][l]=0;
quadForm[j]=quadForm1[j];
loggauss[j]=loggausnew1[j];
}
} else {
for (j=0;j<p;j++){
Gamnew[j][l]=Gam[j][l];
}
}
free(quadForm1);free(loggausnew1);
rhonew[l]=rho[l];
}
/* Within move*/
if (rho[l]==1){
double logpostold=0;
double * quadForm1=malloc(p*sizeof(double));
double * loggausnew1=malloc(p*sizeof(double));
double quad1,quad2;
double loggaussold,loggaussnew;
for (j=0;j<p;j++){
logpostold+=loggauss[j];
logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad1,Gam[j],&loggaussold);
Gamnew[j][l]=1-Gam[j][l];
logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad2,Gamnew[j],&loggaussnew);
double rat=-loggaussold+loggaussnew;
double uni=gsl_ran_flat (rr, 0, 1);
if (Gam[j][l]==0){
quad1=quad2;
loggaussold=loggaussnew;
rat+=log(q1[l])-log(1-q1[l]);
} else {
quad2=quad1;
loggaussnew=loggaussold;
rat=-rat+log(q1[l])-log(1-q1[l]);
}
if (log(uni/(1-uni))<rat){
Gam[j][l]=1;
quadForm[j]=quad1;
loggauss[j]=loggaussold;
} else {
Gam[j][l]=0;
quadForm[j]=quad2;
loggauss[j]=loggaussnew;
}
Gamnew[j][l]=Gam[j][l];
}
free(quadForm1);free(loggausnew1);
}
}
}
free(rhonew); 
}


void sigma2(gsl_rng * rr, int p,int n,double a0, double b0, double* quadForm,double * s2){
int j;
for (j=0;j<p;j++){
double inv=1/(b0+0.5*s2[j]*quadForm[j]);
s2[j]=1/gsl_ran_gamma (rr, a0+n/2.0, inv);
}
}

void EffectZero(int l,gsl_rng * rr,int K,int p,_Bool rho, _Bool * R,double* B,double *B0,_Bool ** Path,double alphab0, double betab0,double alpha, double * lambda2,double *AcceptB0,_Bool ** Gam){
double alphasb0=1;double betasb0=1;
int j,k1;
double B0new=gsl_ran_gamma (rr, alphasb0, 1/betasb0);;
double lograt=0;
double kl=0;
for (j=0;j<p;j++){
double kb=0;double kb1=0;
if (Gam[j][l]==1){
for (k1=0;k1<K;k1++){
kb+=Path[j][k1]*B[k1];
}
kl+=lambda2[j];
kb=kb+*B0;
kb1=kb+B0new;
lograt+=alpha*(log(kb1)-log(kb));
}
}
if (rho==1){
lograt+=(betasb0-betab0-kl)*(B0new-*B0)+(-alphasb0+alphab0)*(log(B0new)-log(*B0));
double uni=gsl_ran_flat (rr, 0, 1);
if (log(uni)<lograt){
*B0=B0new;
*AcceptB0+=1;
} 
} else {
*B0=gsl_ran_gamma (rr, alphab0, 1/betab0);
} 
/*} 
else {
*B0=gsl_ran_gamma (rr, p*alpha+alphab0, 1/(kl+betab0));
*AcceptB0+=1;
//printf("Bo=%lf",kl+betab0);
}
*/

}


void TauLambda(int l,gsl_rng * rr,int K,int p,double *A,double* B,double B0,double * Tau,_Bool ** Path,double alpha, double * lambda2,double *s2,_Bool ** Gam){
/* K is the number of pathways
 *  * R is the binary vector for pathway effect
 *   * Path is the grouping information; it is a binary matrix
 *    * B is the pathway effect
 *     * alphab and betab are priors for B effects
 *      * l is the component index
 *       */
int j,k1;
for (j=0;j<p;j++){
//if (fabs(A[j])>0){
if ((Gam[j][l]==1)&& (fabs(A[j])>0)){
double mu=sqrt(2*lambda2[j]*s2[j]/pow(A[j],2.0));
Tau[j]=1.0/inverseGaussian(rr,  mu, 2*lambda2[j]);
if (Tau[j]<0)
//printf("TauT=%lf \n",Tau[j]);
if (((Tau[j]-Tau[j])!=0)|(Tau[j]<0)){
Tau[j]=1/mu+1/(2*lambda2[j]);
//printf("J=%d, MU=%lf, lamb=%lf, S2=%lf, Beta=%lf B0=%lf\n", j,mu,lambda2[j],s2[j],A[j],B0);
}
double kb=0;
for (k1=0;k1<K;k1++){
kb+=Path[j][k1]*B[k1];
}
lambda2[j]=gsl_ran_gamma (rr, alpha+1, 1/(B0+kb+Tau[j]));
} else {
Tau[j]=gsl_ran_exponential (rr, 1.0/lambda2[j]);
//printf("Tau=%lf ",Tau[j]);
//lambda2[j]=gsl_ran_gamma (rr, alpha+1, 1/(Tau[j]));
lambda2[j]=gsl_ran_gamma (rr, alpha, 1/B0);
}
}
}


void GroupEffect(int l,gsl_rng * rr,_Bool rho,int K,int p, _Bool * R,double *A,double* B,double B0,double * Tau,_Bool ** Path,double alphab, double betab,double w,double alpha, double * lambda2,double *s2,double * AcceptR,_Bool **Gam){
/* K is the number of pathways
 * R is the binary vector for pathway effect
 * Path is the grouping information; it is a binary matrix
 * B is the pathway effect
 * alphab and betab are priors for B effects
 * l is the component index
 */
int j, k,k1;

if (rho==0){
for (k=0;k<K;k++){
B[k]=0;R[k]=0;
}
} else {

double alphas=2;double betas=2; // proposal parameters 
for (k=0;k<K;k++){
_Bool Rknew=0;double Bknew=0;
/* Between model move*/
if (R[k]==1){
Rknew=0; Bknew=0;
double lograt=0;
double kl=0;
for (j=0;j<p;j++){
if (Gam[j][l]==1){
double kb=0;double kb1=0;
for (k1=0;k1<K;k1++){
kb+=Path[j][k1]*B[k1];
}
kb1=kb-Path[j][k]*B[k];
lograt+=alpha*(log(B0+kb1)-log(B0+kb));
kl+=Path[j][k]*lambda2[j];
}
}
lograt+=(betab+kl-betas)*B[k]+(alphas-alphab)*log(B[k])+gsl_sf_lngamma(alphab)-gsl_sf_lngamma(alphas)+alphas*log(betas)-alphab*log(betab)+log(1-w)-log(w);
double uni=gsl_ran_flat (rr, 0, 1);
if (log(uni)<lograt){
B[k]=Bknew; R[k]=Rknew;
AcceptR[k]+=1;
}
}
else if (R[k]==0){
Rknew=1; Bknew=gsl_ran_gamma (rr, alphas, 1/betas);;
double lograt=0;
double kl=0;
for (j=0;j<p;j++){
if (Gam[j][l]==1){
double kb=0;double kb1=0;
for (k1=0;k1<K;k1++){
kb+=Path[j][k1]*B[k1];
}
kb1=kb-Path[j][k]*B[k];
kb=kb1+Path[j][k]*Bknew;
lograt+=alpha*(log(B0+kb)-log(B0+kb1));
kl+=Path[j][k]*lambda2[j];
}
}
lograt+=(betas-betab-kl)*Bknew+(-alphas+alphab)*log(Bknew)-gsl_sf_lngamma(alphab)+gsl_sf_lngamma(alphas)-alphas*log(betas)+alphab*log(betab)-log(1-w)+log(w);
double uni=gsl_ran_flat (rr, 0, 1);
if (log(uni)<lograt){
B[k]=Bknew; R[k]=Rknew;
AcceptR[k]+=1;
}
}
/* Within model move*/

if (R[k]==1){
Bknew=gsl_ran_gamma (rr, alphas, 1/betas);;
double lograt=0;
double kl=0;
for (j=0;j<p;j++){ 
if (Gam[j][l]==1){
double kb=0;double kb1=0;
for (k1=0;k1<K;k1++){
kb+=Path[j][k1]*B[k1];
}
kb1=kb-Path[j][k]*B[k]+Path[j][k]*Bknew;
lograt+=alpha*(log(B0+kb1)-log(B0+kb));
kl+=Path[j][k]*lambda2[j];
}
}
lograt+=(betas-betab-kl)*(Bknew-B[k])+(-alphas+alphab)*(log(Bknew)-log(B[k]));
double uni=gsl_ran_flat (rr, 0, 1);
if (log(uni)<lograt){
B[k]=Bknew;
}
}
}// close for k
} // End for else if (rho==0)
}//close for the function



void LoadAOther(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,_Bool ** Gam){
int s,s1,j,i;
for (j=0;j<p;j++){
for (s=0;s<r;s++){
A[s][j]=0;
}
int NZ1[r];
int nx1=0;
findc(r,Gam[j],0,NZ1, &nx1);
if (nx1>0){
double * SigmaInv=malloc(nx1*nx1*sizeof(double));
for (s=0;s<nx1;s++){
for (s1=0;s1<=s;s1++){
double a=0;
for (i=0;i<n;i++){
        a+=U[i][NZ1[s]]*U[i][NZ1[s1]];
}
  SigmaInv[s*nx1+s1]=SigmaInv[s1*nx1+s]=a/s2[j];
}
SigmaInv[s*nx1+s]+=1/(Tau[s][j]*s2[j]);
}
gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, nx1,nx1);
gsl_linalg_cholesky_decomp (&m.matrix);
gsl_vector *mu =gsl_vector_alloc (nx1);
double xy[nx1];
for (s=0; s<nx1;s++){
xy[s]=0;
for (i=0; i<n;i++){
xy[s]+=U[i][NZ1[s]]*X[i][j]/s2[j];}
}
gsl_vector_view b= gsl_vector_view_array (xy, nx1);
gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
gsl_vector *result=gsl_vector_alloc (nx1);
multivariate_gaussian (rr, mu, &m.matrix,result);
for (s=0; s<nx1;s++){
A[NZ1[s]][j]=gsl_vector_get(result,s);
}
free(SigmaInv); gsl_vector_free (mu);gsl_vector_free (result);
}
}
}




void EstimateLoad(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,_Bool ** Gam){
int s,s1,j,i;
for (j=0;j<p;j++){
for (s=0;s<r;s++){
A[s][j]=0;
}
int NZ1[r];
int nx1=0;
findc(r,Gam[j],0,NZ1, &nx1);
if (nx1>0){
double * SigmaInv=malloc(nx1*nx1*sizeof(double));
for (s=0;s<nx1;s++){
for (s1=0;s1<=s;s1++){
double a=0;
for (i=0;i<n;i++){
        a+=U[i][NZ1[s]]*U[i][NZ1[s1]];
}
  SigmaInv[s*nx1+s1]=SigmaInv[s1*nx1+s]=a/s2[j];
}
SigmaInv[s*nx1+s]+=1/(Tau[s][j]*s2[j]);
}
gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, nx1,nx1);
TRY
   {
gsl_linalg_cholesky_decomp (&m.matrix);
}
CATCH
   {
printf("\nError on the sampling of A.");
}
ETRY;
gsl_vector *mu =gsl_vector_alloc (nx1);
double xy[nx1];
for (s=0; s<nx1;s++){
xy[s]=0;
for (i=0; i<n;i++){
xy[s]+=U[i][NZ1[s]]*X[i][j]/s2[j];}
}
gsl_vector_view b= gsl_vector_view_array (xy, nx1);
gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
//gsl_vector *result=gsl_vector_alloc (nx1);
//multivariate_gaussian (rr, mu, &m.matrix,result);
for (s=0; s<nx1;s++){
A[NZ1[s]][j]=gsl_vector_get(mu,s);
}
free(SigmaInv); gsl_vector_free (mu);
//gsl_vector_free (result);
}
}
}



void proposal(int n,_Bool R[n],_Bool R1[n],float phi, gsl_rng * r)
{
int i;
for (i=0; i<n;i++) 
R1[i]=R[i];
int n1=0;//number of 1
int n0=0;//number of zeros
int v1[n];
findc(n,R,1,v1,&n0);//find indices different from 1 i.e ==0;
int v2[n-n0];
findc(n,R,0,v2,&n1);//find indices different of zeros i.e ==1
double u=gsl_ran_flat (r, 0, 1);
//printf("No==%d, N1==%d\n",n0,n1);
if ((u < phi) || (n0 == 0) || (n1 == 0)) {
int l= gsl_rng_uniform_int(r,n);

   R1[l] = 1 - R[l];
  } else {
int l1=gsl_rng_uniform_int(r,n0);
int l2=gsl_rng_uniform_int(r,n1);

    R1[v1[l1]] = R[v2[l2]];
    R1[v2[l2]] = R[v1[l1]];
  }
}

void SampleUU(gsl_rng * rr,int r, int n,int Np,int * P,double *** A, double ** U,double *** X,double** s2){

//void SampleUU(gsl_rng * rr,int r, int n,int p0,int p1,int p2,double *** A, double ** U,double *** X,double** s2){	
	
int s,s1,j,i;
int sumMark=0;
for (i=0;i<Np;i++) sumMark+=P[i]; 

double * SigmaInv=malloc(r*r*sizeof(double));
for (s=0;s<r;s++){
for (s1=0;s1<=s;s1++){
double a=0;
for (j=0;j<sumMark;j++){
	int k=0;
	for (i=0;i<Np;i++){
	if ((k<=j) && (j<k+P[i])) a+=A[i][s][j-k]*A[i][s1][j-k]/s2[i][j-k]; 
	k+=P[i];
	}
}
SigmaInv[s*r+s1]=SigmaInv[s1*r+s]=a;
}
SigmaInv[s*r+s]+=1;
}

gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, r,r);
//gsl_linalg_cholesky_decomp (&m.matrix);
TRY
   {
//printf("\nError on the sampling of U.");
gsl_linalg_cholesky_decomp (&m.matrix);
//THROW;
}
CATCH
   {
printf("\nError on the sampling of U.");
                      //exit(1);
}
ETRY;
int l;
for (i=0;i<n;i++){
gsl_vector *mu =gsl_vector_alloc (r);
double Ax[r];
for (s=0; s<r;s++){
Ax[s]=0;
for (j=0;j<sumMark;j++){
	        int k=0;
		        for (l=0;l<Np;l++){
     if ((k<=j) && (j<k+P[l])) Ax[s]+=X[l][i][j-k]*A[l][s][j-k]/s2[l][j-k];//a+=A[i][s][j-k]*A[i][s1][j-k]/s2[i][j-k];
					        k+=P[l];
						        }
}
//for (j=0; j<p0+p1+p2;j++){
//if (j<p0)
//Ax[s]+=X[0][i][j]*A[0][s][j]/s2[0][j];
//else if (j<p0+p1)
//Ax[s]+=X[1][i][j-p0]*A[1][s][j-p0]/s2[1][j-p0];
//else Ax[s]+=X[2][i][j-p1-p0]*A[2][s][j-p1-p0]/s2[2][j-p1-p0];
//}
}
gsl_vector_view b= gsl_vector_view_array (Ax, r);
//printf("Size 1 is %zu ",(&m.matrix)->size2); printf("Size 2 is %zu ",mu->size);
gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
gsl_vector *result=gsl_vector_alloc (r);
multivariate_gaussian (rr, mu, &m.matrix,result);
for (s=0; s<r;s++){
U[i][s]=gsl_vector_get(result,s);
//if (i<20)
//printf("%.2lf ",U[i][s]);
}
 gsl_vector_free (mu);gsl_vector_free (result);
}
free(SigmaInv);
}


void SampleU(gsl_rng * rr,int r, int n,int p0,int p1,int p2,double *** A, double ** U,double *** X,double** s2){
int s,s1,j,i;
double * SigmaInv=malloc(r*r*sizeof(double));
for (s=0;s<r;s++){
for (s1=0;s1<=s;s1++){
double a=0;
for (j=0;j<p0+p1+p2;j++){
if (j<p0) a+=A[0][s][j]*A[0][s1][j]/s2[0][j];
else if (j<p0+p1) a+=A[1][s][j-p0]*A[1][s1][j-p0]/s2[1][j-p0];
else a+=A[2][s][j-p1-p0]*A[2][s1][j-p1-p0]/s2[2][j-p1-p0];
}
SigmaInv[s*r+s1]=SigmaInv[s1*r+s]=a;
}
SigmaInv[s*r+s]+=1;
}
gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, r,r);
//gsl_linalg_cholesky_decomp (&m.matrix);
TRY
   {
//printf("\nError on the sampling of U.");
gsl_linalg_cholesky_decomp (&m.matrix);
//THROW;
}
CATCH
   {
printf("\nError on the sampling of U.");
                      //exit(1);
}
ETRY;

for (i=0;i<n;i++){
gsl_vector *mu =gsl_vector_alloc (r);
double Ax[r];
for (s=0; s<r;s++){
Ax[s]=0;
for (j=0; j<p0+p1+p2;j++){
if (j<p0)
Ax[s]+=X[0][i][j]*A[0][s][j]/s2[0][j];
else if (j<p0+p1)
Ax[s]+=X[1][i][j-p0]*A[1][s][j-p0]/s2[1][j-p0];
else Ax[s]+=X[2][i][j-p1-p0]*A[2][s][j-p1-p0]/s2[2][j-p1-p0];
}
}
gsl_vector_view b= gsl_vector_view_array (Ax, r);
//printf("Size 1 is %zu ",(&m.matrix)->size2); printf("Size 2 is %zu ",mu->size);
gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
gsl_vector *result=gsl_vector_alloc (r);
multivariate_gaussian (rr, mu, &m.matrix,result);
for (s=0; s<r;s++){
U[i][s]=gsl_vector_get(result,s);
//if (i<20)
//printf("%.2lf ",U[i][s]);
}
 gsl_vector_free (mu);gsl_vector_free (result);
}
free(SigmaInv);
}


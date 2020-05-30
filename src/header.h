#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) > (b) ? b : a)
#include <setjmp.h>
#define TRY do{ jmp_buf ex_buf__; if( !setjmp(ex_buf__) ){
#define CATCH } else {
#define ETRY } }while(0)
#define THROW longjmp(ex_buf__, 1)

void mainfunction(int *Method1,int * n1,int *P,int *r1,int *Np1,double *datasets,int * IndVar,int *K,int * Paths,int* maxmodel1, int *nbrsample1,int *burninsample1,double *CompoSelMean,double *VarSelMean,double * VarSelMeanGlobal,double * GrpSelMean,double *GrpEffectMean,double * IntGrpMean,double * EstU, double * EstSig2,double *InterceptMean,double *EstLoadMod,double *EstLoad,int *nbrmodel1,double *postgam,double* priorcompsel,double * priorcompselo,double* priorb0,double* priorb,double *priorgrpsel,double *probvarsel);
void SampleUU(gsl_rng * rr,int r, int n,int Np,int * P,double *** A, double ** U,double *** X,double** s2);
double mean(int n,double * x);
void SampleIntercept(gsl_rng * rr,int n, int r, double * intercept, double* sigma2, double sigma20,double ** U, double ** A, double **y);
_Bool compar(int n,_Bool *u,_Bool *v);
void logPostGam(double *logpo,int * IndVar,int Np,int r, int n,int* P,double *** Tau, double ** U,double *** X,double **s2,_Bool ** rho,_Bool *** Gam,double** qv,double* q);
_Bool **  UniqueModel(int nbrsample, int p, _Bool ** rhosample,int * modelidx,int *countd1);
void SamplerhoGamma(gsl_rng * rr,int r, int n,int IndVar,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double* q1,double q2,double* s2,double* quadForm,_Bool** Gam,double *loggauss);
void save1d(char *filename,int n,double *data);
void sampleGam(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X,double *s2,double* quadForm,_Bool ** Gam,double * loggauss,double q);
void logGausQuadForm(int j,int r, int n,int p, double ** Tau, double ** U,double ** X,double s2,double* quadForm,_Bool * Gam,double * loggauss);
void save2d(char *filename,int n,int p,double ** data);
void SampleU(gsl_rng * rr,int r, int n,int p0,int p1,int p2,double *** A, double ** U,double *** X,double** s2);
void EffectZero(int l,gsl_rng * rr,int K,int p,_Bool rho, _Bool * R,double* B,double *B0,_Bool ** Path,double alphab0, double betab0,double alpha, double * lambda2,double *AcceptB0,_Bool ** Gam);
void TauLambda(int l,gsl_rng * rr,int K,int p,double *A,double* B,double B0,double * Tau,_Bool ** Path,double alpha, double * lambda2,double *s2,_Bool ** Gam);
void GroupEffect(int l,gsl_rng * rr,_Bool rho,int K,int p, _Bool * R,double *A,double* B,double B0,double * Tau,_Bool ** Path,double alphab, double betab,double w,double alpha, double * lambda2,double *s2,double * AcceptR,_Bool **Gam);
void sigma2(gsl_rng * rr, int p,int n,double a0, double b0, double* quadForm,double * s2);
void LoadAOther(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,_Bool ** Gam);
void EstimateLoad(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,_Bool ** Gam);
int multivariate_gaussian (const gsl_rng * r,
                               const gsl_vector * mu,
                               const gsl_matrix * L,
                               gsl_vector * result);

int logmultigaussian(const gsl_vector * x,
                                       const gsl_vector * mu,
                                       const gsl_matrix * L,
                                       double * result,double *quadF,
                                       gsl_vector * work);

int logmultigaussianT(const gsl_vector * x, const gsl_vector * y,
                                       const gsl_matrix * L,
                                       double * result,double *quadF,
                                       gsl_vector * work);
double  logpost(int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,_Bool ** Gam,double * loggauss);
void rho1(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,_Bool** Gam,double *loggauss);
void proposal(int n,_Bool R[n],_Bool R1[n],float phi, gsl_rng * r);
void findc(int n,_Bool R[n],int a,int * IDX, int *nx);
double inverseGaussian(gsl_rng * r, double mu, double lambda);
void readBoolVector(char *filename, int nRows,  _Bool * data );
//void rho(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** A, double ** X,double a0, double b0, double q,double s2,double lambda2);
void sort(int n,double *x,int *idx);
double auc(int n, double * esti,_Bool class[n]);
void readDoubleArray(char *filename, int nRows, int nCols, double ** data );
void readBoolArray(char *filename, int nRows, int nCols, _Bool ** data );
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
_Bool **bmatrix(int nrl, int nrh, int ncl, int nch);
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch);

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
#include <time.h>
#include <gsl/gsl_errno.h>
#include "header.h"
// R CMD SHLIB MainAlGo.c  utils.c myfunction.c -lgsl -lgslcblas -o BayesianFA.so
// void mainfunction(const char *Method,int * n1,int *P,int *r1,double *dataset1, double *dataset2, double *outcome,int *K,int *path0, int *path1,int *nbrsample1,int *burninsample1,double *CompoSelMean,double *VarSelMean,double * VarSelMeanGlobal,double * GrpSelMean,double *GrpEffectMean,double * IntGrpMean){

void mainfunction(int *Method1, int *n1, int *P, int *r1, int *Np1, double *datasets, int *IndVar, int *K, int *Paths,
				  int *maxmodel1, int *nbrsample1, int *burninsample1, double *CompoSelMean, double *VarSelMean,
				  double *VarSelMeanGlobal, double *GrpSelMean, double *GrpEffectMean, double *IntGrpMean,
				  double *EstU, double *EstSig2, double *InterceptMean, double *EstLoadMod, double *EstLoad,
				  int *nbrmodel1, double *postgam, double *priorcompsel, double *priorcompselo,
				  double *priorb0, double *priorb, double *priorgrpsel, double *probvarsel, int *chainNbr,
				  double *AcceptB0Rate)
{

	/* GSL errors (e.g. a matrix that is not positive definite) are returned as status codes and
	 * handled by the callers, instead of calling abort() and killing the R session */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
	long seed = chainNbr[0];
	gsl_rng *rr = gsl_rng_alloc(gsl_rng_rand48);
	gsl_rng_set(rr, seed);

	setvbuf(stdout, NULL, _IONBF, 0);
	int Method = Method1[0];
	if (Method == 1)
		printf("\nThe Method is BIP with network information\n");
	else if (Method == 0)
	{
		printf("The Method is BIP without network information\n");
	}
	int i, l, j, k;
	int m;
	clock_t t1 = clock();
	int Np = Np1[0]; // Nber of platforms
	int nbrsample = nbrsample1[0];
	printf("Number of MCMC samples after burn-in is %d\n", nbrsample);
	int burninsample = burninsample1[0];
	printf("Number of burn-in is %d\n", burninsample);
	if ((nbrsample < 1) || (burninsample < 0))
	{ // the posterior means, UniqueModel and the model averaging need at least one sample after burn-in
		printf("Error: the number of MCMC samples after burn-in must be at least 1 and burn-in cannot be negative.\n");
		gsl_rng_free(rr);
		gsl_set_error_handler(old_handler);
		return;
	}
	int n = n1[0];
	printf("Number of samples is %d\n", n);
	// int p0=P[0];
	for (m = 0; m < Np; m++)
	{
		printf("Number of markers in platform %d is %d\n", m, P[m]);
	}
	int r = r1[0];
	printf("Number of components is %d\n", r);
	// int P[Np];P[0]=p0;P[1]=p1;;P[2]=p2;
	double ***X = malloc(Np * sizeof(double **));
	double ***X1 = malloc(Np * sizeof(double **));
	// double dat=0;
	k = 0;
	for (m = 0; m < Np; m++)
	{
		X[m] = dmatrix(0, n - 1, 0, P[m] - 1);
		X1[m] = dmatrix(0, n - 1, 0, P[m] - 1);
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < P[m]; j++)
			{
				X[m][i][j] = X1[m][i][j] = datasets[k];
				k++;
			}
		}
	}
	// int K0=K[0];int K1=K[1];int K2=K[2];
	/*
	for (m=0;m<Np;m++){
	if (IndVar[m]!=1)
	printf("Number of groups for platform %d is %d\n",m,K[m]);
	else printf("Number of groups for outcome is 1");
	}
	*/
	_Bool ***Path = malloc(Np * sizeof(_Bool **));
	//int m1 = 0; //_Bool pp=0;
	int kk = 0;
	for (m = 0; m < Np; m++)
	{
		Path[m] = bmatrix(0, P[m] - 1, 0, K[m] - 1);
		if (IndVar[m] == 1)
		{
			Path[m][0][0] = 1;
		}
		else if (IndVar[m] == 0)
		{
			for (j = 0; j < P[m]; j++)
			{
				for (k = 0; k < K[m]; k++)
				{
					Path[m][j][k] = Paths[kk];
					kk++;
					// if (m==0) pp=path0[j*K[m]+k]; else if (m==1) pp=path1[j*K[m]+k];
					// Path[m][j][k]=pp;
				}
			}
		//	m1 += 1;
		}
		else if (IndVar[m] == 2)
		{ // covariates
			for (j = 0; j < P[m]; j++)
			{
				for (k = 0; k < K[m]; k++)
				{
					Path[m][j][k] = 1;
				}
			}
		}
	}

	double **U = dmatrix(0, n - 1, 0, r - 1);
	double **meanU = dmatrix(0, n - 1, 0, r - 1);
	for (i = 0; i < n; i++)
	{
		for (l = 0; l < r; l++)
		{
			U[i][l] = gsl_ran_ugaussian(rr);
			meanU[i][l] = 0;
		}
	}

	double ***A = malloc(Np * sizeof(double **));
	for (m = 0; m < Np; m++)
	{
		A[m] = dmatrix(0, r - 1, 0, P[m] - 1);
	}

	/* Hyperparameter*/
	double *q = malloc(Np * sizeof(double));
	double **qv = malloc(Np * sizeof(double *));
	double **wg = malloc(Np * sizeof(double *));
	for (m = 0; m < Np; m++)
	{
		q[m] = 0.5;
		qv[m] = malloc(r * sizeof(double)); //
		wg[m] = malloc(r * sizeof(double)); // proba. for group selection
	}
	// double q=0.5;
	// double qv=0.5;
	// double a0=0.01;double b0=0.01;
	double a0 = 1;
	double b0 = 1;
	// double alphab=1;double betab=1;
	double alphab = priorb[0];
	double betab = priorb[1];
	// double al=1;double bl=1;// Hyper for q
	double al = priorcompsel[0];
	double bl = priorcompsel[1]; // Hyper for q
	// double al0=1;double bl0=1;// Hyper for q in the outcome
	double al0 = priorcompselo[0];
	double bl0 = priorcompselo[1]; // Hyper for q in the outcome
	// double alv=1;double blv=1;// Hyper for qv
	// double alphab0=2; double betab0=2;
	// double alv=priorselv[0];double blv=priorselv[1];
	double alphab0 = priorb0[0];
	double betab0 = priorb0[1]; // Hyper for b0
	// double w=0.5;//priors for group selection
	// double alg=1;double blg=1;
	double alg = priorgrpsel[0];
	double blg = priorgrpsel[1];
	double alpha = 1;
	// Initialization

	_Bool **rhoest = malloc(Np * sizeof(_Bool *));
	double **rhomean = malloc(Np * sizeof(double *));
	double **Gamvs = malloc(Np * sizeof(double *));
	//_Bool ** GamvsK=malloc(Np*sizeof(_Bool*));
	_Bool ***R = malloc(Np * sizeof(_Bool **));
	_Bool ***Gam = malloc(Np * sizeof(_Bool **));
	double ***Gammean = malloc(Np * sizeof(double **));
	double ***AcceptR = malloc(Np * sizeof(double **));
	double ***Bmean = malloc(Np * sizeof(double **));
	double ***Rmean = malloc(Np * sizeof(double **));
	double ***B = malloc(Np * sizeof(double **));
	double ***lambda2 = malloc(Np * sizeof(double **));
	double ***Tau = malloc(Np * sizeof(double **));
	double ***Taumean = malloc(Np * sizeof(double **));
	double **B0 = malloc(Np * sizeof(double *));
	double **B0mean = malloc(Np * sizeof(double *));
	double **AcceptB0 = malloc(Np * sizeof(double *));
	double **TryB0 = malloc(Np * sizeof(double *)); // number of MH attempts for B0 after burn-in
	double **sdB0 = malloc(Np * sizeof(double *));  // step size of the log(B0) random walk, tuned during burn-in
	double **quadForm = malloc(Np * sizeof(double *));
	double **loggauss = malloc(Np * sizeof(double *));
	double **s2 = malloc(Np * sizeof(double *));
	double **s2Mean = malloc(Np * sizeof(double *));
	//int mj = 0;
	for (m = 0; m < Np; m++)
	{
		Gamvs[m] = malloc(P[m] * sizeof(double));
		// GamvsK[m]=malloc(P[m]*sizeof(_Bool));
		Gam[m] = bmatrix(0, P[m] - 1, 0, r - 1);
		Gammean[m] = dmatrix(0, P[m] - 1, 0, r - 1);
		rhoest[m] = malloc(r * sizeof(_Bool));
		rhomean[m] = malloc(r * sizeof(double));
		R[m] = bmatrix(0, r - 1, 0, K[m] - 1);
		AcceptR[m] = dmatrix(0, r - 1, 0, K[m] - 1);
		Bmean[m] = dmatrix(0, r - 1, 0, K[m] - 1);
		Rmean[m] = dmatrix(0, r - 1, 0, K[m] - 1);
		B[m] = dmatrix(0, r - 1, 0, K[m] - 1);
		lambda2[m] = dmatrix(0, r - 1, 0, P[m] - 1);
		B0[m] = malloc(r * sizeof(double));
		B0mean[m] = malloc(r * sizeof(double));
		AcceptB0[m] = malloc(r * sizeof(double));
		TryB0[m] = malloc(r * sizeof(double));
		sdB0[m] = malloc(r * sizeof(double));
		Tau[m] = dmatrix(0, r - 1, 0, P[m] - 1);
		Taumean[m] = dmatrix(0, r - 1, 0, P[m] - 1);
		for (l = 0; l < r; l++)
		{
			rhomean[m][l] = 0;
			B0mean[m][l] = 0;
			B0[m][l] = 0.1;
			AcceptB0[m][l] = 0;
			TryB0[m][l] = 0;
			sdB0[m][l] = 1;
			for (j = 0; j < P[m]; j++)
				Gammean[m][j][l] = 0;
			double uni = gsl_ran_flat(rr, 0, 1);
			qv[m][l] = probvarsel[m];
			wg[m][l] = 0.5;		// Prior prob group selection
			if (IndVar[m] == 1) // Response
				qv[m][l] = 0.5;
			else if (IndVar[m] == 2)
				qv[m][l] = 1;
			if (uni < 0.5)
			{
				rhoest[m][l] = 1;
				for (j = 0; j < P[m]; j++)
				{
					if (gsl_ran_flat(rr, 0, 1) < qv[m][l])
					{
						Gam[m][j][l] = 1;
						A[m][l][j] = 0.01;
					}
					else
					{
						Gam[m][j][l] = 0;
						A[m][l][j] = 0;
					}
				}
			}
			else
			{
				rhoest[m][l] = 0;
				for (j = 0; j < P[m]; j++)
				{
					Gam[m][j][l] = 0;
					A[m][l][j] = 0;
				}
			}
			if (IndVar[m] == 1)
			{ // Response
				Gam[m][0][l] = rhoest[m][l];
			}
			else if (IndVar[m] == 2)
			{ // Clinical factors
				rhoest[m][l] = 1;
				for (j = 0; j < P[m]; j++)
				{
					Gam[m][j][l] = 1;
					A[m][l][j] = 0.01;
				}
			}
			for (k = 0; k < K[m]; k++)
			{
				Rmean[m][l][k] = 0;
				Bmean[m][l][k] = 0;
				AcceptR[m][l][k] = 0;
				R[m][l][k] = 0;
				B[m][l][k] = 0;
				// if ((strcmp(Method,"GroupInfo")==0)&& (P[m]!=1)){
				if ((Method == 1) && (IndVar[m] == 0))
				{
					if (rhoest[m][l] == 1)
					{
						double ui = gsl_ran_flat(rr, 0, 1);
						if (ui < wg[m][l])
						{
							R[m][l][k] = 1;
							B[m][l][k] = 0.1;
						}
					}
				}
			}
			for (j = 0; j < P[m]; j++)
			{
				Tau[m][l][j] = 1;
				Taumean[m][l][j] = 0;
				lambda2[m][l][j] = 1;
				if (IndVar[m] == 2)
				{
					Tau[m][l][j] = 100;
				}
			}
		}
		s2[m] = malloc(P[m] * sizeof(double));
		s2Mean[m] = malloc(P[m] * sizeof(double));
		quadForm[m] = malloc(P[m] * sizeof(double));
		loggauss[m] = malloc(P[m] * sizeof(double));
		for (j = 0; j < P[m]; j++)
		{
			loggauss[m][j] = 0;
			s2[m][j] = 0.1;
			s2Mean[m][j] = 0;
			quadForm[m][j] = Gamvs[m][j] = 0; // GamvsK[m][j]=0;
			//mj += 1;
		}
	}

	// Dim of a model
	int dim = 0;
	// Gam entries for all platforms, plus rho entries for non-response platforms (must match the rhomodel fill below)
	for (m = 0; m < Np; m++)
	{
		dim += P[m] * r;
		if (IndVar[m] != 1)
			dim += r;
	}

	int t;
	int N = nbrsample + burninsample;
	double intercp=0;
	// Cross-products reused by every likelihood evaluation: U'U, U'X1[m] and ||x_j||^2
	double **UtU = dmatrix(0, r - 1, 0, r - 1);
	double ***UtX = malloc(Np * sizeof(double **));
	double **xx = malloc(Np * sizeof(double *));
	for (m = 0; m < Np; m++)
	{
		UtX[m] = dmatrix(0, r - 1, 0, P[m] - 1);
		xx[m] = malloc(P[m] * sizeof(double));
	}
	_Bool **rhomodel = malloc(nbrsample * sizeof(_Bool *));
	for (t = 0; t < nbrsample; t++)
	{
		rhomodel[t] = malloc(dim * sizeof(_Bool));
	}
	for (t = 0; t < N; t++)
	{
		if (t == burninsample)
		{ // B0 acceptance rate is computed after burn-in only
			for (m = 0; m < Np; m++)
				for (l = 0; l < r; l++)
					AcceptB0[m][l] = 0;
		}
		crossprod(n, r, r, U, U, UtU);
		for (m = 0; m < Np; m++)
		{
			crossprod(n, r, P[m], U, X1[m], UtX[m]);
			colsumsq(n, P[m], X1[m], xx[m]);
			
			double sumrho = 0;
			for (l = 0; l < r; l++)
			{
				//double sumeta = 0;
				sumrho += rhoest[m][l];
				//for (j = 0; j < P[m]; j++)
				//{
				//	sumeta += Gam[m][j][l];
			//	}
				// qv[m][l]=gsl_ran_beta (rr, alv+rhoest[m][l]*sumeta, blv+rhoest[m][l]*(P[m]-sumeta));
				if (IndVar[m] == 2)
					qv[m][l] = 1;
			}
			if (IndVar[m] == 1)
			{
				q[m] = gsl_ran_beta(rr, al0 + sumrho, bl0 + r - sumrho);
				for (l = 0; l < r; l++)
				{
					qv[m][l] = q[m];
				}
			}
			else if (IndVar[m] == 0)
			{

				q[m] = gsl_ran_beta(rr, al + sumrho, bl + r - sumrho);
			}

			for (j = 0; j < P[m]; j++)
			{
				logGausQuadForm(j, r, n, Tau[m], UtU, UtX[m], xx[m][j], s2[m][j], &quadForm[m][j], Gam[m][j], &loggauss[m][j]);
				if (t == 0)
				{
					loggauss[m][j] = -DBL_MAX;
				}
			}

			if (IndVar[m] != 2)
			{
				SamplerhoGamma(rr, r, n, IndVar[m], P[m], rhoest[m], Tau[m], UtU, UtX[m], xx[m], qv[m], q[m], s2[m], quadForm[m], Gam[m], loggauss[m]);
			}

			sigma2(rr, P[m], n, a0, b0, quadForm[m], s2[m]);
			LoadAOther(rr, r, n, P[m], rhoest[m], Tau[m], A[m], UtU, UtX[m], s2[m], Gam[m]);

			for (l = 0; l < r; l++)
			{
				if (IndVar[m] != 2)
				{
					TauLambda(l, rr, K[m], P[m], A[m][l], B[m][l], B0[m][l], Tau[m][l], Path[m], alpha, lambda2[m][l], s2[m], Gam[m]);
				}
				if ((Method == 1) && (IndVar[m] == 0))
				{
					sumrho = 0;
					for (k = 0; k < K[m]; k++)
						sumrho += R[m][l][k];
					wg[m][l] = gsl_ran_beta(rr, alg + rhoest[m][l] * sumrho, blg + rhoest[m][l] * (K[m] - sumrho));
					GroupEffect(l, rr, rhoest[m][l], K[m], P[m], R[m][l], A[m][l], B[m][l], B0[m][l], Tau[m][l], Path[m], alphab, betab, wg[m][l], alpha, lambda2[m][l], s2[m], AcceptR[m][l], Gam[m]);
				}
				if ((t >= burninsample) && (rhoest[m][l] == 1))
					TryB0[m][l] += 1; // EffectZero makes an MH move only when rho == 1
				double acc0 = AcceptB0[m][l];
				EffectZero(l, rr, K[m], P[m], rhoest[m][l], R[m][l], B[m][l], &B0[m][l], Path[m], alphab0, betab0, alpha, lambda2[m][l], &AcceptB0[m][l], Gam[m], sdB0[m][l]);
				if ((t < burninsample) && (rhoest[m][l] == 1))
				{ // Robbins-Monro tuning towards an acceptance rate of 0.44; the step size is fixed after burn-in
					sdB0[m][l] *= exp((AcceptB0[m][l] - acc0 - 0.44) / sqrt(t + 1.0));
					sdB0[m][l] = MIN(MAX(sdB0[m][l], 0.01), 10);
				}
			}
			if (IndVar[m] == 1)
			{
				SampleIntercept(rr, n, r, &intercp, s2[m], 100.0, U, A[m], X[m]);
				for (i = 0; i < n; i++)
				{
					for (j = 0; j < P[m]; j++)
					{
						X1[m][i][j] = X[m][i][j] - intercp;
					}
				}
			}
		}
		SampleUU(rr, r, n, Np, P, A, Gam, U, X1, s2);

		if (t >= burninsample)
		{
			*InterceptMean += intercp / nbrsample;
			int rm = 0;
			for (m = 0; m < Np; m++)
			{
				for (j = 0; j < P[m]; j++)
				{
					for (l = 0; l < r; l++)
					{
						rhomodel[t - burninsample][rm] = Gam[m][j][l];
						rm++;
					}
				}
			}
			for (m = 0; m < Np; m++)
			{
				if (IndVar[m] != 1)
				{
					for (l = 0; l < r; l++)
					{
						rhomodel[t - burninsample][rm] = rhoest[m][l];
						rm++;
					}
				}
			}

			for (l = 0; l < r; l++)
			{
				for (i = 0; i < n; i++)
				{
					meanU[i][l] += U[i][l] / nbrsample;
				}
			}
			for (m = 0; m < Np; m++)
			{
				for (j = 0; j < P[m]; j++)
				{
					s2Mean[m][j] += s2[m][j] / nbrsample;
					int xx = 1;
					for (l = 0; l < r; l++)
					{
						Gammean[m][j][l] += (double)Gam[m][j][l] / nbrsample;
						Taumean[m][l][j] += Tau[m][l][j] / nbrsample;
						xx *= 1 - Gam[m][j][l];
					}
					if (xx == 0)
						Gamvs[m][j] += 1.0 / nbrsample;
				}

				for (l = 0; l < r; l++)
				{
					B0mean[m][l] += B0[m][l] / nbrsample;
					rhomean[m][l] += (double)rhoest[m][l] / nbrsample;
					for (k = 0; k < K[m]; k++)
					{
						Rmean[m][l][k] += (double)R[m][l][k] / nbrsample;
						Bmean[m][l][k] += B[m][l][k] / nbrsample;
					}
				}
			}

			if (t % ((N + 5) / 5) == 1)
			{
				printf("\n");
				printf("The number of iterations is  %d\n", t);
			}
		}
	}
	printf("\n");
	
	int sumP = 0;
	int sumK = 0;
	for (m = 0; m < Np; m++)
	{
		for (l = 0; l < r; l++)
		{
			CompoSelMean[m * r + l] = rhomean[m][l];
		}
		for (j = 0; j < P[m]; j++)
		{
			for (l = 0; l < r; l++)
			{
				VarSelMean[sumP * r + j * r + l] = Gammean[m][j][l];
			}
		}
		for (l = 0; l < r; l++)
		{
			for (k = 0; k < K[m]; k++)
			{
				GrpSelMean[sumK * r + l * K[m] + k] = Rmean[m][l][k];
				GrpEffectMean[sumK * r + l * K[m] + k] = Bmean[m][l][k];
			}
		}
		for (j = 0; j < P[m]; j++)
		{
			VarSelMeanGlobal[sumP + j] = Gamvs[m][j];
		}
		for (l = 0; l < r; l++)
		{
			IntGrpMean[m * r + l] = B0mean[m][l];
			AcceptB0Rate[m * r + l] = (TryB0[m][l] > 0) ? AcceptB0[m][l] / TryB0[m][l] : -1; // -1: never active
		}
		sumP += P[m];
		sumK += K[m];
	}

	/* Loading estimate for prediction using multiple models*/
	int u = 0;
	for (i = 0; i < n; i++)
	{
		for (l = 0; l < r; l++)
		{
			EstU[u] = meanU[i][l];
			u += 1;
		}
	}

	int countmodel = 0;
	int *modelidx = malloc(nbrsample * sizeof(double));
	_Bool **UniqModel = UniqueModel(nbrsample, dim, rhomodel, modelidx, &countmodel);
	free_bmatrix(rhomodel, 0, nbrsample - 1, 0, dim - 1);
	double *logpo = malloc(countmodel * sizeof(double));
	crossprod(n, r, r, meanU, meanU, UtU);
	for (m = 0; m < Np; m++)
	{
		if (IndVar[m] == 1)
		{
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < P[m]; j++)
				{
					X1[m][i][j] = X[m][i][j] - *InterceptMean;
				}
			}
		}
		crossprod(n, r, P[m], meanU, X1[m], UtX[m]);
		colsumsq(n, P[m], X1[m], xx[m]);
	}
	for (t = 0; t < countmodel; t++)
	{
		int rm = 0;
		//int sumg = 0;
		for (m = 0; m < Np; m++)
		{
			for (j = 0; j < P[m]; j++)
			{
				for (l = 0; l < r; l++)
				{
					Gam[m][j][l] = UniqModel[t][rm];
					rm++;
					//sumg += Gam[m][j][l];
				}
			}
		}
		for (m = 0; m < Np; m++)
		{
			if (IndVar[m] != 1)
			{
				for (l = 0; l < r; l++)
				{
					rhoest[m][l] = UniqModel[t][rm];
					rm++;
				}
			}
			else
			{
				for (l = 0; l < r; l++)
				{
					rhoest[m][l] = Gam[m][0][l];
				}
			}
		}

		logPostGam(&logpo[t], IndVar, Np, r, n, P, Tau, UtU, UtX, xx, s2Mean, rhoest, Gam, qv, q);
	}
	int *highmodelidx = malloc(countmodel * sizeof(int));
	sort(countmodel, logpo, highmodelidx);
	double maxlogpost = logpo[0];
	int maxmodel = maxmodel1[0];
	int nbrmax = MIN(maxmodel, countmodel);
	*nbrmodel1 = nbrmax;
	for (l = 0; l < nbrmax; l++)
	{
		logpo[l] = exp(logpo[l] - maxlogpost);
	}
	double sumpost = nbrmax * mean(nbrmax, logpo);

	for (l = 0; l < nbrmax; l++)
	{
		logpo[l] = logpo[l] / sumpost;
		postgam[l] = logpo[l];
	}
	int ll = 0;
	for (t = 0; t < nbrmax; t++)
	{
		int rm = 0;
		for (m = 0; m < Np; m++)
		{
			for (j = 0; j < P[m]; j++)
			{
				for (l = 0; l < r; l++)
				{
					Gam[m][j][l] = UniqModel[highmodelidx[t]][rm];
					rm++;
				}
			}
		}
		for (m = 0; m < Np; m++)
		{
			if (IndVar[m] != 1)
			{
				for (l = 0; l < r; l++)
				{
					rhoest[m][l] = UniqModel[highmodelidx[t]][rm];
					rm++;
				}
			}
			else
			{
				for (l = 0; l < r; l++)
				{
					rhoest[m][l] = Gam[m][0][l];
				}
			}
		}
		for (m = 0; m < Np; m++)
		{
			EstimateLoad(rr, r, n, P[m], rhoest[m], Taumean[m], A[m], UtU, UtX[m], s2Mean[m], Gam[m]);
			for (l = 0; l < r; l++)
			{
				for (j = 0; j < P[m]; j++)
				{
					EstLoadMod[ll] = A[m][l][j];
					ll++;
				}
			}
		}
	}
	free(modelidx);
	free(logpo);
	free(highmodelidx);
	free_bmatrix(UniqModel, 0, countmodel - 1, 0, dim - 1);

	double thres = 0.5;
	ll = 0;
	int ls = 0;
	for (m = 0; m < Np; m++)
	{
		for (l = 0; l < r; l++)
		{
			if (rhomean[m][l] >= thres)
				rhoest[m][l] = 1;
			else
				rhoest[m][l] = 0;
			for (j = 0; j < P[m]; j++)
			{
				if (rhoest[m][l] == 1)
				{
					if (Gammean[m][j][l] >= thres)
						Gam[m][j][l] = 1;
					else
						Gam[m][j][l] = 0;
				}
				else
				{
					Gam[m][j][l] = 0;
				}
			}
		}
		// X1: the response centred by InterceptMean, as for EstLoadMod; X1[m] == X[m] for the other data types
		crossprod(n, r, P[m], meanU, X1[m], UtX[m]);
		EstimateLoad(rr, r, n, P[m], rhoest[m], Taumean[m], A[m], UtU, UtX[m], s2Mean[m], Gam[m]);

		for (l = 0; l < r; l++)
		{
			for (j = 0; j < P[m]; j++)
			{
				EstLoad[ll] = A[m][l][j];
				ll += 1;
			}
		}
	}
	free_dmatrix(rhomean, 0, Np - 1, 0, r - 1);

	for (m = 0; m < Np; m++)
	{
		for (j = 0; j < P[m]; j++)
		{
			EstSig2[ls] = s2Mean[m][j];
			ls += 1;
		}
	}
	free_dmatrix(B0, 0, Np - 1, 0, r - 1);
	free_dmatrix(B0mean, 0, Np - 1, 0, r - 1);

	for (m = 0; m < Np; m++)
	{
		free_bmatrix(Gam[m], 0, P[m] - 1, 0, r - 1);
		free_dmatrix(Gammean[m], 0, P[m] - 1, 0, r - 1);
		free_bmatrix(R[m], 0, r - 1, 0, K[m] - 1);
		free_bmatrix(Path[m], 0, P[m] - 1, 0, K[m] - 1);
		free_dmatrix(Rmean[m], 0, r - 1, 0, K[m] - 1);
		free_dmatrix(Bmean[m], 0, r - 1, 0, K[m] - 1);
		free_dmatrix(B[m], 0, r - 1, 0, K[m] - 1);
		free_dmatrix(lambda2[m], 0, r - 1, 0, P[m] - 1);
		free_dmatrix(AcceptR[m], 0, r - 1, 0, K[m] - 1);
		free_dmatrix(Tau[m], 0, r - 1, 0, P[m] - 1);
		free_dmatrix(Taumean[m], 0, r - 1, 0, P[m] - 1);
		free_dmatrix(A[m], 0, r - 1, 0, P[m] - 1);
		free_dmatrix(X[m], 0, n - 1, 0, P[m] - 1);
		free_dmatrix(X1[m], 0, n - 1, 0, P[m] - 1);
		free(s2[m]);
		free(s2Mean[m]);
		free(rhoest[m]);
		free(Gamvs[m]); // free(GamvsK[m]);
		free(quadForm[m]);
		free(loggauss[m]);
		free(AcceptB0[m]);
		free(TryB0[m]);
		free(sdB0[m]);
		free(qv[m]);
		free(wg[m]);
		free_dmatrix(UtX[m], 0, r - 1, 0, P[m] - 1);
		free(xx[m]);
	}
	free_dmatrix(UtU, 0, r - 1, 0, r - 1);
	free(UtX);
	free(xx);
	free(Gam);
	free(Gamvs); // free(GamvsK);
	free(Tau);
	free(Taumean);
	free(A);
	free(X);
	free(X1);
	free(Gammean);
	free(R);
	free_dmatrix(U, 0, n - 1, 0, r - 1);
	free_dmatrix(meanU, 0, n - 1, 0, r - 1);
	free(s2);
	free(s2Mean);
	free(rhoest);
	free(q);
	free(qv);
	free(wg);
	free(quadForm);
	free(loggauss);
	free(B);
	free(lambda2);
	gsl_rng_free(rr);
	free(AcceptB0);
	free(TryB0);
	free(sdB0);
	free(Path);
	free(AcceptR);
	free(Rmean);
	free(Bmean);

	t1 = clock() - t1;
	double time_taken = ((double)t1) / CLOCKS_PER_SEC; // in seconds
	printf("\n\nTime taken in seconds is %f\n", time_taken);
	printf("\nTime taken in minutes is %f\n", time_taken / 60);
	gsl_set_error_handler(old_handler);
	
}

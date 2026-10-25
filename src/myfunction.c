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
#include <gsl/gsl_errno.h>
#include "header.h"
void SampleIntercept(gsl_rng *rr, int n, int r, double *intercept, double *sigma2, double sigma20, double **U, double **A, double **y)
{
   int i, l;
   double meany = 0;
   for (i = 0; i < n; i++)
   {
      double ua = 0;
      for (l = 0; l < r; l++)
      {
         ua += U[i][l] * A[l][0];
      }
      meany += y[i][0] - ua;
   }
   meany = meany / n;
   double invsig2 = n / sigma2[0] + 1 / sigma20;
   *intercept = (n / (invsig2 * sigma2[0])) * meany + sqrt(1 / invsig2) * gsl_ran_ugaussian(rr);
}
void logPostGam(double *logpo, int *IndVar, int Np, int r, int n, int *P, double ***Tau, double **UtU, double ***UtX, double **xx, double **s2, _Bool **rho, _Bool ***Gam, double **qv, double *q)
{
   double logpost = 0;
   int m, j, l;
   for (m = 0; m < Np; m++)
   {
      for (j = 0; j < P[m]; j++)
      {
         double logp;
         double quad;
         logGausQuadForm(j, r, n, Tau[m], UtU, UtX[m], xx[m][j], s2[m][j], &quad, Gam[m][j], &logp);
         logpost += logp;
         for (l = 0; l < r; l++)
         {
            /* For omics data, Gam[.][l] is fixed at 0 when rho[l] == 0, so its prior only
             * contributes when the component is selected. For the response, Gam is rho. */
            if ((IndVar[m] == 1) || ((IndVar[m] == 0) && (rho[m][l] == 1)))
               logpost += Gam[m][j][l] * log(qv[m][l]) + (1 - Gam[m][j][l]) * log(1 - qv[m][l]);
         }
      }
      // if (m<Np-1){
      if (IndVar[m] == 0)
      {
         for (l = 0; l < r; l++)
         {
            logpost += rho[m][l] * log(q[m]) + (1 - rho[m][l]) * log(1 - q[m]);
         }
      }
   }
   *logpo = logpost;
}

/* UtU = U'U (r x r), UtX = U'X (r x p) and xx = ||x_j||^2 are precomputed by the caller,
 * so each call costs O(nx1^3) instead of O(n nx1^2). */
void logGausQuadForm(int j, int r, int n, double **Tau, double **UtU, double **UtX, double xx, double s2, double *quadForm, _Bool *Gam, double *loggauss)
{
   int s, s1;
   int NZ1[r];
   int nx1 = 0;
   findc(r, Gam, 0, NZ1, &nx1);

   double result = 0;
   double quadF = xx / s2;
   if (nx1 > 0)
   {
      double Sigma1[nx1 * nx1];
      double y[nx1];
      double sqs2 = sqrt(s2);
      double sumT = 0;
      double logSqrtDet = 0;
      for (s = 0; s < nx1; s++)
      {
         for (s1 = 0; s1 <= s; s1++)
         {
            Sigma1[s * nx1 + s1] = Sigma1[s1 * nx1 + s] = UtU[NZ1[s]][NZ1[s1]];
         }
         Sigma1[s * nx1 + s] += (1.0 / Tau[NZ1[s]][j]);
         sumT += log(Tau[NZ1[s]][j]);
         y[s] = UtX[NZ1[s]][j] / sqs2;
      }
      if (safe_cholesky(Sigma1, nx1, "logGausQuadForm") != GSL_SUCCESS)
      { // zero likelihood: an MH move into this state is rejected; xx/s2 keeps sigma2() finite
         *quadForm = quadF;
         *loggauss = -INFINITY;
         return;
      }
      gsl_matrix_view m1 = gsl_matrix_view_array(Sigma1, nx1, nx1);
      /* y'(U'U+T^-1)^-1 y = ||L^-1 y||^2 */
      gsl_vector_view yv = gsl_vector_view_array(y, nx1);
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &m1.matrix, &yv.vector);
      double qy = 0;
      gsl_blas_ddot(&yv.vector, &yv.vector, &qy);
      for (s = 0; s < nx1; s++)
      {
         logSqrtDet += log(Sigma1[s * nx1 + s]);
      }
      quadF -= qy;
      result = -0.5 * quadF - logSqrtDet - 0.5 * n * log(2.0 * M_PI) - 0.5 * sumT - (n / 2.0) * log(s2);
   }
   else
   {
      result = -(n / 2.0) * log(s2) - 0.5 * n * log(2.0 * M_PI) - 0.5 * quadF;
   }
   *quadForm = quadF;
   *loggauss = result;
}


double logpost(int r, int n, int p, _Bool *rho, double **Tau, double **UtU, double **UtX, double *xx, double q, double *s2, double *quadForm, _Bool **Gam, double *loggauss)
{
   double logpostrho = 0;
   int l, j;
   for (j = 0; j < p; j++)
   {
      logGausQuadForm(j, r, n, Tau, UtU, UtX, xx[j], s2[j], &quadForm[j], Gam[j], &loggauss[j]);
      logpostrho += loggauss[j];
   }
   for (l = 0; l < r; l++)
   {
      logpostrho += rho[l] * log(q) + (1 - rho[l]) * log(1 - q);
   }
   return logpostrho;
}

void rho1(gsl_rng *rr, int r, int n, int p, _Bool *rho, double **Tau, double **UtU, double **UtX, double *xx, double q, double *s2, double *quadForm, _Bool **Gam, double *loggauss)
{
   int l, j;
   double logpostold = 0;
   double logprior = 0;
   for (l = 0; l < r; l++)
   {
      logprior += rho[l] * log(q) + (1 - rho[l]) * log(1 - q);
   }
   for (j = 0; j < p; j++)
   {
      logpostold += loggauss[j];
   }
   logpostold += logprior;
   double *quadForm1 = malloc(p * sizeof(double));
   double *loggauss1 = malloc(p * sizeof(double));
   double phi = 0.5;
   _Bool rhonew[r];
   proposal(r, rho, rhonew, phi, rr);
   double logpostnew = logpost(r, n, p, rhonew, Tau, UtU, UtX, xx, q, s2, quadForm1, Gam, loggauss1);
   double u = gsl_ran_flat(rr, 0, 1);
   if (log(u) < logpostnew - logpostold)
   {
      for (l = 0; l < r; l++)
      {
         rho[l] = rhonew[l];
      }
      for (j = 0; j < p; j++)
      {
         quadForm[j] = quadForm1[j];
         loggauss[j] = loggauss1[j];
      }
   }
}

void SamplerhoGamma(gsl_rng *rr, int r, int n, int IndVar, int p, _Bool *rho, double **Tau, double **UtU, double **UtX, double *xx, double *q1, double q2, double *s2, double *quadForm, _Bool **Gam, double *loggauss)
{
   int l, j;
   _Bool *rhonew = malloc(r * sizeof(_Bool));
   _Bool (*Gamnew)[r] = malloc(p * sizeof(*Gamnew)); /* p x r, on the heap: p can be large */
   for (l = 0; l < r; l++)
   {
      for (j = 0; j < p; j++)
      {
         Gamnew[j][l] = Gam[j][l];
      }
      rhonew[l] = rho[l];
   }

   if (IndVar == 1)
   {
      for (l = 0; l < r; l++)
      {
         double quad1 = 0;
         double quad2 = 0;
         double loggaussold = 0;
         double loggaussnew = 0;
         Gam[0][l] = 0;
         logGausQuadForm(0, r, n, Tau, UtU, UtX, xx[0], s2[0], &quad1, Gam[0], &loggaussold);
         Gamnew[0][l] = 1;
         logGausQuadForm(0, r, n, Tau, UtU, UtX, xx[0], s2[0], &quad2, Gamnew[0], &loggaussnew);
         double rat = loggaussnew + log(q2) - loggaussold - log(1 - q2);
         double uni = gsl_ran_flat(rr, 0, 1);
         // printf("RAT=%.2lf",rat);
         if (log(uni / (1 - uni)) < rat)
         {
            Gamnew[0][l] = Gam[0][l] = rho[l] = 1;
            loggauss[0] = loggaussnew;
            quadForm[0] = quad2;
         }
         else
         {
            Gamnew[0][l] = Gam[0][l] = rho[l] = 0;
            loggauss[0] = loggaussold;
            quadForm[0] = quad1;
         }
      }
   }
   else
   {

      for (l = 0; l < r; l++)
      {
         double logq = 0;
         if (rho[l] == 0)
         {
            rhonew[l] = 1;
            double logpostnew = 0;
            double logpostold = 0;
            double *quadForm1 = malloc(p * sizeof(double));
            double *loggausnew1 = malloc(p * sizeof(double));
            double quad1 = 0;
            double quad2 = 0;
            double loggaussold = 0;
            double loggaussnew = 0;
            for (j = 0; j < p; j++)
            {
               double logqj = 0;
               double logmqj = 0;
               logpostold += loggauss[j];
               Gam[j][l] = 0;
               logGausQuadForm(j, r, n, Tau, UtU, UtX, xx[j], s2[j], &quad1, Gam[j], &loggaussold);
               Gamnew[j][l] = 1;
               logGausQuadForm(j, r, n, Tau, UtU, UtX, xx[j], s2[j], &quad2, Gamnew[j], &loggaussnew);
               double rat = -loggaussold + loggaussnew + log(q1[l]) - log(1 - q1[l]);
               // printf("%lf ",rat);
               double x1 = loggaussnew + log(q1[l]);
               double x2 = loggaussold + log(1 - q1[l]);
               double maxq = MAX(x1, x2);
               double den = maxq + log(exp(x1 - maxq) + exp(x2 - maxq));
               logqj = x1 - den;
               logmqj = x2 - den;

               double uni = gsl_ran_flat(rr, 0, 1);
               if ((log(uni / (1 - uni)) < rat))
               {
                  Gam[j][l] = 1;
                  Gamnew[j][l] = 1;
                  logpostnew += loggaussnew + log(q1[l]);
                  loggausnew1[j] = loggaussnew;
                  logq += logqj;
                  quadForm1[j] = quad2;
               }
               else
               {
                  Gam[j][l] = 0;
                  Gamnew[j][l] = 0;
                  logq += logmqj;
                  logpostnew += loggaussold + log(1 - q1[l]);
                  quadForm1[j] = quad1;
                  loggausnew1[j] = loggaussold;
               }
            }
            logpostnew += log(q2);
            logpostold += log(1 - q2);
            double un = gsl_ran_flat(rr, 0, 1);
            double rat1 = logpostnew - logpostold - logq;
            if (log(un) < rat1)
            {
               rho[l] = 1;
               for (j = 0; j < p; j++)
               {
                  quadForm[j] = quadForm1[j];
                  loggauss[j] = loggausnew1[j];
               }
            }
            else
            {
               rho[l] = 0;
               for (j = 0; j < p; j++)
                  Gam[j][l] = Gamnew[j][l] = 0;
            }
            rhonew[l] = rho[l];
            free(quadForm1);
            free(loggausnew1);
         }
         else
         {
            rhonew[l] = 0;
            double logpostnew = 0;
            double logpostold = 0;
            double *quadForm1 = malloc(p * sizeof(double));
            double *loggausnew1 = malloc(p * sizeof(double));
            double quad1;
            double loggaussnew = 0;
            double loggaussn = 0;
            double logpq = 0;
            for (j = 0; j < p; j++)
            {
               logpostold += loggauss[j];
               Gamnew[j][l] = 1;
               logGausQuadForm(j, r, n, Tau, UtU, UtX, xx[j], s2[j], &quad1, Gamnew[j], &loggaussn);
               Gamnew[j][l] = 0;
               logGausQuadForm(j, r, n, Tau, UtU, UtX, xx[j], s2[j], &quad1, Gamnew[j], &loggaussnew);
               if (p != 1)
               {
                  double x1 = loggaussn + log(q1[l]);
                  double x2 = loggaussnew + log(1 - q1[l]);
                  double maxq = MAX(x1, x2);
                  double den = maxq + log(exp(x1 - maxq) + exp(x2 - maxq));
                  double logqj = x1 - den;
                  double logmqj = x2 - den;
                  logq += Gam[j][l] * log(q1[l]) + (1 - Gam[j][l]) * log(1 - q1[l]);
                  logpq += Gam[j][l] * logqj + (1 - Gam[j][l]) * logmqj;
               }
               logpostnew += loggaussnew;
               loggausnew1[j] = loggaussnew;
               quadForm1[j] = quad1;
               Gamnew[j][l] = Gam[j][l];
            }
            logpostnew += log(1 - q2) + logpq;
            logpostold += log(q2) + logq;
            double un = gsl_ran_flat(rr, 0, 1);
            double rat1 = logpostnew - logpostold;
            if (log(un) < rat1)
            {
               rho[l] = 0;
               for (j = 0; j < p; j++)
               {
                  Gamnew[j][l] = Gam[j][l] = 0;
                  quadForm[j] = quadForm1[j];
                  loggauss[j] = loggausnew1[j];
               }
            }
            else
            {
               for (j = 0; j < p; j++)
               {
                  Gamnew[j][l] = Gam[j][l];
               }
            }
            free(quadForm1);
            free(loggausnew1);
            rhonew[l] = rho[l];
         }
         /* Within move*/
         if (rho[l] == 1)
         {
           // double logpostold = 0;
            double *quadForm1 = malloc(p * sizeof(double));
            double *loggausnew1 = malloc(p * sizeof(double));
            double quad1, quad2;
            double loggaussold, loggaussnew;
            for (j = 0; j < p; j++)
            {
            //   logpostold += loggauss[j];
               logGausQuadForm(j, r, n, Tau, UtU, UtX, xx[j], s2[j], &quad1, Gam[j], &loggaussold);
               Gamnew[j][l] = 1 - Gam[j][l];
               logGausQuadForm(j, r, n, Tau, UtU, UtX, xx[j], s2[j], &quad2, Gamnew[j], &loggaussnew);
               double rat = -loggaussold + loggaussnew;
               double uni = gsl_ran_flat(rr, 0, 1);
               /* After this block quad1/loggaussold hold the Gam[j][l]=1 values and
                * quad2/loggaussnew the Gam[j][l]=0 values */
               if (Gam[j][l] == 0)
               {
                  double tq = quad1;
                  quad1 = quad2;
                  quad2 = tq;
                  double tl = loggaussold;
                  loggaussold = loggaussnew;
                  loggaussnew = tl;
                  rat += log(q1[l]) - log(1 - q1[l]);
               }
               else
               {
                  rat = -rat + log(q1[l]) - log(1 - q1[l]);
               }
               if (log(uni / (1 - uni)) < rat)
               {
                  Gam[j][l] = 1;
                  quadForm[j] = quad1;
                  loggauss[j] = loggaussold;
               }
               else
               {
                  Gam[j][l] = 0;
                  quadForm[j] = quad2;
                  loggauss[j] = loggaussnew;
               }
               Gamnew[j][l] = Gam[j][l];
            }
            free(quadForm1);
            free(loggausnew1);
         }
      }
   }
   free(rhonew);
   free(Gamnew);
}

void sigma2(gsl_rng *rr, int p, int n, double a0, double b0, double *quadForm, double *s2)
{
   int j;
   for (j = 0; j < p; j++)
   {
      double inv = 1 / (b0 + 0.5 * s2[j] * quadForm[j]);
      s2[j] = 1 / gsl_ran_gamma(rr, a0 + n / 2.0, inv);
   }
}

void EffectZero(int l, gsl_rng *rr, int K, int p, _Bool rho, _Bool *R, double *B, double *B0, _Bool **Path, double alphab0, double betab0, double alpha, double *lambda2, double *AcceptB0, _Bool **Gam, double sdB0)
{
   /* Random walk on log(B0): B0new = B0 exp(sdB0 * N(0,1)). The Jacobian of the log
    * transform adds log(B0new/B0) to the acceptance ratio. */
   int j, k1;
   double B0new = *B0 * exp(sdB0 * gsl_ran_ugaussian(rr));
   double lograt = 0;
   double kl = 0;
   for (j = 0; j < p; j++)
   {
      double kb = 0;
      double kb1 = 0;
      if (Gam[j][l] == 1)
      {
         for (k1 = 0; k1 < K; k1++)
         {
            kb += Path[j][k1] * B[k1];
         }
         kl += lambda2[j];
         kb1 = kb + B0new;
         kb = kb + *B0;
         lograt += alpha * (log(kb1) - log(kb));
      }
   }
   if (rho == 1)
   {
      /* Gamma(alphab0, betab0) prior, lambda2 likelihood and Jacobian */
      lograt += -(betab0 + kl) * (B0new - *B0) + alphab0 * (log(B0new) - log(*B0));
      double uni = gsl_ran_flat(rr, 0, 1);
      if (log(uni) < lograt)
      {
         *B0 = B0new;
         *AcceptB0 += 1;
      }
   }
   else
   {
      *B0 = gsl_ran_gamma(rr, alphab0, 1 / betab0);
   }
}

void TauLambda(int l, gsl_rng *rr, int K, int p, double *A, double *B, double B0, double *Tau, _Bool **Path, double alpha, double *lambda2, double *s2, _Bool **Gam)
{
   /* K is the number of pathways
    *  * R is the binary vector for pathway effect
    *   * Path is the grouping information; it is a binary matrix
    *    * B is the pathway effect
    *     * alphab and betab are priors for B effects
    *      * l is the component index
    *       */
   int j, k1;
   for (j = 0; j < p; j++)
   {
      // if (fabs(A[j])>0){
      if ((Gam[j][l] == 1) && (fabs(A[j]) > 0))
      {
         double mu = sqrt(2 * lambda2[j] * s2[j] / pow(A[j], 2.0));
         Tau[j] = 1.0 / inverseGaussian(rr, mu, 2 * lambda2[j]);
         /* Fall back to E[Tau] when the draw is NaN, infinite or not positive */
         if (!isfinite(Tau[j]) || (Tau[j] <= 0))
         {
            Tau[j] = 1 / mu + 1 / (2 * lambda2[j]);
            // printf("J=%d, MU=%lf, lamb=%lf, S2=%lf, Beta=%lf B0=%lf\n", j,mu,lambda2[j],s2[j],A[j],B0);
         }
         double kb = 0;
         for (k1 = 0; k1 < K; k1++)
         {
            kb += Path[j][k1] * B[k1];
         }
         lambda2[j] = gsl_ran_gamma(rr, alpha + 1, 1 / (B0 + kb + Tau[j]));
      }
      else
      {
         Tau[j] = gsl_ran_exponential(rr, 1.0 / lambda2[j]);
         // printf("Tau=%lf ",Tau[j]);
         // lambda2[j]=gsl_ran_gamma (rr, alpha+1, 1/(Tau[j]));
         lambda2[j] = gsl_ran_gamma(rr, alpha, 1 / B0);
      }
   }
}

/* Path[j]'B for variable j, summed in the same order as the original per-j loop */
static double groupsum(int K, _Bool *Pathj, double *B)
{
   int k1;
   double kb = 0;
   for (k1 = 0; k1 < K; k1++)
   {
      kb += Pathj[k1] * B[k1];
   }
   return kb;
}

void GroupEffect(int l, gsl_rng *rr, _Bool rho, int K, int p, _Bool *R, double *A, double *B, double B0, double *Tau, _Bool **Path, double alphab, double betab, double w, double alpha, double *lambda2, double *s2, double *AcceptR, _Bool **Gam)
{
   /* K is the number of pathways
    * R is the binary vector for pathway effect
    * Path is the grouping information; it is a binary matrix
    * B is the pathway effect
    * alphab and betab are priors for B effects
    * l is the component index
    */
   int j, k, a;

   if (rho == 0)
   {
      for (k = 0; k < K; k++)
      {
         B[k] = 0;
         R[k] = 0;
      }
   }
   else
   {

      double alphas = 2;
      double betas = 2; // proposal parameters
      double lgab = gsl_sf_lngamma(alphab);
      double lgas = gsl_sf_lngamma(alphas);
      double lbs = log(betas);
      double lbb = log(betab);
      double l1w = log(1 - w);
      double lw = log(w);
      /* Only variables with Gam[j][l]==1 that belong to group k change the acceptance ratio:
       * for the others the original terms were exactly zero. kb[a] caches Path[j]'B for the
       * a-th active variable and is refreshed whenever B[k] changes. */
      int *act = malloc(p * sizeof(int));
      double *kb = malloc(p * sizeof(double));
      int na = 0;
      for (j = 0; j < p; j++)
      {
         if (Gam[j][l] == 1)
         {
            act[na] = j;
            kb[na] = groupsum(K, Path[j], B);
            na++;
         }
      }
      for (k = 0; k < K; k++)
      {
         _Bool Rknew = 0;
         double Bknew = 0;
         _Bool changed = 0;
         /* Between model move*/
         if (R[k] == 1)
         {
            Rknew = 0;
            Bknew = 0;
            double lograt = 0;
            double kl = 0;
            for (a = 0; a < na; a++)
            {
               j = act[a];
               if (Path[j][k])
               {
                  double kb1 = kb[a] - B[k];
                  lograt += alpha * (log(B0 + kb1) - log(B0 + kb[a]));
                  kl += lambda2[j];
               }
            }
            lograt += (betab + kl - betas) * B[k] + (alphas - alphab) * log(B[k]) + lgab - lgas + alphas * lbs - alphab * lbb + l1w - lw;
            double uni = gsl_ran_flat(rr, 0, 1);
            if (log(uni) < lograt)
            {
               B[k] = Bknew;
               R[k] = Rknew;
               AcceptR[k] += 1;
               changed = 1;
            }
         }
         else if (R[k] == 0)
         {
            Rknew = 1;
            Bknew = gsl_ran_gamma(rr, alphas, 1 / betas);
            double lograt = 0;
            double kl = 0;
            for (a = 0; a < na; a++)
            {
               j = act[a];
               if (Path[j][k])
               {
                  double kb1 = kb[a] - B[k];
                  double kb2 = kb1 + Bknew;
                  lograt += alpha * (log(B0 + kb2) - log(B0 + kb1));
                  kl += lambda2[j];
               }
            }
            lograt += (betas - betab - kl) * Bknew + (-alphas + alphab) * log(Bknew) - lgab + lgas - alphas * lbs + alphab * lbb - l1w + lw;
            double uni = gsl_ran_flat(rr, 0, 1);
            if (log(uni) < lograt)
            {
               B[k] = Bknew;
               R[k] = Rknew;
               AcceptR[k] += 1;
               changed = 1;
            }
         }
         if (changed)
         {
            for (a = 0; a < na; a++)
            {
               if (Path[act[a]][k])
                  kb[a] = groupsum(K, Path[act[a]], B);
            }
            changed = 0;
         }
         /* Within model move*/

         if (R[k] == 1)
         {
            Bknew = gsl_ran_gamma(rr, alphas, 1 / betas);
            double lograt = 0;
            double kl = 0;
            for (a = 0; a < na; a++)
            {
               j = act[a];
               if (Path[j][k])
               {
                  double kb1 = kb[a] - B[k] + Bknew;
                  lograt += alpha * (log(B0 + kb1) - log(B0 + kb[a]));
                  kl += lambda2[j];
               }
            }
            lograt += (betas - betab - kl) * (Bknew - B[k]) + (-alphas + alphab) * (log(Bknew) - log(B[k]));
            double uni = gsl_ran_flat(rr, 0, 1);
            if (log(uni) < lograt)
            {
               B[k] = Bknew;
               for (a = 0; a < na; a++)
               {
                  if (Path[act[a]][k])
                     kb[a] = groupsum(K, Path[act[a]], B);
               }
            }
         }
      } // close for k
      free(act);
      free(kb);
   } // End for else if (rho==0)
} // close for the function

void LoadAOther(gsl_rng *rr, int r, int n, int p, _Bool *rho, double **Tau, double **A, double **UtU, double **UtX, double *s2, _Bool **Gam)
{
   int s, s1, j;
   double Aold[r];
   for (j = 0; j < p; j++)
   {
      for (s = 0; s < r; s++)
      {
         Aold[s] = A[s][j];
         A[s][j] = 0;
      }
      int NZ1[r];
      int nx1 = 0;
      findc(r, Gam[j], 0, NZ1, &nx1);
      if (nx1 > 0)
      {
         double *SigmaInv = malloc(nx1 * nx1 * sizeof(double));
         for (s = 0; s < nx1; s++)
         {
            for (s1 = 0; s1 <= s; s1++)
            {
               SigmaInv[s * nx1 + s1] = SigmaInv[s1 * nx1 + s] = UtU[NZ1[s]][NZ1[s1]] / s2[j];
            }
            SigmaInv[s * nx1 + s] += 1 / (Tau[NZ1[s]][j] * s2[j]);
         }
         if (safe_cholesky(SigmaInv, nx1, "sampling of the loadings A") != GSL_SUCCESS)
         { // keep the previous loadings of variable j, zero where Gam[j][s] == 0
            for (s = 0; s < nx1; s++)
               A[NZ1[s]][j] = Aold[NZ1[s]];
            free(SigmaInv);
            continue;
         }
         gsl_matrix_view m = gsl_matrix_view_array(SigmaInv, nx1, nx1);
         gsl_vector *mu = gsl_vector_alloc(nx1);
         double xy[nx1];
         for (s = 0; s < nx1; s++)
         {
            xy[s] = UtX[NZ1[s]][j] / s2[j];
         }
         gsl_vector_view b = gsl_vector_view_array(xy, nx1);
         gsl_linalg_cholesky_solve(&m.matrix, &b.vector, mu);
         gsl_vector *result = gsl_vector_alloc(nx1);
         multivariate_gaussian(rr, mu, &m.matrix, result);
         for (s = 0; s < nx1; s++)
         {
            A[NZ1[s]][j] = gsl_vector_get(result, s);
         }
         free(SigmaInv);
         gsl_vector_free(mu);
         gsl_vector_free(result);
      }
   }
}

void EstimateLoad(gsl_rng *rr, int r, int n, int p, _Bool *rho, double **Tau, double **A, double **UtU, double **UtX, double *s2, _Bool **Gam)
{
   int s, s1, j;
   for (j = 0; j < p; j++)
   {
      for (s = 0; s < r; s++)
      {
         A[s][j] = 0;
      }
      int NZ1[r];
      int nx1 = 0;
      findc(r, Gam[j], 0, NZ1, &nx1);
      if (nx1 > 0)
      {
         double *SigmaInv = malloc(nx1 * nx1 * sizeof(double));
         for (s = 0; s < nx1; s++)
         {
            for (s1 = 0; s1 <= s; s1++)
            {
               SigmaInv[s * nx1 + s1] = SigmaInv[s1 * nx1 + s] = UtU[NZ1[s]][NZ1[s1]] / s2[j];
            }
            SigmaInv[s * nx1 + s] += 1 / (Tau[NZ1[s]][j] * s2[j]);
         }
         if (safe_cholesky(SigmaInv, nx1, "estimation of the loadings A") != GSL_SUCCESS)
         { // loadings of variable j are left at 0
            free(SigmaInv);
            continue;
         }
         gsl_matrix_view m = gsl_matrix_view_array(SigmaInv, nx1, nx1);
         gsl_vector *mu = gsl_vector_alloc(nx1);
         double xy[nx1];
         for (s = 0; s < nx1; s++)
         {
            xy[s] = UtX[NZ1[s]][j] / s2[j];
         }
         gsl_vector_view b = gsl_vector_view_array(xy, nx1);
         gsl_linalg_cholesky_solve(&m.matrix, &b.vector, mu);
         for (s = 0; s < nx1; s++)
         {
            A[NZ1[s]][j] = gsl_vector_get(mu, s);
         }
         free(SigmaInv);
         gsl_vector_free(mu);
      }
   }
}

void proposal(int n, _Bool R[n], _Bool R1[n], float phi, gsl_rng *r)
{
   int i;
   for (i = 0; i < n; i++)
      R1[i] = R[i];
   int n1 = 0; // number of 1
   int n0 = 0; // number of zeros
   int v1[n];
   findc(n, R, 1, v1, &n0); // find indices different from 1 i.e ==0;
   int v2[n - n0];
   findc(n, R, 0, v2, &n1); // find indices different of zeros i.e ==1
   double u = gsl_ran_flat(r, 0, 1);
   if ((u < phi) || (n0 == 0) || (n1 == 0))
   {
      int l = gsl_rng_uniform_int(r, n);

      R1[l] = 1 - R[l];
   }
   else
   {
      int l1 = gsl_rng_uniform_int(r, n0);
      int l2 = gsl_rng_uniform_int(r, n1);

      R1[v1[l1]] = R[v2[l2]];
      R1[v2[l2]] = R[v1[l1]];
   }
}

void SampleUU(gsl_rng *rr, int r, int n, int Np, int *P, double ***A, _Bool ***Gam, double **U, double ***X, double **s2)
{
   int s, s1, j, i, m, a1;

   /* A[m][s][j] is zero exactly when Gam[m][j][s] == 0 (see LoadAOther), and zero loadings
    * add exactly 0 to SigmaInv and Ax. So for each platform m and component s, only the
    * markers with Gam[m][j][s] == 1 are visited, in the original order. */
   int ***act = malloc(Np * sizeof(int **));
   int **na = malloc(Np * sizeof(int *));
   for (m = 0; m < Np; m++)
   {
      act[m] = malloc(r * sizeof(int *));
      na[m] = malloc(r * sizeof(int));
      for (s = 0; s < r; s++)
      {
         act[m][s] = malloc(P[m] * sizeof(int));
         na[m][s] = 0;
         for (j = 0; j < P[m]; j++)
         {
            if (Gam[m][j][s])
               act[m][s][na[m][s]++] = j;
         }
      }
   }

   double *SigmaInv = malloc(r * r * sizeof(double));
   for (s = 0; s < r; s++)
   {
      for (s1 = 0; s1 <= s; s1++)
      {
         double a = 0;
         for (m = 0; m < Np; m++)
         {
            for (a1 = 0; a1 < na[m][s]; a1++)
            {
               j = act[m][s][a1];
               if (Gam[m][j][s1])
                  a += A[m][s][j] * A[m][s1][j] / s2[m][j];
            }
         }
         SigmaInv[s * r + s1] = SigmaInv[s1 * r + s] = a;
      }
      SigmaInv[s * r + s] += 1;
   }

   /* If the factorisation fails, U keeps its current value for this iteration */
   int ok = (safe_cholesky(SigmaInv, r, "sampling of U") == GSL_SUCCESS);
   gsl_matrix_view mat = gsl_matrix_view_array(SigmaInv, r, r);
   gsl_vector *mu = gsl_vector_alloc(r);
   gsl_vector *result = gsl_vector_alloc(r);
   double Ax[r];
   for (i = 0; ok && (i < n); i++)
   {
      for (s = 0; s < r; s++)
      {
         Ax[s] = 0;
         for (m = 0; m < Np; m++)
         {
            double *Xi = X[m][i];
            double *As = A[m][s];
            double *s2m = s2[m];
            int *acts = act[m][s];
            for (a1 = 0; a1 < na[m][s]; a1++)
            {
               j = acts[a1];
               Ax[s] += Xi[j] * As[j] / s2m[j];
            }
         }
      }
      gsl_vector_view b = gsl_vector_view_array(Ax, r);
      gsl_linalg_cholesky_solve(&mat.matrix, &b.vector, mu);
      multivariate_gaussian(rr, mu, &mat.matrix, result);
      for (s = 0; s < r; s++)
      {
         U[i][s] = gsl_vector_get(result, s);
      }
   }
   gsl_vector_free(mu);
   gsl_vector_free(result);
   free(SigmaInv);
   for (m = 0; m < Np; m++)
   {
      for (s = 0; s < r; s++)
         free(act[m][s]);
      free(act[m]);
      free(na[m]);
   }
   free(act);
   free(na);
}

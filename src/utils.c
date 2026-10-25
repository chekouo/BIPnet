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
#include <gsl/gsl_errno.h>
#include "header.h"

double mean(int n, double *x)
{
  int i;
  double me = 0;
  for (i = 0; i < n; i++)
    me += x[i];
  return me / n;
}
/* UtX = U'X, where U is n x r and X is n x p; UtX is r x p */
void crossprod(int n, int r, int p, double **U, double **X, double **UtX)
{
  int i, s, j;
  for (s = 0; s < r; s++)
    for (j = 0; j < p; j++)
      UtX[s][j] = 0;
  for (i = 0; i < n; i++)
    for (s = 0; s < r; s++)
    {
      double u = U[i][s];
      for (j = 0; j < p; j++)
        UtX[s][j] += u * X[i][j];
    }
}
/* xx[j] = ||x_j||^2 for each column of the n x p matrix X */
void colsumsq(int n, int p, double **X, double *xx)
{
  int i, j;
  for (j = 0; j < p; j++)
    xx[j] = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < p; j++)
      xx[j] += X[i][j] * X[i][j];
}
_Bool compar(int n, _Bool *u, _Bool *v)
{
  int i;
  for (i = 0; i < n; i++)
  {
    if (u[i] != v[i])
      return 0;
  }
  return 1;
}

_Bool **UniqueModel(int nbrsample, int p, _Bool **rhosample, int *modelidx, int *countd1)
{

  int i;
  int j = 0;
  int countd = 0;
  modelidx[countd] = 0;
  for (i = 0; i < nbrsample; i++)
  {
    for (j = 0; j < countd; j++)
    {
      _Bool resultX = compar(p, rhosample[i], rhosample[modelidx[j]]);
      if (resultX == 1)
        break;
    }
    if (j == countd)
    {
      modelidx[countd] = i;
      countd++;
    }
  }
  *countd1 = countd;

  _Bool **UniqModel = malloc(countd * sizeof(_Bool *));
  for (i = 0; i < countd; i++)
  {
    UniqModel[i] = malloc(p * sizeof(_Bool));
    for (j = 0; j < p; j++)
    {
      UniqModel[i][j] = rhosample[modelidx[i]][j];
    }
  }
  Rprintf("\n\n\nNumber of models is= %d", countd);
  Rprintf("\n");
  return UniqModel;
}

int multivariate_gaussian(const gsl_rng *r,
                          const gsl_vector *mu,
                          const gsl_matrix *L,
                          gsl_vector *result)
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
    // gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, L, result);
    gsl_vector_add(result, mu);

    return GSL_SUCCESS;
  }
}

double inverseGaussian(gsl_rng *r, double mu, double lambda)
{
  double v = gsl_ran_gaussian(r, 1); // sample from a normal distribution with a mean of 0 and 1 standard deviation
  double y = v * v;
  /* x = mu + mu^2 y/(2 lambda) - (mu/(2 lambda)) sqrt(4 mu lambda y + mu^2 y^2), written as
   * mu^2 / (mu + mu^2 y/(2 lambda) + ...) to avoid cancellation when mu is large */
  double b = (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu * mu * y * y);
  double x = mu * mu / (mu + (mu * mu * y) / (2 * lambda) + b);
  double test = gsl_ran_flat(r, 0, 1); // sample from a uniform distribution between 0 and 1
  if (test <= (mu) / (mu + x))
    return x;
  else
    return (mu * mu) / x;
}

void sort(int n, double *x, int *idx)
{
  int i, j;
  double a;
  int id;
  for (i = 0; i < n; i++)
    idx[i] = i;
  for (i = 0; i < n; ++i)
  {
    for (j = i + 1; j < n; ++j)
    {
      if (x[i] <= x[j])
      {
        a = x[i];
        id = idx[i];
        idx[i] = idx[j];
        x[i] = x[j];
        idx[j] = id;
        x[j] = a;
      }
    }
  }
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  double **m;

  m = (double **)malloc((unsigned)(nrh - nrl + 1) * sizeof(double *));
  if (!m)
    nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (double *)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
    if (!m[i])
      nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void nrerror(char error_text[])
{
  Rf_error("BIPnet: %s", error_text); // returns to R with an error instead of ending the R session
}
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free((char *)(m[i] + ncl));

  free((char *)(m + nrl));
}
void findc(int n, _Bool R[n], int a, int *IDX, int *nx)
{
  int ii_data[n];
  int idx = 0;
  int ii = 0;
  _Bool exitg2 = 0;
  _Bool guard2 = 0;
  while ((exitg2 == 0) && (ii < n))
  {
    guard2 = 0;
    if (R[ii] != a)
    {
      idx++;
      ii_data[idx - 1] = ii;
      if (idx >= n)
      {
        exitg2 = 1;
      }
      else
      {
        guard2 = 1;
      }
    }
    else
    {
      guard2 = 1;
    }

    if (guard2 == 1)
    {
      ii++;
    }
  }

  int loop_ub = idx;
  for (idx = 0; idx < loop_ub; idx++)
  {
    IDX[idx] = ii_data[idx];
  }
  *nx = loop_ub;
}

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  _Bool **m;

  m = (_Bool **)malloc((nrh - nrl + 1) * sizeof(_Bool *));
  if (!m)
    nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (_Bool *)malloc((nch - ncl + 1) * sizeof(_Bool));
    if (!m[i])
      nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free((_Bool *)(m[i] + ncl));

  free((_Bool *)(m + nrl));
}

/* In-place Cholesky factorisation of the d x d symmetric matrix M (row-major) that never aborts
 * (mainfunction switches the GSL error handler off). If M is not numerically positive definite,
 * a ridge of 1e-10, 1e-8, 1e-6 and then 1e-4 times its mean diagonal is added and the
 * factorisation retried. When the first attempt succeeds, the result is exactly that of
 * gsl_linalg_cholesky_decomp. Returns GSL_SUCCESS, or GSL_EDOM when M has non-finite entries or
 * stays indefinite; a warning naming `where` is printed (the first 10 times only). */
int safe_cholesky(double *M, int d, const char *where)
{
  static int nwarn = 0;
  double M0[d * d];
  int i, k, status, nonfinite = 0;
  double meandiag = 0;
  for (i = 0; i < d * d; i++)
  {
    if (!isfinite(M[i]))
    {
      nonfinite = 1;
      goto fail;
    }
    M0[i] = M[i];
  }
  gsl_matrix_view mv = gsl_matrix_view_array(M, d, d);
  status = gsl_linalg_cholesky_decomp(&mv.matrix);
  if (status == GSL_SUCCESS)
    return GSL_SUCCESS;
  for (i = 0; i < d; i++)
    meandiag += fabs(M0[i * d + i]) / d;
  if (meandiag == 0)
    meandiag = 1;
  double ridge[4] = {1e-10, 1e-8, 1e-6, 1e-4};
  for (k = 0; k < 4; k++)
  {
    for (i = 0; i < d * d; i++)
      M[i] = M0[i];
    for (i = 0; i < d; i++)
      M[i * d + i] += ridge[k] * meandiag;
    status = gsl_linalg_cholesky_decomp(&mv.matrix);
    if (status == GSL_SUCCESS)
    {
      if (nwarn < 10)
        Rprintf("\nWarning (%s): matrix not positive definite, a ridge of %g x mean diagonal was added.\n", where, ridge[k]);
      nwarn++;
      return GSL_SUCCESS;
    }
  }
fail:
  if (nwarn < 10)
    Rprintf("\nWarning (%s): Cholesky factorisation failed (%s); this update is skipped.\n", where,
           nonfinite ? "non-finite entries" : "not positive definite even after adding a ridge");
  if (nwarn == 9)
    Rprintf("Further warnings of this kind are not printed.\n");
  nwarn++;
  return GSL_EDOM;
}

/* Assignment maximising sum_l w[l][perm[l]] over the permutations of 0..r-1 (Hungarian
 * algorithm, O(r^3)). Row l of w is a reference component, column k a current component. */
void max_assignment(int r, double **w, int *perm)
{
  double u[r + 1], v[r + 1], minv[r + 1];
  int p[r + 1], way[r + 1], i, j;
  _Bool used[r + 1];
  for (j = 0; j <= r; j++)
  {
    u[j] = v[j] = 0;
    p[j] = way[j] = 0;
  }
  for (i = 1; i <= r; i++)
  {
    int j0 = 0;
    p[0] = i;
    for (j = 0; j <= r; j++)
    {
      minv[j] = INFINITY;
      used[j] = 0;
    }
    do
    {
      int i0 = p[j0], j1 = 0;
      double delta = INFINITY;
      used[j0] = 1;
      for (j = 1; j <= r; j++)
        if (!used[j])
        {
          double cur = -w[i0 - 1][j - 1] - u[i0] - v[j];
          if (cur < minv[j])
          {
            minv[j] = cur;
            way[j] = j0;
          }
          if (minv[j] < delta)
          {
            delta = minv[j];
            j1 = j;
          }
        }
      for (j = 0; j <= r; j++)
      {
        if (used[j])
        {
          u[p[j]] += delta;
          v[j] -= delta;
        }
        else
          minv[j] -= delta;
      }
      j0 = j1;
    } while (p[j0] != 0);
    do
    {
      int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }
  for (j = 1; j <= r; j++)
    perm[p[j] - 1] = j - 1;
}

/* Signed permutation that best matches the components of U (n x r) to those of Uref:
 * reference component l corresponds to current component perm[l] with sign sgn[l], chosen to
 * maximise sum_l |corr(Uref[,l], U[,perm[l]])|. C (r x r) is workspace. */
void align_components(int n, int r, double **Uref, double **U, int *perm, int *sgn, double **C)
{
  int i, l, k;
  double nref[r], ncur[r];
  for (l = 0; l < r; l++)
  {
    nref[l] = ncur[l] = 0;
    for (k = 0; k < r; k++)
      C[l][k] = 0;
  }
  for (i = 0; i < n; i++)
    for (l = 0; l < r; l++)
    {
      nref[l] += Uref[i][l] * Uref[i][l];
      ncur[l] += U[i][l] * U[i][l];
      for (k = 0; k < r; k++)
        C[l][k] += Uref[i][l] * U[i][k];
    }
  for (l = 0; l < r; l++)
    for (k = 0; k < r; k++)
      C[l][k] = (nref[l] > 0 && ncur[k] > 0) ? C[l][k] / sqrt(nref[l] * ncur[k]) : 0;
  double **W = dmatrix(0, r - 1, 0, r - 1);
  for (l = 0; l < r; l++)
    for (k = 0; k < r; k++)
      W[l][k] = fabs(C[l][k]);
  max_assignment(r, W, perm);
  free_dmatrix(W, 0, r - 1, 0, r - 1);
  for (l = 0; l < r; l++)
    sgn[l] = (C[l][perm[l]] < 0) ? -1 : 1;
}

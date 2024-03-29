#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../tests.h"
#include "../../include/globals.h"
#include "../../core/covariance.matrix.h"
#include "../../math/linear.algebra.h"

#define TRUNCATE_LAMBDA(lambda) \
  if (lambda > 1) \
    lambda = 1; \
  if (lambda < 0) \
    lambda = 0;

/* shrinked mutual information, to be used in C code. */
double c_shmi(int *xx, int llx, int *yy, int lly, int num, int scale) {

int i = 0, j = 0, k = 0, ncomplete = 0;
double **n = NULL, *ni = NULL, *nj = NULL;
double lambda = 0, target = 1/(double)(llx * lly);
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = (double **) Calloc2D(llx, lly, sizeof(double));
  ni = Calloc1D(llx, sizeof(double));
  nj = Calloc1D(lly, sizeof(double));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    if ((xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER)) {

      n[xx[k] - 1][yy[k] - 1]++;
      ncomplete++;

    }/*THEN*/

  /* if there are no complete data points, or if there is just a single complete
   * observation, return independence. */
  if (ncomplete <= 1)
    goto free_and_return;

  /* estimate the optimal lambda for the data. */
  mi_lambda((double *)n, &lambda, target, ncomplete, llx, lly, 0);

  /* switch to the probability scale and shrink the estimates. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
        n[i][j] = lambda * target + (1 - lambda) * n[i][j] / ncomplete;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++) {

    ni[i] += n[i][j];
    nj[j] += n[i][j];

  }/*FOR*/

  /* compute the mutual information from the joint and marginal frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      if (n[i][j] != 0)
        res += n[i][j] * log(n[i][j] / (ni[i] * nj[j]));

  if (scale)
    res *= 2 * ncomplete;

free_and_return:

  Free1D(ni);
  Free1D(nj);
  Free2D(n, llx);

  return res;

}/*C_SHMI*/

/* shrinked conditional mutual information, to be used in C code. */
double c_shcmi(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, int scale) {

int i = 0, j = 0, k = 0, ncomplete = 0;
double ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
double lambda = 0, target = 1/(double)(llx * lly * llz);
double res = 0;

  /* compute the degrees of freedom. */
  *df = (double)(llx - 1) * (double)(lly - 1) * (double)(llz);

  /* initialize the contingency table and the marginal frequencies. */
  n = (double ***) Calloc3D(llx, lly, llz, sizeof(double));
  ni = (double **) Calloc2D(llx, llz, sizeof(double));
  nj = (double **) Calloc2D(lly, llz, sizeof(double));
  nk = Calloc1D(llz, sizeof(double));

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < num; k++)
    if ((xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER) && (zz[k] != NA_INTEGER)) {

      n[xx[k] - 1][yy[k] - 1][zz[k] - 1]++;
      ncomplete++;

    }/*THEN*/

  /* if there are no complete data points, or if there is just a single complete
   * observation, return independence. */
  if (ncomplete <= 1)
    goto free_and_return;

  /* estimate the optimal lambda for the data. */
  mi_lambda((double *)n, &lambda, target, ncomplete, llx, lly, llz);

  /* switch to the probability scale and shrink the estimates. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++)
        n[i][j][k] = lambda * target + (1 - lambda) * n[i][j][k] / ncomplete;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++) {

        ni[i][k] += n[i][j][k];
        nj[j][k] += n[i][j][k];
        nk[k] += n[i][j][k];

      }/*FOR*/

  for (k = 0; k < llz; k++) {

    /* check each level of the conditioning variable to avoid (again)
     * "divide by zero" errors. */
    if (nk[k] == 0)
      continue;

    for (j = 0; j < lly; j++) {

      for (i = 0; i < llx; i++) {

        if (n[i][j][k] > 0)
          res += n[i][j][k] * log( (n[i][j][k] * nk[k]) / (ni[i][k] * nj[j][k]) );

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

free_and_return:

  Free1D(nk);
  Free2D(ni, llx);
  Free2D(nj, lly);
  Free3D(n, llx, lly);

  if (scale)
    res *= 2 * ncomplete;

  return res;

}/*C_SHCMI*/

/* compute the shrinkage intensity lambda for the mutual information. */
void mi_lambda(double *n, double *lambda, double target, int num, int llx,
    int lly, int llz) {

double lden = 0, lnum = 0, temp = 0;

  /* compute the numerator and the denominator of the shrinkage intensity;
   * if the third dimension is a NULL pointer it's a 2-dimensional table. */
  if (llz == 0) {

    for (int i = 0; i < llx; i++)
      for (int j = 0; j < lly; j++) {

        temp = ((double **)n)[i][j] / (double)(num);
        lnum += temp * temp;
        temp = target - ((double **)n)[i][j] / (double)(num);
        lden += temp * temp;

      }/*FOR*/

  }/*THEN*/
  else {

    for (int i = 0; i < llx; i++)
      for (int j = 0; j < lly; j++)
        for (int k = 0; k < llz; k++) {

          temp = ((double ***)n)[i][j][k] / (double)(num);
          lnum += temp * temp;
          temp = target - ((double ***)n)[i][j][k] / (double)(num);
          lden += temp * temp;

      }/*FOR*/

  }/*ELSE*/

   /* compute the shrinkage intensity (avoiding "divide by zero" errors). */
  if (lden == 0)
    *lambda = 1;
  else
    *lambda = (1 - lnum) / ((double)(num - 1) * lden);

  /* bound the shrinkage intensity in the [0,1] interval. */
  TRUNCATE_LAMBDA(*lambda);

}/*MI_LAMBDA*/


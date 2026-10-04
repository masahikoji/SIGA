/* Regeneration kernel for the fixed-score reference randomisation test under the
 * stochastic Pocock-Simon range method.  Compiled with  R CMD SHLIB rt_kernel.c .
 *
 * Draws the same uniforms, in the same order, as the pure-R implementation
 * (one unif_rand() per path per participant), so that with a common seed the
 * regenerated statistics are bit-identical to rt_corrected_pvalues_R().
 *
 * .Call(C_rt_regenerate, X, scores, B, pbc, weights)
 *   X       integer n x K matrix of 0/1 factor levels (participant order)
 *   scores  double  n x L matrix of retained residual vectors
 *   B       integer number of regenerated paths
 *   pbc     double  biased-coin probability
 *   weights double  length K+1 (overall, factors)
 * returns   double  B x L matrix of regenerated statistics T*_j = 0.5 * sum_i z_ij r_il
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <math.h>

SEXP C_rt_regenerate(SEXP X_, SEXP scores_, SEXP B_, SEXP pbc_, SEXP weights_)
{
    const int n = Rf_nrows(X_);
    const int K = Rf_ncols(X_);
    const int L = Rf_ncols(scores_);
    const int B = Rf_asInteger(B_);
    const double pbc = Rf_asReal(pbc_);
    const int *X = INTEGER(X_);
    const double *scores = REAL(scores_);
    const double *w = REAL(weights_);

    if (Rf_nrows(scores_) != n) Rf_error("scores must have n rows");
    if (Rf_length(weights_) != K + 1) Rf_error("weights must have length K + 1");
    if (B < 1) Rf_error("B must be positive");

    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, B, L));
    double *t_rand = REAL(out);
    for (long long q = 0; q < (long long)B * L; q++) t_rand[q] = 0.0;

    int *overall = (int *) R_alloc(B, sizeof(int));
    int *marginal = (int *) R_alloc((size_t)B * 2 * K, sizeof(int)); /* column-major B x 2K */
    for (int b = 0; b < B; b++) overall[b] = 0;
    for (long long q = 0; q < (long long)B * 2 * K; q++) marginal[q] = 0;

    GetRNGstate();
    for (int i = 0; i < n; i++) {
        for (int b = 0; b < B; b++) {
            /* score_plus - score_minus, using |d+1| - |d-1| = 2 sgn(d) for integer d */
            double diff = w[0] * (fabs((double)(overall[b] + 1)) - fabs((double)(overall[b] - 1)));
            for (int j = 0; j < K; j++) {
                int col = 2 * j + X[i + (long long)j * n];
                int cur = marginal[b + (long long)col * B];
                diff += w[j + 1] * (fabs((double)(cur + 1)) - fabs((double)(cur - 1)));
            }
            double u = unif_rand();
            int z;
            if (diff < 0.0)      z = (u < pbc) ? 1 : -1;        /* plus preferred */
            else if (diff > 0.0) z = (u < 1.0 - pbc) ? 1 : -1;  /* minus preferred */
            else                 z = (u < 0.5) ? 1 : -1;        /* tie */
            overall[b] += z;
            for (int j = 0; j < K; j++) {
                int col = 2 * j + X[i + (long long)j * n];
                marginal[b + (long long)col * B] += z;
            }
            for (int l = 0; l < L; l++) {
                t_rand[b + (long long)l * B] += 0.5 * z * scores[i + (long long)l * n];
            }
        }
    }
    PutRNGstate();
    UNPROTECT(1);
    return out;
}

static const R_CallMethodDef callMethods[] = {
    {"C_rt_regenerate", (DL_FUNC) &C_rt_regenerate, 5},
    {NULL, NULL, 0}
};

void R_init_rt_kernel(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

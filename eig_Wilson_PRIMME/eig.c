// -----------------------------------------------------------------
// Eigenvalue computation and helper functions
#include "eig_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Apply Wilson D D^dag with EVENANDODD parity hard-coded
void Matrix_Vec_mult(wilson_vector *src, wilson_vector *res) {
  register site *s;
  register int i;
  Real norm = 0.25 / (kappa * kappa);

  FORALLSITES(i, s)   // Just what copy_wvec does
    s->chi = src[i];

  dslash_w_site(F_OFFSET(chi), F_OFFSET(mp), MINUS, EVENANDODD);
  FORALLSITES(i, s)
    scalar_mult_add_wvec(&(s->chi), &(s->mp), -kappa, &(s->tmp));

  dslash_w_site(F_OFFSET(tmp), F_OFFSET(mp), PLUS, EVENANDODD);
  FORALLSITES(i, s) {
    scalar_mult_add_wvec(&(s->tmp), &(s->mp), -kappa, &(s->psi));
    scalar_mult_wvec(&(s->psi), norm, &(res[i]));
  }
}

// Apply hermitian-Wilson D^dag with EVENANDODD parity hard-coded
void MatVec(wilson_vector *src, wilson_vector *res) {
  register site *s;
  register int i;
  Real norm = 0.5 / kappa;

  FORALLSITES(i, s)
    mult_by_gamma(&(src[i]), &(s->chi), GAMMAFIVE);

  dslash_w_site(F_OFFSET(chi), F_OFFSET(mp), MINUS, EVENANDODD);
  FORALLSITES(i, s) {
    scalar_mult_add_wvec(&(s->chi), &(s->mp), -kappa, &(s->psi));
    scalar_mult_wvec(&(s->psi), norm, &(res[i]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function av_ov gets a number of vectors (stored consecutively)
// which need to be multiplied by the matrix
void av_ov (void *x, void *y, int *Nvecs, primme_params *primme) {
  register site* s;
  int i, j, ivec;
  double *xx;
  Real *yy;
  wilson_vector src[sites_on_node], res[sites_on_node];

  // Copy double precision complex vector x to Real wilson_vector src
  for (ivec = 0; ivec < *Nvecs; ivec++) {
    xx = ((double*) x) + 24 * ivec * sites_on_node;
    FORALLSITES(i, s) {
      yy = &(src[i].d[0].c[0].real);
      for (j = 0; j < 24; j++)
        *(yy++) = *(xx++);
    }

#ifdef DDdag
    Matrix_Vec_mult(src, res);    // D D^dag
#else
    MatVec(src, res);             // H^dag
#endif

    // Copy the resulting wilson_vector res back to complex vector y
    xx = ((double*) y) + 24 * ivec * sites_on_node;
    FORALLSITES(i, s) {
      yy = &(res[i].d[0].c[0].real);
      for (j = 0; j < 24; j++)
        *(xx++) = *(yy++);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function par_GlobalSumDouble is set as primme.globalSumDouble
void par_GlobalSumDouble(void *sendBuf, void *recvBuf,
                         int *count, primme_params *primme) {

  int i;
  for (i = 0; i < *count; i++)
    *((double*)recvBuf + i) = *((double*)sendBuf + i);

  g_vecdoublesum((double*)recvBuf, *count);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function make_evs computes eigenvectors of H5=gamma_5 D through PRIMME
// First calculate eigenvectors of Q^2 then multiply by Q to find signs
int make_evs(int Nvecs, wilson_vector **eigVecs, double *eigVals) {
  register site* s;
  int i, k, l, ivec, iter, ret, maxn = sites_on_node * 12;
  double *rnorms;
  double_complex *workVecs;
  static primme_params primme;

  // Allocate memory for double-precision-only eigenvalue finder
  workVecs = (double_complex *)malloc(Nvecs * maxn * sizeof(double_complex));
  if (workVecs == NULL)
    exit(1);
  rnorms = (double *)malloc(Nvecs * sizeof(double_complex));
  if (rnorms == NULL)
    exit(1);

  // Initialize all the eigenvectors to half-random vectors
  // and copy into double-precision temporary fields
  for (ivec = 0; ivec < Nvecs; ivec++) {
    eigVals[ivec] = 1e16;
    for (k = 0; k < 4; k++) {
      for (l = 0; l < 3; l++) {
        FORALLSITES(i, s) {
          eigVecs[ivec][i].d[k].c[l] = cmplx(((Real)rand()) / RAND_MAX, 0.);
          iter = 12 * (ivec * sites_on_node + i);
          workVecs[iter].real = (double)eigVecs[ivec][i].d[k].c[l].real;
          workVecs[iter].imag = (double)eigVecs[ivec][i].d[k].c[l].imag;
        }
      }
    }
  }

  // Set the parameters of the EV finder
  primme_initialize(&primme);
  primme.n = maxn * number_of_nodes;              // Global size of matrix
  primme.nLocal = maxn;                           // Local volume
  primme.maxOuterIterations = maxIter;
  primme.maxMatvecs = maxIter + 5;
  primme.numProcs = number_of_nodes;
  primme.procID = this_node;
  primme.globalSumDouble = par_GlobalSumDouble;   // Wrapped function
  primme.matrixMatvec = av_ov;                    // Mat-vec product

  primme_set_method(DEFAULT_MIN_MATVECS, &primme);
  primme.printLevel = 1;
  primme.eps = eig_tol;                   // Maximum residual
  primme.numEvals = Nvecs;
  primme.initSize = 0;                    // Number of initial guesses
#ifdef DDdag
  primme.target = primme_smallest;
#else
  primme.target = primme_closest_abs;
  primme.numTargetShifts = 1;
  double target[1];
  target[0] = 0;
  primme.targetShifts = target;
#endif
//  primme_display_params(primme);

  // Call the actual EV finder and check return value
  ret = zprimme(eigVals, (Complex_Z *)workVecs, rnorms, &primme);
  if (ret != 0) {
    node0_printf("Failed with return value %i\n", ret);
    fflush(stdout);
//    exit(1);
  }

  // Print results
  node0_printf("Mass=%.4g --> kappa=%.4g\n", mass, kappa);
  for (ivec = 0; ivec < Nvecs; ivec++)
    node0_printf("EIGENVALUE %d %.8g\n", ivec, eigVals[ivec]);
  fflush(stdout);

  // Clean up
  free(workVecs);
  free(rnorms);
  primme_Free(&primme);
  return primme.stats.numOuterIterations;
}
// -----------------------------------------------------------------

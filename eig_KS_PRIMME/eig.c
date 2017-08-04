// -----------------------------------------------------------------
// Eigenvalue computation and helper functions
#include "eig_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar multiply an SU3 vector in the lattice
// Hard-code EVENANDODD parity
void scalar_mult_latvec(field_offset src, Real scalar, field_offset dest) {
  register int i;
  register site *s;
  register vector *spt, *dpt;

  FORALLSITES(i, s) {
    spt = (vector *)F_PT(s, src);
    dpt = (vector *)F_PT(s, dest);
    scalar_mult_vector(spt, scalar, dpt);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Apply staggered D^dag.D = -D.D, with EVEN parity hard-coded
void Matrix_Vec_mult(vector *src, vector *res) {
  register site *s;
  register int i;

  FOREVENSITES(i, s)
    s->chi = src[i];
  FORODDSITES(i, s)
    clearvec(&(s->chi));

  dslash(F_OFFSET(chi), F_OFFSET(temp), ODD);
  dslash(F_OFFSET(temp), F_OFFSET(psi), EVEN);

  scalar_mult_latvec(F_OFFSET(psi), -1, F_OFFSET(psi));
  FOREVENSITES(i, s)
    res[i] = s->psi;
}

// Apply just staggered D^dag, with EVEN parity hard-coded
void MatVec(vector *src, vector *res) {
  register site *s;
  register int i;

  FOREVENSITES(i, s)
    s->chi = src[i];
  FORODDSITES(i, s)
    clearvec(&(s->chi));

  dslash(F_OFFSET(chi), F_OFFSET(psi), EVEN);
  scalar_mult_latvec(F_OFFSET(psi), -1, F_OFFSET(psi));
  FOREVENSITES(i, s)
    res[i] = s->psi;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function av_ov gets a number of vectors (stored consecutively)
// which need to be multiplied by the matrix
// Note that even sites are stored first
void av_ov (void *x, void *y, int *Nvecs, primme_params *primme) {
  register site *s;
  int i, j, ivec;
  double *xx;
  Real *yy;
  vector src[even_sites_on_node], res[even_sites_on_node];

  // Copy double precision complex vector x to Real vector src
  for (ivec = 0; ivec < *Nvecs; ivec++) {
    xx = ((double*) x) + 6 * ivec * even_sites_on_node;
    FOREVENSITES(i, s) {
      yy = &(src[i].c[0].real);
      for (j = 0; j < 6; j++)
        *(yy++) = *(xx++);
    }

    Matrix_Vec_mult(src, res);

    // Copy the resulting vector res back to complex vector y
    xx = ((double*) y) + 6 * ivec * even_sites_on_node;
    FOREVENSITES(i, s) {
      yy = &(res[i].c[0].real);
      for (j = 0; j < 6; j++)
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
// Function make_evs computes eigenvectors through PRIMME
// Post-processing hits resulting eigenvectors with D^dag
// to explore their helicity
int make_evs(double *start, int Nvecs, vector **eigVecs, double *eigVals) {
  register site *s;
  int i, ivec, iter, ret, maxn = even_sites_on_node * 3;
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

  // Initialize all the eigenvectors to random vectors
  // We will only use the even sites, which come first in the arrays
  for (ivec = 0; ivec < Nvecs; ivec++) {
    eigVals[ivec] = 1e16;
    grsource_imp();
    FOREVENSITES(i, s)
      vec_copy(&(s->g_rand), &(eigVecs[ivec][i]));
  }

  // Copy initial guesses into double-precision temporary fields
  for (ivec = 0; ivec < Nvecs; ivec++) {
    FOREVENSITES(i, s) {
      iter = 3 * (ivec * even_sites_on_node + i);
      workVecs[iter].real = (double)eigVecs[ivec][i].c[0].real;
      workVecs[iter].imag = (double)eigVecs[ivec][i].c[0].imag;
      iter++;
      workVecs[iter].real = (double)eigVecs[ivec][i].c[1].real;
      workVecs[iter].imag = (double)eigVecs[ivec][i].c[1].imag;
      iter++;
      workVecs[iter].real = (double)eigVecs[ivec][i].c[2].real;
      workVecs[iter].imag = (double)eigVecs[ivec][i].c[2].imag;
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
  if (*start <= 0) {
    primme.initSize = 0;                  // Number of initial guesses
    primme.target = primme_smallest;
  }
  else {
    primme.initSize = Nvecs / 2;          // Number of initial guesses
    primme.target = primme_closest_geq;
    primme.numTargetShifts = 1;
    primme.targetShifts = start;
  }
//  primme_display_params(primme);

  // Call the actual EV finder and check return value
  ret = zprimme(eigVals, (Complex_Z *)workVecs, rnorms, &primme);
  if (ret != 0) {
    node0_printf("Failed with return value %i\n", ret);
    fflush(stdout);
//    exit(1);
  }

  // Clean up
  free(workVecs);
  free(rnorms);
  primme_Free(&primme);
  return primme.stats.numOuterIterations;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure the chirality of a unit-normalized fermion state
// Chirality is real since gamma_5 is hermitian
// Hard-code EVENANDODD parity
void measure_chirality(vector *src, double *chirality) {
  register int i;
  register site *s;
  register double cc = 0;
  complex tmp;

  FORALLSITES(i, s)
    vec_copy(&src[i], &(s->tempvec[3]));

  mult_spin_pseudoscalar(F_OFFSET(tempvec[3]), F_OFFSET(temp));

  FORALLSITES(i, s) {
    tmp = su3_dot(&(s->tempvec[3]), &(s->temp));
    cc += tmp.real;
  }
  *chirality = cc;
  g_doublesum(chirality);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print the density and chiral density of a normalized fermion state
// Hard-code EVEN parity
void print_densities(vector *src, char *tag, int y, int z, int t) {

  register int i;
  register site *s;
  complex tmp1, tmp2;

  FOREVENSITES(i, s)
    vec_copy(&src[i], &(s->tempvec[3]));

  mult_spin_pseudoscalar(F_OFFSET(tempvec[3]), F_OFFSET(temp));

  FOREVENSITES(i, s) {
    if (s->y == y && s->z == z && s->t == t) {
      tmp1 = su3_dot(&(s->tempvec[3]), &(s->temp));
      tmp2 = su3_dot(&(s->tempvec[3]), &(s->tempvec[3]));
      node0_printf("%s: %i %.4g %.4g %.4g\n", tag, s->x,
                   tmp2.real, tmp1.real, tmp1.imag);
    }
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Compute eigenvalues and eigenvectors of the Kogut--Susskind dslash^2
// Includes and definitions
#include "block_includes.h"

#define JACOBI_TOL 1.110223e-16
#define MINITER 5
//#define EIG_DEBUG
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Apply the matrix whose eigenvalues we are computing
// Specifically, staggered D^dag D = -D^2, with EVEN parity hard-coded
void Matrix_Vec_mult(su3_vector *src, su3_vector *res) {
  register site *s;
  register int i;

  FORALLSITES(i, s)
    s->chi = src[i];

  dslash(F_OFFSET(chi), F_OFFSET(temp), ODD);
  dslash(F_OFFSET(temp), F_OFFSET(psi), EVEN);
  scalar_mult_latvec(F_OFFSET(psi), -1, F_OFFSET(psi), EVEN);
  FOREVENSITES(i, s)
    res[i] = s->psi;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the norm of the given vector, with EVEN parity hard-coded
void norm(su3_vector *vec, double *norm) {
  register int i;
  register double n = 0;
  register site *s;

  FOREVENSITES(i, s)
    n += magsq_su3vec(&(vec[i]));

  *norm = n;
  g_doublesum(norm);
  *norm = sqrt(*norm);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Normalize the given vector, with EVEN parity hard-coded
void normalize(su3_vector *vec) {
  register int i;
  register site *s;
  double N;

  norm(vec, &N);
  N = 1.0 / N;
  FOREVENSITES(i, s)
    scalar_mult_su3_vector(&(vec[i]), (Real)N, &(vec[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the dot product of the two given vectors
// Hard-code EVEN parity
void dot_product(su3_vector *vec1, su3_vector *vec2, double_complex *dot) {
  register int i;
  register double re = 0, im = 0;
  register site *s;
  complex cc;

  FOREVENSITES(i, s) {
    cc = su3_dot(&(vec1[i]), &(vec2[i]));
    re += cc.real;
    im += cc.imag;
  }
  dot->real = re;
  dot->imag = im;
  g_dcomplexsum(dot);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute vec2 = vec2 - cc * vec1, with EVEN parity hard-coded
void complex_vec_mult_sub(double_complex *cc, su3_vector *vec1,
                          su3_vector *vec2) {

  register int i;
  register site *s;
  complex sc;

  sc.real = (Real)(cc->real);
  sc.imag = (Real)(cc->imag);
  FOREVENSITES(i, s)
    c_scalar_mult_dif_su3vec(&(vec1[i]), &sc, &(vec2[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Project out the N given **vectors from the given *vec
// The vectors are assumed to be orthonormal, EVEN parity is hard-coded
void project_out(su3_vector *vec, su3_vector **vector, int N) {
  register int i;
  double_complex cc;

  for (i = N - 1; i > -1; i--) {
    dot_product(vector[i], vec, &cc);
    complex_vec_mult_sub(&cc, vector[i], vec);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy src to res
void copy_Vector(su3_vector *src, su3_vector *res) {
  memcpy((void *)res, (void *)src, sites_on_node*sizeof(su3_vector));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute vec = vec - rr * vec1, with EVEN parity hard-coded
void double_vec_mult_sub(double *rr, su3_vector *vec1, su3_vector *vec) {
  register  int i;
  register site *s;

  FOREVENSITES(i, s)
    scalar_mult_sub_su3_vector(&(vec[i]), &(vec1[i]), (Real)*rr, &(vec[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute vec = a * vec + b * vec2, with EVEN parity hard-coded
void dax_p_by(double *a, su3_vector *vec, double *b, su3_vector *vec2) {
  register int i;
  register site *s;

  FOREVENSITES(i, s) {
    scalar_mult_su3_vector(&(vec[i]), (Real)(*a), &(vec[i]));
    scalar_mult_add_su3_vector(&(vec[i]), &(vec2[i]), (Real)(*b), &(vec[i]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute vec2 = vec1 + a * vec2, with EVEN parity hard-coded
void vec_plus_double_vec_mult(su3_vector *vec1, double *a,
                              su3_vector *vec2) {

  register int i;
  register site *s;

  FOREVENSITES(i, s) {
    scalar_mult_su3_vector(&(vec2[i]), *a, &(vec2[i]));
    add_su3_vector(&(vec1[i]), &(vec2[i]), &(vec2[i]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute vec2 = vec2 + cc * vec1, with EVEN parity hard-coded
void complex_vec_mult_add(double_complex *cc, su3_vector *vec1,
                          su3_vector *vec2) {

  register int i;
  register site *s;
  complex sc;

  sc.real = (Real)(cc->real);
  sc.imag = (Real)(cc->imag);
  FOREVENSITES(i, s)
    c_scalar_mult_sum_su3vec(&(vec1[i]), &sc, &(vec2[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Hard-code EVEN parity
int Rayleigh_min(su3_vector *vec, su3_vector **eigVec, Real Tolerance,
                 Real RelTol, int Nvecs, int maxIter, int restart) {

  int iter;
  double beta = 1e16, cos_theta, sin_theta;
  double g_norm, old_g_norm, start_g_norm;
  double quot, P_norm, theta, real_vecMp, pMp;
#ifdef EIG_DEBUG
  double vec_norm;
#endif
  double_complex cc;
  su3_vector *Mvec, *grad, *P, *MP;

  Mvec = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
  grad = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
  P = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
  MP = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));

  project_out(vec, eigVec, Nvecs);
  normalize(vec);
  Matrix_Vec_mult(vec, Mvec);
  project_out(Mvec, eigVec, Nvecs);

  // Compute the quotient vec.M.vec, real since M is hermitian
  dot_product(vec, Mvec, &cc);
  quot = cc.real;
  // Compute the grad = M * vec - quot * vec
  copy_Vector(Mvec, grad);
  double_vec_mult_sub(&quot, vec, grad);
  // Set P (the search direction) equal to grad
  copy_Vector(grad, P);
  // Compute the norms of P and grad
  norm(P, &P_norm);
  norm(grad, &g_norm);
  start_g_norm = g_norm;
#ifdef EIG_DEBUG
  node0_printf("Rayleigh_min 0: quot = (%.4g, %.4g), ", quot, cc.imag);
  node0_printf("g = %.4g, b = %.4g, P = %.4g\n", g_norm, beta, P_norm);
#endif

  iter = 0;
  while (g_norm > Tolerance
         && ((iter < maxIter && g_norm / start_g_norm > RelTol)
             || iter < MINITER)) {
    iter++;
    Matrix_Vec_mult(P, MP);
    dot_product(vec, MP, &cc);
    real_vecMp = cc.real;
    dot_product(P, MP, &cc);
    pMp = cc.real;                // p.M.p is real
    theta = 0.5 * atan(2.0 * real_vecMp / (quot * P_norm - pMp / P_norm));
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    if (sin_theta * cos_theta * real_vecMp > 0) {
      theta = theta - 0.5 * PI;   // Choose the minimum, not the maximum
      sin_theta = sin(theta);
      cos_theta = cos(theta);
    }
    // vec = cos(theta) * vec + P * sin(theta) / P_norm
    // Mvec = cos(theta) * Mvec + MP * sin(theta) / P_norm
    sin_theta /= P_norm;
    dax_p_by(&cos_theta, vec, &sin_theta, P);
    dax_p_by(&cos_theta, Mvec, &sin_theta, MP);

    // Renormalize vec
    if (iter % restart == 0) {
#ifdef EIG_DEBUG
      norm(vec, &vec_norm);
      node0_printf("Renormalizing: norm = %g\n", vec_norm);
#endif
      // Project vec on the orthogonal complement of eigVec
      project_out(vec, eigVec, Nvecs);
      normalize(vec);
      Matrix_Vec_mult(vec, Mvec);

      // Recompute the quotient vec.M.vec, real since M is hermitian
      dot_product(vec, Mvec, &cc);
      quot = cc.real;

      // Recompute the grad and g_norm
      copy_Vector(Mvec, grad);
      double_vec_mult_sub(&quot, vec, grad);
      norm(grad, &g_norm);

      // Project P on the orthogonal complement of eigVec
      project_out(P, eigVec, Nvecs);
      // Make P orthogonal to vec
      dot_product(vec, P, &cc);
      complex_vec_mult_sub(&cc, vec, P);
      // Make P orthogonal to grad and recompute the P_norm
      dot_product(grad, P, &cc);
      complex_vec_mult_sub(&cc, grad, P);
      norm(P, &P_norm);
    }

    dot_product(vec, Mvec, &cc);
    // quot = vec.M.vec is real since M is hermitian
    quot = cc.real;
#ifdef EIG_DEBUG
    node0_printf("Rayleigh_min %i: quot = (%.4g, %.4g), ",
                 iter, quot, cc.imag);
    node0_printf("g = %.4g, b = %.4g, P = %.4g\n", g_norm, beta, P_norm);
#endif
    old_g_norm = g_norm;

    copy_Vector(Mvec, grad);
    double_vec_mult_sub(&quot, vec, grad);

    norm(grad, &g_norm);
    beta = cos_theta * g_norm * g_norm / (old_g_norm * old_g_norm);
    // Cut off beta
    if (beta > 2.0)
      beta = 2;

    dot_product(vec, P, &cc);
    cc.real *= beta;
    cc.imag *= beta;
    vec_plus_double_vec_mult(grad, &beta, P); // P = grad + beta * P
    complex_vec_mult_sub(&cc, vec, P);        // P = P - cc * vec
    norm(P, &P_norm);
  }
  project_out(vec, eigVec, Nvecs);
  normalize(vec);

  // Clean up
  free(MP);
  free(P);
  free(grad);
  free(Mvec);
  iter++;
  return iter;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct the projected matrix A and the error of each eigenvector
// Hard-code EVEN parity
void constructArray(su3_vector **eigVec, Matrix *A, double *err) {
  int i, j, Nvecs;
  su3_vector *res, *grad;
  double_complex cc, Aij, Aji;

  Nvecs = A->N;
  res = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
  grad = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
  for (i = 0; i < Nvecs; i++) {
    Matrix_Vec_mult(eigVec[i], res);
    dot_product(res, eigVec[i], &cc);
    A->M[i][i].real = cc.real;
    A->M[i][i].imag = 0;
    copy_Vector(res, grad);
    double_vec_mult_sub(&cc.real, eigVec[i], grad);
    norm(grad, &err[i]);
    for (j = i + 1; j < Nvecs; j++) {
      dot_product(res, eigVec[j], &cc);
      Aij=cc;
      CONJG(cc, Aji);
      dot_product(eigVec[j], res, &cc);
      CSUM(Aji, cc);
      CONJG(cc, cc);
      CSUM(Aij, cc);
      CMULREAL(Aij, 0.5, A->M[i][j]);
      CMULREAL(Aji, 0.5, A->M[j][i]);
    }
  }
  free(grad);
  free(res);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Hard-code EVEN parity
void RotateBasis(su3_vector **eigVec, Matrix *V) {
  register int i, j, k, N;
  register site *s;
  su3_vector **Tmp;

  N = V->N;
  // Allocate the temporary vectors needed
  Tmp = (su3_vector **)malloc(N * sizeof(su3_vector *));
  for (j = 0; j < N; j++)
    Tmp[j] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  for (j = 0; j < N; j++) {
    FOREVENSITES(i, s)
      clearvec(&Tmp[j][i]);

    for (k = 0; k < N; k++)
      complex_vec_mult_add(&V->M[k][j], eigVec[k], Tmp[j]);
  }

  // Copy rotated basis to the eigVec and free temporaries
  for (i = 0; i < N; i++) {
    copy_Vector(Tmp[i], eigVec[i]);
    free(Tmp[i]);
  }
  free(Tmp);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Hard-code EVEN parity
int Kalkreuter(su3_vector **eigVec, double *eigVal, Real Tolerance,
               Real RelTol, int Nvecs, int maxIter,
               int restart, int kiters) {

  register int i, j;
  register site *s;
  int iter = 0, total_iters = 0;
  double max_error = 1e10;
  double *grad, *err;
  Matrix Array, V;
  su3_vector *vec;

  // Allocate the array and eigenvector matrix
  Array = AllocateMatrix(Nvecs);
  V = AllocateMatrix(Nvecs);

  vec = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
  grad = (double *)malloc(Nvecs * sizeof(double));
  err = (double *)malloc(Nvecs * sizeof(double));

  // Initialize all the eigenvectors to random vectors
  for (j = 0; j < Nvecs; j++) {
    eigVal[j] = 1e16;
    grad[j] = 1e10;
    grsource_eig();
    FOREVENSITES(i, s)
      su3vec_copy(&(s->g_rand), &(eigVec[j][i]));
  }

  while (max_error > Tolerance && iter < kiters) {
    iter++;
    for (j = 0; j < Nvecs; j++) {
      if (grad[j] > Tolerance) {
        copy_Vector(eigVec[j], vec);
        total_iters += Rayleigh_min(vec, eigVec, Tolerance, RelTol,
                                    j, maxIter, restart);
        copy_Vector(vec, eigVec[j]);
      }
    }
    constructArray(eigVec, &Array, grad);

#ifdef EIG_DEBUG
    node0_printf("Eigenvalues before diagonalization\n");
    for (i = 0; i < Nvecs; i++)
      node0_printf("  quot(%i) = %g |grad| = %g\n",
                   i, Array.M[i][i].real, grad[i]);
#endif

    Jacobi(&Array, &V, JACOBI_TOL);
    sort_eigenvectors(&Array, &V);
    RotateBasis(eigVec, &V);
    constructArray(eigVec, &Array, grad);

    // Find the maximum error
    max_error = 0;
    for (i = 0; i < Nvecs; i++) {
      err[i] = eigVal[i];
      eigVal[i] = Array.M[i][i].real;
      err[i] = fabs(err[i] - eigVal[i]) / (1.0 - RelTol * RelTol);

      // Stop when g = |M * v - l * v| becomes less than the eigenvalue tolerance
      if (err[i] > max_error)
        max_error = err[i];

      // This results in some overkill
      // It is also possible to stop when the estimated error of the eigenvalue
      // gets smaller than the eigenvalue tolerance
      // This would exploit the quadratic convergence of the algorithm
      // (See Kalkreuter's paper for details)
//      if (grad[i] > max_error)
//        max_error = grad[i];
    }

    node0_printf("\nEigenvalues after diagonalization ");
    node0_printf("at Kalkreuter iteration %i\n", iter);
    for (i = 0; i < Nvecs; i++)
      node0_printf("quot(%d) = %.4g +/- %.4g, |grad| = %.4g\n",
                   i, eigVal[i], err[i], grad[i]);

  }

  node0_printf("\n");
  if (iter == kiters) {
    node0_printf("WARNING: %d Kalkreuter iterations saturated\n", iter);
    node0_printf("         Algorithm may not have converged\n");
  }
  for (i = 0; i < Nvecs; i++)
    node0_printf("EIGENVALUE %d %.6g +/- %.4g\n", i, eigVal[i], err[i]);

  // Deallocate arrays
  deAllocate(&V);
  deAllocate(&Array);
  free(err);
  free(grad);
  free(vec);
  return total_iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure the chirality of a unit-normalized fermion state
// Chirality is real since gamma_5 is hermitian
// Hard-code EVENANDODD parity
void measure_chirality(su3_vector *src, double *chirality) {
  register int i;
  register site *s;
  register double cc = 0;
  complex tmp;

  FORALLSITES(i, s)
    su3vec_copy(&src[i], &(s->tempvec[3]));

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
void print_densities(su3_vector *src, char *tag, int y, int z, int t) {

  register int i;
  register site *s;
  complex tmp1, tmp2;

  FOREVENSITES(i, s)
    su3vec_copy(&src[i], &(s->tempvec[3]));

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

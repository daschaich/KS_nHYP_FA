// -----------------------------------------------------------------
// Runge--Kutta integrator for infinitesimal smearing
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate U = exp(A).U
// Goes to eighth order in the exponential:
//   exp(A) * U = ( 1 + A + A^2/2 + A^3/6 ...) * U
//              = U + A*(U + (A/2)*(U + (A/3)*( ... )))
void exp_mult(int dir, double eps, anti_hermitmat *A) {
  register int i;
  register site *s;
  matrix *link, tmat, tmat2, htemp;
  register Real t2, t3, t4, t5, t6, t7, t8;

  // Take divisions out of site loop (can't be done by compiler)
  t2 = eps / 2.0;
  t3 = eps / 3.0;
  t4 = eps / 4.0;
  t5 = eps / 5.0;
  t6 = eps / 6.0;
  t7 = eps / 7.0;
  t8 = eps / 8.0;

  FORALLSITES(i, s) {
    uncompress_anti_hermitian(&(A[i]), &htemp);
    link = &(s->link[dir]);

    mult_nn(&htemp, link, &tmat);
    scalar_mult_add_matrix(link, &tmat, t8, &tmat2);

    mult_nn(&htemp, &tmat2, &tmat);
    scalar_mult_add_matrix(link, &tmat, t7, &tmat2);

    mult_nn(&htemp, &tmat2, &tmat);
    scalar_mult_add_matrix(link, &tmat, t6, &tmat2);

    mult_nn(&htemp, &tmat2, &tmat);
    scalar_mult_add_matrix(link, &tmat, t5, &tmat2);

    mult_nn(&htemp, &tmat2, &tmat);
    scalar_mult_add_matrix(link, &tmat, t4, &tmat2);

    mult_nn(&htemp, &tmat2, &tmat);
    scalar_mult_add_matrix(link, &tmat, t3, &tmat2);

    mult_nn(&htemp, &tmat2, &tmat);
    scalar_mult_add_matrix(link, &tmat, t2, &tmat2);

    mult_nn(&htemp, &tmat2, &tmat);
    scalar_mult_add_matrix(link, &tmat, eps, &tmat2);

    mat_copy(&tmat2, link);    // This step updates the link U[dir]
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clear A
void clear_antiH(anti_hermitmat *a) {
  a->m01.real = 0.0;
  a->m01.imag = 0.0;
  a->m02.real = 0.0;
  a->m02.imag = 0.0;
  a->m12.real = 0.0;
  a->m12.imag = 0.0;
  a->m00im = 0.0;
  a->m11im = 0.0;
  a->m22im = 0.0;
}

// c <-- c + s * b (output is always last)
void scalar_mult_sum_antiH(anti_hermitmat *b, Real s, anti_hermitmat *c) {
  c->m01.real += s * b->m01.real;
  c->m01.imag += s * b->m01.imag;
  c->m02.real += s * b->m02.real;
  c->m02.imag += s * b->m02.imag;
  c->m12.real += s * b->m12.real;
  c->m12.imag += s * b->m12.imag;
  c->m00im += s * b->m00im;
  c->m11im += s * b->m11im;
  c->m22im += s * b->m22im;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate A = A + f1 * Project_antihermitian_traceless(U.Sdag)
//           U = exp(f2 * A).U
// S is the Lie derivative of the action being used to flow
void update_flow(double f1, double f2) {
  register int i, dir;
  register site *s;
  matrix tmat;
  anti_hermitmat tantiH;

  // Lie derivative of (Wilson) action
  // Need to compute all four before we start updating U
  // Uses tempmat for temporary storage
  staple(S);

  FORALLUPDIR(dir) {
    FORALLSITES(i, s) {
      // Project_antihermitian_traceless(U.Sdag)
      mult_na(&(s->link[dir]), &(S[dir][i]), &tmat);
      make_anti_hermitian(&tmat, &tantiH);
      // A += f1 * U.S
      scalar_mult_sum_antiH(&tantiH, (Real)f1, &(A[dir][i]));
    }
    exp_mult(dir, f2, A[dir]);                  // U = exp(f2 * A).U
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void stout_step_rk() {
  register int i, dir;
  register site *s;

  // Clear A, just in case
  FORALLSITES(i, s) {
    FORALLUPDIR(dir)
      clear_antiH(&A[dir][i]);
  }

  // Runge--Kutta coefficients computed from Eq. 2.4 of arXiv:1203.4469
  update_flow(17.0 * epsilon / 36.0, -9.0 / 17.0);
  update_flow(-8.0 * epsilon / 9.0, 1.0);
  update_flow(3.0 * epsilon / 4.0, -1.0);
}
// -----------------------------------------------------------------

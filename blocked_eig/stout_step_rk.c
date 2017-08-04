// -----------------------------------------------------------------
// Runge--Kutta integrator for infinitesimal smearing
#include "block_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate U = exp(A).U
// Goes to eighth order in the exponential:
//   exp(A) * U = ( 1 + A + A^2/2 + A^3/3 ...) * U
//              = U + A*(U + (A/2)*(U + (A/3)*( ... )))
void exp_mult(int dir, double eps, anti_hermitmat *A) {
  register int i;
  register site *s;
  matrix *link, temp1, temp2, htemp;
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

    mult_su3_nn(&htemp, link, &temp1);
    scalar_mult_add_matrix(link, &temp1, t8, &temp2);

    mult_su3_nn(&htemp, &temp2, &temp1);
    scalar_mult_add_matrix(link, &temp1, t7, &temp2);

    mult_su3_nn(&htemp, &temp2, &temp1);
    scalar_mult_add_matrix(link, &temp1, t6, &temp2);

    mult_su3_nn(&htemp, &temp2, &temp1);
    scalar_mult_add_matrix(link, &temp1, t5, &temp2);

    mult_su3_nn(&htemp, &temp2, &temp1);
    scalar_mult_add_matrix(link, &temp1, t4, &temp2);

    mult_su3_nn(&htemp, &temp2, &temp1);
    scalar_mult_add_matrix(link, &temp1, t3, &temp2);

    mult_su3_nn(&htemp, &temp2, &temp1);
    scalar_mult_add_matrix(link, &temp1, t2, &temp2);

    mult_su3_nn(&htemp, &temp2, &temp1);
    scalar_mult_add_matrix(link, &temp1, eps, &temp2);

    mat_copy(&temp2, link);    // This step updates the link U[dir]
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Project matrix to traceless anti-hermitian matrix
void proj_traceless_antiH(matrix *in, anti_hermitmat *out) {
  int i;
  complex tc;
  matrix tmat;

  su3_adjoint(in, &tmat);
  sub_matrix(in, &tmat, in);
  tc = trace_su3(in);
  CMULREAL(tc, 1 / 3.0, tc);
  for (i = 0; i < 3; i++)
    CSUB(in->e[i][i], tc, in->e[i][i]);

  make_anti_hermitian(in, out);    // Should just amount to copying data
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clear A
void clear_antiH(anti_hermitmat *a) {
  a->m01.real = 0;
  a->m01.imag = 0;
  a->m02.real = 0;
  a->m02.imag = 0;
  a->m12.real = 0;
  a->m12.imag = 0;
  a->m00im = 0;
  a->m11im = 0;
  a->m22im = 0;
}

// C <- A + s * B
void scalar_mult_add_antiH(anti_hermitmat *a, anti_hermitmat *b,
                           Real s, anti_hermitmat *c) {

  c->m01.real = a->m01.real + s * b->m01.real;
  c->m01.imag = a->m01.imag + s * b->m01.imag;
  c->m02.real = a->m02.real + s * b->m02.real;
  c->m02.imag = a->m02.imag + s * b->m02.imag;
  c->m12.real = a->m12.real + s * b->m12.real;
  c->m12.imag = a->m12.imag + s * b->m12.imag;
  c->m00im = a->m00im + s * b->m00im;
  c->m11im = a->m11im + s * b->m11im;
  c->m22im = a->m22im + s * b->m22im;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate A = A + f1 * Project_antihermitian_traceless(U.S)
//           U = exp(f2 * A).U
void update_flow(int dir, anti_hermitmat *A, matrix *S,
                 double f1, double f2) {

  register int i;
  register site *s;
//  double td = 0;
//  complex tc;
  matrix tmat;
  anti_hermitmat tantiH;

  FORALLSITES(i, s) {
    mult_su3_na(&(S[i]), &(s->link[dir]), &tmat); // S.Udag
//    tc = trace_su3(&tmat);
//    td += (double)tc.real;
    proj_traceless_antiH(&tmat, &tantiH);         // Projection
    scalar_mult_add_antiH(&tantiH, &tantiH, -0.5, &tantiH);
    scalar_mult_add_antiH(&(A[i]), &tantiH,
                          (Real)f1, &(A[i]));     // A += f1 * U.S
  }
  exp_mult(dir, -f2, A);                           // U = exp(f2 * A).U
//  node0_printf("Tr(Udag.S)[dir=%d] = %g\n", dir, td / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void stout_step_rk(matrix *S[4], anti_hermitmat *A[4]) {
  register int i, dir;
  register site *s;

  // Clear A, just in case
  FORALLSITES(i, s) {
    for (dir = 0; dir < 4; dir++)
      clear_antiH(&A[dir][i]);
  }

  // Copied from Szabolcs Borsanyi's wilson_flow.c
  staple(S);
  for (dir = 0; dir < 4; dir++)
    update_flow(dir, A[dir], S[dir], 17 * epsilon / 36.0, -9 / 17.0);
  staple(S);
  for (dir = 0; dir < 4; dir++)
    update_flow(dir, A[dir], S[dir], -8 * epsilon / 9.0, 1);
  staple(S);
  for (dir = 0; dir < 4; dir++)
    update_flow(dir, A[dir], S[dir], 3 * epsilon / 4.0, -1);
}
// -----------------------------------------------------------------

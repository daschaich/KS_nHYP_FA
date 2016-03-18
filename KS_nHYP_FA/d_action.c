// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the conjugate gradient should already
// have been run on both even and odd sites, depending on Nf
// The vectors psi[j] contain (M^dag.M)^(-1) chi

#include "ks_dyn_includes.h"
Real ahmat_mag_sq(anti_hermitmat *pt);
double d_hmom_action();
double d_fermion_action();
// -----------------------------------------------------------------



// -----------------------------------------------------------------
double d_action() {
  double ssplaq, stplaq, g_action, h_action, f_action;

  rephase(OFF);
  d_plaquette_a(&ssplaq, &stplaq);
  rephase(ON);
  g_action = -beta * volume * (ssplaq + stplaq);
  h_action = d_hmom_action();
  f_action = d_fermion_action();
  node0_printf("D_ACTION: g, h, f, tot = %.8g %.8g %.8g %.8g\n",
               g_action, h_action, f_action,
               g_action + h_action + f_action);
  return g_action + h_action + f_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to the action, following the clover code
// Checked that order of arguments to dot product doesn't matter
// Checked that su3_rdot works just as well as su3_dot.real
// Level 0 is the slow outer invert and 1 is the fast inner invert
double d_fermion_action() {
  register int i;
  register site *s;
  int j;
  double sum = 0;
  double MsqDiff_x4 = 4.0 * (MH * MH - mass * mass);

  if (num_masses == 1) {
    for (j = 0; j < full_fields; j++) {
      FORALLSITES(i, s)
        sum += (double)su3_rdot(&(s->psi[j][0]), &(s->chi[j][0]));
    }
    if (half_fields == 1) {
      j = full_fields;
      FOREVENSITES(i, s)
        sum += (double)su3_rdot(&(s->psi[j][0]), &(s->chi[j][0]));
    }
  }

  else {  // num_masses = 2
    for (j = 0; j < full_fields; j++) {
      FORALLSITES(i, s) {
        sum += (double)su3_rdot(&(s->chi[j][0]), &(s->chi[j][0]));
        sum += MsqDiff_x4 * (double)su3_rdot(&(s->psi[j][0]), &(s->chi[j][0]));
        sum += (double)su3_rdot(&(s->psi[j][1]), &(s->chi[j][1]));
      }
    }
    if (half_fields == 1) {
      j = full_fields;
      FOREVENSITES(i, s) {
        sum += (double)su3_rdot(&(s->chi[j][0]), &(s->chi[j][0]));
        sum += MsqDiff_x4 * (double)su3_rdot(&(s->psi[j][0]), &(s->chi[j][0]));
        sum += (double)su3_rdot(&(s->psi[j][1]), &(s->chi[j][1]));
      }
    }
  }
  g_doublesum(&sum);
  return(sum);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge momentum contribution to the action
double d_hmom_action() {
  register int i, dir;
  register site *s;
  double sum = 0;

  for(dir = XUP; dir <= TUP; dir++) {
   FORALLSITES(i, s)
    sum += (double)ahmat_mag_sq(&(s->mom[dir]));
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Magnitude squared of an antihermition matrix
Real ahmat_mag_sq(anti_hermitmat *pt) {
  register Real x, sum;
  x = pt->m00im;      sum  = 0.5 * x * x;
  x = pt->m11im;      sum += 0.5 * x * x;
  x = pt->m01.real;   sum += x * x;
  x = pt->m01.imag;   sum += x * x;
  x = pt->m22im;      sum += 0.5 * x * x;
  x = pt->m02.real;   sum += x * x;
  x = pt->m02.imag;   sum += x * x;
  x = pt->m12.real;   sum += x * x;
  x = pt->m12.imag;   sum += x * x;
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a gauge field as an array of four su3_matrices
void gauge_field_copy(field_offset src, field_offset dest) {
register int i, dir, src2, dest2;
register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    for (dir = XUP; dir <= TUP; dir++) {
      su3mat_copy((su3_matrix *)F_PT(s, src2),
                  (su3_matrix *)F_PT(s, dest2));
      src2 += sizeof(su3_matrix);
      dest2 += sizeof(su3_matrix);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Adds adjoint plaquette term
void d_plaquette_a(double *ss_plaq, double *st_plaq) {
  register int i, dir1, dir2;
  register site *s;
  register su3_matrix *m1, *m4;
  su3_matrix mtmp, *su3mat = malloc(sites_on_node * sizeof(*su3mat));
  double ss_sum = 0.0, st_sum = 0.0;
  complex ctt;
  msg_tag *mtag0, *mtag1;

  if (su3mat == NULL) {
    printf("d_plaquette_a: can't malloc su3mat\n");
    fflush(stdout);
    terminate(1);
  }

  for(dir1 = YUP; dir1 <= TUP; dir1++) {
    for(dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(su3_matrix),
                                dir1, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir1]), sizeof(su3_matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->link[dir1]);
        m4 = &(s->link[dir2]);
        mult_su3_an(m4, m1, &su3mat[i]);
      }
      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mult_su3_nn(&su3mat[i], (su3_matrix *)(gen_pt[0][i]), &mtmp);

        if (dir1 == TUP) {
          ctt = complextrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
          st_sum += (ctt.real / 3)
                  + beta_a * (ctt.real * ctt.real + ctt.imag * ctt.imag) / 9;
        }
        else {
          ctt = complextrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
          ss_sum += (ctt.real / 3)
                  + beta_a * (ctt.real * ctt.real + ctt.imag * ctt.imag) / 9;
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  *ss_plaq = ss_sum / ((double)volume);
  *st_plaq = st_sum / ((double)volume);
  free(su3mat);
}
// -----------------------------------------------------------------

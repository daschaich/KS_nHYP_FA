// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the conjugate gradient should already
// have been run on both even and odd sites, depending on Nf
// The vectors psi[j] contain (M^dag.M)^(-1) chi

#include "ks_dyn_includes.h"
Real ahmat_mag_sq(anti_hermitmat *pt);
double hmom_action();
double fermion_action();
// -----------------------------------------------------------------



// -----------------------------------------------------------------
double action() {
  double ssplaq, stplaq, g_action, h_action, f_action;

  rephase(OFF);
  plaquette_a(&ssplaq, &stplaq);
  rephase(ON);
  g_action = -beta * volume * (ssplaq + stplaq);
  h_action = hmom_action();
  f_action = fermion_action();
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
double fermion_action() {
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
double hmom_action() {
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
// Use tempmat for temporary storage
void plaquette_a(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  register su3_matrix *m1, *m4;
  double ss_sum = 0.0, st_sum = 0.0, td;
  complex tc;
  msg_tag *mtag, *mtag2;
  su3_matrix tmat;

  for (dir = YUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      mtag = start_gather_site(F_OFFSET(link[dir2]), sizeof(su3_matrix),
                               dir, EVENANDODD, gen_pt[0]);
      mtag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(su3_matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->link[dir]);
        m4 = &(s->link[dir2]);
        mult_su3_an(m4, m1, &tempmat[i]);
      }

      wait_gather(mtag);
      wait_gather(mtag2);
      FORALLSITES(i, s) {
        m1 = (su3_matrix *)(gen_pt[0][i]);
        m4 = (su3_matrix *)(gen_pt[1][i]);
        mult_su3_nn(&tempmat[i], m1, &tmat);
        tc = complextrace_su3(m4, &tmat);
        td = (tc.real / 3.0);
        td += beta_a * (tc.real * tc.real + tc.imag * tc.imag) / 9.0;
        if (dir == TUP)
          st_sum += td;
        else
          ss_sum += td;
      }
      cleanup_gather(mtag);
      cleanup_gather(mtag2);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  *ss_plaq = ss_sum / ((double)volume);
  *st_plaq = st_sum / ((double)volume);
}
// -----------------------------------------------------------------
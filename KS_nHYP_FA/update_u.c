// -----------------------------------------------------------------
// Update the link matrices
// Go to eighth order in the exponential of the momentum matrices,
// since higher order integrators use large steps
// Evaluation is done as:
//   exp(H) * U = (1 + H + H^2/2 + H^3/6 ...) * U
//              = U + H*(U + (H/2)*(U + (H/3)*( ... )))
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_u(Real eps) {
  register int i, dir;
  register Real t2, t3, t4, t5, t6, t7, t8;
  register site *s;
  matrix *link, temp1, temp2, htemp;

  // Take divisions out of site loop (can't be done by compiler)
  t2 = eps / 2.0;
  t3 = eps / 3.0;
  t4 = eps / 4.0;
  t5 = eps / 5.0;
  t6 = eps / 6.0;
  t7 = eps / 7.0;
  t8 = eps / 8.0;

  FORALLSITES(i, s) {
    FORALLUPDIR(dir) {
      uncompress_anti_hermitian(&(s->mom[dir]), &htemp);
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

      mat_copy(&temp2, link);
    }
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Put these here to keep libraries closer to vanilla MILC

// void scalar_add_diag_su3(su3_matrix *a, Real s)
// A <- A + s * I

// void c_scalar_add_diag_su3(su3_matrix *a, complex *c)
// A <- A + c * I
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void scalar_add_diag_su3(su3_matrix *a, Real s) {
  register int i;

  for(i = 0; i < 3; i++)
    a->e[i][i].real += s;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void c_scalar_add_diag_su3(su3_matrix *a, complex *c) {
  register int i;

  for(i = 0; i < 3; i++){
    a->e[i][i].real += c->real;
    a->e[i][i].imag += c->imag;
  }
}
// -----------------------------------------------------------------

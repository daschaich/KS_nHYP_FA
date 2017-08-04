// -----------------------------------------------------------------
// Add scalar multiplication on matrix
// c <-- c + s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_sum_su3_matrix(su3_matrix *b, Real s, su3_matrix *c) {
  register int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      c->e[i][j].real += s * b->e[i][j].real;
      c->e[i][j].imag += s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_add_su3_matrix(su3_matrix *a, su3_matrix *b, Real s,
                                su3_matrix *c) {

  register int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------

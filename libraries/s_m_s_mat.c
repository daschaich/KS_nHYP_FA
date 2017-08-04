// -----------------------------------------------------------------
// Subtract scalar multiplication on matrix
// c <-- a - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_sub_matrix(matrix *a, matrix *b, Real s,
                                matrix *c) {

  register int i,j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      c->e[i][j].real = a->e[i][j].real - s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag - s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------

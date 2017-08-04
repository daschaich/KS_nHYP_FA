// -----------------------------------------------------------------
// Scalar multiplication on vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_vector(vector *a, Real s, vector *b) {
  register int i;
  for (i = 0; i < 3; i++) {
    b->c[i].real = s * a->c[i].real;
    b->c[i].imag = s * a->c[i].imag;
  }
}
// -----------------------------------------------------------------

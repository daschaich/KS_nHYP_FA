// -----------------------------------------------------------------
// Add scalar multiplication on vector
// a <-- a + s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_sum_vector(vector *a, vector *b, Real s) {
  register int i;
  for (i = 0; i < 3; i++) {
    a->c[i].real += s * b->c[i].real;
    a->c[i].imag += s * b->c[i].imag;
  }
}

void scalar_mult_add_vector(vector *a, vector *b, Real s,
                                vector *c) {

  register int i;
  for (i = 0; i < 3; i++) {
    c->c[i].real = a->c[i].real + s * b->c[i].real;
    c->c[i].imag = a->c[i].imag + s * b->c[i].imag;
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Subtract scalar multiplication on vector
// c <- a - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_sub_su3_vector(su3_vector *a, su3_vector *b, Real s,
                                su3_vector *c) {

  register int i;
  for (i = 0; i < 3; i++) {
    c->c[i].real = a->c[i].real - s * b->c[i].real;
    c->c[i].imag = a->c[i].imag - s * b->c[i].imag;
  }
}
// -----------------------------------------------------------------

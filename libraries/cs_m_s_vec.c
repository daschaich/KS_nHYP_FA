// -----------------------------------------------------------------
// Subtract complex scalar multiplication on vector
// Note that the order of arguments is the opposite of the usual
// c <-- c - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_dif_vec(vector *c, complex *s, vector *b) {
  register int i;
  for (i = 0; i < 3; i++) {
    c->c[i].real -= b->c[i].real * s->real - b->c[i].imag * s->imag;
    c->c[i].imag -= b->c[i].imag * s->real + b->c[i].real * s->imag;
  }
}
// -----------------------------------------------------------------

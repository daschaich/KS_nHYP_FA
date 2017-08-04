// -----------------------------------------------------------------
// Subtract two vectors -- unlike CDIF, output is always last
// c <-- c - b
// c <-- a - b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void dif_vector(vector *b, vector *c) {
  register int i;
  for (i = 0; i < 3; i++) {
    c->c[i].real -= b->c[i].real;
    c->c[i].imag -= b->c[i].imag;
  }
}

void sub_vector(vector *a, vector *b, vector *c) {
  register int i;
  for (i = 0; i < 3; i++) {
    c->c[i].real = a->c[i].real - b->c[i].real;
    c->c[i].imag = a->c[i].imag - b->c[i].imag;
  }
}
// -----------------------------------------------------------------

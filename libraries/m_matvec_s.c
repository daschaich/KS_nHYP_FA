// -----------------------------------------------------------------
// Accumulate matrix--vector multiplication
// c <-- c + a.b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su3_mat_vec_sum(su3_matrix *a, su3_vector *b, su3_vector *c) {
  register int i;
#ifndef FAST
  register int j;
  register complex x;

  for (i = 0; i < 3; i++){
    CMUL(a->e[i][0], b->c[0], x)
    for (j = 1; j < 3; j++)
      CMULSUM(a->e[i][j], b->c[j], x)

    c->c[i].real += x.real;
    c->c[i].imag += x.imag;
  }

#else
  register Real t, ar, ai, br, bi, cr, ci;

  for (i = 0; i < 3; i++) {
    ar = a->e[i][0].real;
    ai = a->e[i][0].imag;
    br = b->c[0].real;
    bi = b->c[0].imag;
    cr = ar * br;
    t = ai * bi; cr -= t;
    ci = ar * bi;
    t = ai * br; ci += t;

    ar = a->e[i][1].real;
    ai = a->e[i][1].imag;
    br = b->c[1].real;
    bi = b->c[1].imag;
    t = ar * br; cr += t;
    t = ai * bi; cr -= t;
    t = ar * bi; ci += t;
    t = ai * br; ci += t;

    ar = a->e[i][2].real;
    ai = a->e[i][2].imag;
    br = b->c[2].real;
    bi = b->c[2].imag;
    t = ar * br; cr += t;
    t = ai * bi; cr -= t;
    t = ar * bi; ci += t;
    t = ai * br; ci += t;

    c->c[i].real += cr;
    c->c[i].imag += ci;
  }
#endif
}
// -----------------------------------------------------------------

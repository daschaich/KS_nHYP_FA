// -----------------------------------------------------------------
// Four adjoint-matrix--vector multiplications, with a single vector
// and four matrices in single array
// c[i] <- adag[i].b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_adj_su3_mat_vec_4dir(su3_matrix *a, su3_vector *b, su3_vector *c) {
#ifndef FAST
  mult_adj_su3_mat_vec(a, b, c);
  mult_adj_su3_mat_vec(a + 1, b, c + 1);
  mult_adj_su3_mat_vec(a + 2, b, c + 2);
  mult_adj_su3_mat_vec(a + 3, b, c + 3);

#else
  register int n;
  register Real c0r, c0i, c1r, c1i, c2r, c2i, br, bi, a0, a1, a2;
  register su3_matrix *mat;
  register su3_vector *vec, *dest;

  mat = a;
  vec = b;
  dest = c;
  for (n = 0; n < 4; n++) {
    br = vec->c[0].real;
    bi = vec->c[0].imag;
    a0 = mat->e[0][0].real;
    a1 = mat->e[0][1].real;
    a2 = mat->e[0][2].real;
    c0r = a0 * br;
    c1r = a1 * br;
    c2r = a2 * br;
    c0i = a0 * bi;
    c1i = a1 * bi;
    c2i = a2 * bi;

    a0 = mat->e[0][0].imag;
    a1 = mat->e[0][1].imag;
    a2 = mat->e[0][2].imag;
    c0r += a0*bi;
    c1r += a1*bi;
    c2r += a2*bi;
    c0i -= a0*br;
    c1i -= a1*br;
    c2i -= a2*br;

    br = vec->c[1].real;
    bi = vec->c[1].imag;
    a0 = mat->e[1][0].real;
    a1 = mat->e[1][1].real;
    a2 = mat->e[1][2].real;
    c0r += a0 * br;
    c1r += a1 * br;
    c2r += a2 * br;
    c0i += a0 * bi;
    c1i += a1 * bi;
    c2i += a2 * bi;

    a0 = mat->e[1][0].imag;
    a1 = mat->e[1][1].imag;
    a2 = mat->e[1][2].imag;
    c0r += a0 * bi;
    c1r += a1 * bi;
    c2r += a2 * bi;
    c0i -= a0 * br;
    c1i -= a1 * br;
    c2i -= a2 * br;

    br = vec->c[2].real;
    bi = vec->c[2].imag;
    a0 = mat->e[2][0].real;
    a1 = mat->e[2][1].real;
    a2 = mat->e[2][2].real;
    c0r += a0 * br;
    c1r += a1 * br;
    c2r += a2 * br;
    c0i += a0 * bi;
    c1i += a1 * bi;
    c2i += a2 * bi;

    a0 = mat->e[2][0].imag;
    a1 = mat->e[2][1].imag;
    a2 = mat->e[2][2].imag;
    c0r += a0 * bi;
    c1r += a1 * bi;
    c2r += a2 * bi;
    c0i -= a0 * br;
    c1i -= a1 * br;
    c2i -= a2 * br;

    dest->c[0].real = c0r;
    dest->c[0].imag = c0i;
    dest->c[1].real = c1r;
    dest->c[1].imag = c1i;
    dest->c[2].real = c2r;
    dest->c[2].imag = c2i;

    mat++;
    dest++;
  }
#endif
}
// -----------------------------------------------------------------

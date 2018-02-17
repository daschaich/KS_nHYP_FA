// -----------------------------------------------------------------
// Compute and compress the traceless anti-hermitian part of a matrix
// dest = 0.5 * (src - src^dag) - Tr[0.5 * (src - src^dag)]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void make_anti_hermitian(matrix *src, anti_hermitmat *dest) {
  Real tr = src->e[0][0].imag + src->e[1][1].imag + src->e[2][2].imag;
  tr *= 0.33333333333333333;

  dest->m00im = src->e[0][0].imag - tr;
  dest->m11im = src->e[1][1].imag - tr;
  dest->m22im = src->e[2][2].imag - tr;
  dest->m01.real = 0.5 * (src->e[0][1].real - src->e[1][0].real);
  dest->m02.real = 0.5 * (src->e[0][2].real - src->e[2][0].real);
  dest->m12.real = 0.5 * (src->e[1][2].real - src->e[2][1].real);
  dest->m01.imag = 0.5 * (src->e[0][1].imag + src->e[1][0].imag);
  dest->m02.imag = 0.5 * (src->e[0][2].imag + src->e[2][0].imag);
  dest->m12.imag = 0.5 * (src->e[1][2].imag + src->e[2][1].imag);
}
// -----------------------------------------------------------------

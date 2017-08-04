// -----------------------------------------------------------------
// Return complex trace of matrix products adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex complextrace_su3( matrix *a, matrix *b ) {
  register int i, j;
  register Real sumr = 0.0, sumi = 0.0;
  complex sum;
  for (i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
    sumr += a->e[i][j].real * b->e[i][j].real
          + a->e[i][j].imag * b->e[i][j].imag;
    sumi += a->e[i][j].real * b->e[i][j].imag
          - a->e[i][j].imag * b->e[i][j].real;
    }
  }
  sum.real = sumr;
  sum.imag = sumi;
  return sum;
}
// -----------------------------------------------------------------

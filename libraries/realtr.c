// -----------------------------------------------------------------
// Return real trace of matrix products a * b and adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real realtrace_nn(matrix *a, matrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      sum += a->e[i][j].real * b->e[j][i].real
           - a->e[i][j].imag * b->e[j][i].imag;
    }
  }
  return sum;
}

Real realtrace(matrix *a, matrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      sum += a->e[i][j].real * b->e[i][j].real
           + a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}
// -----------------------------------------------------------------

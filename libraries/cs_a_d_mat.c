// -----------------------------------------------------------------
// Add complex scalar to each diagonal element of matrix
// a <-- a + s * I
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_add_diag(matrix *a, complex *s) {
  register int i;

  for (i = 0; i < 3; i++) {
    a->e[i][i].real += s->real;
    a->e[i][i].imag += s->imag;
  }
}
// -----------------------------------------------------------------

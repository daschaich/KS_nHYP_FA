// -----------------------------------------------------------------
// Clear the given matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clear_mat(matrix *dest) {
  register int i,j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dest->e[i][j].real = 0.0;
      dest->e[i][j].imag = 0.0;
    }
  }
}
// -----------------------------------------------------------------

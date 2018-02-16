// -----------------------------------------------------------------
// Print the given matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumpmat(matrix *m) {
  int i, j;
  for (i = 0; i < 3; i++){
    for (j = 0; j < 3; j++)
      printf("  (%.4g, %.4g)", m->e[i][j].real, m->e[i][j].imag);
    printf("\n");
  }
  printf("\n");
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Return complex trace of matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex trace_su3(su3_matrix *a) {
  register complex tc;
  CADD(a->e[0][0], a->e[1][1], tc);
  CSUM(tc, a->e[2][2]);
  return tc;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Complex trace of the given matrix
// Return Tr[m]
// c <-- c + Tr[m]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex trace(matrix *m) {
  register complex tc;
  CADD(m->e[0][0], m->e[1][1], tc);
  CSUM(tc, m->e[2][2]);
  return tc;
}

void trace_sum(matrix *m, complex *c) {
  CSUM(*c, m->e[0][0]);
  CSUM(*c, m->e[1][1]);
  CSUM(*c, m->e[2][2]);
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Return real part of dot product of two vectors
// ReTr[adag.b]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real su3_rdot(vector *a, vector *b) {
  register Real tr, sum;
  sum = a->c[0].real * b->c[0].real;
  tr = a->c[0].imag * b->c[0].imag; sum += tr;
  tr = a->c[1].real * b->c[1].real; sum += tr;
  tr = a->c[1].imag * b->c[1].imag; sum += tr;
  tr = a->c[2].real * b->c[2].real; sum += tr;
  tr = a->c[2].imag * b->c[2].imag; sum += tr;
  return sum;
}
// -----------------------------------------------------------------

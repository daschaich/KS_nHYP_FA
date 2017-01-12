// -----------------------------------------------------------------
// Return squared magnitude of SUN vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real magsq_su3vec(su3_vector *a) {
  register Real sum = 0.0;
#ifndef FAST
  register int i;
  for (i = 0; i < DIMF; i++)
    sum += a->c[i].real * a->c[i].real + a->c[i].imag * a->c[i].imag;

#else
  register Real temp;
  temp = a->c[0].real*a->c[0].real; sum += temp;
  temp = a->c[0].imag*a->c[0].imag; sum += temp;
  temp = a->c[1].real*a->c[1].real; sum += temp;
  temp = a->c[1].imag*a->c[1].imag; sum += temp;
  temp = a->c[2].real*a->c[2].real; sum += temp;
  temp = a->c[2].imag*a->c[2].imag; sum += temp;
#endif
  return sum;
}
// -----------------------------------------------------------------

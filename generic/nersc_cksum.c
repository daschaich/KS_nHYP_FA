// -----------------------------------------------------------------
// Utility used in the NERSC archive format
// Compute the low order 32 bits of the unsigned integer sum of the
// float precision real and complex parts of the elements of the gauge matrices
#include "generic_includes.h"

u_int32type nersc_cksum() {
  int i, mu, a, b;
  site *s;
  u_int32type p32, chksum = 0;
  union {
    float       flt;
    u_int32type p32;
  } tmp;

  FORALLSITES(i, s) {
    for (mu = 0; mu < 4; ++mu) {
      for (a = 0; a < 2; a++) {
        for (b = 0; b < 3; b++) {
          tmp.flt = s->link[mu].e[a][b].real;
          p32 = tmp.p32;
          chksum += p32;
          tmp.flt = s->link[mu].e[a][b].imag;
          p32 = tmp.p32;
          chksum += p32;
        }
      }
    }
  }
  g_uint32sum(&chksum);
  return chksum;
}
// -----------------------------------------------------------------

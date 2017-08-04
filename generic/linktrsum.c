// -----------------------------------------------------------------
// Compute the mean global sum of the trace of the gauge links
// Use to check lattice file integrity
#include "generic_includes.h"

void linktrsum(double_complex *linktr) {
  int i, dir;
  site *s;
  matrix *a;

  linktr->real = 0.0;
  linktr->imag = 0.0;
  FORALLSITES(i, s) {
    FORALLUPDIR(dir) {
      a = &s->link[dir];
      CSUM(*linktr, a->e[0][0]);
      CSUM(*linktr, a->e[1][1]);
      CSUM(*linktr, a->e[2][2]);
    }
  }
  g_dcomplexsum(linktr);
  CDIVREAL(*linktr, (4.0 * volume), *linktr);
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Produce gaussian random momenta for the gauge fields
#include "generic_includes.h"
#include <defines.h>          // For SITERAND

void ranmom() {
  register int i, dir;
  register site *s;
  FORALLSITES(i,s) {
    FORALLUPDIR(dir) {
#ifdef SITERAND
      random_anti_hermitian((anti_hermitmat *)&(s->mom[dir]), &(s->site_prn));
#else
      random_anti_hermitian((anti_hermitmat *)&(s->mom[dir]), &node_prn);
#endif
    }
  }
}
// -----------------------------------------------------------------

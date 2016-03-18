// -----------------------------------------------------------------
// Construct a gaussian random vector g_rand, computed at all sites
// Stripped everything we don't use (ferm_links_t, z2rsource stuff)

#include "eig_includes.h"

void grsource_imp() {
  register int i, j;
  register site *s;
  FORALLSITES(i, s) {
    for (j = 0; j < 3; j++) {
#ifdef SITERAND
      s->g_rand.c[j].real = gaussian_rand_no(&(s->site_prn));
      s->g_rand.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
      s->g_rand.c[j].real = gaussian_rand_no(&node_prn);
      s->g_rand.c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }
  }
}
// -----------------------------------------------------------------

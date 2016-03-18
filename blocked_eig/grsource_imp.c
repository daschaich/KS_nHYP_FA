// -----------------------------------------------------------------
// Gaussian random vector routines -- one for eig, one for pbp
#include "block_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// g_rand must always be computed at all sites
// It is never used directly, so we can overwrite it for each level
// The parity is where dest is computed: EVEN, ODD or EVENANDODD
// For single pseudofermion
//   dest = M^dag g_rand
void grsource_imp(field_offset dest, Real M, int parity) {
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
  // Hit g_rand with M^dag
  dslash(F_OFFSET(g_rand), dest, parity);
  scalar_mult_latvec(dest, -1.0, dest, parity);
  scalar_mult_add_latvec(dest, F_OFFSET(g_rand), 2.0 * M, dest, parity);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct a gaussian random vector g_rand, computed at all sites
// Stripped everything we don't use (ferm_links_t, z2rsource stuff)
void grsource_eig() {
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

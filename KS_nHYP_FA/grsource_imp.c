// -----------------------------------------------------------------
// Construct a gaussian random vector g_rand, and either
//   dest = Mdag g_rand
//   dest = (M_1^dag)^{-1} M_0^dag g_rand
// depending on which Hasenbusch level we are dealing with
// Stripped everything we don't use (ferm_links_t, z2rsource stuff)

// g_rand must always be computed at all sites
// It is never used directly, so we can overwrite it for each level
// The parity is where dest is computed: EVEN, ODD or EVENANDODD
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// For single pseudofermion or inner Hasenbusch pseudofermion
// dest = M^dag g_rand
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
// For outer pseudofermion with Hasenbusch preconditioning
// dest_0 = (M_1^dag)^{-1} M_0^dag g_rand
// EVENANDODD seems required except for final M_1 application
void grsource_Hasen(field_offset dest, int *iters, int parity) {
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
  // Hit g_rand with M_0^dag (small mass), store in ttt[0][0]
  dslash(F_OFFSET(g_rand), F_OFFSET(ttt[0][0]), EVENANDODD);
  scalar_mult_latvec(F_OFFSET(ttt[0][0]), -1.0,
                     F_OFFSET(ttt[0][0]), EVENANDODD);
  scalar_mult_add_latvec(F_OFFSET(ttt[0][0]), F_OFFSET(g_rand), 2.0 * mass,
                         F_OFFSET(ttt[0][0]), EVENANDODD);

  // Hit ttt[0][0] with (M_1^dag M_1)^{-1} (large mass), store in ttt[0][1]
  clear_latvec(F_OFFSET(ttt[0][1]), EVENANDODD);
  *iters += ks_congrad(F_OFFSET(ttt[0][0]),
                       F_OFFSET(ttt[0][1]), MH, EVENANDODD);

  // Hit ttt[0][1] with M_1 (large mass)
  // so that only (M_1^dag){-1} ttt[0][1] remains
  dslash(F_OFFSET(ttt[0][1]), dest, parity);
  scalar_mult_add_latvec(dest, F_OFFSET(ttt[0][1]), 2.0 * MH, dest, parity);
  clear_latvec(F_OFFSET(ttt[0][0]), EVENANDODD);
  clear_latvec(F_OFFSET(ttt[0][1]), EVENANDODD);
}
// -----------------------------------------------------------------

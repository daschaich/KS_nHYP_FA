// -----------------------------------------------------------------
// Matrix inversion for general even plus odd quark actions
// dslash and ks_congrad made more v6-like
// ferm_links_t, relresid and unused functions stripped out

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"
#include "../KS_nHYP_FA/ks_dyn_includes.h"    // For ks_congrad
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest <-- M^(-1) src
int mat_invert_uml(field_offset src, field_offset dest,
                   field_offset temp, Real mass) {

  register int i;
  register site *s;
  int cgn;
  Real norm = 1.0 / (2.0 * mass);

  if (src == temp) {
    printf("MAT_INVERT_UML: src = temp\n");
    exit(0);
  }

  // "Precondition" both even and odd sites
  // temp <-- M^dag src
  dslash(src, F_OFFSET(ttt), EVENANDODD);
  scalar_mult_add_latvec(F_OFFSET(ttt), src,
                         -2.0 * mass, temp, EVENANDODD);
  scalar_mult_latvec(temp, -1.0, temp, EVENANDODD);

  // dest_e <-- (M^dag M)^-1 temp_e  (even sites only)
  cgn = ks_congrad(temp, dest, mass, EVEN);

  // Reconstruct odd site solution
  // dest_o <-- 1/2m (Dslash_oe * dest_e + src_o)
  dslash(dest, F_OFFSET(ttt), ODD);
  FORODDSITES(i, s) {
    sub_su3_vector((su3_vector *)F_PT(s, src), &(s->ttt),
                   (su3_vector *)F_PT(s, dest));
    scalar_mult_su3_vector((su3_vector *)F_PT(s, dest), norm,
                           (su3_vector *)F_PT(s, dest));
  }

  // Polish off odd sites to correct for possible roundoff error
  // dest_o <-- (M^dag M)^-1 temp_o (odd sites only)
  cgn += ks_congrad(temp, dest, mass, ODD);

  return cgn;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Matrix inversion for general even plus odd quark actions
// dslash and ks_congrad made more v6-like
// ferm_links_t, relresid and unused functions stripped out

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"
#include "../KS_nHYP_FA/ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int mat_invert_uml(field_offset src, field_offset dest,
                   field_offset temp, Real mass) {

  int cgn;
  register int i;
  register site *s;

  if (src == temp) {
    printf("MAT_INVERT_UML: src = temp\n");
    exit(0);
  }

  // "Precondition" both even and odd sites
  // temp <- M^dag * src
  dslash(src, F_OFFSET(ttt1[0]), EVENANDODD);
  scalar_mult_add_latvec(F_OFFSET(ttt1[0]), src,
                         -2.0 * mass, temp, EVENANDODD);
  scalar_mult_latvec(temp, -1.0, temp, EVENANDODD);

  // dest_e <- (M_adj M)^-1 temp_e  (even sites only)
  cgn = ks_congrad(temp, dest, mass, EVEN);

  // Reconstruct odd site solution
  // dest_o <- 1/2m (Dslash_oe * dest_e + src_o)
  dslash(dest, F_OFFSET(ttt1[0]), ODD);
  FORODDSITES(i, s) {
    sub_su3_vector((su3_vector *)F_PT(s, src), &(s->ttt1[0]),
                   (su3_vector *)F_PT(s, dest) );
    scalar_mult_su3_vector((su3_vector *)F_PT(s, dest),
                           1.0 / (2.0 * mass),
                           (su3_vector *)F_PT(s, dest));
  }

  // Polish off odd sites to correct for possible roundoff error
  // dest_o <- (M_adj M)^-1 temp_o  (odd sites only)
  cgn += ks_congrad(temp, dest, mass, ODD);

  return cgn;
}
// -----------------------------------------------------------------

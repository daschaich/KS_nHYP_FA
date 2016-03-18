// -----------------------------------------------------------------
// Helper functions copied from generic_ks/d_congrad5.c and
// generic_ks/dslash.c

#include "block_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND   // For loopend.h
#include "../include/loopend.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clear an su3_vector in the lattice
void clear_latvec(field_offset v, int parity) {
  register int i, j;
  register site *s;
  register su3_vector *vv;

  switch(parity) {
    case EVEN:
      FOREVENSITES(i, s) {
        vv = (su3_vector *)F_PT(s, v);
        for (j = 0; j < 3; j++) {
          vv->c[j].real = 0;
          vv->c[j].imag = 0;
        }
      }
      break;

    case ODD:
      FORODDSITES(i, s) {
        vv = (su3_vector *)F_PT(s, v);
        for(j = 0; j < 3; j++) {
          vv->c[j].real = 0;
          vv->c[j].imag = 0;
        }
      }
      break;

    case EVENANDODD:
      FORALLSITES(i, s) {
        vv = (su3_vector *)F_PT(s, v);
        for(j = 0; j < 3; j++) {
          vv->c[j].real = 0;
          vv->c[j].imag = 0;
        }
      }
      break;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy an su3_vector in the lattice
void copy_latvec(field_offset src, field_offset dest, int parity) {
  register int i;
  register site *s;
  register su3_vector *spt, *dpt;

  switch(parity) {
    case EVEN:
      FOREVENSITESDOMAIN(i, s) {
        s = &(lattice[i]);
        spt = (su3_vector *)F_PT(s, src);
        dpt = (su3_vector *)F_PT(s, dest);
        *dpt = *spt;
      }
      break;
    case ODD:
      FORODDSITESDOMAIN(i, s) {
        s = &(lattice[i]);
        spt = (su3_vector *)F_PT(s, src);
        dpt = (su3_vector *)F_PT(s, dest);
        *dpt = *spt;
      }
      break;
    case EVENANDODD:
      FORALLSITESDOMAIN(i, s) {
        s = &(lattice[i]);
        spt = (su3_vector *)F_PT(s, src);
        dpt = (su3_vector *)F_PT(s, dest);
        *dpt = *spt;
      }
      break;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar multiply an SU3 vector in the lattice
void scalar_mult_latvec(field_offset src, Real scalar,
                        field_offset dest, int parity) {

  register int i;
  register site *s;
  register su3_vector *spt, *dpt;

  switch(parity) {
    case EVEN:
      FOREVENSITES(i, s) {
        spt = (su3_vector *)F_PT(s, src);
        dpt = (su3_vector *)F_PT(s, dest);
        scalar_mult_su3_vector(spt, scalar, dpt);
      }
      break;
    case ODD:
      FORODDSITES(i, s) {
        spt = (su3_vector *)F_PT(s, src);
        dpt = (su3_vector *)F_PT(s, dest);
        scalar_mult_su3_vector(spt, scalar, dpt);
      }
      break;
    case EVENANDODD:
      FORALLSITES(i, s) {
        spt = (su3_vector *)F_PT(s, src);
        dpt = (su3_vector *)F_PT(s, dest);
        scalar_mult_su3_vector(spt, scalar, dpt);
      }
      break;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar multiply and add an SU3 vector in the lattice
// dest = src1 + scalar * src2
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
                            Real scalar, field_offset dest, int parity) {

  register int i;
  register site *s;
  register su3_vector *spt1, *spt2, *dpt;

  FORSOMEPARITY(i, s, parity) {
    spt1 = (su3_vector *)F_PT(s, src1);
    spt2 = (su3_vector *)F_PT(s, src2);
    dpt = (su3_vector *)F_PT(s, dest);
    if (i < loopend - FETCH_UP)
      prefetch_VVV((su3_vector *)F_PT((s+FETCH_UP), src1),
                   (su3_vector *)F_PT((s+FETCH_UP), src2),
                   (su3_vector *)F_PT((s+FETCH_UP), dest));

    scalar_mult_add_su3_vector(spt1, spt2, scalar, dpt);
  } END_LOOP
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Convenience functions to clean up all gathers
static void cleanup_one_gather_set(msg_tag *tags[]) {
  int i;
  for(i = XUP; i <= TUP; i++) {
    cleanup_gather(tags[i]);
    cleanup_gather(tags[OPP_DIR(i)]);
  }
}

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]) {
  cleanup_one_gather_set(tags1);
  cleanup_one_gather_set(tags2);
}
// -----------------------------------------------------------------

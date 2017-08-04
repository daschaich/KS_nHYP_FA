// -----------------------------------------------------------------
// Helper functions copied from generic_ks/d_congrad5.c and
// generic_ks/dslash.c

#include "ks_dyn_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND   // For loopend.h
#include "../include/loopend.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clear an vector in the lattice
void clear_latvec(field_offset v, int parity) {
  register int i, j;
  register site *s;
  register vector *vv;

  switch(parity) {
    case EVEN:
      FOREVENSITES(i, s) {
        vv = (vector *)F_PT(s, v);
        for (j = 0; j < 3; j++) {
          vv->c[j].real = 0;
          vv->c[j].imag = 0;
        }
      }
      break;

    case ODD:
      FORODDSITES(i, s) {
        vv = (vector *)F_PT(s, v);
        for(j = 0; j < 3; j++) {
          vv->c[j].real = 0;
          vv->c[j].imag = 0;
        }
      }
      break;

    case EVENANDODD:
      FORALLSITES(i, s) {
        vv = (vector *)F_PT(s, v);
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
// Copy an vector in the lattice
void copy_latvec(field_offset src, field_offset dest, int parity) {
  register int i;
  register site *s;
  register vector *spt, *dpt;

  switch(parity) {
    case EVEN:
      FOREVENSITESDOMAIN(i, s) {
        s = &(lattice[i]);
        spt = (vector *)F_PT(s, src);
        dpt = (vector *)F_PT(s, dest);
        *dpt = *spt;
      }
      break;
    case ODD:
      FORODDSITESDOMAIN(i, s) {
        s = &(lattice[i]);
        spt = (vector *)F_PT(s, src);
        dpt = (vector *)F_PT(s, dest);
        *dpt = *spt;
      }
      break;
    case EVENANDODD:
      FORALLSITESDOMAIN(i, s) {
        s = &(lattice[i]);
        spt = (vector *)F_PT(s, src);
        dpt = (vector *)F_PT(s, dest);
        *dpt = *spt;
      }
      break;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar multiply a vector in the lattice
void scalar_mult_latvec(field_offset src, Real scalar,
                        field_offset dest, int parity) {

  register int i;
  register site *s;
  register vector *spt, *dpt;

  switch(parity) {
    case EVEN:
      FOREVENSITES(i, s) {
        spt = (vector *)F_PT(s, src);
        dpt = (vector *)F_PT(s, dest);
        scalar_mult_vector(spt, scalar, dpt);
      }
      break;
    case ODD:
      FORODDSITES(i, s) {
        spt = (vector *)F_PT(s, src);
        dpt = (vector *)F_PT(s, dest);
        scalar_mult_vector(spt, scalar, dpt);
      }
      break;
    case EVENANDODD:
      FORALLSITES(i, s) {
        spt = (vector *)F_PT(s, src);
        dpt = (vector *)F_PT(s, dest);
        scalar_mult_vector(spt, scalar, dpt);
      }
      break;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar multiply and add vectors in the lattice
// c <-- a + scalar * b
void scalar_mult_add_latvec(field_offset a, field_offset b,
                            Real scalar, field_offset c, int parity) {

  register int i;
  register site *s;
  register vector *apt, *bpt, *cpt;

  FORSOMEPARITY(i, s, parity) {
    apt = (vector *)F_PT(s, a);
    bpt = (vector *)F_PT(s, b);
    cpt = (vector *)F_PT(s, c);
    if (i < loopend - FETCH_UP)
      prefetch_VVV((vector *)F_PT((s + FETCH_UP), a),
                   (vector *)F_PT((s + FETCH_UP), b),
                   (vector *)F_PT((s + FETCH_UP), c));

    scalar_mult_add_vector(apt, bpt, scalar, cpt);
  } END_LOOP
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Convenience functions to clean up all gathers
static void cleanup_one_gather_set(msg_tag *tags[]) {
  int i;
  FORALLUPDIR(i) {
    cleanup_gather(tags[i]);
    cleanup_gather(tags[OPP_DIR(i)]);
  }
}

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]) {
  cleanup_one_gather_set(tags1);
  cleanup_one_gather_set(tags2);
}
// -----------------------------------------------------------------

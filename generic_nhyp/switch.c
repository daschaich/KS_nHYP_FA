// -----------------------------------------------------------------
// Exchange the thin and smeared gauge fields.
// Only the s->link ever carries the staggered phases
// This is the only structure the fermions are supposed to talk to
#include "nhyp_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy original s->link to gauge_field_thin, block and then
// put gauge_field into s->link
void block_and_fatten() {
  register int i, dir;
  register site* s;

  // Take out the staggered phases before copying to gauge_field_thin
  rephase(OFF);
  for (dir = 0; dir < 4; dir++) {
    FORALLSITES(i, s)
      gauge_field_thin[dir][i] = s->link[dir];
  }

  block_nhyp();

  for (dir = 0; dir < 4; dir++) {
    FORALLSITES(i, s)
      s->link[dir] = gauge_field[dir][i];
  }

  // Now the s->link has the fat field
  // Restore the staggered phases
  rephase(ON);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Put pre-existing gauge_field into s->link, applying staggered phases
void just_fatten() {
  register int i, dir, j, k;
  register site* s;

  for (dir = 0; dir < 4; dir++) {
    FORALLSITES(i, s) {
      s->link[dir] = gauge_field[dir][i];
      for (j = 0; j < 3; j++) {
        for(k = 0; k < 3; k++) {
          s->link[dir].e[j][k].real *= s->phase[dir];
          s->link[dir].e[j][k].imag *= s->phase[dir];
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Put gauge_field_thin into s->link
// For consistency we take the phases out first,
// but that doesn't matter since it is overwritten
void leanlinks() {
  register int i, dir;
  register site* s;

  rephase(OFF);
  for (dir = 0; dir < 4; dir++) {
    FORALLSITES(i, s)
      s->link[dir] = gauge_field_thin[dir][i];
  }

  rephase(ON);
}
// -----------------------------------------------------------------

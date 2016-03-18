// -----------------------------------------------------------------
// Exchange the thin and smeared gauge fields.
// Only the s->link ever carries the staggered phases
// This is the only structure the fermions are supposed to talk to
#include "nhyp_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy original s->link to gauge_field_thin, block and then
// put gauge_field into s->link
// Don't include staggered phases
void unphased_block_and_fatten() {
  register int i,dir;
  register site* s;

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
}
// -----------------------------------------------------------------

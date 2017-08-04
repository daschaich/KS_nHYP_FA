// -----------------------------------------------------------------
// nHYP smear N times, saving the original thin links
#include "ks_dyn_includes.h"
void block_N_fatten(int N) {
  register int i, dir;
  register site* s;
  int j;

  block_and_fatten();
  if (N > 1) {
    // Save unsmeared links in gauge_field_save
    for (dir = 0; dir < 4; dir++) {
      FORALLSITES(i, s)
        mat_copy(gauge_field_thin[dir] + i, gauge_field_save[dir] + i);
    }

    // block_and_fatten starts by copying s->link into gauge_field_thin
    for (j = 1; j < N; j++)
      block_and_fatten();

    // Reset gauge_field_thin to hold unsmeared links
    for (dir = 0; dir < 4; dir++) {
      FORALLSITES(i, s)
        mat_copy(gauge_field_save[dir] + i, gauge_field_thin[dir] + i);
    }
  }
}
// -----------------------------------------------------------------

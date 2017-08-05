// -----------------------------------------------------------------
// Main procedure for HYP-smeared static potential
// Does nHYP smearing, transforms to axial gauge
// and computes (time-like) Wilson loops

// Includes and definitions
#define CONTROL
#include "hvy_qpot_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  register int i;
  register site *s;
  int prompt, dir;
  double dtime, ssplaq, stplaq;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Smear and check resulting plaquettes
  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      gauge_field_thin[dir][i] = s->link[dir];
  }
  block_nhyp();   // Smears gauge_field_thin into gauge_field
  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      s->link[dir] = gauge_field[dir][i];
  }
  plaquette(&ssplaq, &stplaq);
  node0_printf("Plaquettes after smearing: %.8g %.8g\n", ssplaq, stplaq);

  // Fix to axial gauge and check resulting plaquettes
  if (startflag != CONTINUE)
    ax_gauge();

  plaquette(&ssplaq, &stplaq);
  node0_printf("Plaquettes in axial gauge: %.8g %.8g\n", ssplaq, stplaq);

  // Compute on-axis time-like Wilson loops
  // The argument is officially the number of smearing steps, labels output
  w_loop1(1);

  // Compute off-axis time-like Wilson loops, if desired
  if (off_axis_flag == 1)
    w_loop2(1);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  return 0;
}
// -----------------------------------------------------------------

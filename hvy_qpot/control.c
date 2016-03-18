// -----------------------------------------------------------------
// Main procedure for SU(3) HYP-smeared static potential
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

  // Setup
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc,&argv);
  g_sync();
  prompt = setup();

  // Load input and run (loop removed)
  if (readin(prompt) == 0) {
    dtime = -dclock();

    // Smear and check resulting plaquettes
    for (dir = 0; dir < 4; dir++) {
      FORALLSITES(i, s)
        gauge_field_thin[dir][i] = s->link[dir];
    }
    block_nhyp();   // Smears gauge_field_thin into gauge_field
    for (dir = 0; dir < 4; dir++) {
      FORALLSITES(i, s)
        s->link[dir] = gauge_field[dir][i];
    }
    d_plaquette(&ssplaq, &stplaq);
    node0_printf("Plaquettes after smearing: %.8g %.8g\n", ssplaq, stplaq);

    // Fix to axial gauge and check resulting plaquettes
    if (startflag != CONTINUE)
      ax_gauge();

    d_plaquette(&ssplaq, &stplaq);
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
  } // readin(prompt) == 0
  return 0;
}
// -----------------------------------------------------------------

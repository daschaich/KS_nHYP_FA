// -----------------------------------------------------------------
// Main procedure for smearing
// Includes and definitions
#define CONTROL
#include "smearing_includes.h"

int main(int argc, char *argv[]) {
  int ismear, prompt;
  double ssplaq, stplaq, dtime;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Smear however many times requested, monitoring plaquette
  plaquette(&ssplaq, &stplaq);
  node0_printf("Plaquettes before smearing:  %.8g %.8g\n", ssplaq, stplaq);
  for (ismear = 0; ismear < nsmear; ismear++) {
    unphased_block_and_fatten();
    plaquette(&ssplaq, &stplaq);
    node0_printf("Plaquettes after smearing %d: %.8g %.8g\n",
                 ismear + 1, ssplaq, stplaq);
  }

  // Fix to Landau gauge if requested
  if (fixflag == LANDAU_GAUGE_FIX) {
    node0_printf("Fixing to Landau gauge...\n");
    double gtime = -dclock();
    gaugefix(8, 1.5, 5000, 1.0e-7, -1, -1, 0, NULL, NULL, 0, NULL, NULL);
    gtime += dclock();
    node0_printf("GFIX time = %.4g seconds\n", gtime);

    plaquette(&ssplaq, &stplaq);
    node0_printf("Plaquettes after gauge fix:  %.8g %.8g\n", ssplaq, stplaq);
  }
  else if (fixflag == NO_GAUGE_FIX) { // Braces suppress compiler warning
    node0_printf("Gauge fixing skipped\n");
  }
  else {
    node0_printf("ERROR: only LANDAU_GAUGE_FIX ");
    node0_printf("and NO_GAUGE_FIX supported\n");
    terminate(1);
  }

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  // Save lattice if requested -- no phases in links
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);

  normal_exit(0);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Main procedure for SU(3) S4b order parameters
#define CONTROL
#include "spectrum_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char **argv) {
  int iters, prompt;
  double dtime;

  // Setup
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc, &argv);
  g_sync();
  prompt = setup();

  // Load input and run
  // readin does a rephase(ON) on the read-in lattice
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Gauge fix
  // First argument picks Coulomb gauge
  // Next are relax_boost, max_iter, tolerance (set in defines.h),
  // diffmat, sumvec, nvector, vector_offset, vector_parity,
  // nantiherm, antiherm_offset and antiherm_parity
//  node0_printf("GFIX SKIPPED\n");
  rephase(OFF);
  gaugefix(TUP, 1.8, 500, GAUGE_FIX_TOL,
           F_OFFSET(tempmat1), F_OFFSET(tempvec[0]), 0, NULL, NULL,
           0, NULL, NULL);

  rephase(ON);

#ifndef PBOUND
  node0_printf("\nANTIPERIODIC BC in time direction\n");
#else
  node0_printf("\nPERIODIC BC in time direction\n");
#endif

  block_and_fatten();
  iters = fpi_2(&mass);
  iters += spectrum2(mass, F_OFFSET(chi), F_OFFSET(psi));
  iters += nl_spectrum(mass, F_OFFSET(chi), F_OFFSET(psi),
                       F_OFFSET(tempmat1), F_OFFSET(tempmat2));
  leanlinks();
  fflush(stdout);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n\n", iters);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------

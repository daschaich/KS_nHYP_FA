// -----------------------------------------------------------------
// Main procedure for SU(3) S4b order parameters
#define CONTROL
#include "S4b_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char **argv) {
  int avm_iters, prompt;
  double ssplaq, stplaq, act, dtime;
  complex plp = cmplx(99, 99), xplp = cmplx(99, 99);

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
  if (readin(prompt) == 0) {
    dtime = -dclock();

    // Take out the staggered phases to measure gauge observables
    rephase(OFF);
    d_plaquette(&ssplaq, &stplaq);
    plp = ploop(TUP);
    xplp = ploop(XUP);
    node0_printf("GMES %.8g %.8g 0 %.8g %.8g\n",    // No update step
                 plp.real, plp.imag, ssplaq, stplaq);
    node0_printf("POLYA %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, xplp.real, xplp.imag);

    // Action adds adjoint plaquette term
    d_plaquette_a(&ssplaq, &stplaq);
    act = (ssplaq + stplaq) / 2;
    node0_printf("ACT %.8g\n", act);
    rephase(ON);

#ifndef PBOUND
    node0_printf("\nANTIPERIODIC BC in time direction\n");
#else
    node0_printf("\nPERIODIC BC in time direction\n");
#endif

    // Even and odd plaquettes and links in various directions
    block_and_fatten();
    rephase(OFF);
    d_plaquette(&ssplaq, &stplaq);    // To check meas_plaq()
    node0_printf("Plaquettes after smearing: %.8g %.8g\n", ssplaq, stplaq);
    meas_plaq();
    rephase(ON);
    avm_iters = meas_link(F_OFFSET(chi), F_OFFSET(psi), mass);
    leanlinks();
    fflush(stdout);

    node0_printf("RUNNING COMPLETED\n");
    dtime += dclock();
    node0_printf("Time = %.4g seconds\n", dtime);
    node0_printf("total_iters = %d\n\n", avm_iters);
    fflush(stdout);
  } // readin(prompt) == 0
  return 0;
}
// -----------------------------------------------------------------

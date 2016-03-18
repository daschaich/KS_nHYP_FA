// -----------------------------------------------------------------
// Main procedure for nHYP-smeared staggered SU(3) evolution
// Includes and definitions
// For benchmarking, define CGTIME in defines.h
#define CONTROL
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int meascount = 0, traj_done, prompt;
  int s_iters = 0, avm_iters = 0, avs_iters = 0;
  Real f_eps0, f_eps1, g_eps;
  double ssplaq, stplaq, act, dtime;
  complex plp = cmplx(99, 99), xplp = cmplx(99, 99);

  // Set up
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc, &argv);
  g_sync();
  prompt = setup();

  // Load input and run (loop removed)
  // readin does a rephase(ON) on the read-in lattice
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  if (beta < 0)
    exit(1);

  node0_printf("\nFUNDAMENTAL--ADJOINT GAUGE ACTION ");
  node0_printf("WITH beta_a/beta_f COEFFICIENT %2.3f\n", beta_a);

  f_eps0 = traj_length / (Real)nsteps[0];
  f_eps1 = f_eps0 / (2 * (Real)nsteps[1]);
  g_eps = f_eps1 / (2 * (Real)nsteps[MAX_MASSES]);
  node0_printf("f_eps0 %.4g f_eps1 %.4g g_eps %.4g\n", f_eps0, f_eps1, g_eps);

  // Perform warmup trajectories
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();

  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories, reunitarizations and measurements
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    // Do the trajectories
    s_iters = update();
    avs_iters += s_iters;

    // Take out the staggered phases to measure gauge observables
    rephase(OFF);
    plaquette(&ssplaq, &stplaq);
    plp = ploop(TUP);
    xplp = ploop(XUP);
    node0_printf("GMES %.8g %.8g %d %.8g %.8g\n",
                 plp.real, plp.imag, s_iters, ssplaq, stplaq);
    node0_printf("POLYA %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, xplp.real, xplp.imag);

    // Action adds adjoint plaquette term
    plaquette_a(&ssplaq, &stplaq);
    act = (ssplaq + stplaq) / 2;
    node0_printf("ACT %.8g\n", act);
    rephase(ON);

    // Measure every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
#ifndef PBOUND
      node0_printf("\nANTIPERIODIC BC in time direction\n");
#else
      node0_printf("\nPERIODIC BC in time direction\n");
#endif

      meascount++;
      block_N_fatten(Nsmear);
      rephase(OFF);
      plaquette(&ssplaq, &stplaq);    // To check meas_plaq()
      node0_printf("Plaquettes after smearing: %.8g %.8g\n",
                   ssplaq, stplaq);
      meas_plaq();
      rephase(ON);
      avm_iters += meas_link(F_OFFSET(chi[0][0]),
                             F_OFFSET(psi[0][0]), mass);
      leanlinks();
      fflush(stdout);
    }
  } // End loop over trajectories

  node0_printf("RUNNING COMPLETED\n");
  node0_printf("Average CG iters for steps: %.4g\n",
               (double)avs_iters / trajecs);
  if (meascount > 0) {
    node0_printf("Average CG iters for measurements: %.4g\n",
                 (double)avm_iters / meascount);
  }
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n\n", total_iters);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET) {
    rephase(OFF);
    save_lattice(saveflag, savefile, stringLFN);
    rephase(ON);
  }
  normal_exit(0);
  return 0;
}
// -----------------------------------------------------------------

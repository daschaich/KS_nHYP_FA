// -----------------------------------------------------------------
// Main procedure for nHYP-smeared staggered SU(3) reversibility test
// Includes and definitions
// For benchmarking, define CGTIME in defines.h
#define CONTROL
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// A copy of update(), reversing instead of refreshing the momenta
// and not generating a new pseudofermion configuration
int reverse() {
  register int i, dir;
  register site *s;
  int iters = 0;
  Real cg_time[2], old_cg_time[2], next_cg_time[2];
  double startaction = 0, endaction, change;
  double gnorm = 0;
  double fnorm[2];
  fnorm[0] = 0;
  fnorm[1] = 0;
  old_cg_time[0] = -1e6;
  old_cg_time[1] = -1e6;
  cg_time[0] = -1e6;
  cg_time[1] = -1e6;

  // Reverse the momenta (anti_hermitmat defined in include/su3.h)
  FORALLSITES(i, s) {
    for (dir = XUP; dir <= TUP; dir++) {
      s->mom[dir].m01.real *= -1.0;
      s->mom[dir].m02.real *= -1.0;
      s->mom[dir].m12.real *= -1.0;
      s->mom[dir].m01.imag *= -1.0;
      s->mom[dir].m02.imag *= -1.0;
      s->mom[dir].m12.imag *= -1.0;
      s->mom[dir].m00im *= -1.0;
      s->mom[dir].m11im *= -1.0;
      s->mom[dir].m22im *= -1.0;
    }
  }

  // Do conjugate gradient to get (Madj M)^{-1} * chi
#ifdef CG_DEBUG
  node0_printf("UPDATE: calling congrad for level 0 (startaction)\n");
#endif
  doCG(0, mass, &iters);
  if (num_masses == 2) {
#ifdef CG_DEBUG
    node0_printf("UPDATE: calling congrad for level 1 (startaction)\n");
#endif
    doCG(1, MH, &iters);
  }

  cg_time[0] = 0;
  cg_time[1] = 0;

  // Save action and do microcanonical updating
  // No need to save links, since will automatically accept
  startaction = action();
  if (num_masses == 2)
    iters += update_step_three(old_cg_time, cg_time, next_cg_time, fnorm, &gnorm);
  else
    iters += update_step_two(old_cg_time, cg_time, next_cg_time, fnorm, &gnorm);

  // Do conjugate gradient to get (Madj M)^{-1} * chi
#ifdef CG_DEBUG
  node0_printf("UPDATE: calling congrad for level 0 (endaction)\n");
#endif
  doCG(0, mass, &iters);
  if (num_masses == 2) {
#ifdef CG_DEBUG
    node0_printf("UPDATE: calling congrad for level 1 (endaction)\n");
#endif
    doCG(1, MH, &iters);
  }

  // Find new action, always accept
  endaction = action();
  change = endaction - startaction;

  // Warn about overflow
  if (fabs((double)change) > 1e20) {
    node0_printf("WARNING: Correcting apparent overflow: Delta S = %.4g\n",
                 change);
    change = 1e20;
  }
  node0_printf("delta S = %.4g, start S = %.12g, end S = %.12g\n",
               change, startaction, endaction);

  node0_printf("IT_PER_TRAJ %d\n", iters);
  node0_printf("MONITOR_FORCE_GAUGE %.4g\n",
               gnorm / (double)(4 * nsteps[0] * nsteps[1]));
  node0_printf("MONITOR_FORCE_FERMION0 %.4g\n",
               fnorm[0] / (double)(2 * nsteps[0]));
  node0_printf("MONITOR_FORCE_FERMION1 %.4g\n",
               fnorm[1] / (double)(4 * nsteps[0] * nsteps[1]));

  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]){
  int prompt, s_iters = 0;
  double ssplaq, stplaq, startAct, endAct, dtime;
  complex plp = cmplx(99, 99), xplp = cmplx(99, 99);

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if(remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Always accept
#ifdef HMC_ALGORITHM
#undef HMC_ALGORITHM
#endif

  // Load input and run (NB loop removed)
  // readin does a rephase(ON) on the read-in lattice
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  node0_printf("\nFUNDAMENTAL--ADJOINT GAUGE ACTION ");
  node0_printf("WITH beta_a/beta_f COEFFICIENT %2.3f\n", beta_a);

  // Gauge observables at start of trajectory
  plaquette(&ssplaq, &stplaq);
  plp = ploop(TUP);
  xplp = ploop(XUP);
  // Flip signs to account for KS phases
  ssplaq *= -1;       stplaq *= -1;
  plp.real *= -1;     plp.imag *= -1;
  xplp.real *= -1;    xplp.imag *= -1;
  node0_printf("GMES %.8g %.8g %d %.8g %.8g\n",
               plp.real, plp.imag, s_iters, ssplaq, stplaq);
  node0_printf("POLYA %.8g %.8g %.8g %.8g\n",
               plp.real, plp.imag, xplp.real, xplp.imag);

  // Action at start of trajectory
  plaquette_a(&ssplaq, &stplaq);
  startAct = 0.5 * (ssplaq + stplaq);
  node0_printf("ACT %.8g\n", startAct);

  // Evolve forward for one trajectory
  s_iters = update();
  total_iters = s_iters;

  // Gauge observables at end of trajectory
  plaquette(&ssplaq, &stplaq);
  plp = ploop(TUP);
  xplp = ploop(XUP);
  // Flip signs to account for KS phases
  ssplaq *= -1;       stplaq *= -1;
  plp.real *= -1;     plp.imag *= -1;
  xplp.real *= -1;    xplp.imag *= -1;
  node0_printf("GMES %.8g %.8g %d %.8g %.8g\n",
               plp.real, plp.imag, s_iters, ssplaq, stplaq);
  node0_printf("POLYA %.8g %.8g %.8g %.8g\n",
               plp.real, plp.imag, xplp.real, xplp.imag);

  // Action at end of trajectory
  plaquette_a(&ssplaq, &stplaq);
  endAct = (ssplaq + stplaq) / 2;
  node0_printf("ACT %.8g\n", endAct);

  // Reverse momenta and evolve backwards for one trajectory
  s_iters = reverse();
  total_iters += s_iters;

  // Gauge observables hopefully back at the start of the trajectory
  plaquette(&ssplaq, &stplaq);
  plp = ploop(TUP);
  xplp = ploop(XUP);
  // Flip signs to account for KS phases
  ssplaq *= -1;       stplaq *= -1;
  plp.real *= -1;     plp.imag *= -1;
  xplp.real *= -1;    xplp.imag *= -1;
  node0_printf("GMES %.8g %.8g %d %.8g %.8g\n",
               plp.real, plp.imag, s_iters, ssplaq, stplaq);
  node0_printf("POLYA %.8g %.8g %.8g %.8g\n",
               plp.real, plp.imag, xplp.real, xplp.imag);

  // Action hopefully back at the start of the trajectory
  plaquette_a(&ssplaq, &stplaq);
  endAct = 0.5 * (ssplaq + stplaq);
  node0_printf("ACT %.8g\n", endAct);

  // Done
  node0_printf("RUNNING COMPLETED\n");
  node0_printf("Initial action: %.8g\n", startAct);
  node0_printf("Final action:   %.8g\n", endAct);
  node0_printf("Difference:     %.4g\n", fabs(endAct - startAct));
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------

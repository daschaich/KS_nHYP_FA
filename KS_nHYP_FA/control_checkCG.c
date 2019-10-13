// -----------------------------------------------------------------
// Main procedure for checking CG algorithm
// Corresponding check in v6/generic_ks/grsource_imp() suggests
// that parity must be EVENANDODD for this to work
// Includes and definitions
#define CONTROL
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Return maximum difference between any component of end and start
double check_latvec(field_offset end, field_offset start) {
  register int i;
  register site *s;
  int j;
  double maxDiff = 0, diff;
  vector *eptr, *sptr;

  FORALLSITES(i, s) {
    eptr = (vector *)F_PT(s, end);
    sptr = (vector *)F_PT(s, start);
    for (j = 0; j < 3; j++) {
      diff = fabs(eptr->c[j].real - sptr->c[j].real);
      if (diff > maxDiff)
        maxDiff = diff;

      diff = fabs(eptr->c[j].imag - sptr->c[j].imag);
      if (diff > maxDiff)
        maxDiff = diff;
    }
  }
  return maxDiff;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  register int i;
  register site *s;
  int j, prompt, total_iters = 0;
  double dtime, diff;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if(remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run
  // readin does a rephase(ON) on the read-in lattice
  if (readin(prompt) == 0) {
    dtime = -dclock();
    node0_printf("\nFUNDAMENTAL--ADJOINT GAUGE ACTION ");
    node0_printf("WITH beta_a/beta_f COEFFICIENT %2.3f\n", beta_a);
    clear_latvec(F_OFFSET(ttt[0][0]), EVENANDODD);
    clear_latvec(F_OFFSET(chi[0][0]), EVENANDODD);
    clear_latvec(F_OFFSET(psi[0][0]), EVENANDODD);
    block_N_fatten(Nsmear);

    // chi = M^dag ttt = (-Dslash + 2m) ttt
    // Gaussian random source
    grsource_imp(F_OFFSET(chi[0][0]), mass, EVENANDODD);
    copy_latvec(F_OFFSET(g_rand), F_OFFSET(ttt[0][0]), EVENANDODD);

//    // Point source, adding M^dag code from grsource_imp.c
//    FORALLSITES(i, s) {
//      if (i == 0)
//        s->ttt[0][0].c[0].real = 1;
//    }
//    dslash(F_OFFSET(ttt[0][0]), F_OFFSET(chi[0][0]), EVENANDODD);
//    scalar_mult_latvec(F_OFFSET(chi[0][0]), -1.0, F_OFFSET(chi[0][0]),
//                       EVENANDODD);
//    scalar_mult_add_latvec(F_OFFSET(chi[0][0]), F_OFFSET(ttt[0][0]),
//                           2.0 * mass, F_OFFSET(chi[0][0]), EVENANDODD);

    // Calculate psi = (M^dag.M)^(-1) M^dag ttt = M^(-1) ttt
    // M^dag M = -D^2 + 4m^2
    total_iters = ks_congrad(F_OFFSET(chi[0][0]), F_OFFSET(psi[0][0]),
                             mass, EVENANDODD);

    // Now hit M psi = (Dslash + 2m) psi = ttt
    dslash(F_OFFSET(psi[0][0]), F_OFFSET(chi[0][0]), EVENANDODD);
    scalar_mult_add_latvec(F_OFFSET(chi[0][0]),
                           F_OFFSET(psi[0][0]), 2.0 * mass,
                           F_OFFSET(chi[0][0]), EVENANDODD);

    // Check result against original source
    diff = check_latvec(F_OFFSET(chi[0][0]), F_OFFSET(ttt[0][0]));
    node0_printf("Max diff: %.4g\n", diff);
    node0_printf("ALL: \n");
    FORALLSITES(i, s) {
      for (j = 0; j < 3; j++) {\
        node0_printf("site %d color %d: (%.4g, %.4g)\n", i, j,
                     fabs(s->chi[0][0].c[j].real - s->ttt[0][0].c[j].real),
                     fabs(s->chi[0][0].c[j].imag - s->ttt[0][0].c[j].imag));
      }
    }

    // Done
    node0_printf("RUNNING COMPLETED\n");
    dtime += dclock();
    node0_printf("Time = %.4g seconds\n", dtime);
    node0_printf("total_iters = %d\n", total_iters);
    fflush(stdout);
  } // readin(prompt) == 0
  return 0;
}
// -----------------------------------------------------------------

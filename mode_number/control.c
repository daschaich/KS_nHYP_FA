// -----------------------------------------------------------------
// Main procedure for stochastic mode number calculation
// Includes and definitions
#define CONTROL
#include "mode_includes.h"
// -----------------------------------------------------------------


// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  register int i, dir, ismear, ipbp;
  register site *s;
  int prompt, count;
  Real ssplaq, stplaq;
  double dtime, nu, ave, err;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run
  // readin does a rephase(ON) on the read-in lattice
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Smear however many times requested
  // Print out plaquette each time (unsmeared plaquette printed by setup)
  for (ismear = 0; ismear < Nsmear; ismear++) {
    block_and_fatten();   // Smears gauge_field_thin into gauge_field
    FORALLUPDIR(dir) {
      FORALLSITES(i, s)
        gauge_field_thin[dir][i] = gauge_field[dir][i];
    }

    rephase(OFF);
    plaquette(&ssplaq, &stplaq);
    rephase(ON);
    node0_printf("Plaquettes after smearing %d: %.8g %.8g\n",
                 ismear + 1, ssplaq, stplaq);
  }
  node0_printf("\n");

  // Load coefficients of the polynomial
  // "star" is actually (Omega / Omega_*)^2, so take the square root
  coefficients();
  starSq = star;
  star = sqrt(starSq);
  node0_printf("Step function Norder %d epsilon %.4g delta %.4g star %.4g\n",
               Norder, epsilon, delta, star);

  // Generate and save stochastic sources
  for (ipbp = 0; ipbp < npbp; ipbp++) {
    grsource_plain();
    FORALLSITES(i, s)
      vec_copy(&(s->g_rand), &(source[ipbp][i]));
  }

  // NB: The Dirac operator normalization means that
  // the input lambda is twice the mass in the Dirac operator
  // So we use M = lambda/2 and print 2M=lambda
  for (count = 0; count < Npts; count++) {
    ave = 0;
    err = 0;
    M += spacing / 2;         // !!! Note factor of 2

    for (ipbp = 0; ipbp < npbp; ipbp++) {
      // Hit gaussian random vector twice with step function
      FORALLSITES(i, s)
        vec_copy(&(source[ipbp][i]), &(s->g_rand));

      step(F_OFFSET(g_rand), F_OFFSET(chi));
      step(F_OFFSET(chi), F_OFFSET(R2));

      // Take the norm to obtain nu(omega)
      nu = 0;
      FOREVENSITES(i, s)
        nu += (double)magsq_vec(&(s->R2));
      g_doublesum(&nu);
      ave += nu;
      err += nu * nu;
      node0_printf("NU(%.4g): %.4g ( %d of %d )\n",
                   2 * M, nu, ipbp + 1, npbp);          // Note factor of 2
    }

    // Average measurements and standard error
    ave /= npbp;
    err /= npbp;
    err -= ave * ave;
    err = sqrt(err / (npbp - 1));
    node0_printf("nu(%.4g): %.6g %.4g ( ave over %d )\n",
                 2 * M, ave, err, npbp);                // Note factor of 2
  }
  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------

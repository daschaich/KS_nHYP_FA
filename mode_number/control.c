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

  // Smear however many times requested
  // Print out plaquette each time (unsmeared plaquette printed by setup)
  for (ismear = 0; ismear < nsmear; ismear++) {
    block_and_fatten();   // Smears gauge_field_thin into gauge_field
    for (dir = XUP; dir <= TUP; dir++) {
      FORALLSITES(i, s)
        gauge_field_thin[dir][i] = gauge_field[dir][i];
    }

    rephase(OFF);
    d_plaquette(&ssplaq, &stplaq);
    rephase(ON);
    node0_printf("Plaquettes after smearing %d: %.8g %.8g\n",
        ismear + 1, ssplaq, stplaq);
  }
  node0_printf("\n");

  // Load Norder and coefficients of the polynomial
  // We copy in "star" from Agostino's code
  // This is actually (Omega / Omega_*)^2, so we take the square root
  coefficients();
  starSq = star;
  star = sqrt(starSq);
  node0_printf("Step function Norder %d, epsilon %.4g, delta %.4g, ",
               Norder, epsilon, delta);
  node0_printf("star %.4g\n", star);

  // Generate and save stochastic sources
  source = (su3_vector **)malloc(npbp * sizeof(su3_vector *));
  for (ipbp = 0; ipbp < npbp; ipbp++) {
    grsource_plain();
    FIELD_ALLOC(source[ipbp], su3_vector);
    FORALLSITES(i, s)
      su3vec_copy(&(s->g_rand), &(source[ipbp][i]));
  }

  // NB: The Dirac operator normalization means that
  // the input lambda is twice the mass in the Dirac operator
  // So we use M = lambda/2 and print 2M=lambda
  for (count = 0; count < Npts; count++) {
    ave = 0;
    err = 0;
    M += spacing / 2;

    for (ipbp = 0; ipbp < npbp; ipbp++) {
      //grsource_plain();    // Construct gaussian random vector g_rand

      // Hit gaussian random vector twice with step function
      //step(F_OFFSET(rand[ipbp]), F_OFFSET(chi));
      FORALLSITES(i, s)
        su3vec_copy(&(source[ipbp][i]), &(s->g_rand));

      step(F_OFFSET(g_rand), F_OFFSET(chi));
      step(F_OFFSET(chi), F_OFFSET(R2));

      // Take the norm to produce nu(omega)
      nu = 0;
      FOREVENSITES(i, s)
        nu += (double)magsq_su3vec(&(s->R2));
      g_doublesum(&nu);
      ave += nu;
      err += nu * nu;
      node0_printf("NU(%.4g): %.4g ( %d of %d )\n", 2 * M, nu, ipbp + 1, npbp);
    }

    // Average measurements and standard error
    ave /= npbp;
    err /= npbp;
    err -= ave * ave;
    err = sqrt(err / (npbp - 1));
    node0_printf("nu(%.4g): %.6g %.4g ( ave over %d )\n", 2 * M, ave, err, npbp);
  }
  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
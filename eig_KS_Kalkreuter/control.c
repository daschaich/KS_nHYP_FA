// -----------------------------------------------------------------
// Main procedure for calculating staggered eigenvalues
// Includes and definitions
#define CONTROL
#include "eig_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char **argv) {
  register int i, dir, ismear, ivec;
  register site *s;
  int prompt, total_iters;
  Real ssplaq, stplaq;
  double dtime, chirality;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run
  // readin does a rephase(ON) on the read-in lattice
  if (readin(prompt) == 0) {
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

    // Allocate eigenvectors
    eigVal = (double *)malloc(Nvecs * sizeof(double));
    eigVec = (su3_vector **)malloc(Nvecs * sizeof(su3_vector *));
    for (i = 0; i < Nvecs; i++)
      eigVec[i] = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));

    // Hard-code EVEN parity in Kalkreuter
    total_iters = Kalkreuter(eigVec, eigVal, eig_tol, error_decr,
                             Nvecs, maxIter, restart, kiters);

    for (ivec = 0; ivec < Nvecs; ivec++) {
      // Construct the odd part of the vector
      //   i Dslash Psi / sqrt(eigVal)
      // Ignore the i factor, which is irrelevant for the chirality
      // Halve EVENANDODD result since the EVEN vector is normalized to 1
      FORALLSITES(i, s)
        s->chi = eigVec[ivec][i];

      dslash(F_OFFSET(chi), F_OFFSET(psi), ODD);
      FORODDSITES(i, s)
        scalar_mult_su3_vector(&(s->psi), 1.0 / sqrt(eigVal[ivec]),
                               &(eigVec[ivec][i]));

      // Hard-code EVENANDODD parity in measure_chirality
      measure_chirality(eigVec[ivec], &chirality);
      node0_printf("CHIRALITY %d %.4g\n", ivec, chirality / 2);
    }

    // Perhaps in the future
    // Hard-code EVEN parity in print_densities
//    for (i = 0; i < Nvecs; i++) {
//      sprintf(label, "DENSITY(%i)", i);
//      print_densities(eigVec[i], label, ny / 2, nz / 2, nt / 2);
//    }

    for (i = 0; i < Nvecs; i++)
      free(eigVec[i]);
    free(eigVec);
    free(eigVal);

    node0_printf("RUNNING COMPLETED\n");
    dtime += dclock();
    node0_printf("Time = %.4g seconds\n", dtime);
    node0_printf("total_iters = %d\n", total_iters);
    fflush(stdout);
  } // readin(prompt) == 0
  return 0;
}
// -----------------------------------------------------------------

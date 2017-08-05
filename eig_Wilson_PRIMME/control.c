// -----------------------------------------------------------------
// Main procedure for calculating hermitian-Wilson eigenvalues
#define CONTROL
#include "eig_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  register int i, dir, ismear, ivec;
  register site *s;
  int prompt;
  Real ssplaq, stplaq;
  double dtime;

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

  // Smear however many times requested
  // Print out plaquette each time (unsmeared plaquette printed by setup)
  for (ismear = 0; ismear < nsmear; ismear++) {
    unphased_block_and_fatten();   // Smears gauge_field_thin into gauge_field
    FORALLUPDIR(dir) {
      FORALLSITES(i, s)
        gauge_field_thin[dir][i] = gauge_field[dir][i];
    }

    plaquette(&ssplaq, &stplaq);
    node0_printf("Plaquettes after smearing %d: %.8g %.8g\n",
        ismear + 1, ssplaq, stplaq);
  }
  node0_printf("\n");

  // Turn on anti-periodic boundary conditions
  boundary_flip_nhyp(MINUS);

  // Allocate eigenvectors now that we know Nvecs
  eigVal = malloc(Nvecs * sizeof(double));
  eigVec = malloc(Nvecs * sizeof(wilson_vector *));
  for (ivec = 0; ivec < Nvecs; ivec++)
    eigVec[ivec] = malloc(sites_on_node * sizeof(wilson_vector));
  tmpvec = malloc(sites_on_node * sizeof(wilson_vector));

  // Compute the eigenvectors
  total_iters = make_evs(Nvecs, eigVec, eigVal);

  // Clean up and close down
  for (i = 0; i < Nvecs; i++)
    free(eigVec[i]);
  free(eigVec);
  free(eigVal);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------

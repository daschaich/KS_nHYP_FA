// -----------------------------------------------------------------
// Main procedure for staggered eigenvalues
#define CONTROL
#include "eig_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  register int i, j, dir, ismear, ivec;
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

  // Allocate eigenvectors now that we know Nvecs
  eigVal = malloc((block + 1) * sizeof(double));
  eigVec = malloc((block + 1) * sizeof(vector *));
  for (ivec = 0; ivec < block + 1; ivec++)
    eigVec[ivec] = malloc(sites_on_node * sizeof(vector));

  // Cycle over blocks, after special case first
  // Calculate eigenvectors -- I haven't checked them,
  // so just print out the eigenvalues for now
  total_iters = make_evs(&start, block, eigVec, eigVal);
  for (ivec = 0; ivec < block; ivec++)
    node0_printf("EIGENVALUE %d %.8g\n", ivec, eigVal[ivec]);
  fflush(stdout);

  // Set next start to (a little under) the last eigenvector
  start = (eigVal[block - 1] + eigVal[block - 2]) / 2;

  for (i = block; i < Nvecs; i += block) {
    // Clear eigenvalues and try using eigenvectors as initial guesses
    // Number of initial guesses to use is primme.initSize in eig.c
    for (ivec = 0; ivec < block + 1; ivec++) {
      eigVal[ivec] = 0;
      FORALLSITES(j, s)
        vec_copy(&eigVec[block - ivec][j], &(eigVec[ivec][j]));
    }
    total_iters += make_evs(&start, block + 1, eigVec, eigVal);

    // Check last overlapping eigenvalue
    node0_printf("CHECK: start=%.8g --> %.8g\n", start, eigVal[0]);
    for (ivec = 1; ivec < block + 1; ivec++)
      node0_printf("EIGENVALUE %d %.8g\n", i + ivec - 1, eigVal[ivec]);
    fflush(stdout);

    // Set next start to (a little under) the last eigenvector
    start = (eigVal[block] + eigVal[block - 1]) / 2;
  }

  // Clean up and close down
  for (i = 0; i < block + 1; i++)
    free(eigVec[i]);
  free(eigVec);
  free(eigVal);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);

  normal_exit(0);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------

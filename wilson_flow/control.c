// -----------------------------------------------------------------
// Main procedure for SU(3) Wilson flow
#define CONTROL
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[])  {
  register int dir;
  int prompt;
  double dtime;

  su3_matrix *S[4];
  anti_hermitmat *A[4];

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Allocate fields used by integrator
  FORALLUPDIR(dir) {
    S[dir] = malloc(sites_on_node * sizeof(su3_matrix));
    A[dir] = malloc(sites_on_node * sizeof(anti_hermitmat));
  }

  // Run Wilson flow!
  // Adaptive step size stuff all handled in this routine
  wflow(S, A);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  FORALLUPDIR(dir) {
    free(S[dir]);
    free(A[dir]);
  }

  // Save lattice if requested -- no phases in links
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);
  normal_exit(0);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Main procedure for SU(3) MCRG-blocked measurements
#define CONTROL
#include "mcrg_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[])  {
  register int i, dir;
  register site *s;
  int j, k, prompt, num, bl, blmax;
  double dtime, rr[10], rro[10];
  complex plp, xplp;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run (loop removed)
  if (readin(prompt) == 0) {
    dtime = -dclock();
    // Set up loop tables
    make_loop_table2();

    // Maximum number of blockings determined by smallest dimension
    if (nx < nt)
      j = nx;
    else
      j = nt;

    // Works even if we can only block down to odd j > 4
    blmax = 0;
    while (j % 2 == 0 && j > 2) {    // While j is even
      j /= 2;
      blmax++;
    }

    node0_printf("Simple blocking from twice-nHYP smeared links\n");
    blocked_gauge_loops(0, rr);
    for (j = 0; j < nloop; j++) {
      for (k = 0; k < nreps; k++)
       node0_printf("LOOPS %d %d 0 0 %.8g\n",
                    j, k, rr[j + k * nloop] / volume / loop_num[j]);
    }

    // Save the original links
    FORALLSITES(i, s) {
      for (dir = XUP; dir <= TUP; dir++)
        su3mat_copy(&(s->link[dir]), &(s->link[dir + 8]));
    }

    // Checked that this reproduces the original ploop() result
    plp = blocked_ploop(0, TUP);
    xplp = blocked_ploop(0, XUP);
    node0_printf("POLYA ORIG %.6g %.6g %.6g %.6g\n",
                 (double)plp.real, (double)plp.imag,
                 (double)xplp.real, (double)xplp.imag);

    // Loop over different outer smearing parameters (MAX 100)
    // Number and values of outer smearing parameters from infile
    for (num = 0; num < num_alpha; num++) {
      if (alpha_mcrg[num] <= 0) {
        node0_printf("WARNING: Skipping alpha_mcrg[%d] = %.4g <= 0\n",
                     num, alpha_mcrg[num]);
        continue;
      }

      // Restore the original links
      FORALLSITES(i, s) {
        for (dir = XUP; dir <= TUP; dir++)
          su3mat_copy(&(s->link[dir + 8]), &(s->link[dir]));
      }

      // Loop over blocking levels (automatically determined)
      for (bl = 1; bl <= blmax; bl++) {
        for (j = 0; j < nloop; j++) {
          for (k = 0; k < nreps; k++)
            rro[j + k * nloop] = rr[j + k * nloop];
        }

        // nHYP smear twice with appropriate outer alpha
        alpha_smear[0] = alpha_mcrg[num];
        block_nhyp_mcrg(num, bl);
        block_nhyp_mcrg(num, bl);
        // Finally block
        // The first argument 0 only does the center staple
        block_mcrg(0, bl);

        blocked_gauge_loops(bl, rr);
        for (j = 0; j < nloop; j++) {
          for (k = 0; k < nreps; k++)
            rr[j + k * nloop] /= (volume * loop_num[j]);
        }
        for (j = 0; j < nloop; j++) {
          for (k = 0; k < nreps; k++)
            node0_printf("LOOPS %d %d %d %.4g %.8g\n",
                         j, k, bl, alpha_mcrg[num], rr[j + k * nloop]);
        }
//        for (j = 0; j < nloop; j++) {
//          node0_printf("CC %d %d %.4g ", j, bl, alpha_mcrg[num]);
//          for (k = 0; k < nloop; k++)
//            printf(" %.8g ", rr[j] * rr[k]);
//
//          node0_printf("\n");
//          node0_printf("CCO %d %d %.4g ", j, bl, alpha_mcrg[num]);
//          for (k = 0; k < nloop; k++)
//            printf(" %.8g ", rr[j] * rro[k]);
//
//          node0_printf("\n");
//        }

        // Calculate and print blocked Polyakov loops
        plp = blocked_ploop(bl, TUP);
        xplp = blocked_ploop(bl, XUP);
        node0_printf("POLYA NHYP %d %.6g %.6g %.6g %.6g %.6g\n",
                     bl, alpha_mcrg[num],
                     (double)plp.real, (double)plp.imag,
                     (double)xplp.real, (double)xplp.imag);
      } // End loop over blocking levels
    } // End loop over outer smearing parameters

    node0_printf("RUNNING COMPLETED\n");
    dtime += dclock();
    node0_printf("Time = %.4g seconds\n", dtime);
    fflush(stdout);
  } // readin(prompt) == 0
  return 0;
}
// -----------------------------------------------------------------

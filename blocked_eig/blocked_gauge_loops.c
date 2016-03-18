// -----------------------------------------------------------------
// Use ordinary gathers, accumulate in staple.e[0][0]
#include "block_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void blocked_gauge_loops(int block, double *result) {
  register int i, k;
  register site *s;
  int rep, dirs[10], sign[10], length;
  double g_action = 0, action, act2, total_action;
  complex trace;

  // For loop_table
  int ln, iloop;
  su3_matrix *tmat;
  tmat = (su3_matrix *)malloc(sites_on_node * sizeof(su3_matrix));
  if (tmat == NULL)
    exit(1);

  // Gauge action
  for (iloop = 0; iloop < nloop; iloop++) {
    for (rep = 0; rep < nreps; rep++)
      result[iloop + nloop * rep] = 0;
  }

  for (iloop = 0; iloop < nloop; iloop++) {
    length = loop_length[iloop];
    // Loop over rotations and reflections
    for (ln = 0; ln < loop_num[iloop]; ln++) {
      // Set up dirs and sign
      for (k = 0; k < length; k++) {
        if (loop_table[iloop][ln][k] < 4) {
          sign[k] = 1;
          dirs[k] = (loop_table[iloop][ln][k]) % 4;
        }
        else {
          sign[k] = -1;
          dirs[k] = (7 - loop_table[iloop][ln][k]) % 4;
        }
      }

      // Unfortunately, I don't seem able to combine these
      if (block > 0)
        blocked_path(block, dirs, sign, length, tmat);
      else
        path(dirs, sign, length, tmat);

      FORALLSITES(i, s) {
        trace = trace_su3(&tmat[i]);
        // Avoid excessively large numbers
        // Measure relative to some reasonable value of the loops
        action = (double)trace.real;
        total_action = (double)loop_coeff[iloop][0] * action;
        act2 = action;
        result[iloop] += act2;
        for (rep = 1; rep < nreps; rep++) {
          act2 *= action;
          result[iloop + nloop * rep] += act2;
          total_action += (double)loop_coeff[iloop][rep] * act2;
        }
        g_action  += total_action;
      }
    } // End loop over rotations and reflections
  } // End loop over iloop
  g_doublesum(&g_action);
  for (iloop = 0; iloop < nloop; iloop++) {
    for (rep = 0; rep < nreps; rep++) {
      g_doublesum(&result[iloop + nloop * rep]);
//      node0_printf("LOOP %d %d %.8g \n",
//                   iloop, loop_num[iloop],
//                   result[iloop] / volume / loop_num[iloop]);
    }
  }
  free(tmat);
}
// -----------------------------------------------------------------

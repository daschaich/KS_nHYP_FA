// -----------------------------------------------------------------
// Construct staples for Wilson and Symanzik flow
// Use tempmat, tempmat2 and tempsym for temporary storage
#include "flow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute 1-hop staples defined by directions dir and dir2
// Add to stp with given coefficient
void staple_1hop(int dir, int dir2, field_offset lnk1, field_offset lnk2,
                 Real coeff, matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  matrix tmat, tmat2;

  // Get blocked_link[dir2] from direction dir
  tag0 = start_gather_site(lnk2, sizeof(matrix), dir,
                           EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(matrix), dir2,
                           EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat,
  // then gathered to x
  FORALLSITES(i, s)
    mult_an((matrix*)F_PT(s, lnk2), (matrix*)F_PT(s, lnk1), &(tempmat[i]));

  // Finish lower staple
  wait_gather(tag0);
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mult_nn(&(tempmat[i]), (matrix *)gen_pt[0][i], &tmat);
    mat_copy(&tmat, &(tempmat[i]));       // Overwrite tempmat
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat, sizeof(matrix),
                            OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_nn((matrix *)F_PT(s, lnk2), (matrix *)gen_pt[1][i], &tmat);
    mult_na(&tmat, (matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_sum_matrix(&tmat2, coeff, &(stp[i]));
  }

  // Finally add the lower staple
  wait_gather(tag2);
  FORALLSITES(i, s)
    scalar_mult_sum_matrix((matrix *)gen_pt[2][i], coeff, &(stp[i]));

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute 2-hop staples defined by directions dir and dir2
// Add to stp with given coefficient
//
void staple_2hop(int dir, int dir2, field_offset lnk1, field_offset lnk2,
                 Real coeff, matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4, *tag5, *tag6, *tag7, *tag8;
  matrix tmat, tmat2;

  // First compute top of upper staple
  // Get blocked_link[dir2] from direction dir
  tag0 = start_gather_site(lnk2, sizeof(matrix), dir,
                           EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(matrix), dir2,
                           EVENANDODD, gen_pt[1]);

  // Store product in tempmat
  wait_gather(tag0);
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mult_nn((matrix*)F_PT(s, lnk2), (matrix *)gen_pt[1][i], &tmat);
    mult_na(&tmat, (matrix *)gen_pt[0][i], &(tempmat[i]));
  }

//TODO: MUCH MORE STILL TO DO...
  terminate(1);

  cleanup_gather(tag0);
  cleanup_gather(tag1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Sum staples for direction dir over all other directions
void staple(matrix *stp[NDIMS]) {
  register int i;
  register site *s;
  int dir, dir2;
  Real coeff_1hop = 0.0, coeff_2hop = 0.0;

  // Set staple coefficients
  if (flowflag == WILSON) {
    coeff_1hop = 1.0;
    coeff_2hop = 0.0;   // Should never be used, but let's make sure
  }         // Braces suppress compiler complaint
  else if (flowflag == SYMANZIK) {
    coeff_1hop =  5.0 /  3.0;
    coeff_2hop = -1.0 / 12.0;
  }         // Braces suppress compiler complaint
  else {    // This should have been caught by setup.c
    node0_printf("Error: Unrecognized flow type %d, aborting\n", flowflag);
    terminate(1);
  }

  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      clear_mat(&(stp[dir][i]));    // Clear the staple collectors

    FORALLUPDIR(dir2) {
      if (dir == dir2)
        continue;

      // Add 1-hop contributions to staples
      staple_1hop(dir, dir2, F_OFFSET(link[dir]), F_OFFSET(link[dir2]),
                  coeff_1hop, stp[dir]);

      // Add 2-hop contributions to staples
      if (flowflag != WILSON) {
        staple_2hop(dir, dir2, F_OFFSET(link[dir]), F_OFFSET(link[dir2]),
                    coeff_2hop, stp[dir]);
      }
    }
  }
}
// -----------------------------------------------------------------

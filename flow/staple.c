// -----------------------------------------------------------------
// Construct staples for Wilson and Symanzik flow
// Use tempmat, tempmat2 and tempsym for temporary storage
#include "flow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute 3-link staples defined by directions dir and dir2
// Add to stp with given coefficient
void staple_3link(int dir, int dir2, field_offset lnk1, field_offset lnk2,
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
// Compute 5-link staples defined by directions dir and dir2
// Include all six ways of removing one link from the 1x2 rectangle
// Add to stp with given coefficient
void staple_5link(int dir, int dir2, field_offset lnk1, field_offset lnk2,
                 Real coeff, matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4, *tag5, *tag6, *tag7, *tag8;
  matrix tmat, tmat2;

  // Get link[dir2] from direction dir and link[dir] from direction dir2
  // These will be used in several terms
  tag0 = start_gather_site(lnk2, sizeof(matrix), dir,
                           EVENANDODD, gen_pt[0]);
  tag1 = start_gather_site(lnk1, sizeof(matrix), dir2,
                           EVENANDODD, gen_pt[1]);

  // (1) Construct top of upper staple in direction dir2
  // Store it in tempmat and gather from direction dir2
  wait_gather(tag0);
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mult_nn((matrix*)F_PT(s, lnk2), (matrix *)gen_pt[1][i], &tmat);
    mult_na(&tmat, (matrix *)gen_pt[0][i], &(tempmat[i]));
  }
  tag2 = start_gather_field(tempmat, sizeof(matrix), dir2,
                            EVENANDODD, gen_pt[2]);

  // (2) Construct bottom of lower staple in direction -dir2
  // Store it in tempmat2 and gather from direction -dir2
  FORALLSITES(i, s) {
    mult_an((matrix*)F_PT(s, lnk2), (matrix*)F_PT(s, lnk1), &tmat);
    mult_nn(&tmat, (matrix *)gen_pt[0][i], &(tempmat2[i]));
  }
  tag3 = start_gather_field(tempmat2, sizeof(matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[3]);

  // (3--4) Construct end of staples protruding in direction dir
  // Store it in tempsym[0] and gather from direction dir
  FORALLSITES(i, s) {
    mult_nn((matrix*)F_PT(s, lnk1), (matrix *)gen_pt[0][i], &tmat);
    mult_na(&tmat, (matrix *)gen_pt[1][i], &(tempsym[0][i]));
  }
  tag4 = start_gather_field(tempsym[0], sizeof(matrix), dir,
                            EVENANDODD, gen_pt[4]);

  // (5--6) Construct end of staples protruding in direction -dir
  // Store it in tempsym[1] and gather from direction -dir
  FORALLSITES(i, s) {
    mult_nn((matrix*)F_PT(s, lnk2), (matrix *)gen_pt[1][i], &tmat);
    mult_an((matrix*)F_PT(s, lnk1), &tmat, &(tempsym[1][i]));
  }
  tag5 = start_gather_field(tempsym[1], sizeof(matrix), OPP_DIR(dir),
                            EVENANDODD, gen_pt[5]);

  // (1) Finish constructing upper staple in direction dir2
  // Add it to stp and clean up
  wait_gather(tag2);
  FORALLSITES(i, s) {
    mult_nn((matrix*)F_PT(s, lnk2), (matrix *)gen_pt[2][i], &tmat);
    mult_na(&tmat, (matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_sum_matrix(&tmat2, coeff, &(stp[i]));
  }
  cleanup_gather(tag2);     // Done with this

  // (2) Finish constructing lower staple in direction -dir2
  // Store it in tempsym[2] and gather from direction -dir2
  wait_gather(tag3);
  FORALLSITES(i, s) {
    mult_an((matrix*)F_PT(s, lnk2), (matrix *)gen_pt[3][i], &tmat);
    mult_nn(&tmat, (matrix *)gen_pt[0][i], &(tempsym[2][i]));
  }
  tag6 = start_gather_field(tempsym[2], sizeof(matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[6]);
  cleanup_gather(tag3);     // Done with this

  // (4) Construct lower staple protruding in direction dir
  // Store it in tempsym[3] and gather from direction -dir2
  wait_gather(tag4);
  FORALLSITES(i, s) {
    mult_nn((matrix*)F_PT(s, lnk1), (matrix *)gen_pt[4][i], &tmat);
    mult_an((matrix*)F_PT(s, lnk2), &tmat, &(tempsym[3][i]));
  }
  tag7 = start_gather_field(tempsym[3], sizeof(matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[7]);

  // (6) Construct lower staple protruding in direction -dir
  // Store it in tempsym[4] and gather from direction -dir2
  wait_gather(tag5);
  FORALLSITES(i, s) {
    mult_an((matrix *)gen_pt[5][i], (matrix*)F_PT(s, lnk1), &tmat);
    mult_nn(&tmat, (matrix *)gen_pt[0][i], &(tempsym[4][i]));
  }
  tag8 = start_gather_field(tempsym[4], sizeof(matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[8]);

  // (3) Construct and add upper staple protruding in direction dir
  FORALLSITES(i, s) {
    mult_nn((matrix*)F_PT(s, lnk2), (matrix *)gen_pt[1][i], &tmat);
    mult_na(&tmat, (matrix *)gen_pt[4][i], &tmat2);
    scalar_mult_sum_matrix(&tmat2, coeff, &(stp[i]));
  }
  cleanup_gather(tag4);     // Done with this

  // (5) Construct and add upper staple protruding in direction -dir
  FORALLSITES(i, s) {
    mult_nn((matrix *)gen_pt[5][i], (matrix *)gen_pt[1][i], &tmat);
    mult_na(&tmat, (matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_sum_matrix(&tmat2, coeff, &(stp[i]));
  }
  cleanup_gather(tag0);     // Done with these
  cleanup_gather(tag1);
  cleanup_gather(tag5);

  // (2) Add lower staple in direction -dir2
  wait_gather(tag6);
  FORALLSITES(i, s)
    scalar_mult_sum_matrix((matrix *)gen_pt[6][i], coeff, &(stp[i]));
  cleanup_gather(tag6);     // Done with this

  // (4) Add lower staple protruding in direction dir
  wait_gather(tag7);
  FORALLSITES(i, s)
    scalar_mult_sum_matrix((matrix *)gen_pt[7][i], coeff, &(stp[i]));
  cleanup_gather(tag7);     // Done with this

  // (6) Add lower staple protruding in direction dir
  wait_gather(tag8);
  FORALLSITES(i, s)
    scalar_mult_sum_matrix((matrix *)gen_pt[8][i], coeff, &(stp[i]));
  cleanup_gather(tag8);     // Done with this
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Sum staples for direction dir over all other directions
void staple(matrix *stp[NDIMS]) {
  register int i;
  register site *s;
  int dir, dir2;
  Real coeff_3link = 0.0, coeff_5link = 0.0;

  // Set staple coefficients
  if (flowflag == WILSON) {
    coeff_3link = 1.0;
    coeff_5link = 0.0;   // Should never be used, but let's make sure
  }         // Braces suppress compiler complaint
  else if (flowflag == SYMANZIK) {
    coeff_3link =  5.0 /  3.0;
    coeff_5link = -1.0 / 12.0;
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

      // Add 3-link contributions to staples
      staple_3link(dir, dir2, F_OFFSET(link[dir]), F_OFFSET(link[dir2]),
                   coeff_3link, stp[dir]);

      // Add 5-link contributions to staples
      if (flowflag != WILSON) {
        staple_5link(dir, dir2, F_OFFSET(link[dir]), F_OFFSET(link[dir2]),
                     coeff_5link, stp[dir]);
      }
    }
  }
}
// -----------------------------------------------------------------

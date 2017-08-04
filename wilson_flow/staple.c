// -----------------------------------------------------------------
// Construct staples
// Use tempmat for temporary storage
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Add staples defined by directions dir and dir2
void directional_staple(int dir, int dir2, field_offset lnk1,
                        field_offset lnk2, su3_matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix tmat, tmat2;

  // Get blocked_link[dir2] from direction dir
  tag0 = start_gather_site(lnk2, sizeof(su3_matrix), dir,
                           EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(su3_matrix), dir2,
                           EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat,
  // then gathered to x
  FORALLSITES(i, s) {
    mult_su3_an((su3_matrix*)F_PT(s, lnk2), (su3_matrix*)F_PT(s, lnk1),
                &(tempmat[i]));
  }

  // Finish lower staple
  wait_gather(tag0);
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mult_su3_nn(&(tempmat[i]), (su3_matrix *)gen_pt[0][i], &tmat);
    su3mat_copy(&tmat, &(tempmat[i]));       // Initialize tempmat
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat, sizeof(su3_matrix),
                            OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_su3_nn((su3_matrix *)F_PT(s, lnk2), (su3_matrix *)gen_pt[1][i],
                &tmat);
    mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &tmat2);
    sum_su3_matrix(&tmat2, &(stp[i]));
  }

  // Finally add the lower staple
  wait_gather(tag2);
  FORALLSITES(i, s)
    sum_su3_matrix((su3_matrix *)gen_pt[2][i], &(stp[i]));

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Sums staples for direction dir over all other directions
void staple(su3_matrix *stp[NDIMS]) {
  register int i;
  register site *s;
  int dir, dir2;

  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      clear_su3mat(&(stp[dir][i]));

    FORALLUPDIR(dir2) {
      if (dir != dir2)
        directional_staple(dir, dir2, F_OFFSET(link[dir]),
                           F_OFFSET(link[dir2]), stp[dir]);
    }
  }
}
// -----------------------------------------------------------------

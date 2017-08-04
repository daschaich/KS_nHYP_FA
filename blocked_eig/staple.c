// -----------------------------------------------------------------
// Construct staples
#include "block_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void directional_staple(int dir1, int dir2, field_offset lnk1,
                        field_offset lnk2, matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  matrix tmat1, tmat2;

  // Get blocked_link[dir2] from direction dir1
  tag0 = start_gather_site(lnk2, sizeof(matrix), dir1,
                      EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir1] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(matrix), dir2,
                      EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat1,
  // then gathered to x
  FORALLSITES(i, s)
    mult_su3_an((matrix*)F_PT(s,lnk2), (matrix*)F_PT(s,lnk1),
                tempmat1 + i);

   wait_gather(tag0);
   wait_gather(tag1);

  // Finish lower staple
  FORALLSITES(i, s) {
    mult_su3_nn(tempmat1 + i, (matrix *)gen_pt[0][i], &tmat1);
    mat_copy(&tmat1, tempmat1 + i);
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat1, sizeof(matrix),
                      OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_su3_nn((matrix*)F_PT(s,lnk2), (matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (matrix *)gen_pt[0][i], &tmat2);
    add_matrix(stp + i, &tmat2, stp + i);
  }

  // Finally add the lower staple.
  wait_gather(tag2);
  FORALLSITES(i, s)
    add_matrix(stp+i, (matrix *)gen_pt[2][i], stp+i);

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void staple(matrix *stp[4]) {
  register int i;
  register site *s;
  int dir1, dir2;

  for (dir1 = XUP; dir1 <= TUP; dir1++) {
    FORALLSITES(i, s)
      clear_mat(&(stp[dir1][i]));

    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir1 != dir2)
        directional_staple(dir1, dir2, F_OFFSET(link[dir1]),
                           F_OFFSET(link[dir2]), stp[dir1]);
    }
  }
}
// -----------------------------------------------------------------

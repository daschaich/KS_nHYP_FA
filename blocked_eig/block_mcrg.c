// -----------------------------------------------------------------
// Blocking function for MCRG-blocked measurements
// Handles "neighboring" sites separated by 2^block links
#include "block_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void staple_mcrg(int dir1, int block) {
  register int i, dir2;
  register site *s;
  int j, bl, start, disp[4];    // Displacement vector for general gather
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4, *tag5, *tag6;
  matrix tmat1, tmat2;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  start = 1;                    // Indicates staple sum not initialized
  // Loop over other directions
  for (dir2 = XUP; dir2 <= TUP; dir2++) {
    if (dir2 != dir1) {
      // Get link[dir2] from direction 2 * dir1
      clear_disp(disp);
      disp[dir1] = 2 * bl;
      tag0 = start_general_gather_site(F_OFFSET(link[dir2]),
                                       sizeof(matrix), disp,
                                       EVENANDODD, gen_pt[0]);
      wait_general_gather(tag0);

      // Get link[dir1] from direction dir2
      clear_disp(disp);
      disp[dir2] = bl;
      tag1 = start_general_gather_site(F_OFFSET(link[dir1]),
                                       sizeof(matrix), disp,
                                       EVENANDODD, gen_pt[1]);
      wait_general_gather(tag1);

      // Get link[dir1] from direction dir1 + dir2
      clear_disp(disp);
      disp[dir1] = bl;
      disp[dir2] = bl;
      tag2 = start_general_gather_site(F_OFFSET(link[dir1]),
                                       sizeof(matrix), disp,
                                       EVENANDODD, gen_pt[2]);
      wait_general_gather(tag2);

      // Get link[dir2] from direction -dir2
      clear_disp(disp);
      disp[dir2] = -bl;
      tag3 = start_general_gather_site(F_OFFSET(link[dir2]),
                                       sizeof(matrix), disp,
                                       EVENANDODD, gen_pt[3]);
      wait_general_gather(tag3);

      // Get link[dir1] from direction -dir2
      clear_disp(disp);
      disp[dir2] = -bl;
      tag4 = start_general_gather_site(F_OFFSET(link[dir1]),
                                       sizeof(matrix), disp,
                                       EVENANDODD, gen_pt[4]);
      wait_general_gather(tag4);

      // Get link[dir1] from direction dir1 - dir2
      clear_disp(disp);
      disp[dir1] = bl;
      disp[dir2]= -bl;
      tag5 = start_general_gather_site(F_OFFSET(link[dir1]),
                                       sizeof(matrix), disp,
                                       EVENANDODD, gen_pt[5]);
      wait_general_gather(tag5);

      // Get link[dir2] from direction 2 * dir1 - dir2
      clear_disp(disp);
      disp[dir1] = 2 * bl;
      disp[dir2] = -bl;

      tag6 = start_general_gather_site(F_OFFSET(link[dir2]),
                                       sizeof(matrix), disp,
                                       EVENANDODD, gen_pt[6]);
      wait_general_gather(tag6);

      // Upper staple
      if (start) {          // The first contribution to the staple
        FORALLSITES(i, s) {
          mult_su3_nn(&(s->link[dir2]), (matrix *)(gen_pt[1][i]), &tmat1);
          mult_su3_nn(&tmat1, (matrix *)(gen_pt[2][i]), &tmat2);
          mult_su3_na(&tmat2, (matrix *)(gen_pt[0][i]), &(s->tempmat2));
        }
        start = 0;
      }
      else {
        FORALLSITES(i, s) {
          mult_su3_nn(&(s->link[dir2]), (matrix *)(gen_pt[1][i]), &tmat1);
          mult_su3_nn(&tmat1, (matrix *)(gen_pt[2][i]), &tmat2);
          mult_su3_na(&tmat2, (matrix *)(gen_pt[0][i]), &tmat1);
          add_matrix(&(s->tempmat2), &tmat1, &(s->tempmat2));
        }
      }
      cleanup_general_gather(tag0);
      cleanup_general_gather(tag1);
      cleanup_general_gather(tag2);

      // Lower staple
      FORALLSITES(i, s) {
        mult_su3_an((matrix *)(gen_pt[3][i]),
                    (matrix *)(gen_pt[4][i]), &tmat1);
        mult_su3_nn(&tmat1, (matrix *)(gen_pt[5][i]), &tmat2);
        mult_su3_nn(&tmat2, (matrix *)(gen_pt[6][i]), &tmat1);
        add_matrix(&(s->tempmat2), &tmat1, &(s->tempmat2));
      }
      cleanup_general_gather(tag3);
      cleanup_general_gather(tag4);
      cleanup_general_gather(tag5);
      cleanup_general_gather(tag6);
    }
  } // End loop over dir2 != dir1
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Closely resembles ../generic_nhyp/block_nhyp.c:block_nhyp3()
// If num == 0, do only center link, not the full staple
void block_mcrg(int num, int block) {
  register int dir, i;
  register site *s;
  int disp[4], j;
  Real f[3];
  complex ctmp;
  msg_tag *tag0;
  matrix tmat, Omega, eQ, Id, Q, Q2;

  for (dir = XUP; dir <= TUP; dir++) {
    // First the central spine of the blocking
    clear_disp(disp);
    disp[dir] = 1;
    for (j = 1; j < block; j++)
      disp[dir] *= 2;             // Block size

    tag0 = start_general_gather_site(F_OFFSET(link[dir]),
                                     sizeof(matrix), disp,
                                     EVENANDODD, gen_pt[0]);
    wait_general_gather(tag0);

    FORALLSITES(i, s)
      mult_su3_nn(&(s->link[dir]), (matrix *)(gen_pt[0][i]),
                  &(s->staple));

    if (num == 0) {               // Do only center link
      FORALLSITES(i, s)
        mat_copy(&(s->staple), &(s->link[dir + 4]));
    }
    else {                         // Do the full staple
      staple_mcrg(dir, block);

      FORALLSITES(i, s) {
        // Make Omega
        scalar_mult_matrix(&(s->staple), 1 - alpha_smear[num], &Q);
        scalar_mult_add_matrix(&Q, &(s->tempmat2),
                                   alpha_smear[num] / 6.0, &Omega);
        mult_su3_an(&Omega, &Omega, &Q);

        // IR stabilization regulator set in defines.h
        scalar_add_diag_su3(&Q, IR_STAB);
#ifndef NHYP_DEBUG
        compute_fhb(&Q, f, NULL, 0);
#else
        compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

        // Compute Q**2
        mult_su3_nn(&Q, &Q, &Q2);

        // Compute Q^(-1/2) via Eq. (3.8)
        ctmp = cmplx(f[0], 0);
        diag_su3(&Id, &ctmp);
        scalar_mult_add_matrix(&Id, &Q, f[1], &tmat);
        scalar_mult_add_matrix(&tmat, &Q2, f[2], &eQ);

        // Multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2)
        mult_su3_nn(&Omega, &eQ, &tmat);
        mat_copy(&tmat, &(s->link[dir + 4]));
      }
    } // End of staple
  } // End loop over dir

  FORALLSITES(i, s) {
    for (dir = XUP; dir <= TUP; dir++)
      mat_copy(&(s->link[dir + 4]), &(s->link[dir]));
  }
}
// -----------------------------------------------------------------

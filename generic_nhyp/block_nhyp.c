// -----------------------------------------------------------------
// Construct nHYP-smeared links
// Use tempmat for temporary storage
// Reference:
//   Anna Hasenfratz, Roland Hoffmann and Stefan Schaefer
//   ``Hypercubic smeared links for dynamical fermions''
//   JHEP 0705:029,2007. [hep-lat/0702028]
#include "nhyp_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Do three smearing levels to construct nHYP links
void block_nhyp() {
#ifdef TIMING
  TIC(3)
#endif
  block_nhyp1();
  block_nhyp2();
  block_nhyp3();
#ifdef TIMING
  TOC(3, time_block_nhyp)
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp3() {
  register int dir, dir2, i;
  register site *s;
  Real f[3];
  Real ftmp1, ftmp2;
  matrix Omega, Q;
  matrix tmat, eQ, Q2;

  ftmp1 = alpha_smear[0] / (6 * (1 - alpha_smear[0]));
  ftmp2 = 1 - alpha_smear[0];

  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, s)
      clear_mat(&Staple3[dir][i]);

    // Compute the staple
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir2 != dir)
        staple_nhyp(dir, dir2, hyplink2[dir2][dir],
                    hyplink2[dir][dir2], Staple3[dir]);
    }

    FORALLSITES(i, s) {
      // Make Omega
      scalar_mult_add_matrix(gauge_field_thin[dir] + i,
                                 Staple3[dir] + i, ftmp1, &Q);

      scalar_mult_matrix(&Q, ftmp2, &Omega);
      Staple3[dir][i] = Omega;     // Save for force
      mult_an(&Omega, &Omega, &Q);

      // IR stabilization regulator set in defines.h
      scalar_add_diag(&Q, IR_STAB);
#ifndef NHYP_DEBUG
      compute_fhb(&Q, f, NULL, 0);
#else
      compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

      // Compute Q^2
      mult_nn(&Q, &Q, &Q2);

      // Compute Q^(-1/2) via Eq. (3.8)
      scalar_mult_matrix(&Q, f[1], &tmat);
      scalar_mult_add_matrix(&tmat, &Q2, f[2], &eQ);
      scalar_add_diag(&eQ, f[0]);

      // Multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2)
      mult_nn(&Omega, &eQ, gauge_field[dir] + i);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp2() {
  register int dir, dir2, dir3, dir4, i;
  register site *st;
  Real f[3], ftmp1, ftmp2;
  matrix Omega, Q;
  matrix tmat, eQ, Q2;

  ftmp1 = alpha_smear[1] / (4 * (1 - alpha_smear[1]));
  ftmp2 = (1 - alpha_smear[1]);

  for (dir = XUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir2 != dir) {
        FORALLSITES(i, st)
          clear_mat(Staple2[dir2][dir] + i);

        // Compute the staple
        for (dir3 = 0; dir3 < 4; dir3++) {
          if (dir3 != dir && dir3 != dir2) {
            for (dir4 = XUP; dir4 <= TUP; dir4++) {
              if (dir4 != dir && dir4 != dir2 && dir4 != dir3)
                break;
            }
            staple_nhyp(dir, dir3, hyplink1[dir4][dir],
                        hyplink1[dir4][dir3], Staple2[dir2][dir]);
          }
        }

        FORALLSITES(i, st) {
          // Make Omega
          scalar_mult_add_matrix(gauge_field_thin[dir] + i,
                                     Staple2[dir2][dir] + i, ftmp1, &Q);

          scalar_mult_matrix(&Q, ftmp2, &Omega);
          Staple2[dir2][dir][i] = Omega;     // Save for force

          mult_an(&Omega, &Omega, &Q);
          scalar_add_diag(&Q, IR_STAB);
#ifndef NHYP_DEBUG
          compute_fhb(&Q, f, NULL, 0);
#else
          compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

          // Compute Q^2
          mult_nn(&Q, &Q, &Q2);

          // Compute Q^(-1/2) via Eq. (3.8)
          scalar_mult_matrix(&Q, f[1], &tmat);
          scalar_mult_add_matrix(&tmat, &Q2, f[2], &eQ);
          scalar_add_diag(&eQ, f[0]);

          // Multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2) 
          mult_nn(&Omega, &eQ, hyplink2[dir2][dir] + i);
        }
      }
    } // Loop over dir2
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp1() {
  register int dir1, dir2, i;
  register site *s;
  Real f[3], ftmp1, ftmp2;
  matrix Omega, Q, tmat, eQ, Q2;

  ftmp1 = alpha_smear[2] / (2 * (1 - alpha_smear[2]));
  ftmp2 = (1 - alpha_smear[2]);

  // dir1 is the direction of the original link
  // dir2 is the other direction that defines the staple
  for (dir1 = XUP; dir1 <= TUP; dir1++) {
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir1 != dir2) {
       FORALLSITES(i, s)
         clear_mat(Staple1[dir2][dir1] + i);

        // Compute the staple
        staple_nhyp(dir1, dir2, gauge_field_thin[dir1],
                    gauge_field_thin[dir2], Staple1[dir2][dir1]);

        FORALLSITES(i, s) {
          // Make Omega
          scalar_mult_add_matrix(gauge_field_thin[dir1] + i,
                                     Staple1[dir2][dir1] + i, ftmp1, &Q);

          scalar_mult_matrix(&Q, ftmp2, &Omega);
          Staple1[dir2][dir1][i] = Omega;     // Save for force

          mult_an(&Omega, &Omega, &Q);
          scalar_add_diag(&Q, IR_STAB);
#ifndef NHYP_DEBUG
          compute_fhb(&Q, f, NULL, 0);
#else
          compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

          // Compute Q^2
          mult_nn(&Q, &Q, &Q2);

          // Compute Q^(-1/2) via Eq. (3.8)
          scalar_mult_matrix(&Q, f[1], &tmat);
          scalar_mult_add_matrix(&tmat, &Q2, f[2], &eQ);
          scalar_add_diag(&eQ, f[0]);

          // Multiply Omega by eQ = (Omega^dag Omega)^(-1/2)
          mult_nn(&Omega, &eQ, hyplink1[dir2][dir1] + i);
        }
      }
    } // Loop over dir2
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dir1 is the direction of the original link
// dir2 is the other direction that defines the staple
// stp must be cleared before being this function is called!
void staple_nhyp(int dir1, int dir2, matrix *lnk1,
                 matrix *lnk2, matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  matrix tmat1, tmat2;

  // Get blocked_link[dir2] from direction dir1
  tag0 = start_gather_field(lnk2, sizeof(matrix), dir1,
                            EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir1] from direction dir2
  tag1 = start_gather_field(lnk1, sizeof(matrix), dir2,
                            EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat,
  // then gathered to x
  FORALLSITES(i, s)
    mult_an(lnk2 + i, lnk1 + i, tempmat + i);

  wait_gather(tag0);
  wait_gather(tag1);

  // Finish lower staple
  FORALLSITES(i, s) {
    mult_nn(tempmat + i, (matrix *)gen_pt[0][i], &tmat1);
    mat_copy(&tmat1, tempmat + i);
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat, sizeof(matrix),
                            OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_nn(lnk2 + i, (matrix *)gen_pt[1][i], &tmat1);
    mult_na(&tmat1, (matrix *)gen_pt[0][i], &tmat2);
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

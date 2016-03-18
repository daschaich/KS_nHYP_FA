// -----------------------------------------------------------------
// Construct nHYP-smeared links
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
  su3_matrix Omega, Q;
  su3_matrix tmat, eQ, Q2;

  ftmp1 = alpha_smear[0] / (6 * (1 - alpha_smear[0]));
  ftmp2 = 1 - alpha_smear[0];

  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, s)
      clear_su3mat(&Staple3[dir][i]);

    // Compute the staple
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir2 != dir)
        staple_nhyp(dir, dir2, hyplink2[dir2][dir],
                    hyplink2[dir][dir2], Staple3[dir]);
    }

    FORALLSITES(i, s) {
      // Make Omega
      scalar_mult_add_su3_matrix(gauge_field_thin[dir] + i,
                                 Staple3[dir] + i, ftmp1, &Q);

      scalar_mult_su3_matrix(&Q, ftmp2, &Omega);
      Staple3[dir][i] = Omega;     // Save for force
      mult_su3_an(&Omega, &Omega, &Q);

      // IR stabilization regulator set in defines.h
      scalar_add_diag_su3(&Q, IR_STAB);
#ifndef NHYP_DEBUG
      compute_fhb(&Q, f, NULL, 0);
#else
      compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

      // Compute Q^2
      mult_su3_nn(&Q, &Q, &Q2);

      // Compute Q^(-1/2) via Eq. (3.8)
      scalar_mult_su3_matrix(&Q, f[1], &tmat);
      scalar_mult_add_su3_matrix(&tmat, &Q2, f[2], &eQ);
      scalar_add_diag_su3(&eQ, f[0]);

      // Multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2)
      mult_su3_nn(&Omega, &eQ, gauge_field[dir] + i);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp2() {
  register int dir, dir2, dir3, dir4, i;
  register site *st;
  Real f[3], ftmp1, ftmp2;
  su3_matrix Omega, Q;
  su3_matrix tmat, eQ, Q2;

  ftmp1 = alpha_smear[1] / (4 * (1 - alpha_smear[1]));
  ftmp2 = (1 - alpha_smear[1]);

  for (dir = XUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir2 != dir) {
        FORALLSITES(i, st)
          clear_su3mat(Staple2[dir2][dir] + i);

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
          scalar_mult_add_su3_matrix(gauge_field_thin[dir] + i,
                                     Staple2[dir2][dir] + i, ftmp1, &Q);

          scalar_mult_su3_matrix(&Q, ftmp2, &Omega);
          Staple2[dir2][dir][i] = Omega;     // Save for force

          mult_su3_an(&Omega, &Omega, &Q);
          scalar_add_diag_su3(&Q, IR_STAB);
#ifndef NHYP_DEBUG
          compute_fhb(&Q, f, NULL, 0);
#else
          compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

          // Compute Q^2
          mult_su3_nn(&Q, &Q, &Q2);

          // Compute Q^(-1/2) via Eq. (3.8)
          scalar_mult_su3_matrix(&Q, f[1], &tmat);
          scalar_mult_add_su3_matrix(&tmat, &Q2, f[2], &eQ);
          scalar_add_diag_su3(&eQ, f[0]);

          // Multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2) 
          mult_su3_nn(&Omega, &eQ, hyplink2[dir2][dir] + i);
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
  su3_matrix Omega, Q, tmat, eQ, Q2;

  ftmp1 = alpha_smear[2] / (2 * (1 - alpha_smear[2]));
  ftmp2 = (1 - alpha_smear[2]);

  // dir1 is the direction of the original link
  // dir2 is the other direction that defines the staple
  for (dir1 = XUP; dir1 <= TUP; dir1++) {
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir1 != dir2) {
       FORALLSITES(i, s)
         clear_su3mat(Staple1[dir2][dir1] + i);

        // Compute the staple
        staple_nhyp(dir1, dir2, gauge_field_thin[dir1],
                    gauge_field_thin[dir2], Staple1[dir2][dir1]);

        FORALLSITES(i, s) {
          // Make Omega
          scalar_mult_add_su3_matrix(gauge_field_thin[dir1] + i,
                                     Staple1[dir2][dir1] + i, ftmp1, &Q);

          scalar_mult_su3_matrix(&Q, ftmp2, &Omega);
          Staple1[dir2][dir1][i] = Omega;     // Save for force

          mult_su3_an(&Omega, &Omega, &Q);
          scalar_add_diag_su3(&Q, IR_STAB);
#ifndef NHYP_DEBUG
          compute_fhb(&Q, f, NULL, 0);
#else
          compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

          // Compute Q^2
          mult_su3_nn(&Q, &Q, &Q2);

          // Compute Q^(-1/2) via Eq. (3.8)
          scalar_mult_su3_matrix(&Q, f[1], &tmat);
          scalar_mult_add_su3_matrix(&tmat, &Q2, f[2], &eQ);
          scalar_add_diag_su3(&eQ, f[0]);

          // Multiply Omega by eQ = (Omega^dag Omega)^(-1/2)
          mult_su3_nn(&Omega, &eQ, hyplink1[dir2][dir1] + i);
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
void staple_nhyp(int dir1, int dir2, su3_matrix *lnk1,
                 su3_matrix *lnk2, su3_matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix tmat1, tmat2;

  // Get blocked_link[dir2] from direction dir1
  tag0 = start_gather_field(lnk2, sizeof(su3_matrix), dir1,
                            EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir1] from direction dir2
  tag1 = start_gather_field(lnk1, sizeof(su3_matrix), dir2,
                            EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat1,
  // then gathered to x
  FORALLSITES(i, s)
    mult_su3_an(lnk2 + i, lnk1 + i, tempmat1 + i);

  wait_gather(tag0);
  wait_gather(tag1);

  // Finish lower staple
  FORALLSITES(i, s) {
    mult_su3_nn(tempmat1 + i, (su3_matrix *)gen_pt[0][i], &tmat1);
    su3mat_copy(&tmat1, tempmat1 + i);
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat1, sizeof(su3_matrix),
                            OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_su3_nn(lnk2 + i, (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], &tmat2);
    add_su3_matrix(stp + i, &tmat2, stp + i);
  }

  // Finally add the lower staple.
  wait_gather(tag2);
  FORALLSITES(i, s)
    add_su3_matrix(stp+i, (su3_matrix *)gen_pt[2][i], stp+i);

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------

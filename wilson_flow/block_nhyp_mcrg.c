// -----------------------------------------------------------------
// Specialized nHYP smearing for MCRG-blocked measurements
// Handles "neighboring" sites separated by 2^block links
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Trivial helper function called several times below
void clear_disp(int *disp) {
  register int i;
  for (i = XUP; i <= TUP; i++)
    disp[i] = 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void diag_su3(su3_matrix* Q, complex *f) {
  Q->e[0][0].real = f->real;
  Q->e[0][0].imag = f->imag;
  Q->e[0][1].real = 0;
  Q->e[0][1].imag = 0;
  Q->e[0][2].real = 0;
  Q->e[0][2].imag = 0;
  Q->e[1][0].real = 0;
  Q->e[1][0].imag = 0;
  Q->e[1][1].real = f->real;
  Q->e[1][1].imag = f->imag;
  Q->e[1][2].real = 0;
  Q->e[1][2].imag = 0;
  Q->e[2][0].real = 0;
  Q->e[2][0].imag = 0;
  Q->e[2][1].real = 0;
  Q->e[2][1].imag = 0;
  Q->e[2][2].real = f->real;
  Q->e[2][2].imag = f->imag;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void staple1_mcrg(int dir, int dir3, int dir4, int block) {
  register int i, dir2;
  register site *s;
  int j, bl, disp[4];           // Displacement vector for general gather
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4;
  su3_matrix tmat, tmat2;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  // Loop over other directions
  for (dir2 = XUP; dir2 <= TUP; dir2++) {
    if (dir2 != dir && dir2 !=dir3 && dir2 != dir4) {
      // Get link[dir2] from direction dir
      clear_disp(disp);
      disp[dir] = bl;
      tag0 = start_general_gather_site(F_OFFSET(link[dir2]),
                                       sizeof(su3_matrix), disp,
                                       EVENANDODD, gen_pt[0]);
      wait_general_gather(tag0);

      // Get link[dir] from direction dir2
      clear_disp(disp);
      disp[dir2] = bl;
      tag1 = start_general_gather_site(F_OFFSET(link[dir]),
                                       sizeof(su3_matrix), disp,
                                       EVENANDODD, gen_pt[1]);
      wait_general_gather(tag1);

      // Get link[dir2] from direction -dir2
      clear_disp(disp);
      disp[dir2] = -bl;
      tag2 = start_general_gather_site(F_OFFSET(link[dir2]),
                                       sizeof(su3_matrix), disp,
                                       EVENANDODD, gen_pt[2]);
      wait_general_gather(tag2);

      // Get link[dir] from direction -dir2
      clear_disp(disp);
      disp[dir2] = -bl;
      tag3 = start_general_gather_site(F_OFFSET(link[dir]),
                                       sizeof(su3_matrix), disp,
                                       EVENANDODD, gen_pt[3]);
      wait_general_gather(tag3);

      // Get link[dir2] from displacement dir - dir2
      clear_disp(disp);
      disp[dir] = bl;
      disp[dir2] = -bl;
      tag4 = start_general_gather_site(F_OFFSET(link[dir2]),
                                       sizeof(su3_matrix), disp,
                                       EVENANDODD, gen_pt[4]);
      wait_general_gather(tag4);

      // Upper staple
      FORALLSITES(i, s) {
        mult_su3_nn(&(s->link[dir2]), (su3_matrix *)gen_pt[1][i], &tmat);
        mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &(s->tempmat2));
      }
      cleanup_general_gather(tag0);
      cleanup_general_gather(tag1);

      // Lower staple
      FORALLSITES(i, s) {
        mult_su3_an((su3_matrix *)gen_pt[2][i],
                    (su3_matrix *)gen_pt[3][i], &tmat);
        mult_su3_nn(&tmat, (su3_matrix *)gen_pt[4][i], &tmat2);
        sum_su3_matrix(&tmat2, &(s->tempmat2));
      }
      cleanup_general_gather(tag2);
      cleanup_general_gather(tag3);
      cleanup_general_gather(tag4);
    }
  } // End loop over other directions
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void staple2_mcrg(int dir, int dir4, int block) {
  register int i, dir2;
  register site *s;
  int start, j, bl, disp[4];    // Displacement vector for general gather
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4;
  su3_matrix tmat, tmat2;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  // Loop over other directions
  start = 1;                    // Indicates staple sum not initialized
  for (dir2 = XUP; dir2 <= TUP; dir2++) {
    if (dir2 != dir && dir2 != dir4) {
      // Get hyplink1[hyp1ind[dir][dir4][dir2]] from direction dir
      clear_disp(disp);
      disp[dir] = bl;
      tag0 = start_general_gather_site(
                   F_OFFSET(hyplink1[hyp1ind[dir][dir4][dir2]]),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[0]);
      wait_general_gather(tag0);

      // Get hyplink1[hyp1ind[dir2][dir4][dir]] from direction dir2
      clear_disp(disp);
      disp[dir2] = bl;
      tag1 = start_general_gather_site(
                   F_OFFSET(hyplink1[hyp1ind[dir2][dir4][dir]]),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[1]);
      wait_general_gather(tag1);

      // Get hyplink1[hyp1ind[dir][dir4][dir2]] from direction -dir2
      clear_disp(disp);
      disp[dir2] = -bl;
      tag2 = start_general_gather_site(
                   F_OFFSET(hyplink1[hyp1ind[dir][dir4][dir2]]),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[2]);
      wait_general_gather(tag2);

      // Get hyplink1[hyp1ind[dir2][dir4][dir]] from direction -dir2
      clear_disp(disp);
      disp[dir2] = -bl;
      tag3 = start_general_gather_site(
                   F_OFFSET(hyplink1[hyp1ind[dir2][dir4][dir]]),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[3]);
      wait_general_gather(tag3);

      // Get hyplink1[hyp1ind[dir][dir4][dir2]] from direction dir - dir2
      clear_disp(disp); disp[dir] = bl;disp[dir2] = -bl;
      tag4 = start_general_gather_site(
                   F_OFFSET(hyplink1[hyp1ind[dir][dir4][dir2]]),
                   sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4]);
      wait_general_gather(tag4);

      // Upper staple
      if (start) {        // The first contribution to the staple
        FORALLSITES(i, s) {
          mult_su3_nn(&(s->hyplink1[hyp1ind[dir][dir4][dir2]]),
                      (su3_matrix *)gen_pt[1][i], &tmat);
          mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &(s->tempmat2));
        }
        start = 0;
      }
      else {
        FORALLSITES(i, s) {
          mult_su3_nn(&(s->hyplink1[hyp1ind[dir][dir4][dir2]]),
                      (su3_matrix *)gen_pt[1][i], &tmat);
          mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &tmat2);
          sum_su3_matrix(&tmat2, &(s->tempmat2));
        }
      }
      cleanup_general_gather(tag0);
      cleanup_general_gather(tag1);

      // Lower staple
      FORALLSITES(i, s) {
        mult_su3_an((su3_matrix *)gen_pt[2][i], (su3_matrix *)gen_pt[3][i], &tmat);
        mult_su3_nn(&tmat, (su3_matrix *)gen_pt[4][i], &tmat2);
        add_su3_matrix(&(s->tempmat2), &tmat2, &(s->tempmat2));
      }
      cleanup_general_gather(tag2);
      cleanup_general_gather(tag3);
      cleanup_general_gather(tag4);
    }
  } // End loop over other directions
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void staple3_mcrg(int dir, int block) {
  register int i, dir2;
  register site *s;
  int start, j, bl, disp[4];      // displacement vector for general gather
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4;
  su3_matrix tmat,tmat2;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  // Loop over other directions
  start = 1;                    // Indicates staple sum not initialized
  FORALLUPDIR(dir2) {
    if (dir2 == dir)
      continue;
    // Get hyplink2[hyp2ind[dir][dir2]] from direction dir
    clear_disp(disp);
    disp[dir] = bl;
    tag0 = start_general_gather_site(F_OFFSET(hyplink2[hyp2ind[dir][dir2]]),
                                     sizeof(su3_matrix), disp,
                                     EVENANDODD, gen_pt[0]);
    wait_general_gather(tag0);

    // Get hyplink2[hyp2ind[dir2][dir]] from direction dir2
    clear_disp(disp);
    disp[dir2] = bl;
    tag1 = start_general_gather_site(F_OFFSET(hyplink2[hyp2ind[dir2][dir]]),
                                     sizeof(su3_matrix), disp,
                                     EVENANDODD, gen_pt[1]);
    wait_general_gather(tag1);

    // Get  from direction -dir2
    clear_disp(disp);
    disp[dir2] = -bl;
    tag2 = start_general_gather_site(F_OFFSET(hyplink2[hyp2ind[dir][dir2]]),
                                     sizeof(su3_matrix), disp,
                                     EVENANDODD, gen_pt[2]);
    wait_general_gather(tag2);

    // Get hyplink2[hyp2ind[dir2][dir]] from direction -dir2
    clear_disp(disp);
    disp[dir2] = -bl;
    tag3 = start_general_gather_site(F_OFFSET(hyplink2[hyp2ind[dir2][dir]]),
                                     sizeof(su3_matrix), disp,
                                     EVENANDODD, gen_pt[3]);
    wait_general_gather(tag3);

    // Get hyplink2[hyp2ind[dir][dir2]] from displacement +dir-dir2
    clear_disp(disp);
    disp[dir] = bl;
    disp[dir2] = -bl;
    tag4 = start_general_gather_site(F_OFFSET(hyplink2[hyp2ind[dir][dir2]]),
                                     sizeof(su3_matrix), disp,
                                     EVENANDODD, gen_pt[4]);
    wait_general_gather(tag4);

    // Upper staple
    if (start) {          // The first contribution to the staple
      FORALLSITES(i, s) {
        mult_su3_nn(&(s->hyplink2[hyp2ind[dir][dir2]]),
                    (su3_matrix *)gen_pt[1][i], &tmat);
        mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &(s->tempmat2));
      }
      start = 0;
    }
    else {
      FORALLSITES(i, s) {
        mult_su3_nn(&(s->hyplink2[hyp2ind[dir][dir2]]),
                    (su3_matrix *)gen_pt[1][i], &tmat);
        mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &tmat2);
        add_su3_matrix(&(s->tempmat2), &tmat2, &(s->tempmat2));
      }
    }
    cleanup_general_gather(tag0);
    cleanup_general_gather(tag1);

    // Lower staple
    FORALLSITES(i, s) {
      mult_su3_an((su3_matrix *)gen_pt[2][i], (su3_matrix *)gen_pt[3][i],
                  &tmat);
      mult_su3_nn(&tmat, (su3_matrix *)gen_pt[4][i], &tmat2);
      sum_su3_matrix(&tmat2, &(s->tempmat2));
    }
    cleanup_general_gather(tag2);
    cleanup_general_gather(tag3);
    cleanup_general_gather(tag4);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp1_mcrg(int num, int block) {
  register int dir, dir2, dir3, i;
  register site *s;
  Real f[3], ftmp1, ftmp2;
  complex tc;
  su3_matrix tmat, Omega, eQ, Id, Q, Q2;

  ftmp1 = alpha_smear[2] / (2.0 * (1.0 - alpha_smear[2]));
  ftmp2 = 1.0 - alpha_smear[2];

  // Loop over link directions
  FORALLUPDIR(dir) {
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;
      for (dir3 = dir2 + 1; dir3 <= TUP; dir3++) {
        if (dir == dir3)
          continue;

        // Compute the staple
        staple1_mcrg(dir, dir2, dir3, block);

        FORALLSITES(i, s) {
          // Make Omega
          scalar_mult_add_su3_matrix(&(s->link[dir]), &(s->tempmat2),
                                     ftmp1, &Q);
          scalar_mult_su3_matrix(&Q, ftmp2, &Omega);
          mult_su3_an(&Omega, &Omega, &Q);
          scalar_add_diag_su3(&Q, IR_STAB);
#ifndef NHYP_DEBUG
          compute_fhb(&Q, f, NULL, 0);
#else
          compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

          // Compute Q**2
          mult_su3_nn(&Q, &Q, &Q2);

          // Compute Q^(-1/2) via Eq. (3.8)
          tc = cmplx(f[0],0.);
          diag_su3(&Id, &tc);
          scalar_mult_add_su3_matrix(&Id, &Q, f[1], &tmat);
          scalar_mult_add_su3_matrix(&tmat, &Q2, f[2], &eQ);

          // Multiply Omega by eQ = (Omega^dag Omega)^(-1/2)
          mult_su3_nn(&Omega, &eQ, &(s-> hyplink1[hyp1ind[dir2][dir3][dir]]));
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp2_mcrg(int num, int block) {
  register int dir, dir2, i;
  register site *s;
  Real f[3], ftmp1, ftmp2;
  complex tc;
  su3_matrix tmat, Omega, eQ, Id, Q, Q2;

  ftmp1 = alpha_smear[1] / (4.0 * (1.0 - alpha_smear[1]));
  ftmp2 = 1.0 - alpha_smear[1];

  FORALLUPDIR(dir) {
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      // Compute the staple
      staple2_mcrg(dir, dir2, block);
      FORALLSITES(i, s) {
        // Make Omega
        scalar_mult_add_su3_matrix(&(s->link[dir]), &(s->tempmat2),
                                   ftmp1, &Q);
        scalar_mult_su3_matrix(&Q, ftmp2, &Omega);
        mult_su3_an(&Omega,&Omega,&Q);
        scalar_add_diag_su3(&Q, IR_STAB);
#ifndef NHYP_DEBUG
        compute_fhb(&Q, f, NULL, 0);
#else
        compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

        // Compute Q**2
        mult_su3_nn(&Q, &Q, &Q2);

        // Compute Q^(-1/2) via Eq. (3.8)
        tc=cmplx(f[0], 0);
        diag_su3(&Id, &tc);
        scalar_mult_add_su3_matrix(&Id, &Q, f[1], &tmat);
        scalar_mult_add_su3_matrix(&tmat, &Q2, f[2], &eQ);

        // Multiply Omega by eQ = (Omega^dag Omega)^(-1/2)
        mult_su3_nn(&Omega, &eQ, &(s->hyplink2[hyp2ind[dir2][dir]]));
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp3_mcrg(int num, int block) {
  register int dir, i;
  register site *s;
  Real f[3], ftmp1, ftmp2;
  complex tc;
  su3_matrix tmat, Omega, eQ, Id, Q, Q2;

  ftmp1 = alpha_smear[0]/ (6.0 * (1.0 - alpha_smear[0]));
  ftmp2 = 1.0 - alpha_smear[0];

  FORALLUPDIR(dir) {
    // Compute the staple
    staple3_mcrg(dir, block);
    FORALLSITES(i, s) {
      // Make Omega
      scalar_mult_add_su3_matrix(&(s->link[dir]),&(s->tempmat2),
                                 ftmp1, &Q);
      scalar_mult_su3_matrix(&Q, ftmp2, &Omega);
      mult_su3_an(&Omega, &Omega, &Q);
      scalar_add_diag_su3(&Q, IR_STAB);
#ifndef NHYP_DEBUG
      compute_fhb(&Q, f, NULL, 0);
#else
      compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

      // Compute Q**2
      mult_su3_nn(&Q, &Q, &Q2);

      // compute Q^(-1/2) via Eq. (3.8)
      tc = cmplx(f[0], 0.0);
      diag_su3(&Id, &tc);
      scalar_mult_add_su3_matrix(&Id, &Q, f[1], &tmat);
      scalar_mult_add_su3_matrix(&tmat, &Q2, f[2], &eQ);

      // Multiply Omega by eQ = (Omega^dag Omega)^(-1/2)
      mult_su3_nn(&Omega, &eQ, &(s->link[dir]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp_mcrg(int num, int block) {
  register int i, j, k, count = 0, dir, dir2, dir3;
  register site *s;

  // First step is to fill the index arrays
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++)
      hyp2ind[i][j] = 999;
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (i != j)
        hyp2ind[i][j] = count++;
    }
  }

  count = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 4; k++)
        hyp1ind[i][j][k] = 999;
    }
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < j; k++) {
        if (j != i && k != i) {
          hyp1ind[j][k][i] = count;
          hyp1ind[k][j][i] = count++;
        }
      }
    }
  }

  FORALLSITES(i, s) {
    FORALLUPDIR(dir) {
      FORALLUPDIR(dir2) {
        if (dir2 != dir) {
          for (dir3 = dir2 + 1; dir3 <= TUP; dir3++) {
            if (dir != dir3) {
              su3mat_copy(&(s->link[dir]),
                          &(s->hyplink1[hyp1ind[dir2][dir3][dir]]));
            }
          }
        }
      }
    }
  }

  FORALLSITES(i, s) {
    FORALLUPDIR(dir) {
      FORALLUPDIR(dir2) {
        if (dir2 != dir) {
          su3mat_copy(&(s->link[dir]),
                      &(s->hyplink2[hyp2ind[dir2][dir]]));
        }
      }
    }
  }

  // Now we are ready to go
  node0_printf("Alpha %.4g %.4g %.4g\n",
               alpha_smear[0], alpha_smear[1], alpha_smear[2]);

  if (alpha_smear[2] != 0)
    block_nhyp1_mcrg(num, block);
  if (alpha_smear[1] != 0)
    block_nhyp2_mcrg(num, block);

  block_nhyp3_mcrg(num, block);
}
// -----------------------------------------------------------------

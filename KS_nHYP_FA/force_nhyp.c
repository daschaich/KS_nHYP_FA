// -----------------------------------------------------------------
// Derivative of the nHYP link w.r.t. the constituent links
// Two contributions:
// 1) The (thin) link itself : added to Sigma in Sigma_update1
// 2) The (fat) staples      : made in compute-sigma23 and put in SimgaH
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Helper routine for Sigma_update1
void make_2hermitian(su3_matrix *A) {
  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = i; j < 3; j++) {
      A->e[i][j].real =  A->e[i][j].real + A->e[j][i].real;
      A->e[i][j].imag =  A->e[i][j].imag - A->e[j][i].imag;
      A->e[j][i].real =  A->e[i][j].real;
      A->e[j][i].imag = -A->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute Gamma
void Sigma_update1(const int dir, su3_matrix *sigma_off, su3_matrix *stp,
                   su3_matrix **lambda, const Real alpha1,
                   const Real alpha2, const int sigfresh) {

  register int i, j;
  register site *s;
  Real f[3], bb[3][3];
  complex traces[3], tc, tc2, tc3;
  su3_matrix Gamma, Qisqrt, Q, Q2, Omega, tmat, SigmaOmega;

  FORALLSITES(i, s) {
    // Make Omega, Q, Q^2
    Omega = stp[i];
    mult_su3_an(&Omega, &Omega, &Q);
    mult_su3_nn(&Q, &Q, &Q2);

    // Compute inverse sqrt
#ifndef NHYP_DEBUG
    compute_fhb(&Q, f, bb, 1);
#else
    compute_fhb(&Omega, &Q, f, bb, 1);
#endif

    scalar_mult_su3_matrix(&Q, f[1], &tmat);
    scalar_mult_add_su3_matrix(&tmat, &Q2, f[2], &Qisqrt);
    scalar_add_diag_su3(&Qisqrt, f[0]);

    // We'll need Sigma * Omega a few times
    mult_su3_nn(sigma_off + i, &Omega, &SigmaOmega);

    // Now the B matrices and their traces with Sigma*Omega
    tc = trace_su3(&SigmaOmega);
    tc2 = complextrace_su3(&Q, &SigmaOmega);
    tc3 = complextrace_su3(&Q2, &SigmaOmega);
    for (j = 0; j < 3; j++) {
      traces[j].real = tc.real * bb[j][0] + tc2.real * bb[j][1]
                                          + tc3.real * bb[j][2];
      traces[j].imag = tc.imag * bb[j][0] + tc2.imag * bb[j][1]
                                          + tc3.imag * bb[j][2];
    }

    // The contributions to A tr(B_i Sigma Omega) Q^(i)
    c_scalar_mult_su3mat(&Q, &traces[1], &Gamma);
    c_scalar_mult_sum_su3mat(&Q2, &traces[2], &Gamma);
    c_scalar_add_diag_su3(&Gamma, &traces[0]);

    // The terms proportional to f_i
    scalar_mult_sum_su3_matrix(&SigmaOmega, f[1], &Gamma);
    mult_su3_nn(&SigmaOmega, &Q, &tmat);
    scalar_mult_sum_su3_matrix(&tmat, f[2], &Gamma);
    mult_su3_nn(&Q, &SigmaOmega, &tmat);
    scalar_mult_sum_su3_matrix(&tmat, f[2], &Gamma);

    // Gamma = (A + Adag)Qdag + Q^{-1/2}Sigma
    make_2hermitian(&Gamma);
    mult_su3_na(&Gamma, &Omega, &tmat);

    mult_su3_nn(&Qisqrt, sigma_off + i, &Gamma);
    add_su3_matrix(&Gamma, &tmat, &Gamma);
    scalar_mult_su3_matrix(&Gamma, alpha2, &tmat);
    su3_adjoint(&tmat, lambda[dir] + i);

    // The derivative which contributes to the new global Sigma
    // If this is the first level, then Sigma has to be initialized.
    // On later levels, we accumulate contributions
    if (sigfresh == 0)
        scalar_mult_su3_matrix(&Gamma, alpha1, Sigma[dir]+i);
    else
        scalar_mult_add_su3_matrix(Sigma[dir] + i, &Gamma, alpha1,
                                   Sigma[dir] + i);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute next-level Sigma (sig)
// lnk is the links the routine operates on
// dir1 is the direction of the link and the 'main' direction of Sigma
// lambda = Gamma^dag
void compute_sigma23(su3_matrix *sig, su3_matrix *lnk1, su3_matrix *lnk2,
                     su3_matrix *lambda1, su3_matrix *lambda2,
                     int dir1, int dir2, int parity) {

  register int i;
  register site *st;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4;
  msg_tag *tag5, *tag6, *tag7, *tag8, *tag9;
  su3_matrix tmat, tmat2;
  int disp[4];    // Displacement vector for general gather

  // Loop over other directions
  // Compute force from plaquettes in the dir1, dir2 plane

  // Displacement vector for bb_link 2 sites away
  FORALLUPDIR(i)
    disp[i] = 0;
  disp[dir1] = 1;
  disp[dir2] = -1;

  // Get link[dir2] from displacement dir1 - dir2
  tag4 = start_general_gather_field(lnk1, sizeof(su3_matrix), disp,
                                    parity, gen_pt[4]);

  // Get link[dir2] from direction dir1
  tag0 = start_gather_field(lnk1, sizeof(su3_matrix), dir1,
                            parity, gen_pt[0]);

  // Get link[dir1] from direction dir2
  tag1 = start_gather_field(lnk2, sizeof(su3_matrix), dir2,
                            parity, gen_pt[1]);

  // Get link[dir2] from direction -dir2
  tag2 = start_gather_field(lnk1, sizeof(su3_matrix), OPP_DIR(dir2),
                            parity, gen_pt[2]);

  // Get link[dir1] from direction -dir2
  tag3 = start_gather_field(lnk2, sizeof(su3_matrix), OPP_DIR(dir2),
                            parity, gen_pt[3]);

  // Get LambdaU[dir2] from direction dir1
  tag5 = start_gather_field(lambda2, sizeof(su3_matrix), dir1,
                            parity, gen_pt[5]);

  // Get LambdaU[dir1] from direction dir2
  tag6 = start_gather_field(lambda1, sizeof(su3_matrix), dir2,
                            parity, gen_pt[6]);

  // Get LambdaU[dir2] from direction -dir2
  tag7 = start_gather_field(lambda2, sizeof(su3_matrix), OPP_DIR(dir2),
                            parity, gen_pt[7]);

  // Get LambdaU[dir1] from direction -dir2
  tag8 = start_gather_field(lambda1, sizeof(su3_matrix), OPP_DIR(dir2),
                            parity, gen_pt[8]);

  wait_general_gather(tag4);

  // Get LambdaU[dir2] from displacement dir1 - dir2
  tag9 = start_general_gather_field(lambda2, sizeof(su3_matrix), disp,
                                    parity, gen_pt[9]);

  wait_gather(tag0);
  wait_gather(tag1);
  wait_gather(tag2);
  wait_gather(tag3);
  wait_gather(tag5);
  wait_gather(tag6);
  wait_gather(tag7);
  wait_gather(tag8);
  wait_general_gather(tag9);

  FORSOMEPARITY(i, st, parity) {
    // Term 1
    mult_su3_nn(lambda2 + i, (su3_matrix *)gen_pt[1][i], &tmat);
    mult_su3_na((su3_matrix *)gen_pt[0][i], &tmat, sig + i);

    // Term 2
    mult_su3_nn((su3_matrix *)gen_pt[8][i],
                (su3_matrix *)gen_pt[4][i], &tmat);
    mult_su3_an(&tmat, (su3_matrix *)gen_pt[2][i], &tmat2);
    add_su3_matrix(sig + i, &tmat2, sig + i);

    // Term 3
    mult_su3_nn((su3_matrix *)gen_pt[3][i],
                (su3_matrix *)gen_pt[9][i], &tmat);
    mult_su3_an(&tmat, (su3_matrix *)gen_pt[2][i], &tmat2);
    add_su3_matrix(sig + i, &tmat2, sig + i);

    // Term 4
    mult_su3_nn((su3_matrix *)gen_pt[3][i],
                (su3_matrix *)gen_pt[4][i], &tmat);
    mult_su3_an(&tmat, (su3_matrix *)gen_pt[7][i], &tmat2);
    add_su3_matrix(sig + i, &tmat2, sig + i);

    // Term 5
    mult_su3_na((su3_matrix *)gen_pt[5][i],
                (su3_matrix *)gen_pt[1][i], &tmat);
    mult_su3_na(&tmat, lnk1+i, &tmat2);
    add_su3_matrix(sig + i, &tmat2, sig + i);

    // Term 6
    mult_su3_na((su3_matrix *)gen_pt[0][i],
                (su3_matrix *)gen_pt[6][i], &tmat);
    mult_su3_na(&tmat, lnk1 + i, &tmat2);
    add_su3_matrix(sig + i, &tmat2, sig + i);
  }

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
  cleanup_gather(tag3);
  cleanup_general_gather(tag4);
  cleanup_gather(tag5);
  cleanup_gather(tag6);
  cleanup_gather(tag7);
  cleanup_gather(tag8);
  cleanup_general_gather(tag9);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Third-level force
void nhyp_force3(int dir3, int dir2) {
  register int i, dir, dir1, dir4;
  register site *s;

  FORALLUPDIR(dir) {
    if (dir == dir2 || dir == dir3)
      continue;
    FORALLUPDIR(dir4) {
      if (dir4 != dir && dir4 != dir3 && dir4 != dir2)
        break;
    }
    Sigma_update1(dir, SigmaH[dir], Staple1[dir4][dir], Lambda2,
                  1.0 - alpha_smear[2], alpha_smear[2] / 2.0, 1);
  }

  FORALLUPDIR(dir) {
    if (dir == dir2 || dir == dir3)
      continue;
    FORALLUPDIR(dir1) {
      if (dir1 == dir || dir1 == dir2 || dir1 == dir3)
        continue;
      compute_sigma23(tempmat, gauge_field_thin[dir1], gauge_field_thin[dir],
                      Lambda2[dir], Lambda2[dir1], dir, dir1, EVENANDODD);

      FORALLSITES(i, s)
        add_su3_matrix(Sigma[dir] + i, tempmat + i, Sigma[dir] + i);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Second-level force
void nhyp_force2(int dir2) {
  register int i;
  register site *s;
  int dir, dir3, dir4;
  int imap[4][4] = {{0, 0, 1, 2}, {0, 0, 0, 3}, {1, 0, 0, 0}, {2, 3, 0, 0}};
  int iimap;

  // dir3 is the main direction of the twice-smeared hyplink2
  // dir2 is the secondary direction of the twice-smeared hyplink2
  // dir is the main direction of the once-smeared link
  FORALLUPDIR(dir3) {
    if (dir3 != dir2) {
      Sigma_update1(dir3, SigmaH[dir3], Staple2[dir2][dir3], Lambda1,
                    1.0 - alpha_smear[1], alpha_smear[1] / 4.0, 1);
    }
  }

  FORALLUPDIR(dir3) {
    if (dir3 == dir2)
      continue;
    FORALLUPDIR(dir) {
      if (dir3 == dir || dir == dir2)
        continue;
      FORALLUPDIR(dir4) {
        if (dir4 != dir && dir4 != dir2 && dir4 != dir3)
          break;
      }
      compute_sigma23(SigmaH[dir], hyplink1[dir4][dir3], hyplink1[dir4][dir],
                      Lambda1[dir], Lambda1[dir3], dir, dir3, EVENANDODD);
    }

    // This part is really awkward
    // nhyp_force3 is symmetric in the arguments, but the input is not
    // So add the dir2, dir3 and the dir3, dir2 terms
    // To do so, one has to sort them...
    // To save memory, store only the upper triangular part of the 4x4 matrix
    // For that, only a 4 'vector' is necessary
    // Which field to use is given by the imap array hard-coded above
    iimap = imap[dir2][dir3];
    if (dir2 < dir3) {
      FORALLSITES(i, s) {
        FORALLUPDIR(dir)
          su3mat_copy(SigmaH[dir] + i, SigmaH2[iimap][dir] + i);
      }
    }
    else {
      FORALLSITES(i, s) {
        FORALLUPDIR(dir) {
          add_su3_matrix(SigmaH[dir] + i, SigmaH2[iimap][dir] + i,
                         SigmaH[dir] + i);
        }
      }
      nhyp_force3(dir2, dir3);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// First-level force
void nhyp_force1() {
  int dir, dir2;

  // Loop over the link directions,
  // compute sigma from the link itself,
  // construct the Lambdas
  FORALLUPDIR(dir) {
    Sigma_update1(dir, Sigma[dir], Staple3[dir], LambdaU,
                  1.0 - alpha_smear[0], alpha_smear[0] / 6.0, 0);
  }

  // Construct field_offsets pointing to the links in dir2 direction
  // where dir1 is excluded.  Here, this makes no sense, however this
  // mechanism makes compute_sigma23 re-useable on each level
  FORALLUPDIR(dir2) {
    FORALLUPDIR(dir) {
      if (dir == dir2)
        continue;
      compute_sigma23(SigmaH[dir], hyplink2[dir][dir2], hyplink2[dir2][dir],
                      LambdaU[dir], LambdaU[dir2], dir, dir2, EVENANDODD);
      }
    nhyp_force2(dir2);
  }
}
// -----------------------------------------------------------------

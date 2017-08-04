// -----------------------------------------------------------------
// Update the momentum matrices
// Use tempmat and tempmat2 for temporary storage
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta with the gauge force
double gauge_force(Real eps) {
  register int i, dir, dir2;
  register site *s;
  register Real eb3 = eps * beta / 3.0;
  msg_tag *tag0, *tag1, *tag2;
  int start;
  matrix tmat, tmat2;
  complex tc;
  double norm = 0;

  // Loop over directions, update mom[dir]
  FORALLUPDIR(dir) {
    start = 1; // Indicates staple sum (in tempmat2) not initialized

    // Loop over other directions
    // Compute force from plaquettes in the dir, dir2 plane
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      // Get link[dir2] from direction dir
      tag0 = start_gather_site(F_OFFSET(link[dir2]),
                               sizeof(matrix),
                               dir, EVENANDODD, gen_pt[0]);

      // Start gather for the "upper staple"
     tag2 = start_gather_site(F_OFFSET(link[dir]),
                              sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[2]);

      // Begin the computation "at the dir2DOWN point"
      // We will later gather the intermediate result "to the home point"
      wait_gather(tag0);
      FORALLSITES(i, s) {
        mult_su3_an(&(s->link[dir2]), &(s->link[dir]), &tmat);
        mult_su3_nn(&tmat, (matrix *)gen_pt[0][i],
                    (matrix *)&(tempmat[i]));
      }

      // Gather lower staple "up to home site"
      tag1 = start_gather_field(tempmat, sizeof(matrix),
                                OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

      // The "upper" staple
      // One of the links has already been gathered,
      // since it was used in computing
      // the "lower" staple of the site above (in dir2)
      wait_gather(tag2);
      if (start) {  // Initialize staple sum in tempmat2
        FORALLSITES(i, s) {
          mult_su3_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
          mult_su3_na(&tmat, (matrix *)gen_pt[0][i], &(tempmat2[i]));
          mult_su3_na(&(s->link[dir]), &(tempmat2[i]), &tmat);
          tc = trace_su3(&tmat);
          tc.real = 1.0 - 2.0 * beta_a / 3.0 * tc.real;
          tc.imag = -2.0 * beta_a / 3.0 * tc.imag;

          c_scalar_mult_mat(&(tempmat2[i]), &tc, &tmat);
          mat_copy(&tmat, &(tempmat2[i]));
        }
        start = 0;
      }
      else {
        FORALLSITES(i, s) {
          mult_su3_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
          mult_su3_na(&tmat, (matrix *)gen_pt[0][i], &tmat2);
          mult_su3_na(&(s->link[dir]), &tmat2, &tmat);
          tc = trace_su3(&tmat);
          tc.real = 1.0 - 2.0 * beta_a / 3.0 * tc.real;
          tc.imag *= -2.0 * beta_a / 3.0;

          c_scalar_mult_sum_mat(&tmat2, &tc, &(tempmat2[i]));
        }
      }

      wait_gather(tag1);
      FORALLSITES(i, s) {
        mult_su3_na(&(s->link[dir]), (matrix *)gen_pt[1][i], &tmat);
        tc = trace_su3(&tmat);
        tc.real = 1.0 - 2.0 * beta_a / 3.0 * tc.real;
        tc.imag *= -2.0 * beta_a / 3.0;

        c_scalar_mult_sum_mat((matrix *)gen_pt[1][i], &tc,
                                 &(tempmat2[i]));
      }
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
    }

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, s) {
      mult_su3_na(&(s->link[dir]), &(tempmat2[i]), &tmat);
      uncompress_anti_hermitian(&(s->mom[dir]), &tmat2);
      scalar_mult_add_matrix(&tmat2, &tmat, eb3, &(tempmat2[i]));
      make_anti_hermitian(&(tempmat2[i]), &(s->mom[dir]));
      norm += (double)realtrace_su3(&tmat, &tmat);
    }
  }

  g_doublesum(&norm);
  return eb3 * sqrt(norm) / (double)volume;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Based on v6 code; ignore diagonal term since killed by Lie derivative
// Assumes that the CG has been run, with the answer in psi
double fermion_force(int level, Real eps) {
  register int i, dir;
  register site *s;
  int j;
  Real ferm_epsilon = 2 * eps;
  double norm = 0;
  msg_tag *tag0, *tag1;
  anti_hermitmat ahtmp;
  matrix tmat, tmat2;
  vector tvec;

  // Zero the force collectors
  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      clear_mat(&(Sigma[dir][i]));
  }

  // Rescale outer level to account for Hasenbusch preconditioning
  if (num_masses == 2 && level == 0)
    ferm_epsilon *= 4.0 * (MH * MH - mass * mass);

  // Hit psi with dslash -- off-diagonal piece of M^dag?
  just_fatten();
  for (j = 0; j < full_fields; j++)
    dslash(F_OFFSET(psi[j][level]), F_OFFSET(ttt[j][level]), EVENANDODD);
  if (half_fields == 1) {
    j = full_fields;
    dslash(F_OFFSET(psi[j][level]), F_OFFSET(ttt[j][level]), ODD);
  }

  // Fermionic contribution to the force
  // su3_projector overwrites its target, so add via temporary tmat2
  for (j = 0; j < full_fields; j++) {
    // For all sites, gather first ttt[j][level] before entering loop
    tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(vector),
                             XUP, EVENANDODD, gen_pt[0]);

    FORALLUPDIR(dir) {
      // For all sites, gather psi[j][level]
      tag1 = start_gather_site(F_OFFSET(psi[j][level]), sizeof(vector),
                               dir, EVENANDODD, gen_pt[1]);

      wait_gather(tag0);
      FORALLSITES(i, s) {
        mult_mat_vec(&(s->link[dir]), (vector *)gen_pt[0][i], &tvec);
        su3_projector(&tvec, &(s->psi[j][level]), &tmat2);
        sum_matrix(&tmat2, &(Sigma[dir][i]));
      }
      cleanup_gather(tag0);

      // For all sites, gather next ttt[j][level]
      if (dir < TUP)
        tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(vector),
                                 dir + 1, EVENANDODD, gen_pt[0]);

      wait_gather(tag1);
      FORALLSITES(i, s) {
        mult_mat_vec(&(s->link[dir]), (vector *)gen_pt[1][i], &tvec);
        su3_projector(&(s->ttt[j][level]), &tvec, &tmat2);
        sum_matrix(&tmat2, &(Sigma[dir][i]));
      }
      cleanup_gather(tag1);
    }
  }
  if (half_fields == 1) {
    j = full_fields;
    // For even sites, gather first ttt[j][level] before entering loop
    tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(vector),
                             XUP, EVEN, gen_pt[0]);

    FORALLUPDIR(dir) {
      // For odd sites, gather psi[j][level]
      tag1 = start_gather_site(F_OFFSET(psi[j][level]), sizeof(vector),
                               dir, ODD, gen_pt[1]);

      wait_gather(tag0);
      FOREVENSITES(i, s) {
        mult_mat_vec(&(s->link[dir]), (vector *)gen_pt[0][i], &tvec);
        su3_projector(&tvec, &(s->psi[j][level]), &tmat2);
        sum_matrix(&tmat2, &(Sigma[dir][i]));
      }
      cleanup_gather(tag0);

      // For even sites, gather next ttt[j][level]
      if (dir < TUP)
        tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(vector),
                                 dir + 1, EVEN, gen_pt[0]);

      wait_gather(tag1);
      FORODDSITES(i, s) {
        mult_mat_vec(&(s->link[dir]), (vector *)gen_pt[1][i], &tvec);
        su3_projector(&(s->ttt[j][level]), &tvec, &tmat2);
        sum_matrix(&tmat2, &(Sigma[dir][i]));
      }
      cleanup_gather(tag1);
    }
  }
  leanlinks();  // Restore the thin links

  // Multiply a HYP field dagger from the left on the force
  FORALLUPDIR(dir) {
    FORALLSITES(i, s) {
      mult_su3_an(&(gauge_field[dir][i]), &(Sigma[dir][i]), &tmat);
      mat_copy(&tmat, &(Sigma[dir][i]));
    }
  }

  // More chain rule: dH/dU^T = dV/dU^T dH/dV^T
  // nhyp_force receives and returns the force in Sigma
  // The order of contributions is crucial!
  for (j = Nsmear; j > 1; j--) {
    // Save j-1-smeared links in gauge_field_temp
    block_N_fatten(j - 1);
    leanlinks();
    FORALLUPDIR(dir) {
      FORALLSITES(i, s)
        mat_copy(&(gauge_field[dir][i]), &(gauge_field_temp[dir][i]));
    }

    block_N_fatten(j);   // gauge_field holds j-smeared links
    // gauge_field_thin receives j-1-smeared links from gauge_field_temp
    FORALLUPDIR(dir) {
      FORALLSITES(i, s)
        mat_copy(&(gauge_field_temp[dir][i]), &(gauge_field_thin[dir][i]));
    }
    leanlinks();          // s->link points to j-1-smeared links
    nhyp_force1();

    // Restore unsmeared links to gauge_field_thin
    FORALLUPDIR(dir) {
      FORALLSITES(i, s)
        mat_copy(&(gauge_field_save[dir][i]), &(gauge_field_thin[dir][i]));
    }
    leanlinks();            // s->link points to unsmeared links
  }
  // Once-smeared contribution
  block_and_fatten();   // Overwrite block_N_fatten above
  leanlinks();
  nhyp_force1();

  if (Nsmear > 1) {   // Get the right number of smearings in gauge_field
    block_N_fatten(Nsmear);   // gauge_field holds j-smeared links
    leanlinks();
  }

  // Finally we can update the momenta
  FORALLUPDIR(dir) {
    FORALLSITES(i, s) {
      uncompress_anti_hermitian(&(s->mom[dir]), &tmat2);
      mult_su3_nn(&(gauge_field_thin[dir][i]), &(Sigma[dir][i]), &tmat);
      make_anti_hermitian(&tmat, &ahtmp);

      uncompress_anti_hermitian(&ahtmp, &tmat);
      // tmat2 += ferm_epsilon * tmat
      scalar_mult_sum_matrix(&tmat, ferm_epsilon, &tmat2);
      make_anti_hermitian(&tmat2, &(s->mom[dir]));
      norm += (double)realtrace_su3(&tmat, &tmat);
    }
  }
  g_doublesum(&norm);
  return 2.0 * ferm_epsilon * sqrt(norm) / (double)volume;
}
// -----------------------------------------------------------------

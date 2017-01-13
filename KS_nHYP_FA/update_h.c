// -----------------------------------------------------------------
// Update the momentum matrices
// Includes and definitions
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta with the gauge force
double gauge_force(Real eps) {
  register int i, dir, dir2;
  register site *st;
  register Real eb3 = eps * beta / 3.0;
  msg_tag *tag0, *tag1, *tag2;
  int start;
  su3_matrix tmat, tmat2;
  complex tc;
  double norm = 0;

  // Loop over directions, update mom[dir]
  FORALLUPDIR(dir) {
    start = 1; // Indicates staple sum not initialized
    FORALLSITES(i, st)
      clear_su3mat(&(st->staple));    // Initialize staple

    // Loop over other directions
    // Compute force from plaquettes in the dir, dir2 plane
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      // Get link[dir2] from direction dir
      tag0 = start_gather_site(F_OFFSET(link[dir2]),
                               sizeof(su3_matrix),
                               dir, EVENANDODD, gen_pt[0]);

      // Start gather for the "upper staple"
     tag2 = start_gather_site(F_OFFSET(link[dir]),
                              sizeof(su3_matrix),
                              dir2, EVENANDODD, gen_pt[2]);

      // Begin the computation "at the dir2DOWN point"
      // We will later gather the intermediate result "to the home point"
      wait_gather(tag0);
      FORALLSITES(i, st) {
        mult_su3_an(&(st->link[dir2]), &(st->link[dir]), &tmat);
        mult_su3_nn(&tmat, (su3_matrix *)gen_pt[0][i],
                    (su3_matrix *)&(st->tempmat1));
      }

      // Gather lower staple "up to home site"
      tag1 = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix),
                               OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

      // The "upper" staple
      // One of the links has already been gathered,
      // since it was used in computing
      // the "lower" staple of the site above (in dir2)
      wait_gather(tag2);
      if (start) {  // This is the first contribution to staple
        FORALLSITES(i, st) {
          mult_su3_nn(&(st->link[dir2]), (su3_matrix *)gen_pt[2][i], &tmat);
          mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &(st->staple));
          mult_su3_na(&(st->link[dir]), &(st->staple), &tmat);
          tc = trace_su3(&tmat);
          tc.real = 1.0 - 2.0 * beta_a / 3.0 * tc.real;
          tc.imag = -2.0 * beta_a / 3.0 * tc.imag;

          c_scalar_mult_su3mat(&(st->staple), &tc, &tmat);
          su3mat_copy(&tmat, &(st->staple));
        }
        start = 0;
      }
      else {
        FORALLSITES(i, st) {
          mult_su3_nn(&(st->link[dir2]), (su3_matrix *)gen_pt[2][i], &tmat);
          mult_su3_na(&tmat, (su3_matrix *)gen_pt[0][i], &tmat2);
          mult_su3_na(&(st->link[dir]), &tmat2, &tmat);
          tc = trace_su3(&tmat);
          tc.real = 1.0 - 2.0 * beta_a / 3.0 * tc.real;
          tc.imag = -2.0 * beta_a / 3.0 * tc.imag;

          c_scalar_mult_add_su3mat(&(st->staple), &tmat2, &tc,
                                   &(st->staple));
        }
      }

      wait_gather(tag1);
      FORALLSITES(i, st) {
        mult_su3_na(&(st->link[dir]), (su3_matrix *)gen_pt[1][i], &tmat);
        tc = trace_su3(&tmat);
        tc.real = 1 - 2 * beta_a / 3.0 * tc.real;
        tc.imag = -2 * beta_a / 3.0 * tc.imag;

        c_scalar_mult_add_su3mat(&(st->staple),
                                 (su3_matrix *)gen_pt[1][i], &tc,
                                 &(st->staple));
      }
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
    }

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, st) {
      mult_su3_na(&(st->link[dir]), &(st->staple), &tmat);
      uncompress_anti_hermitian(&(st->mom[dir]), &tmat2);
      scalar_mult_add_su3_matrix(&tmat2, &tmat, eb3, &(st->staple));
      make_anti_hermitian(&(st->staple), &(st->mom[dir]));
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
  register site *st;
  int j;
  Real ferm_epsilon = 2 * eps;
  double norm = 0;
  msg_tag *tag0, *tag1;
  anti_hermitmat ahtmp;
  su3_matrix tmat, tmat2, tmat3;
  su3_vector tvec;

  // Zero the force collectors
  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, st)
      clear_su3mat(&Sigma[dir][i]);
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
    tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(su3_vector),
                             XUP, EVENANDODD, gen_pt[0]);

    for (dir = XUP; dir <= TUP; dir++) {
      // For all sites, gather psi[j][level]
      tag1 = start_gather_site(F_OFFSET(psi[j][level]), sizeof(su3_vector),
                               dir, EVENANDODD, gen_pt[1]);

      wait_gather(tag0);
      FORALLSITES(i, st) {
        mult_su3_mat_vec(&(st->link[dir]), (su3_vector *)gen_pt[0][i], &tvec);
        su3_projector(&tvec, &(st->psi[j][level]), &tmat2);
        add_su3_matrix(&Sigma[dir][i], &tmat2, &Sigma[dir][i]);
      }
      cleanup_gather(tag0);

      // For all sites, gather next ttt[j][level]
      if (dir < TUP)
        tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(su3_vector),
                                 dir + 1, EVENANDODD, gen_pt[0]);

      wait_gather(tag1);
      FORALLSITES(i, st) {
        mult_su3_mat_vec(&(st->link[dir]), (su3_vector *)gen_pt[1][i], &tvec);
        su3_projector(&(st->ttt[j][level]), &tvec, &tmat2);
        add_su3_matrix(&Sigma[dir][i], &tmat2, &Sigma[dir][i]);
      }
      cleanup_gather(tag1);
    }
  }
  if (half_fields == 1) {
    j = full_fields;
    // For even sites, gather first ttt[j][level] before entering loop
    tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(su3_vector),
                             XUP, EVEN, gen_pt[0]);

    for (dir = XUP; dir <= TUP; dir++) {
      // For odd sites, gather psi[j][level]
      tag1 = start_gather_site(F_OFFSET(psi[j][level]), sizeof(su3_vector),
                               dir, ODD, gen_pt[1]);

      wait_gather(tag0);
      FOREVENSITES(i, st) {
        mult_su3_mat_vec(&(st->link[dir]), (su3_vector *)gen_pt[0][i], &tvec);
        su3_projector(&tvec, &(st->psi[j][level]), &tmat2);
        add_su3_matrix(&Sigma[dir][i], &tmat2, &Sigma[dir][i]);
      }
      cleanup_gather(tag0);

      // For even sites, gather next ttt[j][level]
      if (dir < TUP)
        tag0 = start_gather_site(F_OFFSET(ttt[j][level]), sizeof(su3_vector),
                                 dir + 1, EVEN, gen_pt[0]);

      wait_gather(tag1);
      FORODDSITES(i, st) {
        mult_su3_mat_vec(&(st->link[dir]), (su3_vector *)gen_pt[1][i], &tvec);
        su3_projector(&(st->ttt[j][level]), &tvec, &tmat2);
        add_su3_matrix(&Sigma[dir][i], &tmat2, &Sigma[dir][i]);
      }
      cleanup_gather(tag1);
    }
  }
  leanlinks();  // Restore the thin links

  // Multiply a HYP field dagger from the left on the force
  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, st) {
      mult_su3_an(gauge_field[dir] + i, Sigma[dir] + i, &tmat);
      su3mat_copy(&tmat, Sigma[dir] + i);
    }
  }

  // More chain rule: dH/dU^T = dV/dU^T dH/dV^T
  // nhyp_force receives and returns the force in Sigma
  // The order of contributions is crucial!
  for (j = Nsmear; j > 1; j--) {
    // Save j-1-smeared links in gauge_field_temp
    block_N_fatten(j - 1);
    leanlinks();
    for (dir = XUP; dir <= TUP; dir++) {
      FORALLSITES(i, st)
        su3mat_copy(gauge_field[dir] + i, gauge_field_temp[dir] + i);
    }

    block_N_fatten(j);   // gauge_field holds j-smeared links
    // gauge_field_thin receives j-1-smeared links from gauge_field_temp
    for (dir = XUP; dir <= TUP; dir++) {
      FORALLSITES(i, st)
        su3mat_copy(gauge_field_temp[dir] + i, gauge_field_thin[dir] + i);
    }
    leanlinks();          // s->link points to j-1-smeared links
    nhyp_force1();

    // Restore unsmeared links to gauge_field_thin
    for (dir = XUP; dir <= TUP; dir++) {
      FORALLSITES(i, st)
        su3mat_copy(gauge_field_save[dir] + i, gauge_field_thin[dir] + i);
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
  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, st) {
      uncompress_anti_hermitian(&(st->mom[dir]), &tmat2);
      mult_su3_nn(gauge_field_thin[dir] + i, Sigma[dir] + i, &tmat);
      make_anti_hermitian(&tmat, &ahtmp);

      uncompress_anti_hermitian(&ahtmp, &tmat);
      // tmat3 = tmat2 + ferm_epsilon * tmat
      scalar_mult_add_su3_matrix(&tmat2, &tmat, ferm_epsilon, &tmat3);
      make_anti_hermitian(&tmat3, &(st->mom[dir]));
      norm += (double)realtrace_su3(&tmat, &tmat);
    }
  }
  g_doublesum(&norm);
  return 2.0 * ferm_epsilon * sqrt(norm) / (double)volume;
}
// -----------------------------------------------------------------

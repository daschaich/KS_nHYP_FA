// -----------------------------------------------------------------
// Main procedure for SU(3) Wilson flow
#define CONTROL
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void mcrg_block(Real t, int blmax) {
  register int i, dir;
  register site *s;
  int j, k, bl;
  double dtime, rr[10];
  complex plp, xplp;
  dtime = -dclock();

  // Save the original links
  FORALLSITES(i, s) {
    FORALLUPDIR(dir)
      su3mat_copy(&(s->link[dir]), &(s->link[dir + 8]));
  }

  // Set up loop tables
  make_loop_table2();

  node0_printf("Simple blocking from twice-nHYP smeared links\n");
  blocked_gauge_loops(0, rr);
  for (j = 0; j < nloop; j++) {
    for (k = 0; k < nreps; k++)
     node0_printf("LOOPS %g %d %d 0 0 %.8g\n",
                  t, j, k, rr[j + k * nloop] / volume / loop_num[j]);
  }

  // Checked that this reproduces the original ploop() result
  plp = blocked_ploop(0, TUP);
  xplp = blocked_ploop(0, XUP);
  node0_printf("POLYA ORIG %g %.6g %.6g %.6g %.6g\n",
               t, (double)plp.real, (double)plp.imag,
               (double)xplp.real, (double)xplp.imag);

  // Loop over blocking levels (automatically determined)
  for (bl = 1; bl <= blmax; bl++) {
    // nHYP smear twice with appropriate outer alpha
    block_nhyp_mcrg(0, bl);
    block_nhyp_mcrg(0, bl);
    // Finally block
    // The first argument 0 only does the center staple
    block_mcrg(0, bl);

    blocked_gauge_loops(bl, rr);
    for (j = 0; j < nloop; j++) {
      for (k = 0; k < nreps; k++)
        rr[j + k * nloop] /= (volume * loop_num[j]);
    }
    for (j = 0; j < nloop; j++) {
      for (k = 0; k < nreps; k++)
        node0_printf("LOOPS %g %d %d %d %.4g %.8g\n",
                     t, j, k, bl, alpha_smear[0], rr[j + k * nloop]);
    }
    // Calculate and print blocked Polyakov loops
    plp = blocked_ploop(bl, TUP);
    xplp = blocked_ploop(bl, XUP);
    node0_printf("POLYA NHYP %g %d %.6g %.6g %.6g %.6g %.6g\n",
                 t, bl, alpha_smear[0],
                 (double)plp.real, (double)plp.imag,
                 (double)xplp.real, (double)xplp.imag);
  } // End loop over blocking levels

  // Restore the original links
  FORALLSITES(i, s) {
    FORALLUPDIR(dir)
      su3mat_copy(&(s->link[dir + 8]), &(s->link[dir]));
  }

  node0_printf("BLOCKING COMPLETED\n");
  dtime += dclock();
  node0_printf("Blocking time %.4g seconds\n", dtime);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[])  {
  register int i, dir;
  register site *s;
  int j, prompt, istep, block_count = 0, blmax = 0;
  int max_scale, eps_scale = 1, N_base = 0;     // All three always positive
  double dtime, t = 0.0, E = 0.0, tSqE = 0.0, old_tSqE, der_tSqE;
  double ssplaq, stplaq, plaq = 0, old_plaq, baseline = 0.0;
  double check, old_check, topo, old_topo;
  double slope_tSqE, slope_check, slope_topo, prev_tSqE;
  double interp_t, interp_plaq, interp_E, interp_tSqE;
  double interp_check, interp_topo, interp_der;
  su3_matrix *S[4];
  anti_hermitmat *A[4];

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Allocate fields used by integrator
  FORALLUPDIR(dir) {
    S[dir] = malloc(sites_on_node * sizeof(su3_matrix));
    A[dir] = malloc(sites_on_node * sizeof(anti_hermitmat));
  }

  // Determine maximum number of blockings from smallest dimension
  // Works even if we can only block down to odd j > 4
  if (nx < nt)
    j = nx;
  else
    j = nt;
  while (j % 2 == 0 && j > 2) {    // While j is even
    j /= 2;
    blmax++;
  }

  // Special case: measure MCRG-blocked observables before any flow
  if (num_block > 0 && tblock[0] == 0) {
    mcrg_block(tblock[block_count], blmax);
    block_count++;
  }

  // Wilson flow!
  epsilon = start_eps;
  max_scale = (int)floor(max_eps / start_eps);     // Maximum scaling factor
  for (istep = 0; fabs(t) <  fabs(tmax) - 0.5 * fabs(epsilon); istep++) {
    stout_step_rk(S, A);
    t += epsilon;

    // Find 8F_munu = sum_{clover} (U - Udag)
    // Subtract the (lattice artifact?) trace at each lattice site
    make_field_strength(F_OFFSET(link), F_OFFSET(FS));

    // Save previous data for slope and interpolation
    // old_plaq is only used for adjusting step size below
    old_tSqE = tSqE;
    old_check = check;
    old_topo = topo;
    old_plaq = plaq;

    // Compute t^2 E and its slope
    E = 0.0;
    FORALLSITES(i, s) {
      for (dir = 0; dir < 6; dir++)
        E -= (double)realtrace_su3_nn(&(s->FS[dir]), &(s->FS[dir]));
    }
    g_doublesum(&E);
    E /= (volume * 64.0); // Normalization factor of 1/8 for each F_munu
    tSqE = t * t * E;
    der_tSqE = fabs(t) * (tSqE - old_tSqE) / fabs(epsilon);
    // Any negative signs in t and epsilon should cancel out anyway...

    // Might as well extract topology
    // Normalization is 1/(64*4pi^2) again with 1/8 for each F_munu
    topo = 0.0;
    FORALLSITES(i, s) {
      topo -= (double)realtrace_su3_nn(&(s->FS[0]), &(s->FS[5])); // XYZT
      topo -= (double)realtrace_su3_nn(&(s->FS[3]), &(s->FS[2])); // XTYZ
      topo -= (double)realtrace_su3(&(s->FS[1]), &(s->FS[4]));    // XZ(YT)^dag
    }
    g_doublesum(&topo);
    topo *= 0.000395785873603;

    // Check with plaquette
    plaquette(&ssplaq, &stplaq);
    plaq = 0.5 * (ssplaq + stplaq);
    check = 12.0 * t * t * (3.0 - plaq);

    // If necessary, interpolate from previous t-eps to current t
    // before printing out results computed above
    if (eps_scale > 1) {
      slope_tSqE = (tSqE - old_tSqE) / epsilon;
      slope_check = (check - old_check) / epsilon;
      slope_topo = (topo - old_topo) / epsilon;
      prev_tSqE = old_tSqE;
      interp_t = t - epsilon;
      for (j = 0; j < eps_scale - 1; j++) {
        interp_t += start_eps;
        interp_tSqE = old_tSqE + j * slope_tSqE;
        interp_E = interp_tSqE / (interp_t * interp_t);
        interp_der = interp_t * (interp_tSqE - prev_tSqE) / start_eps;
        prev_tSqE = interp_tSqE;

        interp_check = old_check + j * slope_check;
        interp_plaq = 3.0 - check / (12.0 * interp_t * interp_t);

        interp_topo = old_topo + j * slope_topo;
        node0_printf("WFLOW %g %g %g %g %g %g %g (interp)\n",
                     interp_t, interp_plaq, interp_E, interp_tSqE,
                     interp_der, interp_check, interp_topo);
      }
    }
    node0_printf("WFLOW %g %g %g %g %g %g %g\n",
                 t, plaq, E, tSqE, der_tSqE, check, topo);

    // Do MCRG blocking at specified t
    // Use start_eps rather than epsilon to get more accurate targeting
    if (block_count < num_block
        && fabs(t + 0.5 * start_eps) >= fabs(tblock[block_count])) {
      mcrg_block(tblock[block_count], blmax);
      block_count++;
    }

    // Choose epsilon for the next step
    // For t < 5 keep epsilon = start_eps fixed
    // and set up a baseline Delta_plaq * epsilon
    // This avoids problems interpolating during the initial rise
    // The t < 5 choice may be conservative (t < 2 might suffice)
    if (fabs(t) < 5.0) {
      baseline += (plaq - old_plaq) * epsilon;
      N_base++;
    }

    // For t > 5, set epsilon = eps_scale * start_eps
    // where eps_scale = floor(baseline / Delta_plaq)
    // Any signs should cancel in the Delta_plaq factors
    // Round down if eps_scale exceeds max_scale
    // Also reduce to land on next tblock or final tmax
    else {
      // Finish setting up the baseline if this is the first t > 5
      if (N_base > 0) {
        baseline /= N_base;
        N_base = -99;
      }

      // Basic scaling to keep Delta_plaq * epsilon fixed
      eps_scale = (int)floor(baseline / ((plaq - old_plaq) * start_eps));
      if (fabs(eps_scale) > fabs(max_scale))
        eps_scale = max_scale;
      epsilon = eps_scale * start_eps;

      // Make sure we don't overshoot tmax or next tblock
      // This can probably be made more elegant
      // Need to make sure epsilon never becomes zero
      if (fabs(t + epsilon) > fabs(tmax)) {
        eps_scale = (int)floor((tmax - t) / start_eps);
        if (eps_scale < 1)
          eps_scale = 1;
        epsilon = eps_scale * start_eps;
      }
      else if (block_count < num_block
               && fabs(t + epsilon) > fabs(tblock[block_count])) {
        eps_scale = (int)floor((tblock[block_count] - t) / start_eps);
        if (eps_scale < 1)
          eps_scale = 1;
        epsilon = eps_scale * start_eps;
      }
    }
  }

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  FORALLUPDIR(dir) {
    free(S[dir]);
    free(A[dir]);
  }

  // Save lattice if requested -- no phases in links
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);
  normal_exit(0);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------

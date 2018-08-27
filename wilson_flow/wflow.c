// -----------------------------------------------------------------
// Run Wilson flow with adaptive step sizes
// Just need to pass through fields used by integrator
#include "wflow_includes.h"

void wflow() {
  register int i;
  register site *s;
  int j, istep, block_count = 0, blmax = 0;
  int max_scale, eps_scale = 1, N_base = 0;     // All three always positive
  double t = 0.0, E = 0.0, E_ss = 0.0, E_st = 0.0, tSqE = 0.0, der_tSqE;
  double ssplaq, stplaq, plaq = 0, old_plaq, baseline = 0.0;
  double check, old_check, topo, old_topo;
  double old_Ess, old_Est, old_tSqE, slope_Ess, slope_Est;
  double slope_tSqE, slope_check, slope_topo, prev_tSqE;
  double Delta_t, interp_t, interp_plaq;
  double interp_Ess, interp_Est, interp_E, interp_tSqE;
  double interp_check, interp_topo, interp_der;

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
  if (num_block > 0 && tblock[0] == 0.0) {
    mcrg_block(tblock[0], blmax);
    block_count++;
  }

  // Go
  epsilon = start_eps;
  max_scale = (int)floor(max_eps / start_eps);     // Maximum scaling factor
  for (istep = 0; fabs(t) <  fabs(tmax) - 0.5 * fabs(epsilon); istep++) {
    stout_step_rk();
    t += epsilon;

    // Find 8F_munu = sum_{clover} (U - Udag)
    // Subtract the (lattice artifact?) trace at each lattice site
    make_field_strength(F_OFFSET(link), F_OFFSET(FS));

    // Save previous data for slope and interpolation
    // old_plaq is only used for adjusting step size below
    old_Ess = E_ss;
    old_Est = E_st;
    old_tSqE = tSqE;
    old_check = check;
    old_topo = topo;
    old_plaq = plaq;

    // Compute t^2 E and its slope
    // Accumulate space--space and space--time separately for anisotropy
    // (Just look at usual t^2*E for those, but call them 'E' for short)
    E_ss = 0.0;
    E_st = 0.0;
    FORALLSITES(i, s) {
      E_ss -= (double)realtrace_nn(&(s->FS[0]), &(s->FS[0]));   // XY
      E_ss -= (double)realtrace_nn(&(s->FS[1]), &(s->FS[1]));   // XZ
      E_ss -= (double)realtrace_nn(&(s->FS[2]), &(s->FS[2]));   // YZ
      E_st -= (double)realtrace_nn(&(s->FS[3]), &(s->FS[3]));   // XT
      E_st -= (double)realtrace_nn(&(s->FS[4]), &(s->FS[4]));   // YT
      E_st -= (double)realtrace_nn(&(s->FS[5]), &(s->FS[5]));   // ZT
    }
    g_doublesum(&E_ss);
    g_doublesum(&E_st);

    // Normalization factor of 1/8 for each F_munu
    E = (E_ss + E_st) / (volume * 64.0);
    E_ss *= t * t / (volume * 64.0);
    E_st *= t * t / (volume * 64.0);
    tSqE = E_ss + E_st;
    der_tSqE = fabs(t) * (tSqE - old_tSqE) / fabs(epsilon);
    // Any negative signs in t and epsilon should cancel out anyway...

    // Might as well extract topology
    // Normalization is 1/(64*4pi^2) again with 1/8 for each F_munu
    topo = 0.0;
    FORALLSITES(i, s) {
      topo -= (double)realtrace_nn(&(s->FS[0]), &(s->FS[5])); // XYZT
      topo -= (double)realtrace_nn(&(s->FS[3]), &(s->FS[2])); // XTYZ
      topo -= (double)realtrace(&(s->FS[1]), &(s->FS[4]));    // XZ(YT)^dag
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
      slope_Ess = (E_ss - old_Ess) / epsilon;
      slope_Est = (E_st - old_Est) / epsilon;
      slope_tSqE = (tSqE - old_tSqE) / epsilon;
      slope_check = (check - old_check) / epsilon;
      slope_topo = (topo - old_topo) / epsilon;
      prev_tSqE = old_tSqE;   // For interpolating the derivative
      interp_t = t - epsilon;
      for (j = 0; j < eps_scale - 1; j++) {
        interp_t += start_eps;
        Delta_t = (j + 1) * start_eps;
        interp_Ess = old_Ess + Delta_t * slope_Ess;
        interp_Est = old_Est + Delta_t * slope_Est;
        interp_E = interp_Ess + interp_Est;
        interp_tSqE = old_tSqE + Delta_t * slope_tSqE;
        interp_der = interp_t * (interp_tSqE - prev_tSqE) / start_eps;
        prev_tSqE = interp_tSqE;

        interp_check = old_check + Delta_t * slope_check;
        interp_plaq = 3.0 - check / (12.0 * interp_t * interp_t);

        interp_topo = old_topo + Delta_t * slope_topo;
        node0_printf("WFLOW %g %g %g %g %g %g %g %g %g (interp)\n",
                     interp_t, interp_plaq, interp_E, interp_tSqE,
                     interp_der, interp_check, interp_topo,
                     interp_Ess, interp_Est);
      }
    }
    node0_printf("WFLOW %g %g %g %g %g %g %g %g %g\n",
                 t, plaq, E, tSqE, der_tSqE, check, topo, E_ss, E_st);

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
}
// -----------------------------------------------------------------

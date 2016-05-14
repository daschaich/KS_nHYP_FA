// -----------------------------------------------------------------
// Main procedure for SU(3) Wilson flow
#define CONTROL
#include "wflow_includes.h"
#define TOL 1e-6        // For floating-point comparisons
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
    for (dir = XUP; dir <= TUP; dir++)
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

  // Save the original links
  FORALLSITES(i, s) {
    for (dir = XUP; dir <= TUP; dir++)
      su3mat_copy(&(s->link[dir]), &(s->link[dir + 8]));
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
     for (dir = XUP; dir <= TUP; dir++)
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
  int j, prompt, istep, block_count = 0, blmax = 0, L, count;
  double dtime, t = 0.0, E, old_value, new_value = 0.0, der_value;
  double ssplaq, stplaq, plaq, check, topo, eta = 0.0;
  double oldPlaq, dPlaq = -1.0;
  double magicC = 0.05, magicT;
  complex tc;
  su3_matrix t_mat, *S[4];
  anti_hermitmat *A[4];

  // Setup
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc, &argv);
  g_sync();
  prompt = setup();
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Allocate fields used by integrator
  for (dir = 0; dir < 4; dir++) {
    S[dir] = malloc(sites_on_node * sizeof(su3_matrix));
    A[dir] = malloc(sites_on_node * sizeof(anti_hermitmat));
  }

  // Determine maximum number of blockings from smallest dimension
  // Works even if we can only block down to odd j > 4
  // Simultaneously compute appropriate 0.05L)^2/8 
  if (nx < nt) {
    j = nx;
    L = nx;
  }
  else {
    j = nt;
    L = nt;
  }
  while (j % 2 == 0 && j > 2) {    // While j is even
    j /= 2;
    blmax++;
  }

  // Special case: measure MCRG-blocked observables before any flow
  if (num_block > 0 && tblock[0] == 0) {
    mcrg_block(tblock[block_count], blmax);
    block_count++;
  }

  // Only adjust epsilon if start_eps < max_eps
  if (fabs(start_eps - max_eps) > TOL) {
    // Find the first magicT = (cL)^2 / 8 not within the first ten steps,
    // with L = min(nx, nt) set above
    magicT = magicC * magicC  * L * L / 8.0;
    while (magicT < 10.0 * start_eps) {
      magicC += 0.05;
      magicT = magicC * magicC * L * L / 8.0;
    }

    // Figure out how many steps until the next magicT
    // Then add one more and reduce epsilon to land bang on magicT
    count = ceil(magicT / start_eps) + 1;
    epsilon = magicT / (Real)count;
    node0_printf("EPS %g to reach t = %g for c = %g in %d steps\n",
                 epsilon, magicT, magicC, count);
  }
  else
    epsilon = max_eps;

  // Main flow loop
  plaquette(&ssplaq, &stplaq);
  oldPlaq = 0.5 * (ssplaq + stplaq);
  for (istep = 0;
       tmax == 0 || fabs(t) < fabs(tmax) - 0.5 * fabs(start_eps);
       istep++) {

    stout_step_rk(S, A);
    t += epsilon;

    // Finds 8F_munu = sum_{clover} (U - Udag)
    // Subtracts the (lattice artifact?) trace at each lattice site
    make_field_strength(F_OFFSET(link), F_OFFSET(FS));

    // The variables necessary for scale setting: used if tmax==0 only
    old_value = new_value;
    E = 0;
    FORALLSITES(i, s) {
      for (dir = 0; dir < 6; dir++) {
        mult_su3_nn(&(s->FS[dir]), &(s->FS[dir]), &t_mat);
        tc = trace_su3(&t_mat);
        E -= (double)tc.real;
      }
    }
    g_doublesum(&E);
    E /= (volume * 64.0); // Normalization factor of 1/8 for each F_munu
    new_value = t * t * E;
    der_value = fabs(t) * (new_value - old_value) / fabs(epsilon);
    // Any negative signs in t and epsilon should cancel out anyway...

    // Might as well extract topology
    // Normalization is 1/4pi^2 again with 1/8 for each F_munu
    topo = 0.0;
    FORALLSITES(i, s) {
      mult_su3_nn(&(s->FS[0]), &(s->FS[5]), &t_mat); // XYZT
      tc = trace_su3(&t_mat);
      topo -= (double)tc.real;

      mult_su3_nn(&(s->FS[3]), &(s->FS[2]), &t_mat); // XTYZ
      tc = trace_su3(&t_mat);
      topo -= (double)tc.real;

      mult_su3_na(&(s->FS[1]), &(s->FS[4]), &t_mat); // XZTY;  TY=(YT)^dag
      tc = trace_su3(&t_mat);
      topo -= (double)tc.real;
    }
    g_doublesum(&topo);
    topo *= 0.000395785873603;

    // Check with plaquette
    plaquette(&ssplaq, &stplaq);
    plaq = 0.5 * (ssplaq + stplaq);
    check = 12.0 * t * t * (3.0 - plaq);
    node0_printf("WFLOW %g %g %g %g %g %g %g\n",
                 t, plaq, E, new_value, der_value, check, topo);

    // For the tmax == 0 case, figure out if we can stop
    // Never stop before t=1
    // stepi%20 == 0 suppresses fluctuations in the flow length
    // between different configurations in the same ensemble
    if (tmax == 0 && fabs(t) > 1 && istep % 20 == 0) {
      if (new_value > 0.45 && der_value > 0.35)
        break;
    }

    // Do MCRG blocking at specified t
    // Since epsilon should only increase, this can go before the adjustment
    if (block_count < num_block
        && fabs(t + 0.5 * epsilon) >= fabs(tblock[block_count])) {
      mcrg_block(tblock[block_count], blmax);
      block_count++;
    }

    // Adjust epsilon when appropriate
    count--;
    dPlaq = fabs(plaq - oldPlaq);         // Force positivity

    // Average over the first ten epsilon * dPlaq to define eta
    if (istep < 10)
      eta += 0.1 * epsilon * dPlaq;

    // Only adjust epsilon if start_eps < max_eps
    if (count == 0 && fabs(epsilon - max_eps) > TOL) {
      // We should be at the magicT we wanted
      if (fabs(t - magicT) > TOL)
        node0_printf("WARNING: wanted t=%g but have %g\n", magicT, t);

      // Just end if we're already done
      if (fabs(tmax - t) < 0.5 * fabs(epsilon))
        break;

      epsilon = eta / dPlaq;
      if (epsilon > max_eps)
        epsilon = max_eps;

      // Figure out how many steps until the next magicT
      // Then add one more and reduce epsilon to land bang on magicT
      // Recall L = min(nx, nt)
      magicC += 0.05;
      magicT = magicC * magicC * L * L / 8.0;
      count = ceil((magicT - t) / start_eps) + 1;
      epsilon = (magicT - t) / (Real)count;
      node0_printf("EPS %g to reach t = %g for c = %g in %d steps\n",
                   epsilon, magicT, magicC, count);
    }
    oldPlaq = plaq;
  }

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  for (dir = 0; dir < 4; dir++) {
    free(S[dir]);
    free(A[dir]);
  }

  // Save lattice if requested -- no phases in links
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);
  normal_exit(0);
  return 0;
}
// -----------------------------------------------------------------

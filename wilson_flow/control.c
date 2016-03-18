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
  int j, prompt, istep, block_count = 0, blmax = 0;
  double dtime, t=0, E, old_value, new_value=0, der_value;
  double ssplaq, stplaq, td, check, topo;
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
  dtime =-dclock();

  // Allocate fields used by integrator
  for (dir = 0; dir < 4; dir++) {
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
  for (istep = 0; tmax == 0 || fabs(t) <  fabs(tmax) - fabs(epsilon) / 2; istep++) {
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
    E /= (volume * 64); // Normalization factor of 1/8 for each F_munu
    new_value = t * t * E;
    der_value = fabs(t) * (new_value - old_value) / fabs(epsilon);
    // Any negative signs in t and epsilon should cancel out anyway...

    // Might as well extract topology
    topo = 0;
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
    // Same normalization
    topo /= (volume * 64 * 0.02533029591058444286); // 1 / (volume / 4pi^2)
    // Correct normalization
//    topo *= 0.02533029591058444286 / 64; // 1/4pi^2

    // Check with plaquette
    d_plaquette(&ssplaq, &stplaq);
    td = (ssplaq + stplaq) / 2;
    check = 12 * t * t * (3 - td);
    node0_printf("WFLOW %g %g %g %g %g %g %g\n",
                 t, td, E, new_value, der_value, check, topo);

    // Does MCRG blocking at specified t
    if (block_count < num_block
        && fabs(t + epsilon / 2) >= fabs(tblock[block_count])) {
      mcrg_block(tblock[block_count], blmax);
      block_count++;
    }

    // For the tmax == 0 case, figure out if we can stop
    // Never stop before t=1
    // stepi%20 == 0 suppresses fluctuations in the flow length
    // between different configurations in the same ensemble
    if (tmax == 0 && fabs(t) > 1 && istep % 20 == 0) {
      if (new_value > 0.45 && der_value > 0.35)
        break;
    }
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

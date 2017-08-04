// -----------------------------------------------------------------
// Main procedure for SU(3) MCRG-blocked measuremenst
#define CONTROL
#include "block_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[])  {
  register int i, dir, ismear, ivec;
  register site *s;
  int j, k, prompt, istep, total_iters, avm_iters = 0;
  Real ssplaq, stplaq;
  double rr[10], t = 0, old_value, new_value = 0, der_value, E, chirality;
  double td, check, topo, dtime, tot_time;
  complex plp, xplp, tc;
  matrix **bak, t_mat, *S[4];
  anti_hermitmat *A[4];

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run (loop removed)
  if (readin(prompt) == 0) {
    dtime = -dclock();
    tot_time = -dclock();
    make_loop_table2();                 // Set up loop tables

    // Check initial Wilson and Polyakov loops
    blocked_gauge_loops(0, rr);
    for (j = 0; j < nloop; j++) {
      for (k = 0; k < nreps; k++)
       node0_printf("LOOPS %d %d 0 0 %.8g\n",
                    j, k, rr[j + k * nloop] / (volume * loop_num[j]));
    }
    plp = blocked_ploop(0, TUP);
    xplp = blocked_ploop(0, XUP);
    node0_printf("POLYA %.6g %.6g %.6g %.6g\n",
                 (double)plp.real, (double)plp.imag,
                 (double)xplp.real, (double)xplp.imag);

    // Allocate fields used by integrator -- must not be used after blocking
    FIELD_ALLOC(tempmat1, matrix);
    FIELD_ALLOC(tempmat2, matrix);
    for (dir = 0; dir < 4; dir++) {
      S[dir] = (matrix *)malloc(sites_on_node * sizeof(matrix));
      A[dir] = (anti_hermitmat *)malloc(sites_on_node * sizeof(anti_hermitmat));
    }

    // Wilson flow!
    if (tmax != 0) {  // Braces suppress annoying compiler warning
      node0_printf("Wilson flow with dt = %g and tmax = %g\n", epsilon, tmax);
    }
    else
      node0_printf("No Wilson flow\n");

    for (istep = 0; fabs(t) <  fabs(tmax) - fabs(epsilon) / 2; istep++) {
      stout_step_rk(S, A);
      t += epsilon;

      // Finds 8F_munu = sum_{clover} (U - Udag)
      // Subtracts the (lattice artifact?) trace at each lattice site
      make_field_strength(F_OFFSET(link), F_OFFSET(fieldstrength));

      // The variables necessary for scale setting: used if tmax==0 only
      old_value = new_value;
      E = 0;
      FORALLSITES(i, s) {
        for (dir = 0; dir < 6; dir++) {
          mult_su3_nn(&(s->fieldstrength[dir]),
                      &(s->fieldstrength[dir]), &t_mat);
          tc = trace_su3(&t_mat);
          E -= (double)tc.real;
        }
      }
      g_doublesum(&E);
      E /= (volume * 64); // Normalization factor of 1/8 for each F_munu
      new_value = t * t * E;
      der_value = t * (new_value - old_value) / epsilon;

      // Might as well extract topology
      topo = 0;
      FORALLSITES(i, s) {
        mult_su3_nn(&(s->fieldstrength[0]), // XYZT
                    &(s->fieldstrength[5]), &t_mat);
        tc = trace_su3(&t_mat);
        topo -= (double)tc.real;

        mult_su3_nn(&(s->fieldstrength[3]), // XTYZ
                    &(s->fieldstrength[2]), &t_mat);
        tc = trace_su3(&t_mat);
        topo -= (double)tc.real;

        mult_su3_na(&(s->fieldstrength[1]), // XZTY
                    &(s->fieldstrength[4]), &t_mat);
        tc = trace_su3(&t_mat);
        topo -= (double)tc.real;
      }
      g_doublesum(&topo);
      // Same normalization
      topo /= (volume * 64 * 0.02533029591058444286); // 1/4pi^2

      // Check with plaquette
      d_plaquette(&ssplaq, &stplaq);
      td = (ssplaq + stplaq) / 2;
      check = 12 * t * t * (3 - td);
      node0_printf("WFLOW %g %g %g %g %g %g %g\n",
                   t, td, E, new_value, der_value, check, topo);
    }
    // Free fields used by integrator
    free(tempmat1);
    free(tempmat2);
    for (dir = 0; dir < 4; dir++) {
      free(S[dir]);
      free(A[dir]);
    }

    // Check Wilson and Polyakov loops after Wilson flow
    blocked_gauge_loops(0, rr);
    for (j = 0; j < nloop; j++) {
      for (k = 0; k < nreps; k++)
       node0_printf("WFLOWED LOOPS %d %d 0 0 %.8g\n",
                    j, k, rr[j + k * nloop] / (volume * loop_num[j]));
    }
    plp = blocked_ploop(0, TUP);
    xplp = blocked_ploop(0, XUP);
    node0_printf("WFLOWED POLYA %.6g %.6g %.6g %.6g\n",
                 (double)plp.real, (double)plp.imag,
                 (double)xplp.real, (double)xplp.imag);
    node0_printf("WILSON FLOW COMPLETED\n");
    dtime += dclock();
    node0_printf("Wilson flow time %.4g seconds\n", dtime);
    fflush(stdout);
    dtime = -dclock();

    // Smearing parameters fixed, read from infile
    // Only block once
    // nHYP smear twice
    block_nhyp_mcrg(0, 1);
    block_nhyp_mcrg(0, 1);
    // Finally block, only doing center staple
    block_mcrg(0, 1);

    // Print out blocked loops to check against resized lattice
    blocked_gauge_loops(1, rr);
    for (j = 0; j < nloop; j++) {
      for (k = 0; k < nreps; k++)
        node0_printf("LOOPS %d %d %d %.4g %.8g\n",
                     j, k, 1, alpha_smear[0],
                     rr[j + k * nloop] / (volume * loop_num[j]));
    }
    plp = blocked_ploop(1, TUP);
    xplp = blocked_ploop(1, XUP);
    node0_printf("POLYA NHYP %d %.6g %.6g %.6g %.6g %.6g\n",
                 1, alpha_smear[0],
                 (double)plp.real, (double)plp.imag,
                 (double)xplp.real, (double)xplp.imag);
    fflush(stdout);

    // Allocate backup space and save one set of blocked links
    bak = (matrix **)malloc((volume / 16) * sizeof(matrix *));
    for (j = 0; j < (volume / 16); j++)
      bak[j] = (matrix *)malloc(4 * sizeof(matrix));

    node0_printf("(dx, dy, dz, dt) = (%i, %i, %i, %i)\n", dx, dy, dz, dt);
    FORALLSITES(i, s) {
      if (s->x % 2 == dx && s->y % 2 == dy && s->z % 2 == dz
                                           && s->t % 2 == dt) {
        j = (s->x - dx) * nt * nz * ny / 16 + (s->y - dy) * nt * nz / 8
                     + (s->z - dz) * nt / 4 + (s->t - dt) / 2;
        for (dir = XUP; dir <= TUP; dir++)
          mat_copy(&(s->link[dir]), &(bak[j][dir]));
      }
    }

    // Reset with all lattice dimensions halved
    free_lattice();
    nx /= 2;
    ny /= 2;
    nz /= 2;
    nt /= 2;
    volume /= 16;
    setup_layout();
    make_lattice();
    reset_machine();    // Clear old gathers
    make_nn_gathers();
    make_fields();      // Allocate space for fields on reduced lattice
    phaseset();         // Set up staggered phase vectors, boundary conditions

    // Create dummy lattice of appropriate size
    startflag = FRESH;    // Makes IO ignore startfile
    sprintf(startfile, "NULL"); // Check that this doesn't matter
    startlat_p = reload_lattice(startflag, startfile);

    // Copy saved set of blocked links into new lattice
    FORALLSITES(i, s) {
      j = s->x * nt * nz * ny + s->y * nt * nz
                              + s->z * nt + s->t;
      for (dir = XUP; dir <= TUP; dir++)
        mat_copy(&(bak[j][dir]), &(s->link[dir]));
    }
    for (j = 0; j < volume; j++)
      free(bak[j]);
    free(bak);

    blocked_gauge_loops(0, rr);
    for (j = 0; j < nloop; j++) {
      for (k = 0; k < nreps; k++)
       node0_printf("CHECK LOOPS %d %d 0 0 %.8g\n",
                    j, k, rr[j + k * nloop] / (volume * loop_num[j]));
    }
    plp = blocked_ploop(0, TUP);
    xplp = blocked_ploop(0, XUP);
    node0_printf("POLYA CHECK %.6g %.6g %.6g %.6g\n",
                 (double)plp.real, (double)plp.imag,
                 (double)xplp.real, (double)xplp.imag);
    node0_printf("BLOCKING COMPLETED, sites_on_node = %d\n", sites_on_node);
    dtime += dclock();
    node0_printf("Block time %.4g seconds\n", dtime);
    fflush(stdout);
    dtime = -dclock();

    // Add staggered phases for eigenvalue calculation
    phases_in = OFF;
    rephase(ON);

    // Allow different HYP parameters
    alpha_smear[0] = alpha_store[0];
    alpha_smear[1] = alpha_store[1];
    alpha_smear[2] = alpha_store[2];

    // Save lattice before smearing if requested (no phases)
    if (saveflag != FORGET) {
      rephase(OFF);
      save_lattice(saveflag, savefile, stringLFN);
      rephase(ON);
    }

    // Smear however many times requested
    // Print out plaquette each time (unsmeared plaquette printed by setup)
    for (ismear = 0; ismear < nsmear; ismear++) {
      block_and_fatten();   // Smears gauge_field_thin into gauge_field
      for (dir = XUP; dir <= TUP; dir++) {
        FORALLSITES(i, s)
          gauge_field_thin[dir][i] = gauge_field[dir][i];
      }
      rephase(OFF);
      d_plaquette(&ssplaq, &stplaq);
      node0_printf("Plaquettes after smearing %d: %.8g %.8g\n",
                   ismear + 1, ssplaq, stplaq);
      rephase(ON);
      avm_iters += meas_link(F_OFFSET(chi), F_OFFSET(psi), mass);
      fflush(stdout);
    }

    // Allocate eigenvectors
    eigVal = (double *)malloc(Nvecs * sizeof(double));
    eigVec = (vector **)malloc(Nvecs * sizeof(vector *));
    for (i = 0; i < Nvecs; i++)
      eigVec[i] = (vector *)malloc(sites_on_node * sizeof(vector));

    // Hard-code EVEN parity in Kalkreuter
    total_iters = Kalkreuter(eigVec, eigVal, eig_tol, error_decr,
                             Nvecs, maxIter, restart, kiters);

    for (ivec = 0; ivec < Nvecs; ivec++) {
      // Construct the odd part of the vector
      //   i Dslash Psi / sqrt(eigVal)
      // Ignore the i factor, which is irrelevant for the chirality
      // Halve EVENANDODD result since the EVEN vector is normalized to 1
      FORALLSITES(i, s)
        s->chi = eigVec[ivec][i];

      dslash(F_OFFSET(chi), F_OFFSET(psi), ODD);
      FORODDSITES(i, s)
        scalar_mult_vector(&(s->psi), 1.0 / sqrt(eigVal[ivec]),
                               &(eigVec[ivec][i]));

      // Hard-code EVENANDODD parity in measure_chirality
      measure_chirality(eigVec[ivec], &chirality);
      node0_printf("CHIRALITY %d %.4g\n", ivec, chirality / 2);
    }
    dtime += dclock();    // Done with dtime
    node0_printf("Eigenvalue time %.4g seconds\n", dtime);

    // Perhaps in the future
    // Hard-code EVEN parity in print_densities
//    for (i = 0; i < Nvecs; i++) {
//      sprintf(label, "DENSITY(%i)", i);
//      print_densities(eigVec[i], label, ny / 2, nz / 2, nt / 2);
//    }

    for (i = 0; i < Nvecs; i++)
      free(eigVec[i]);
    free(eigVec);
    free(eigVal);
    free_fields();
    free_lattice();

    node0_printf("RUNNING COMPLETED\n");
    tot_time += dclock();
    node0_printf("Time = %.4g seconds\n", tot_time);
    node0_printf("total_iters = %d\n", total_iters);
    fflush(stdout);
  }
  return 0;
}
// -----------------------------------------------------------------

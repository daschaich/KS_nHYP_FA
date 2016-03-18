// -----------------------------------------------------------------
// SU(3) MCRG-blocked measurements setup
#include "block_includes.h"
#include <string.h>

int initial_set();
#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed, and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU(3) Kogut--Susskind blocked eigenvalues\n");
    printf("MIMD version 7ish with Kostas's Kalkreuter\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    printf("Twice-HYP smeared links, smearing parameters from infile\n");
    printf("  IR_STAB = %.4g\n", (Real)IR_STAB);
    printf("  EPS_SQ = %.4g\n", (Real)EPS_SQ);
#ifdef NHYP_DEBUG
    printf("NHYP_DEBUG turned on\n");
#endif
#ifdef NO_UNIT_CHECK
    printf("NOT checking unitarity when loading lattice\n");
#endif
    time_stamp("start");

    status = get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    if(status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
// Amount reported is pure guesswork/imagination because I don't really care
void make_fields() {
  FIELD_ALLOC_VEC(gauge_field, su3_matrix, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, su3_matrix, 4);

  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, su3_matrix, 4);
  FIELD_ALLOC_VEC(Staple3, su3_matrix, 4);
  FIELD_ALLOC(tempmat1, su3_matrix);
  FIELD_ALLOC(tempmat2, su3_matrix);

  // Check the total number of su3_matrices; this may not be accurate
  node0_printf("Mallocing %.1f MBytes per node for fields\n",
               (double)sites_on_node * 98 * sizeof(su3_matrix) / 1e6);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Free all fields allocated by make_fields
void free_fields() {
  register int dir, dir2;
  for (dir = 0; dir < 4; dir++) {
    free(gauge_field[dir]);
    free(gauge_field_thin[dir]);
    free(Staple3[dir]);

    for (dir2 = 0; dir2 < 4; dir2++) {
      if(dir != dir2) {
        free(hyplink1[dir][dir2]);
        free(hyplink2[dir][dir2]);
        free(Staple1[dir][dir2]);
        free(Staple2[dir][dir2]);
      }
    }
  }
  free(tempmat1);
  free(tempmat2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(3) MCRG-blocked measurements
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_mcrg0", &par_buf.alpha_mcrg0);
    IF_OK status += get_f(stdin, prompt, "alpha_mcrg1", &par_buf.alpha_mcrg1);
    IF_OK status += get_f(stdin, prompt, "alpha_mcrg2", &par_buf.alpha_mcrg2);

    // Where in the hypercube the blocked lattice originates
    IF_OK status += get_i(stdin, prompt, "dx", &par_buf.dx);
    IF_OK status += get_i(stdin, prompt, "dy", &par_buf.dy);
    IF_OK status += get_i(stdin, prompt, "dz", &par_buf.dz);
    IF_OK status += get_i(stdin, prompt, "dt", &par_buf.dt);
    if (par_buf.dx < 0 || par_buf.dx > 1
     || par_buf.dy < 0 || par_buf.dy > 1
     || par_buf.dz < 0 || par_buf.dz > 1
     || par_buf.dt < 0 || par_buf.dt > 1) {
      node0_printf("ERROR: Hypercube origin must be within 2^4 lattice\n");
      node0_printf("ERROR: (dx, dy, dz, dt) = (%i, %i, %i, %i)\n",
                   par_buf.dx, par_buf.dy, par_buf.dz, par_buf.dt);
      status++;
    }

    // Always calculating massless eigenvalues
    IF_OK status += get_i(stdin, prompt, "nsmear", &par_buf.nsmear);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    IF_OK status += get_f(stdin, prompt, "tmax", &par_buf.tmax);
    IF_OK status += get_f(stdin, prompt, "epsilon", &par_buf.epsilon);

    IF_OK status += get_i(stdin, prompt, "Nvecs", &par_buf.Nvecs);
    IF_OK status += get_f(stdin, prompt, "eig_tol", &par_buf.eig_tol);
    IF_OK status += get_i(stdin, prompt, "maxIter", &par_buf.maxIter);
    IF_OK status += get_i(stdin, prompt, "restart", &par_buf.restart);
    IF_OK status += get_i(stdin, prompt, "kiters", &par_buf.kiters);
    IF_OK status += get_f(stdin, prompt, "error_decr", &par_buf.error_decr);

    // Mass for pbp
    IF_OK status += get_f(stdin, prompt, "mass", &par_buf.mass);

    // Number of stochastic sources for pbp
    IF_OK status += get_i(stdin, prompt, "npbp", &par_buf.npbp);

    // Maximum conjugate gradient iterations
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);

    // Maximum conjugate gradient restarts
    IF_OK status += get_i(stdin, prompt, "max_cg_restarts", &par_buf.nrestart);

    // Error per site for conjugate gradient
    IF_OK {
      status += get_f(stdin, prompt, "error_per_site", &x);
      par_buf.rsqmin = x * x;
    }

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin, prompt, &(par_buf.saveflag),
                                       par_buf.savefile);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  alpha_smear[0] = par_buf.alpha_mcrg0;
  alpha_smear[1] = par_buf.alpha_mcrg1;
  alpha_smear[2] = par_buf.alpha_mcrg2;

  dx = par_buf.dx;
  dy = par_buf.dy;
  dz = par_buf.dz;
  dt = par_buf.dt;

  nsmear = par_buf.nsmear;
  alpha_store[0] = par_buf.alpha_hyp0;
  alpha_store[1] = par_buf.alpha_hyp1;
  alpha_store[2] = par_buf.alpha_hyp2;

  tmax = par_buf.tmax;
  epsilon = par_buf.epsilon;

  Nvecs = par_buf.Nvecs;
  eig_tol = par_buf.eig_tol;
  maxIter = par_buf.maxIter;
  restart = par_buf.restart;
  kiters = par_buf.kiters;
  error_decr = par_buf.error_decr;

  mass = par_buf.mass;
  npbp = par_buf.npbp;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.savefile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, volume + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  return prompt;
}
// -----------------------------------------------------------------

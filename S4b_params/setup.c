// -----------------------------------------------------------------
// SU(3) S4b order parameters setup
#include <string.h>
#include "S4b_includes.h"

#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if (mynode() == 0) {
    // Print banner
    printf("SU(3) with Kogut--Susskind fermions\n");
    printf("S4b phase order parameter measurement\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    printf("nHYP links, reading alpha_smear parameters from infile\n");
    printf("  IR_STAB = %.4g\n", (Real)IR_STAB);
    printf("  EPS_SQ = %.4g\n", (Real)EPS_SQ);
#ifdef NHYP_DEBUG
    printf("NHYP_DEBUG turned on\n");
#endif
#ifdef CG_DEBUG
    printf("CG_DEBUG turned on\n");
#endif
#ifdef CGTIME
    printf("CGTIME turned on\n");
#endif
    time_stamp("start");

    status = get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt, "nflavors", &par_buf.nflavors);
    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    if (status > 0)
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
  iseed = par_buf.iseed;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  one_ov_vol = 1.0 / (Real)volume;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
void make_fields() {
  FIELD_ALLOC_VEC(gauge_field, matrix, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, matrix, 4);

  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, matrix, 4);
  FIELD_ALLOC_VEC(Staple3, matrix, 4);
  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);

  node0_printf("Mallocing %.1f MBytes per core for fields\n",
               (double)sites_on_node * 62 * sizeof(matrix) / 1e6);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume, seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, volume + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Allocate space for fields
  make_fields();
  // Set up staggered phase vectors, boundary conditions
  phaseset();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(3) S4b order parameters
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // beta, adjoint ratio and mass
    IF_OK status += get_f(stdin, prompt, "beta", &par_buf.beta);
    IF_OK status += get_f(stdin, prompt, "beta_a", &par_buf.beta_a);
    IF_OK status += get_f(stdin, prompt, "mass", &par_buf.mass);

    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // Number of stochastic sources
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

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  beta = par_buf.beta;
  beta_a = par_buf.beta_a;
  mass = par_buf.mass;
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  npbp = par_buf.npbp;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;

  startflag = par_buf.startflag;
  strcpy(startfile, par_buf.startfile);

  // Do whatever is needed to get lattice
  if (startflag == CONTINUE)
    rephase(OFF);

  startlat_p = reload_lattice(startflag, startfile);
  // If a lattice was read in, put in staggered phase factors
  // and antiperiodic boundary condition
  phases_in = OFF;
  rephase(ON);
  return 0;
}
// -----------------------------------------------------------------

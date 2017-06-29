// -----------------------------------------------------------------
// Staggered eigenvalue setup
#include "eig_includes.h"
#include <string.h>   // For strcpy

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;

#define IF_OK if (status == 0)
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read and distribute lattice size and random number seed
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU3 Kogut--Susskind eigenvalue calculation\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    printf("nHYP links, reading alpha_smear parameters from infile\n");
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
  total_iters = 0;
  return(prompt);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
void make_fields() {
  FIELD_ALLOC_VEC(gauge_field, su3_matrix, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, su3_matrix, 4);

  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, su3_matrix, 4);
  FIELD_ALLOC_VEC(Staple3, su3_matrix, 4);

  FIELD_ALLOC(tempmat1, su3_matrix);

  // Check the total number of su3_matrices; this may not be accurate
  node0_printf("Mallocing %.1f MBytes per core for fields\n",
               (double)sites_on_node * 98 * sizeof(su3_matrix) / 1e6);
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
// Read in parameters for SU(3) eigenvalues
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Always calculating massless eigenvalues
    IF_OK status += get_i(stdin, prompt, "nsmear", &par_buf.nsmear);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    IF_OK status += get_i(stdin, prompt, "Nvecs", &par_buf.Nvecs);
    IF_OK status += get_f(stdin, prompt, "eig_tol", &par_buf.eig_tol);
    IF_OK status += get_i(stdin, prompt, "maxIter", &par_buf.maxIter);
    IF_OK status += get_i(stdin, prompt, "restart", &par_buf.restart);
    IF_OK status += get_i(stdin, prompt, "kiters", &par_buf.kiters);
    IF_OK status += get_f(stdin, prompt, "error_decr", &par_buf.error_decr);

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

  nsmear = par_buf.nsmear;
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  Nvecs = par_buf.Nvecs;
  eig_tol = par_buf.eig_tol;
  maxIter = par_buf.maxIter;
  restart = par_buf.restart;
  kiters = par_buf.kiters;
  error_decr = par_buf.error_decr;

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

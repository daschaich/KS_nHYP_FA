// -----------------------------------------------------------------
// SU(3) stochastic mode number calculation setup
#include "mode_includes.h"
#include <string.h>

#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU(3) with Kogut--Susskind fermions\n");
    printf("Stochastic mode number measurement\n");
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
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
void make_fields() {
  int ipbp;
  Real size = (Real)(8.0 * sizeof(matrix));
  FIELD_ALLOC_VEC(gauge_field, matrix, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, matrix, 4);

  size += (Real)((4.0 + 4.0 * 12.0) * sizeof(matrix));
  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, matrix, 4);
  FIELD_ALLOC_VEC(Staple3, matrix, 4);

  size += (Real)(sizeof(matrix));
  FIELD_ALLOC(tempmat, matrix);

  size *= sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume, seed
  prompt = initial_set();
  // Initialize the node source number generator
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
  return(prompt);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(3) Monte Carlo
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status, ipbp;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Smearing parameters
    IF_OK status += get_i(stdin, prompt, "Nsmear", &par_buf.Nsmear);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // Number of stochastic sources
    IF_OK status += get_i(stdin, prompt, "npbp", &par_buf.npbp);

    // Which order polynomial to use
    IF_OK status += get_i(stdin, prompt, "order", &par_buf.order);

    // Number of Omegas and the interval
    IF_OK status += get_i(stdin, prompt, "npts", &par_buf.npts);
    IF_OK status += get_f(stdin, prompt, "startomega", &par_buf.startomega);
    IF_OK status += get_f(stdin, prompt, "spacing", &par_buf.spacing);

    // Maximum conjugate gradient iterations and restarts
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);
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

  Nsmear = par_buf.Nsmear;
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  // Include some mallocs here (which is called after make_fields)
  npbp = par_buf.npbp;
  source = malloc(npbp * sizeof(vector *));   // Stochastic sources
  for (ipbp = 0; ipbp < npbp; ipbp++)
    FIELD_ALLOC(source[ipbp], vector);

  Norder = par_buf.order;
  coeffs = malloc((Norder + 1) * sizeof(double)); // Polynomial coefficients

  Npts = par_buf.npts;
  spacing = par_buf.spacing;
  M = 0.5 * par_buf.startomega;    // !!! Note factor of 2
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

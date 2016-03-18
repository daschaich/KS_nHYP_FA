// -----------------------------------------------------------------
// SU(3) vacuum polarization setup
#include "vacpol_includes.h"

int initial_set();
void make_fields();
#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
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
// On node zero, read lattice size, seed and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU(3) with Kogut--Susskind fermions\n");
    printf("Vacuum polarization\n");
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

  // For QDP
  ndim = 4;
  latsize[0] = nx;
  latsize[1] = ny;
  latsize[2] = nz;
  latsize[3] = nt;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for vacuum polarization calculation
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Fermion mass and smearing parameters
    IF_OK status += get_f(stdin, prompt, "mass", &par_buf.mass);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // Point source location and maximum QSq
    IF_OK status += get_i(stdin, prompt, "src_x", &par_buf.x_src);
    IF_OK status += get_i(stdin, prompt, "src_y", &par_buf.y_src);
    IF_OK status += get_i(stdin, prompt, "src_z", &par_buf.z_src);
    IF_OK status += get_i(stdin, prompt, "src_t", &par_buf.t_src);

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

    // Get output file tag
    IF_OK status += get_s(stdin, prompt, "outpat", par_buf.outpat);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  mass = par_buf.mass;
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  x_src = par_buf.x_src;
  y_src = par_buf.y_src;
  z_src = par_buf.z_src;
  t_src = par_buf.t_src;

  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;

  startflag = par_buf.startflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(outpat, par_buf.outpat);

  // Do whatever is needed to get lattice
  if (startflag == CONTINUE)
    rephase(OFF);
  else
    startlat_p = reload_lattice(startflag, startfile);

  // If a lattice was read in, put in staggered phase factors
  // and antiperiodic boundary condition
  phases_in = OFF;
  rephase(ON);
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
// Amount Malloced is pure guesswork/imagination
void make_fields() {
  FIELD_ALLOC_VEC(gauge_field, su3_matrix, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, su3_matrix, 4);

  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, su3_matrix, 4);
  FIELD_ALLOC_VEC(Staple3, su3_matrix, 4);
  FIELD_ALLOC(tempmat1, su3_matrix);

  node0_printf("Mallocing %.1f MBytes per node for fields\n",
               (double)sites_on_node * 61 * sizeof(su3_matrix) / 1e6);
}
// -----------------------------------------------------------------

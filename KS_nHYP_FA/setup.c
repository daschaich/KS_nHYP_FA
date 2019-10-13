// -----------------------------------------------------------------
// Dynamical nHYP staggered setup
#include "ks_dyn_includes.h"
#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed, nflavors and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU(3) with Kogut--Susskind fermions\n");
    printf("Microcanonical simulation with refreshing\n");
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
#ifdef CG_DEBUG
    printf("CG_DEBUG turned on\n");
#endif
#ifdef CGTIME
    printf("CGTIME turned on\n");
#endif
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("Phi algorithm\n");
#else  // Quit!
    printf("Only works for phi algorithm!\n");
    exit(1);
#endif
    time_stamp("start");
    status = get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt, "nflavors", &par_buf.nflavors);
#ifdef PHI_ALGORITHM
    IF_OK if ((par_buf.nflavors % 4) != 0) {
      printf("Error: Use phi algorithm only for multiples of four flavors\n");
      status++;
    }
#endif
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

  nflavors = par_buf.nflavors;
  if (nflavors % 8 == 0) {
    half_fields = 0;
    full_fields = nflavors / 8;
  }
  else {
    half_fields = 1;
    full_fields = (nflavors - 4) / 8;
  }
  this_node = mynode();
  node0_printf("full_fields=%d, half_fields=%d\n", full_fields, half_fields);
  // Check that we can handle this nflavors
  if (full_fields > MAX_FIELDS) {
    node0_printf("Error: Can only handle nflavors up to %d\n",
                 8 * MAX_FIELDS);
    node0_printf("       Recompile with different MAX_FIELDS for more\n");
    exit(1);
  }

  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;
  iseed = par_buf.iseed;

  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  one_ov_vol = 1.0 / (Real)volume;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
// Changed to mimic allocation order in v6/KS_nHYP_FA/control.c
void make_fields() {
  FIELD_ALLOC_VEC(gauge_field, matrix, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, matrix, 4);

  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, matrix, 4);
  FIELD_ALLOC_MAT(SigmaH2, matrix, 4, 4);

  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);
  FIELD_ALLOC_VEC(Staple3, matrix, 4);
  FIELD_ALLOC_VEC(LambdaU, matrix, 4);
  FIELD_ALLOC_VEC(Lambda1, matrix, 4);
  FIELD_ALLOC_VEC(Lambda2, matrix, 4);
  FIELD_ALLOC_VEC(SigmaH, matrix, 4);
  FIELD_ALLOC_VEC(Sigma, matrix, 4);

  node0_printf("Mallocing %.1f MBytes per core for fields\n",
               (double)sites_on_node * 98 * sizeof(matrix) / 1e6);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume, nflavors, seed
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
// Read in parameters for SU(3) Monte Carlo
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status, i;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Warms, trajecs
    IF_OK status += get_i(stdin, prompt, "warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt, "trajecs", &par_buf.trajecs);
    IF_OK status += get_f(stdin, prompt, "traj_length", &par_buf.traj_length);

    // Number of pseudo_fermions, note restrictions
    IF_OK status += get_i(stdin, prompt, "number_of_PF", &par_buf.num_masses);
    if (par_buf.num_masses > MAX_MASSES || par_buf.num_masses < 1) {
      printf("Error: num_masses = %d must be <= %d and > 0!\n",
             par_buf.num_masses, MAX_MASSES);
      status++;
    }

    // Number of steps for all possible Hasenbusch masses, and gauge
    // This combined with trajectory length defines step size
    for (i = 0; i < MAX_MASSES; i++)
      IF_OK status += get_i(stdin, prompt, "nstep", &par_buf.nsteps[i]);

    // Gauge step is last in nsteps
    IF_OK status += get_i(stdin, prompt, "nstep_gauge",
                          &par_buf.nsteps[MAX_MASSES]);

    // Trajectories between propagator measurements
    IF_OK status += get_i(stdin, prompt, "traj_between_meas",
                          &par_buf.propinterval);

    // beta, adjoint ratio, mass and Hasenbusch mass
    IF_OK status += get_f(stdin, prompt, "beta", &par_buf.beta);
    if (beta < 0) {
      printf("Error: negative beta %.4g not allowed\n", beta);
      status++;
    }
    IF_OK status += get_f(stdin, prompt, "beta_a", &par_buf.beta_a);
    IF_OK status += get_f(stdin, prompt, "mass", &par_buf.mass);
    IF_OK status += get_f(stdin, prompt, "Hasenbusch_mass", &par_buf.MH);

    // Smearing parameters
    IF_OK status += get_i(stdin, prompt, "Nsmear", &par_buf.Nsmear);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // Number of stochastic sources
    IF_OK status += get_i(stdin, prompt, "npbp", &par_buf.npbp);

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

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  traj_length = par_buf.traj_length;
  num_masses = par_buf.num_masses;
  for (i = 0; i < MAX_MASSES + 1; i++)
    nsteps[i] = par_buf.nsteps[i];
  if (num_masses > 1)
    MH = par_buf.MH;

  propinterval = par_buf.propinterval;
  npbp = par_buf.npbp;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;

  beta = par_buf.beta;
  beta_a = par_buf.beta_a;
  mass = par_buf.mass;
  Nsmear = par_buf.Nsmear;
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.savefile);

  // Do whatever is needed to get lattice
  if (startflag == CONTINUE)
    rephase(OFF);
  else
    startlat_p = reload_lattice(startflag, startfile);
  // If a lattice was read in, put in staggered phase factors
  // and antiperiodic boundary condition
  phases_in = OFF;
  rephase(ON);

  // Only make gauge_field_save and gauge_field_temp if needed
  if (Nsmear > 1) {
    FIELD_ALLOC_VEC(gauge_field_save, matrix, 4);
    FIELD_ALLOC_VEC(gauge_field_temp, matrix, 4);
  }
  else if (Nsmear == 0) {
    node0_printf("Nsmear=0 breaks things; just set alpha=0\n");
    exit(1);
  }
  return 0;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// SU(3) MCRG-blocked measurements setup
#include "mcrg_includes.h"
#include <string.h>

int initial_set();
#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume
  prompt = initial_set();
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed, and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU(3) MCRG-blocked measurements\n");
    printf("MIMD version 7ish\n");
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

    if(status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0 )
    normal_exit(0);

  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(3) MCRG-blocked measurements
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status, i;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // A maximum of 100 outer nHYP smearing parameters should be plenty
    IF_OK status += get_i(stdin, prompt, "num_alpha", &par_buf.num_alpha);
    if (par_buf.num_alpha > 100) {
      node0_printf("ERROR: Need to recompile for num_alpha > 100\n");
      status++;
    }
    for (i = 0; i < par_buf.num_alpha; i++)
      IF_OK status += get_f(stdin, prompt, "alpha_mcrg", &par_buf.alpha_mcrg[i]);

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

  num_alpha = par_buf.num_alpha;
  for (i = 0; i < num_alpha; i++)
    alpha_mcrg[i] = par_buf.alpha_mcrg[i];

  // alpha_smear[0] should be overwritten by same value
  alpha_smear[0] = par_buf.alpha_mcrg[0];
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  startflag = par_buf.startflag;
  strcpy(startfile,par_buf.startfile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice( startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------

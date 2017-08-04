// -----------------------------------------------------------------
// SU(3) Wilson flow setup
#include <string.h>
#include "wflow_includes.h"

#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if (mynode() == 0) {
    // Print banner
    printf("Wilson flow with optional MCRG blocking\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    time_stamp("start");

    status = get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);

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

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  return prompt;
}
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
  // Allocate temporary fields and Wilson flow stuff
  FIELD_ALLOC(tempmat, su3_matrix);
  FIELD_ALLOC(tempmat2, su3_matrix);
  FIELD_ALLOC_VEC(S, su3_matrix, NDIMS);
  FIELD_ALLOC_VEC(A, anti_hermitmat, NDIMS);

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(3) Wilson flow measurements
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status, i;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Wilson flow parameters
    IF_OK status += get_f(stdin, prompt, "start_eps", &par_buf.start_eps);
    IF_OK status += get_f(stdin, prompt, "max_eps", &par_buf.max_eps);
    if (par_buf.start_eps == 0) {
      node0_printf("ERROR: epsilon=%g won't get you anywhere\n",
                   par_buf.start_eps);
      status++;
    }
    if (fabs(par_buf.start_eps) > fabs(par_buf.max_eps)) {
      node0_printf("WARNING: Resetting start_eps to max_eps\n");
      par_buf.start_eps = par_buf.max_eps;
    }
    if (par_buf.start_eps * par_buf.max_eps < 0)
      node0_printf("WARNING: start_eps and max_eps have different signs\n");

    IF_OK status += get_f(stdin, prompt, "tmax", &par_buf.tmax);
    if (par_buf.start_eps * par_buf.tmax < 0)
      node0_printf("WARNING: start_eps and tmax have different signs\n");

    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);
    // A maximum of 100 tvalues to perform blocking should be enough
    IF_OK status += get_i(stdin, prompt, "num_block", &par_buf.num_block);
    if (par_buf.num_block > 100) {
      node0_printf("ERROR: Need to recompile for num_block > 100\n");
      status++;
    }
    for (i = 0; i < par_buf.num_block; i++) {
      IF_OK status += get_f(stdin, prompt, "tblock", &par_buf.tblock[i]);
      // Make sure we're going in the right direction
      if (i > 0 && fabs(par_buf.tblock[i]) <= fabs(par_buf.tblock[i - 1])) {
        node0_printf("ERROR: We require tblock be sorted\n");
        node0_printf("ERROR: tblock[%d]=%g; tblock[%d]=%g\n",
                     i, par_buf.tblock[i], i - 1, par_buf.tblock[i - 1]);
        status++;
        break;
      }
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

  start_eps = par_buf.start_eps;
  max_eps = par_buf.max_eps;
  tmax = par_buf.tmax;
  num_block = par_buf.num_block;
  for (i = 0; i < num_block; i++)
    tblock[i] = par_buf.tblock[i];

  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.savefile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------

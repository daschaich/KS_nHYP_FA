// -----------------------------------------------------------------
// SU(3) lattice loading setup
#include "load_includes.h"
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
    printf("SU(3) lattice loading\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
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
// Read in parameters for SU(3) MCRG-blocked measurements
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

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

  startflag = par_buf.startflag;
  strcpy(startfile,par_buf.startfile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------

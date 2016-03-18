// -----------------------------------------------------------------
// SU(3) Wilson flow setup
#include <string.h>
#include "smearing_includes.h"

int initial_set();
#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
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
  node0_printf("Mallocing %.1f MBytes per node for fields\n",
               (double)sites_on_node * 98 * sizeof(su3_matrix) / 1e6);
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
  // Allocate space for temporary su3_matrix field
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed, and send to others
int initial_set() {
  int prompt,status;
  if (mynode() == 0) {
    // Print banner
    printf("Smearing\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
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
// Read in parameters for SU(3) Wilson flow measurements
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Smearing parameters
    IF_OK status += get_i(stdin, prompt, "nsmear", &par_buf.nsmear);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out whether or not to gauge fix to Landau gauge
    IF_OK status += ask_gauge_fix(stdin, prompt, &par_buf.fixflag);

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

  nsmear = par_buf.nsmear;
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  startflag = par_buf.startflag;
  fixflag = par_buf.fixflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.savefile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice( startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------

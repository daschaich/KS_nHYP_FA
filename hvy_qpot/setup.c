// -----------------------------------------------------------------
// SU(3) HYP-smeared static potential setup
#include "hvy_qpot_includes.h"
#include <string.h>

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

  // Print banner, get volume
  prompt = initial_set();
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Allocate space for fields
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, and send to other nodes
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU(3) HYP-smeared static potential\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    printf("nHYP links, reading alpha_smear parameters from infile\n");
    printf("  IR_STAB = %.4g\n", (Real)IR_STAB);
    printf("  EPS_SQ = %.4g\n", (Real)EPS_SQ);
#ifdef NHYP_DEBUG
    printf("  NHYP_DEBUG turned on\n");
#endif

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
  return(prompt);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(3) HYP-smeared static potential
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // Find out whether or not to calculate off-axis Wilson loops
    IF_OK status += get_i(stdin, prompt, "off_axis_flag",
                          &par_buf.off_axis_flag);

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

  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  off_axis_flag = par_buf.off_axis_flag;
  startflag = par_buf.startflag;
  strcpy(startfile, par_buf.startfile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  return 0;
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

  FIELD_ALLOC(tempmat1, matrix);
  FIELD_ALLOC_VEC(Staple3, matrix, 4);

  node0_printf("Mallocing %.1f MBytes per core for fields\n",
               (double)sites_on_node * 61 * sizeof(matrix) / 1e6);
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;             // 1 if it is time to stop
  int nx, ny, nz, nt;       // Lattice dimensions
  int off_axis_flag;        // Whether to calculate off-axis Wilson loops
  int startflag;            // What to do for beginning lattice
  char startfile[MAXFILENAME];

  // Smearing parameters
  Real alpha_hyp0, alpha_hyp1, alpha_hyp2;
} params;
#endif
// -----------------------------------------------------------------

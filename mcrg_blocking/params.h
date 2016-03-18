// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;             // 1 if it is time to stop
  int nx, ny, nz, nt;       // Lattice dimensions
  int startflag;            // What to do for beginning lattice
  char startfile[MAXFILENAME];

  // Smearing parameters -- 100 alphas is probably overkill
  int num_alpha;
  Real alpha_hyp1, alpha_hyp2;
  Real alpha_mcrg[100];
} params;
#endif
// -----------------------------------------------------------------

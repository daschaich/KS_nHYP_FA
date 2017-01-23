// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;             // 1 if it is time to stop
  int nx, ny, nz, nt;       // Lattice dimensions
  int iseed;                // For random numbers
  int startflag;            // What to do for beginning lattice
  char startfile[MAXFILENAME];

  // Smearing parameters
  int Nsmear;
  Real alpha_hyp0, alpha_hyp1, alpha_hyp2;

  // Inversion parameters
  int npbp;                 // Number of stochastic sources
  int niter;                // Maximum number of CG iterations
  int nrestart;             // Maximum number of CG restarts
  Real rsqmin;              // For deciding on convergence

  // Mode number parameters
  int npts, order;
  Real startomega;
  Real spacing;
} params;
#endif
// -----------------------------------------------------------------

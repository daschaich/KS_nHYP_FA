// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;       // 1 if it is time to stop

  // Initialization parameters
  int nx, ny, nz, nt; // Lattice dimensions
  int iseed;          // For random numbers

  int Nvecs;
  Real eig_tol;
  Real error_decr;    // Error decrease per Rayleigh minimization
  int maxIter;        // Maximum number of Rayleigh iterations
  int restart;        // Restart Rayleigh minimization every restart iterations
  int kiters;         // Kalkreuter iterations

  int startflag;      // What to do for beginning lattice
  char startfile[MAXFILENAME];

  // Smearing parameters
  int nsmear;
  Real alpha_hyp0, alpha_hyp1, alpha_hyp2;
} params;
#endif
// -----------------------------------------------------------------

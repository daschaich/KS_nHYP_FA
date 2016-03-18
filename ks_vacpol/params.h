// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;         // 1 if it is time to stop
  int nx, ny, nz, nt;   // Lattice dimensions
  int iseed;            // For random numbers
  Real mass;            // Fermion mass
  int startflag;        // What to do for beginning lattice
  char startfile[MAXFILENAME], outpat[MAXFILENAME];

  // Smearing parameters
  Real alpha_hyp0, alpha_hyp1, alpha_hyp2;

  // Inversion parameters
  int niter;                        // Maximum number of CG iterations
  int nrestart;                     // Maximum number of CG restarts
  int x_src, y_src, z_src, t_src;   // Propagator source location
  Real rsqmin;                      // For deciding on convergence
} params;
#endif
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;             // 1 if it is time to stop
  int nx, ny, nz, nt;       // Lattice dimensions
  int dx, dy, dz, dt;       // Origin of blocked lattice
  int iseed;          // For random numbers
  Real mass;          // Fermion mass
  Real epsilon;
  Real tmax;          // How far to run Wilson flow
  int Nvecs;
  Real eig_tol;
  Real error_decr;    // Error decrease per Rayleigh minimization
  int maxIter;        // Maximum number of Rayleigh iterations
  int restart;        // Restart Rayleigh minimization every restart iterations
  int kiters;         // Kalkreuter iterations
  int startflag;      // What to do for beginning lattice
  int saveflag;       // What to do with lattice at end
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable

  // Smearing parameters
  int nsmear;
  Real alpha_mcrg0, alpha_mcrg1, alpha_mcrg2;
  Real alpha_hyp0, alpha_hyp1, alpha_hyp2;

  // Inversion parameters
  int npbp;                 // Number of stochastic sources
  int niter;                // Maximum number of CG iterations
  int nrestart;             // Maximum number of CG restarts
  Real rsqmin;              // For deciding on convergence
} params;
#endif
// -----------------------------------------------------------------

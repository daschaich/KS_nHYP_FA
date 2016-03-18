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
  int nflavors;       // The number of flavors

  int warms;                // The number of warmup trajectories
  int trajecs;              // The number of real trajectories
  Real traj_length;         // The length of each trajectory
  // Steps for all possible Hasenbusch masses, plus one for gauge
  int nsteps[MAX_MASSES + 1];
  int propinterval;         // Number of trajectories between measurements
  int startflag;            // What to do for beginning lattice
  int saveflag;             // What to do with lattice at end
  int num_masses;           // Number of Hasenbusch masses
  Real beta, mass;          // Gauge coupling, fermion mass
  Real beta_a;              // Adjoint-to-fundamental gauge coupling

  // Smearing parameters
  int Nsmear;
  Real alpha_hyp0, alpha_hyp1, alpha_hyp2;

  // Inversion parameters
  int npbp;                     // Number of stochastic sources
  int niter;                    // Maximum number of CG iterations
  int nrestart;                 // Maximum number of CG restarts
  Real rsqmin;                  // For deciding on convergence
  Real MH;                      // Hasenbusch mass
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable
} params;
#endif
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"  // For MAXFILENAME
#include "../include/io_lat.h"  // For gauge_file
#include "../include/io_lat.h"    // For gauge_file
#include "../include/su3.h"
#include "../include/random.h"    // For double_prn
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

#ifdef SITERAND
  // The state information for a random number generator
  double_prn site_prn;
  // Align to double word boundary (kludge for Intel compiler)
  int space1;
#endif

  // Gauge field
  su3_matrix link[4];

  // Staggered phases, which have been absorbed into the matrices
  // Also the antiperiodic boundary conditions
  Real phase[4];

  // Program-dependent fields
  // Staggered complex vectors
  su3_vector g_rand;      // Gaussian random vector
  //su3_vector *rand;  // Gaussian random vectors
  su3_vector chi;         // Source vector
  su3_vector psi;         // Solution vector
  su3_vector p;           // CG change vector
  su3_vector mp;          // CG temporary vector
  su3_vector r;           // CG residual vector

  su3_vector R1;          // Passed as destination for intermediate X
  su3_vector R2;          // Passed as destination for step
  su3_vector Xsrc;        // X.src in step.c
  su3_vector bj, bjp1, bjp2;  // b[j], b[j + 1], b[j + 2] for Clenshaw
  su3_vector Zbjp1;           // Z.b[j + 1] in step.c

  // Temporary vectors and matrices
  su3_vector tempvec[4];  // One for each direction
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int volume;          // Volume of lattice
EXTERN int npbp;            // Number of stochastic sources
EXTERN int iseed;           // Random number seed
EXTERN int niter, nrestart;
EXTERN Real beta, beta_a, mass, rsqmin;
EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME];
EXTERN int startflag;   // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int total_iters;
EXTERN int phases_in;   // 1 if KS and BC phases absorbed into matrices

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

EXTERN gauge_file *startlat_p;

// Each node maintains a structure with
// the pseudorandom number generator state
EXTERN double_prn node_prn;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 8   // Needed by ../generic/make_lattice.c
EXTERN char **gen_pt[N_POINTERS];

EXTERN su3_matrix *gauge_field[4];
EXTERN su3_matrix *gauge_field_thin[4];

// nHYP stuff -- hyplinks, Staples, tempmat used by generic_nhyp/block_nhyp.c
EXTERN int nsmear;
EXTERN Real alpha_smear[3];
EXTERN su3_matrix *hyplink1[4][4];
EXTERN su3_matrix *hyplink2[4][4];
EXTERN su3_matrix *Staple1[4][4];
EXTERN su3_matrix *Staple2[4][4];
EXTERN su3_matrix *Staple3[4];

EXTERN su3_matrix *tempmat1;

// Mode number stuff
EXTERN int Npts;
EXTERN Real M, spacing;
EXTERN su3_vector **source;

// Step function stuff
// For now Norder selects between options hard-coded in coeffs.c
EXTERN int Norder;
EXTERN double epsilon;
EXTERN double delta;    // Unused, but may be useful to record in the output
EXTERN double starSq, star;     // Ratio (Omega / Omega_*)^2 and its sqrt
EXTERN double *coeffs;
#endif // _LATTICE_H
// -----------------------------------------------------------------


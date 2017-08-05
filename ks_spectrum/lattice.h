// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    // For MAXFILENAME
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
  matrix link[4];

  // Staggered phases, which have been absorbed into the matrices
  // Also the antiperiodic boundary conditions
  Real phase[4];

  // Staggered complex vectors
  vector g_rand;      // Gaussian random vector
  vector chi;         // Gaussian random source vector
  vector psi;         // Solution vector = Kinverse * chi
  vector p;           // CG change vector
  vector mp;          // CG temporary vector
  vector r;           // CG residual vector
  vector ttt;         // For ../generic_ks/mat_invert.c
                      // and ../generic_ks/flavor_ops.c
  vector M_inv;

  // Spectrum stuff
  vector prop, aprop;       // For ../generic_ks/spectrum_nlpi2.c
  vector propmat[3];        // For three source colors
  vector propmat2[3];       // For three source colors in nl_spectrum.c

  // Temporary vectors (one for each direction) and matrix
  vector tempvec[4];
  matrix gfix_scratch;    // For gauge fixing
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int npbp;            // Number of stochastic sources
EXTERN int niter, nrestart;
EXTERN Real beta, mass, rsqmin;

EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME];
EXTERN int startflag;   // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int total_iters;
EXTERN int phases_in;   // 1 if KS and BC phases absorbed into links

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Each node maintains a structure
// with the pseudorandom number generator state
EXTERN double_prn node_prn;

EXTERN gauge_file *startlat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 16   // Needed by ../generic_ks/nl_spectrum.c
EXTERN char **gen_pt[N_POINTERS];

EXTERN matrix *gauge_field[4];
EXTERN matrix *gauge_field_thin[4];

// Temporary matrices
EXTERN matrix *tempmat, *tempmat2;

// Spectrum stuff -- source timeslice and increment
EXTERN int src_start, src_inc, n_src;

// nHYP stuff
EXTERN double alpha_smear[3];
EXTERN matrix *hyplink1[4][4];
EXTERN matrix *hyplink2[4][4];

EXTERN matrix *Staple1[4][4];
EXTERN matrix *Staple2[4][4];
EXTERN matrix *Staple3[4];
#endif
// -----------------------------------------------------------------

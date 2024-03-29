// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"   // For MAXFILENAME
#include "../include/io_lat.h"   // For gauge_file
#include "../include/su3.h"
#include "../include/random.h"   // For double_prn
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
  vector psi;         // Solution vector
  vector chi;         // Source vector
  vector temp;        // For Matrix_Vec_mult
  vector tempvec[4];  // For dslash
  vector ttt;         // For generic_ks/flavor_ops.c
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int nx, ny, nz, nt;    // Lattice dimensions
EXTERN int volume;            // Volume of lattice
EXTERN int iseed;             // Random number seed
EXTERN Real resid;
EXTERN int total_iters;
EXTERN int phases_in;   // Flag if KS and BC phases absorbed into matrices

EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME];
EXTERN int startflag;   // Beginning lattice: CONTINUE, RELOAD, FRESH

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Each node maintains a structure with
// the pseudorandom number generator state
EXTERN double_prn node_prn;

EXTERN gauge_file *startlat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// We need 8 mostly, but 10 for force_nhyp and 16 for nl_spectrum
#define N_POINTERS 8  // Needed by ../generic/make_lattice.c
EXTERN char **gen_pt[N_POINTERS];

// Eigenvalue stuff
EXTERN int Nvecs, block;
EXTERN double *eigVal, start;
EXTERN vector **eigVec;
EXTERN Real eig_tol;          // Tolerance for the eigenvalue computation
EXTERN int maxIter;           // Maximum iterations

EXTERN matrix *gauge_field[4];
EXTERN matrix *gauge_field_thin[4];

// nHYP stuff -- hyplinks, Staples, tempmat used by generic_nhyp/block_nhyp.c
EXTERN int nsmear;
EXTERN double alpha_smear[3];
EXTERN matrix *hyplink1[4][4];
EXTERN matrix *hyplink2[4][4];
EXTERN matrix *Staple1[4][4];
EXTERN matrix *Staple2[4][4];
EXTERN matrix *Staple3[4];

EXTERN matrix *tempmat;

#endif // _LATTICE_H
// -----------------------------------------------------------------

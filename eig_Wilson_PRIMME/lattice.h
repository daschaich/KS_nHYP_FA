// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"   // For MAXFILENAME
#include "../include/io_lat.h"   // For gauge_file
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

  // Gauge field
  matrix link[4];

  // Wilson complex vectors
  wilson_vector psi;  // Solution vector
  wilson_vector chi;  // Source vector
  wilson_vector mp;   // Another CG vector
  wilson_vector tmp;  // Another temporary CG vector for D D^dag

  // Wilson complex half-vector
  // MAXHTMP defined in defines.h; must be at least 8
  half_wilson_vector htmp[MAXHTMP];    // Needed by dslash_w_site
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
EXTERN Real resid;
EXTERN int total_iters;
EXTERN Real mass;
EXTERN Real kappa;

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
EXTERN int Nvecs;
EXTERN double *eigVal;
EXTERN wilson_vector **eigVec;
EXTERN wilson_vector *tmpvec; // Used in sign check in eig.c
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

EXTERN matrix *tempmat1;

#endif // _LATTICE_H
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"  // For MAXFILENAME
#include "../include/io_lat.h"  // For gauge_file
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

  // No random numbers
  // Doubled gauge field provides workspace
  su3_matrix link[4];
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

// Global variables
EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int volume;          // Volume of lattice

EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME];
EXTERN int startflag;   // Beginning lattice: CONTINUE, RELOAD, FRESH

// Some of these global variables are node dependent
// They are set in "setup_layout()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

EXTERN gauge_file *startlat_p;

// No random numbers

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 8   // Needed by ../generic/make_lattice.c
EXTERN char **gen_pt[N_POINTERS];
#endif // _LATTICE_H
// -----------------------------------------------------------------

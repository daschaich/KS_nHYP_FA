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
  // Gauge field
  matrix link[4];

  // Program-dependent fields -- accumulators for Wilson loops
  matrix s_link, s_link_f, t_link_f, diag, staple;
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
EXTERN int off_axis_flag;   // Whether to calculate off-axis loops
EXTERN int tot_smear;       // Label output with number of smearing steps

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

// No random numbers

EXTERN gauge_file *startlat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 8   // Needed by ../generic/make_lattice.c
EXTERN char **gen_pt[N_POINTERS];

EXTERN matrix *gauge_field[4];
EXTERN matrix *gauge_field_thin[4];

// nHYP stuff
EXTERN double alpha_smear[3];
EXTERN matrix *hyplink1[4][4];
EXTERN matrix *hyplink2[4][4];

EXTERN matrix *Staple1[4][4];
EXTERN matrix *Staple2[4][4];
EXTERN matrix *Staple3[4];

EXTERN matrix *tempmat;
#endif // _LATTICE_H
// -----------------------------------------------------------------

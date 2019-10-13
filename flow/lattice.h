// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/dirs.h"      // For NDIMS
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

  // HYP stuff needed for MCRG measurements
  // Tripled gauge field provides workspace and backup storage
  matrix link[12];
  matrix hyplink1[12], hyplink2[12];
  matrix FS[6];       // Field strength for F^2 and topological charge
} site;

// Defines for index on field_strength
#define FS_XY 0
#define FS_XZ 1
#define FS_YZ 2
#define FS_XT 3
#define FS_YT 4
#define FS_ZT 5
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
EXTERN Real one_ov_vol;

EXTERN double g_ssplaq, g_stplaq;   // Global plaqs for I/O
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];  // ILDG LFN if applicable
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag;     // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int saveflag;      // 1 if we will save the lattice;

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

EXTERN gauge_file *startlat_p;

// Loop stuff
#define nloop 6
#define nreps 1
#define max_num 300

// Global definitions for general action (checked?)
EXTERN int loop_ind[nloop][10], loop_length[nloop];
EXTERN int loop_table[nloop][max_num][10];
EXTERN int loop_num[nloop], loop_char[max_num];
EXTERN Real loop_coeff[nloop][nreps];
EXTERN int ch, loop_ch[nloop][max_num];
EXTERN Real loop_term[max_num][nreps];
EXTERN int hyp1ind[4][4][4];
EXTERN int hyp2ind[4][4];

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 8   // Needed by ../generic/make_lattice.c
EXTERN char **gen_pt[N_POINTERS];

// Temporary fields
// tempsym will only be malloc'd for Symanzik flow 2-hop staples
#define N_TEMPSYM 5   // Needed by ../generic/make_lattice.c
EXTERN matrix *tempmat, *tempmat2, *tempsym[N_TEMPSYM];

// Wilson and Symanzik flow stuff
EXTERN Real tmax, start_eps, max_eps, epsilon;
EXTERN matrix *S[NDIMS];
EXTERN anti_hermitmat *A[NDIMS];
EXTERN int flowflag;    // Specify Wilson vs. Symanzik flow

// MCRG blocking stuff
EXTERN Real alpha_smear[3];
EXTERN int num_block;
EXTERN Real tblock[100];
#endif
// -----------------------------------------------------------------

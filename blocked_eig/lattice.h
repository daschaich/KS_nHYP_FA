// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"  // For MAXFILENAME
#include "../include/generic_quark_types.h"
#include "../include/io_lat.h"  // For gauge_file
#include "../include/su3.h"
#include "../include/random.h"    // For double_prn

// Defines for index on field_strength
#define FS_XY 0
#define FS_XZ 1
#define FS_YZ 2
#define FS_XT 3
#define FS_YT 4
#define FS_ZT 5
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

  // Doubled gauge field provides workspace
  su3_matrix link[8];

  // Staggered phases, which have been absorbed into the matrices
  // Also the antiperiodic boundary conditions
  Real phase[4];

  // Staggered complex vectors
  su3_vector g_rand;            // Gaussian random vector
  su3_vector chi;               // Gaussian random source vector
  su3_vector psi;               // Solution vector = Kinverse * chi
  su3_vector p;                 // CG change vector
  su3_vector mp;                // CG temporary vector
  su3_vector r;                 // CG residual vector
  su3_vector ttt1[1];           // For ../generic_ks/mat_invert.c
  su3_vector M_inv;

  // Temporary vectors and matrices
  su3_vector temp;                  // For Matrix_Vec_mult
  su3_vector tempvec[4];            // For dslash
  su3_matrix fieldstrength[6];      // For Wilson flow

  // Accumulators for Wilson and Polyakov loops
  su3_matrix hyplink1[12], hyplink2[12], tempmat1, tempmat2, staple;
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int dx, dy, dz, dt;  // Hypercube origin
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int npbp;            // Number of stochastic sources
EXTERN int niter, nrestart;
EXTERN Real mass, rsqmin;

EXTERN int total_iters;
EXTERN int phases_in;   // Flag if KS and BC phases absorbed into matrices
EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];  // ILDG LFN if applicable
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag;   // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int saveflag;  // 1 if we will save the lattice;

// Eigenvalue stuff
EXTERN int Nvecs;
EXTERN double *eigVal;
EXTERN su3_vector **eigVec;
EXTERN Real eig_tol;      // Tolerance for the eigenvalue computation
EXTERN Real error_decr;   // Error decrease per Rayleigh minimization
EXTERN int maxIter;       // max  Rayleigh iterations
EXTERN int restart;       // Restart  Rayleigh every so many iterations
EXTERN int kiters;        // Kalkreuter iterations
// Some of these global variables are node dependent
// They are set in "setup_layout()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

// Each node maintains a structure with
// the pseudorandom number generator state
EXTERN double_prn node_prn;

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

EXTERN su3_matrix *gauge_field[4];
EXTERN su3_matrix *gauge_field_thin[4];

// nHYP stuff
EXTERN double alpha_store[3];
EXTERN int nsmear;
EXTERN double alpha_smear[3];
EXTERN su3_matrix *hyplink1[4][4];
EXTERN su3_matrix *hyplink2[4][4];
EXTERN su3_matrix *Staple1[4][4];
EXTERN su3_matrix *Staple2[4][4];
EXTERN su3_matrix *Staple3[4];

// Used and freed in Wilson flow calculation
EXTERN su3_matrix *tempmat1;
EXTERN su3_matrix *tempmat2;    // Used in Polyakov loop calculation

// Wilson flow stuff
EXTERN double epsilon, tmax;
#endif // _LATTICE_H
// -----------------------------------------------------------------

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
  su3_matrix link[4];

  // Program-dependent fields
#ifdef HMC_ALGORITHM
  su3_matrix old_link[4];  // For accept/reject
#endif

  // Antihermitian momentum matrices in each direction
  anti_hermitmat mom[4];

  // Staggered phases, which have been absorbed into the matrices
  // Also the antiperiodic boundary conditions
  Real phase[4];

  // Staggered single-spin three-element complex vectors
  su3_vector g_rand;            // Gaussian random vector
  su3_vector chi[MAX_FIELDS][MAX_MASSES];  // Gaussian random source vector
  su3_vector psi[MAX_FIELDS][MAX_MASSES];  // Solution vector = Kinverse * chi
  su3_vector p;                 // CG change vector
  su3_vector mp;                // CG temporary vector
  su3_vector r;                 // CG residual vector
  su3_vector ttt[MAX_FIELDS][MAX_MASSES];  // Temporary vectors for D*psi
#ifdef PHI_ALGORITHM
  su3_vector old_psi[MAX_FIELDS][MAX_MASSES];  // For predicting next psi
#endif
  su3_vector M_inv;

  // Temporary vectors and matrices
  su3_vector tempvec[4];  // One for each direction
  su3_matrix tempmat1, tempmat2, staple;
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
EXTERN int warms, trajecs, niter, nrestart, propinterval, nflavors;
EXTERN int full_fields, half_fields;    // Derived from nflavors
EXTERN Real traj_length;

// Global Hasenbusch variables
EXTERN int num_masses;    // Maximum number of masses, <= MAX_MASSES
EXTERN Real MH;
EXTERN int nsteps[MAX_MASSES + 1];
EXTERN Real max_gf, max_ff[2];

// Global action and evolution variables
EXTERN Real rsqmin, beta, beta_a, mass;
EXTERN Real epsilon;
EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];   // ILDG LFN if applicable
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag; // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int saveflag;  // 1 if we will save the lattice
EXTERN int total_iters;
EXTERN int phases_in; // 1 if KS and BC phases absorbed into matrices

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// We need 8 mostly, but 10 for force_nhyp and 16 for nl_spectrum
#define N_POINTERS 10
EXTERN char **gen_pt[N_POINTERS];

EXTERN su3_matrix *gauge_field[4];
EXTERN su3_matrix *gauge_field_thin[4];
EXTERN su3_matrix *gauge_field_save[4];
EXTERN su3_matrix *gauge_field_temp[4];

// nHYP stuff
EXTERN int Nsmear;
EXTERN double alpha_smear[3];
EXTERN su3_matrix *hyplink1[4][4];
EXTERN su3_matrix *hyplink2[4][4];
EXTERN su3_matrix *Sigma[4];
EXTERN su3_matrix *SigmaH[4];
EXTERN su3_matrix *SigmaH2[4][4];

EXTERN su3_matrix *Staple1[4][4];
EXTERN su3_matrix *Staple2[4][4];
EXTERN su3_matrix *Staple3[4];

EXTERN su3_matrix *LambdaU[4];
EXTERN su3_matrix *Lambda1[4];
EXTERN su3_matrix *Lambda2[4];

EXTERN su3_matrix *tempmat, *tempmat2;    // Used in Polyakov loop calculation

// Up to 20 concurrent timers for timing
#ifdef TIMING
EXTERN double tmptime[20];
EXTERN double time_dcongrad;      // From congrad.  Use tmptime[0]
EXTERN double time_fermion_force; // From update.  Use tmptime[1]
EXTERN double time_block_nhyp;    // From block_nhyp.  Use tmptime[3]
EXTERN double time_compute_fhb;   // Not currently using tmptime[4]
#endif

#endif
// -----------------------------------------------------------------

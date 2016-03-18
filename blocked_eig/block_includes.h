// -----------------------------------------------------------------
// Include files for SU(3) MCRG-blocked measurements
#include <stdio.h>
#include <stdlib.h>
#include <string.h>             // For strlen and memcpy
#include <math.h>
#include "../include/config.h"  // Keep this first
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
#include "jacobi.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
void make_fields();
void free_fields();
int readin(int prompt);
void stout_step_rk(su3_matrix *S[4], anti_hermitmat *A[4]);
void staple(su3_matrix *stp[4]);

void grsource_imp();
void grsource_eig();
int meas_link(field_offset chi, field_offset psi, Real mass);

// CG stuff
int ks_congrad(field_offset src, field_offset dest, Real M, int parity);
void doCG(int level, Real M, int *iters);

// Dslash stuff and associated helper functions
// (latter were in v6/generic_ks/d_congrad5.c for some reason)
void dslash(field_offset chi, field_offset psi, int parity);
void dslash_special(field_offset chi, field_offset psi, int parity,
                    msg_tag **tag, int start);
void clear_latvec(field_offset v, int parity);
void copy_latvec(field_offset src, field_offset dest, int parity);
void scalar_mult_latvec(field_offset src, Real scalar,
                        field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
                            Real scalar, field_offset dest, int parity);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// nHYP stuff specific to for MCRG-blocked measurements
void clear_disp(int *disp);
void diag_su3(su3_matrix* Q, complex *f);
void block_nhyp_mcrg(int num, int block);
void block_mcrg(int num, int block);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wilson and Polyakov loop stuff
void make_loop_table2();
void blocked_gauge_loops(int block, double *result);
void path(int *dir, int *sign, int length, su3_matrix *resmat);
void blocked_path(int block, int *dir, int *sign,
                  int length, su3_matrix *resmat);

complex blocked_ploop(int block, int dir);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Eigenvalue routines
int Kalkreuter(su3_vector **eigVec, double *eigVal, Real Tolerance,
               Real RelTol, int Nvecs, int maxIter, int restart, int kiters);

void measure_chirality(su3_vector *src, double *chirality);
void print_densities(su3_vector *src, char *tag, int y, int z, int t);
// -----------------------------------------------------------------

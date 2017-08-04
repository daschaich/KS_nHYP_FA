// -----------------------------------------------------------------
// Includes and definitions for nHYP smearing
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/config.h"  // Keep this first
#include <lattice.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/dirs.h"
#include "../include/comdefs.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"    // For rephase
#include "../include/generic_nhyp.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
void unphased_block_and_fatten();
void block_and_fatten();
void just_fatten();
void leanlinks();
void block_nhyp();
void staple_nhyp(int dir1, int dir2, matrix *lnk1,
                 matrix *lnk2, matrix *stp);

void block_nhyp1();
void block_nhyp2();
void block_nhyp3();
#ifndef NHYP_DEBUG
void compute_fhb(matrix *Q, Real *f, Real b[3][3], int compute_b);
#else
void compute_fhb(matrix *Omega, matrix *Q, Real *f,
                 Real b[3][3], int compute_b);
#endif
void scalar_add_diag_su3(matrix *a, Real s);
void c_scalar_add_diag_su3(matrix *a, complex *c);
// -----------------------------------------------------------------

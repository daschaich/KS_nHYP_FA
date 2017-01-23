// -----------------------------------------------------------------
// Includes and definitions for nHYP smearing
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/config.h"  // Keep this first
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/dirs.h"
#include "../include/comdefs.h"
#include "../include/generic.h"
#include <lattice.h>
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
void block_and_fatten();
void unphased_block_and_fatten();
void just_fatten();
void leanlinks();
void block_nhyp();
void boundary_flip_nhyp(int sign);
void staple_nhyp(int dir1, int dir2, su3_matrix *lnk1,
                 su3_matrix *lnk2, su3_matrix *stp);

void block_nhyp1();
void block_nhyp2();
void block_nhyp3();
#ifndef NHYP_DEBUG
void compute_fhb(su3_matrix *Q, Real *f, Real b[3][3], int compute_b);
#else
void compute_fhb(su3_matrix *Omega, su3_matrix *Q, Real *f,
                 Real b[3][3], int compute_b);
#endif
// -----------------------------------------------------------------

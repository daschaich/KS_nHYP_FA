// -----------------------------------------------------------------
// Include files for SU(3) Wilson flow measurements
#include <stdio.h>
#include <stdlib.h>
#include <string.h>             // For strlen
#include <math.h>
#include "../include/config.h"  // Keep this first
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "defines.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);
void wflow();
void stout_step_rk();
void staple(matrix *stp[NDIMS]);

// nHYP stuff specific to for MCRG-blocked measurements
void mcrg_block(Real t, int blmax);
void clear_disp(int *disp);
void diag_su3(matrix *Q, complex *f);
void block_nhyp_mcrg(int num, int block);
void block_mcrg(int num, int block);

// Wilson loop stuff
void make_loop_table2();
void blocked_gauge_loops(int block, double *result);
void path(int *dir, int *sign, int length);
void blocked_path(int block, int *dir, int *sign, int length);

// Polyakov loop stuff
complex blocked_ploop(int block, int dir);
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Include files for staggered eigenvalues
#include <stdio.h>
#include <stdlib.h>
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
int readin(int prompt);
void grsource_imp();
void dslash(field_offset chi, field_offset psi, int parity);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Eigenvalue routines
int Kalkreuter(vector **eigVec, double *eigVal, Real Tolerance,
               Real RelTol, int Nvecs, int maxIter, int restart, int kiters);

void measure_chirality(vector *src, double *chirality);
void print_densities(vector *src, char *tag, int y, int z, int t);
// -----------------------------------------------------------------

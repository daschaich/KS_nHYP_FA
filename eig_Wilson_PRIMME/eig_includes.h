// -----------------------------------------------------------------
// Include files for SU(3) hermitian-Wilson eigenvalues
#include <stdio.h>
#include <stdlib.h>
#include <string.h>             // For strlen
#include <math.h>
#include "../include/config.h"  // Keep this first
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
#include "primme.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);

// Eigenvalue routines
int make_evs(int Nvecs, wilson_vector **eigVec, double *eigVal);
void measure_chirality(wilson_vector *src, double *chirality) ;
void print_densities(wilson_vector *src, char *tag, int y, int z, int t);
// -----------------------------------------------------------------

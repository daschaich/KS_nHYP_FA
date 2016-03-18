// -----------------------------------------------------------------
// Include files for SU(3) HYP-smeared static potential
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
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);

// Currently fixing tot_smear=1, but allows some cheap flexibility
void w_loop1(int tot_smear);
void w_loop2(int tot_smear);
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Include files for staggered eigenvalues
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
#include "../include/generic_ks.h"
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);
void grsource_plain();

// Dslash and CG stuff, also associated helper functions
// (latter were in v6/generic_ks/d_congrad5.c for some reason)
int ks_congrad(field_offset src, field_offset dest, Real M, int parity);
void dslash(field_offset chi, field_offset psi, int parity);
void cleanup_one_gather_set(msg_tag *tags[]);
void clear_latvec(field_offset v, int parity);
void copy_latvec(field_offset src, field_offset dest, int parity);
void scalar_mult_latvec(field_offset src, Real scalar,
                        field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
                            Real scalar, field_offset dest, int parity);

// Mode number and step function stuff
void coefficients();
void step(field_offset src, field_offset dest);
// -----------------------------------------------------------------

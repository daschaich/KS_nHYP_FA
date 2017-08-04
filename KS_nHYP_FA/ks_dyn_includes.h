// -----------------------------------------------------------------
// Include files for dynamical nHYP-smeared staggered HMC
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
#include "../include/generic_ks.h"
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high-level code
int setup();
int readin(int prompt);

void rephase(int flag);
void gauge_field_copy(field_offset src, field_offset dest);

// CG stuff
int ks_congrad(field_offset src, field_offset dest, Real M, int parity);
void doCG(int level, Real M, int *iters);
int update_step(Real *old_cg_time, Real *cg_time, Real *next_cg_time,
                double *fnorm, double *gnorm);

// Dslash stuff and associated helper functions
// (latter were in v6/generic_ks/d_congrad5.c for some reason)
void dslash(field_offset chi, field_offset psi, int parity);
void dslash_special(field_offset chi, field_offset psi, int parity,
                    msg_tag **tag, int start);
void clear_latvec(field_offset v, int parity);
void copy_latvec(field_offset src, field_offset dest, int parity);
void scalar_mult_latvec(field_offset src, Real scalar,
                        field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset a, field_offset b,
                            Real scalar, field_offset c, int parity);

// Evolution stuff
void grsource_imp(field_offset dest, Real mass, int parity);
void grsource_Hasen(field_offset dest, int *iters, int parity);
double action();
void plaquette_a(double *ss_plaq, double *st_plaq);
void meas_plaq();
int meas_link(field_offset chi, field_offset psi, Real mass);
int update();
void update_h(Real eps);
void update_u(Real eps);
double gauge_force(Real eps);
double fermion_force(int level, Real eps);

// nHYP force stuff
void block_N_fatten(int N);
void nhyp_force1();
// -----------------------------------------------------------------

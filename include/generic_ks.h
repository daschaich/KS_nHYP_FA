// -----------------------------------------------------------------
// Macros and declarations for generic KS routines
#ifndef _GENERIC_KS_H
#define _GENERIC_KS_H

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/generic_quark_types.h"
#include "../include/comdefs.h"
#ifdef HAVE_QOP
#include <qop.h>
#endif
#ifdef HAVE_QIO
#include <qio.h>
#endif

#define ALL_T_SLICES -1
#define MAXDESCRP 128
#define MAXSRCLABEL 8
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// rephase.c
void phaseset();
void rephase(int flag);

void prefetch_vector(vector *);
void prefetch_matrix(matrix *);

// f_meas.c
int f_meas_imp(field_offset chi_off, field_offset psi_off, Real mass);

// fpi_2.c
int fpi_2(Real *masses);

/* flavor_ops.c */
void sym_shift(int dir, field_offset src,field_offset dest);
void zeta_shift(int n, int *d, field_offset src, field_offset dest);
void eta_shift(int n, int *d, field_offset src, field_offset dest);

void mult_flavor_vector(int mu, field_offset src, field_offset dest);
void mult_flavor_tensor(int mu, int nu, field_offset src, field_offset dest);
void mult_flavor_pseudovector(int mu, field_offset src, field_offset dest);
void mult_flavor_pseudoscalar(field_offset src, field_offset dest);

void mult_spin_vector(int mu, field_offset src, field_offset dest);
void mult_spin_tensor(int mu, int nu, field_offset src, field_offset dest);
void mult_spin_pseudovector(int mu, field_offset src, field_offset dest);
void mult_spin_pseudoscalar(field_offset src, field_offset dest);

// mat_invert.c
int mat_invert_uml(field_offset src, field_offset dest,
                   field_offset temp, Real mass);

// nl_spectrum.c
int nl_spectrum(Real vmass, field_offset tempvec1, field_offset tempvec2,
                field_offset tempmat1, field_offset tempmat2);

/* show_generic_ks_opts.c */
void show_generic_ks_opts();

// spectrum2.c
int spectrum2(Real vmass, field_offset temp1, field_offset temp2);

// spectrum_nlpi2.c
int spectrum_nlpi2(Real fmass, Real amass);
#endif
// -----------------------------------------------------------------

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
// congrad.c
void dslash(field_offset chi, field_offset psi, int parity);
void dslash_special(field_offset chi, field_offset psi, int parity,
                    msg_tag **tag, int start);
int ks_congrad(field_offset src, field_offset dest, Real M, int parity);

// congrad_helpers.c
void cleanup_one_gather_set(msg_tag *tags[]);
void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]);
void clear_latvec(field_offset v, int parity);
void copy_latvec(field_offset src, field_offset dest, int parity);
void scalar_mult_latvec(field_offset src, Real scalar,
                        field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset a, field_offset b,
                            Real scalar, field_offset c, int parity);

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
int nl_spectrum(Real vmass, field_offset tempvec1, field_offset tempvec2);

// rephase.c
void rephase(int flag);

// spectrum2.c
int spectrum2(Real vmass, field_offset temp1, field_offset temp2);

// spectrum_nlpi2.c
int spectrum_nlpi2(Real fmass, Real amass);
#endif
// -----------------------------------------------------------------

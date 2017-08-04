#ifndef _GENERIC_QOPQDP_H
#define _GENERIC_QOPQDP_H
/******************** generic_qopqdp.h ******************************/

#include <qop.h>
#include "../include/su3.h"

/*********************************************************************/
/* These definitions should come from qopqdp.h 
   Then we can put these (or whatever we want) in generic_qopmilc.h */

typedef struct {
   float real;	
   float imag;
} QOPRAW_F_Complex;

typedef struct {
   double real;	
   double imag;
} QOPRAW_D_Complex;

#define QOPRAW_real(a) (a).real
#define QOPRAW_imag(a) (a).imag
#define QOPRAW_c_eq_r_plus_ir(c,a,b) {QOPRAW_real(c) = (a); QOPRAW_imag(c) = (b);}

#define NS 4

#define rawcvdef(P,C,NC) QOPRAW_##P##_Complex c[NC]             /* ColorVector */
#define rawdfdef(P,C,NC) QOPRAW_##P##_Complex d[NC][NS]         /* DiracFermion */
#define rawgfdef(P,C,NC) QOPRAW_##P##_Complex e[NC][NC]         /* GaugeField */
#define rawfodef(P,C,NC) QOPRAW_##P##_Complex e[NC][NC]         /* Force */

#define QOPRAW_elem_V(a,ic) (a).c[ic]
#define QOPRAW_elem_D(a,ic,is) (a).d[ic][is]
#define QOPRAW_elem_M(a,ic,jc) (a).e[ic][jc]

/* everything below here is determined from the above macros */

/* SU(3) single precision */
typedef struct { rawcvdef(F,3,3); } QOPRAW_F3_ColorVector;
typedef struct { rawdfdef(F,3,3); } QOPRAW_F3_DiracFermion;
typedef struct { rawgfdef(F,3,3); } QOPRAW_F3_GaugeField;
typedef struct { rawfodef(F,3,3); } QOPRAW_F3_Force;

#define QOPRAW_F3_elem_V(a,ic) QOPRAW_elem_V(a,ic)
#define QOPRAW_F3_elem_D(a,ic,is) QOPRAW_elem_D(a,ic,is)
#define QOPRAW_F3_elem_M(a,ic,jc) QOPRAW_elem_M(a,ic,jc)

/* SU(3) double precision */
typedef struct { rawcvdef(D,3,3); } QOPRAW_D3_ColorVector;
typedef struct { rawdfdef(D,3,3); } QOPRAW_D3_DiracFermion;
typedef struct { rawgfdef(D,3,3); } QOPRAW_D3_GaugeField;
typedef struct { rawfodef(D,3,3); } QOPRAW_D3_Force;

#define QOPRAW_D3_elem_V(a,ic) QOPRAW_elem_V(a,ic)
#define QOPRAW_D3_elem_D(a,ic,is) QOPRAW_elem_D(a,ic,is)
#define QOPRAW_D3_elem_M(a,ic,jc) QOPRAW_elem_M(a,ic,jc)

/*********************************************************************/

/* map_milc_to_qopqdp.c */

fmatrix ** create_raw4_F_G (void);
fmatrix ** create_raw4_F_F (void);
fvector * create_raw_F_V(void);
QOPRAW_F3_DiracFermion * create_raw_F_D(void);

void destroy_raw4_F_G (fmatrix *raw[]);
void destroy_raw4_F_F (fmatrix *raw[]);
void destroy_raw_F_V (fvector *raw);
void destroy_raw_F_D (QOPRAW_F3_DiracFermion *raw);

fmatrix ** create_raw4_F_G_from_site(field_offset src, int milc_parity);
fmatrix ** create_raw4_F_F_from_site(field_offset src, int milc_parity);
fvector * create_raw_F_V_from_site(field_offset src, int milc_parity);
QOPRAW_F3_DiracFermion * create_raw_F_D_from_site(field_offset src, int milc_parity);

fmatrix ** create_raw4_F_G_from_field(matrix *src, int milc_parity);
fmatrix ** create_raw4_F_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
fvector * create_raw_F_V_from_field(vector *src, int milc_parity);
QOPRAW_F3_DiracFermion * create_raw_F_D_from_field(wilson_vector *src, int milc_parity);

void unload_raw4_F_G_to_site(field_offset dest, fmatrix *raw[], 
			     int milc_parity);
void unload_raw4_F_F_to_site(field_offset dest, fmatrix *raw[], 
			     int milc_parity);
void unload_raw_F_V_to_site(field_offset dest, fvector *raw, 
			    int milc_parity);
void unload_raw_F_D_to_site(field_offset dest, QOPRAW_F3_DiracFermion *raw, 
			    int milc_parity);

void unload_raw4_F_G_to_field(matrix *dest, fmatrix *raw[], 
			      int milc_parity);
void unload_raw4_F_F_to_field(anti_hermitmat *dest, fmatrix *raw[], 
			      int milc_parity);
void unload_raw_F_V_to_field(vector *dest, fvector *raw, 
			     int milc_parity);
void unload_raw_F_D_to_field(wilson_vector *dest, QOPRAW_F3_DiracFermion *raw, 
			     int milc_parity);

dmatrix ** create_raw4_D_G (void);
void destroy_raw4_D_G (dmatrix *raw[]);
dmatrix ** create_raw4_D_F (void);
void destroy_raw4_D_F (dmatrix *raw[]);
dvector * create_raw_D_V(void);
void destroy_raw_D_V (dvector *raw);
void destroy_raw_D_D (QOPRAW_D3_DiracFermion *raw);
dmatrix ** create_raw4_D_G_from_site(field_offset src, int milc_parity);
dmatrix ** create_raw4_D_G_from_field(matrix *src, int milc_parity);
dmatrix ** create_raw4_D_F_from_site(field_offset src, int milc_parity);
dmatrix ** create_raw4_D_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
dvector * create_raw_D_V_from_site(field_offset src, int milc_parity);
dvector * create_raw_D_V_from_field(vector *src, int milc_parity);
QOPRAW_D3_DiracFermion * create_raw_D_D_from_site(field_offset src, int milc_parity);
QOPRAW_D3_DiracFermion * create_raw_D_D_from_field(wilson_vector *src, int milc_parity);
void unload_raw4_D_G_to_site(field_offset dest, dmatrix *raw[], 
			     int milc_parity);
void unload_raw4_D_G_to_field(matrix *dest, dmatrix *raw[], 
			      int milc_parity);
void unload_raw4_D_F_to_site(field_offset dest, dmatrix *raw[], 
			     int milc_parity);
void unload_raw4_D_F_to_field(anti_hermitmat *dest, dmatrix *raw[], 
			      int milc_parity);
void unload_raw_D_V_to_site(field_offset dest, dvector *raw, 
			    int milc_parity);
void unload_raw_D_V_to_field(vector *dest, dvector *raw, 
			     int milc_parity);
void unload_raw_D_D_to_site(field_offset dest, QOPRAW_D3_DiracFermion *raw, 
			    int milc_parity);
void unload_raw_D_D_to_field(wilson_vector *dest, QOPRAW_D3_DiracFermion *raw, 
			     int milc_parity);

#endif /* GENERIC_QOPQDP_H */

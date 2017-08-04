#ifndef _GENERIC_QOPMILC_H
#define _GENERIC_QOPMILC_H
/******************** generic_qopmilc.h *****************************/

#include <qop.h>
#include "../include/su3.h"

/* map_milc_to_qopmilc.c */

fmatrix ** create_raw4_F_G (void);
fmatrix ** create_raw4_F_F (void);
fvector * create_raw_F_V(void);
fwilson_vector * create_raw_F_D(void);

void destroy_raw4_F_G (fmatrix *raw[]);
void destroy_raw4_F_F (fmatrix *raw[]);
void destroy_raw_F_V (fvector *raw);
void destroy_raw_F_D (fwilson_vector *raw);

fmatrix ** create_raw4_F_G_from_site(field_offset src, int milc_parity);
fmatrix ** create_raw4_F_F_from_site(field_offset src, int milc_parity);
fvector * create_raw_F_V_from_site(field_offset src, int milc_parity);
fwilson_vector * create_raw_F_D_from_site(field_offset src, int milc_parity);

fmatrix ** create_raw4_F_G_from_field(matrix *src, int milc_parity);
fmatrix ** create_raw4_F_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
fvector * create_raw_F_V_from_field(vector *src, int milc_parity);
fwilson_vector * create_raw_F_D_from_field(wilson_vector *src, int milc_parity);

void unload_raw4_F_G_to_site(field_offset dest, fmatrix *raw[], 
			     int milc_parity);
void unload_raw4_F_F_to_site(field_offset dest, fmatrix *raw[], 
			     int milc_parity);
void unload_raw_F_V_to_site(field_offset dest, fvector *raw, 
			    int milc_parity);
void unload_raw_F_D_to_site(field_offset dest, fwilson_vector *raw, 
			    int milc_parity);

void unload_raw4_F_G_to_field(matrix *dest, fmatrix *raw[], 
			      int milc_parity);
void unload_raw4_F_F_to_field(anti_hermitmat *dest, fmatrix *raw[], 
			      int milc_parity);
void unload_raw_F_V_to_field(vector *dest, fvector *raw, 
			     int milc_parity);
void unload_raw_F_D_to_field(wilson_vector *dest, fwilson_vector *raw, 
			     int milc_parity);

dmatrix ** create_raw4_D_G (void);
void destroy_raw4_D_G (dmatrix *raw[]);
dmatrix ** create_raw4_D_F (void);
void destroy_raw4_D_F (dmatrix *raw[]);
dvector * create_raw_D_V(void);
void destroy_raw_D_V (dvector *raw);
void destroy_raw_D_D (dwilson_vector *raw);
dmatrix ** create_raw4_D_G_from_site(field_offset src, int milc_parity);
dmatrix ** create_raw4_D_G_from_field(matrix *src, int milc_parity);
dmatrix ** create_raw4_D_F_from_site(field_offset src, int milc_parity);
dmatrix ** create_raw4_D_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
dvector * create_raw_D_V_from_site(field_offset src, int milc_parity);
dvector * create_raw_D_V_from_field(vector *src, int milc_parity);
dwilson_vector * create_raw_D_D_from_site(field_offset src, int milc_parity);
dwilson_vector * create_raw_D_D_from_field(wilson_vector *src, int milc_parity);
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
void unload_raw_D_D_to_site(field_offset dest, dwilson_vector *raw, 
			    int milc_parity);
void unload_raw_D_D_to_field(wilson_vector *dest, dwilson_vector *raw, 
			     int milc_parity);

#endif /* GENERIC_QOPMILC_H */

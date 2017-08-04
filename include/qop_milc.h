#ifndef QOP_MILC_H
#define QOP_MILC_H

#include "../include/su3.h"

/* For MILC test implementation of SciDAC Level 3 routines */

#if QOP_Precision == 1

#define QOP_ColorVector_struct QOP_F3_ColorVector_struct
#define QOP_DiracFermion_struct QOP_F3_DiracFermion_struct
#define QOP_GaugeField_struct QOP_F3_GaugeField_struct
#define QOP_Force_struct QOP_F3_Force_struct
#define QOP_FermionLinksAsqtad_struct QOP_F3_FermionLinksAsqtad_struct

struct QOP_ColorVector_struct {
  fvector *v;
  int evenodd;
};

struct QOP_DiracFermion_struct {
  fwilson_vector *d;
  int evenodd;
};

struct QOP_GaugeField_struct {
  fmatrix *g;
  int evenodd;
};

struct QOP_Force_struct {
  fmatrix *f;
  int evenodd;
};

#else

#define QOP_ColorVector_struct QOP_D3_ColorVector_struct
#define QOP_DiracFermion_struct QOP_D3_DiracFermion_struct
#define QOP_GaugeField_struct QOP_D3_GaugeField_struct
#define QOP_Force_struct QOP_D3_Force_struct
#define QOP_FermionLinksAsqtad_struct QOP_D3_FermionLinksAsqtad_struct

struct QOP_ColorVector_struct {
  dvector *v;
  int evenodd;
};

struct QOP_DiracFermion_struct {
  dwilson_vector *d;
  int evenodd;
};

struct QOP_GaugeField_struct {
  dmatrix *g;
  int evenodd;
};

struct QOP_Force_struct {
  dmatrix *f;
  int evenodd;
};

#endif

struct QOP_FermionLinksAsqtad_struct {
  struct QOP_GaugeField_struct *fat;
  struct QOP_GaugeField_struct *lng;
  int evenodd;
};

/* Precision conversion routines for unpacking the data members above 
   into flat arrays */
/* qop_milc_utilities.c */

matrix * create_links_from_qop_milc_F(fmatrix *src);
void destroy_links_from_qop_milc_F(matrix *g);
matrix *create_links_from_qop_milc_D(dmatrix *src);
void destroy_links_from_qop_milc_D(matrix *g);
vector *create_latvec_from_qop_milc_F(fvector *src);
void destroy_latvec_from_qop_milc_F(vector *v);
vector *create_latvec_from_qop_milc_D(dvector *src);
void destroy_latvec_from_qop_milc_D(vector *v);
void copy_latvec_to_qop_milc_F( fvector *dest, vector *src);
void copy_latvec_to_qop_milc_D( dvector *dest, vector *src);

#include <qop.h>

#endif /* QOP_MILC_H */

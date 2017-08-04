/*********************** map_milc_to_qop.c **************************/
/* Functions for mapping MILC data layouts to raw QOP layouts       */
/* C. DeTar 10/19/2005                                              */

#include "generic_includes.h"
#include "../include/generic_qop.h"

/* Create empty raw links */

#define make_create_raw4(P, T, MILCTYPE) \
MILCTYPE ** \
create_raw4_##P##_##T (void){ \
  MILCTYPE **raw = NULL; \
  int dir; \
  raw = (MILCTYPE **)malloc(4*sizeof(MILCTYPE *)); \
  FORALLUPDIR(dir){ \
    raw[dir] = (MILCTYPE *)malloc(sites_on_node*sizeof(MILCTYPE)); \
    if(raw[dir] == NULL){ \
      printf("create4_raw: No room\n"); \
      return NULL; \
    } \
  } \
  return raw; \
}

/* Destroy raw links */

#define make_destroy_raw4(P, T, MILCTYPE) \
void \
destroy_raw4_##P##_##T (MILCTYPE *raw[]){ \
  int dir; \
  FORALLUPDIR(dir){ \
    if(raw[dir] != NULL) \
      free(raw[dir]); \
  } \
  free(raw); \
}

/* Create empty raw field */

#define make_create_raw(P, T, MILCTYPE) \
MILCTYPE * \
create_raw_##P##_##T(void){ \
  MILCTYPE *raw = NULL; \
  raw = (MILCTYPE *)malloc(sites_on_node*sizeof(MILCTYPE)); \
  if(raw == NULL){ \
    printf("create_raw: No room\n"); \
    return NULL; \
  } \
  return raw; \
}

/* Destroy raw field */

#define make_destroy_raw(P, T, MILCTYPE) \
void \
destroy_raw_##P##_##T (MILCTYPE *raw){ \
  if(raw != NULL) free(raw); \
}

/* Copy types with possible conversion */

/* Convert (or copy) MILC types between specific and prevailing precision */

#if (PRECISION==1)

static void 
f2p_mat(matrix *dest, fmatrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fmatrix));
}

static void 
p2f_mat(fmatrix *dest, matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fmatrix));
}

static void 
d2p_mat(matrix *dest, dmatrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2d_mat(dmatrix *dest, matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
f2p_vec(vector *dest, fvector *src){
  memcpy((void *)dest, (void *)src, sizeof(fvector));
}

static void 
p2f_vec(fvector *dest, vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fvector));
}

static void 
d2p_vec(vector *dest, dvector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
p2d_vec(dvector *dest, vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

#else

static void 
f2p_mat(matrix *dest, fmatrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

/* Convert (or copy) matrix from prevailing to single precision */
static void 
p2f_mat(fmatrix *dest, matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
d2p_mat(matrix *dest, dmatrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dmatrix));
}

static void 
p2d_mat(dmatrix *dest, matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dmatrix));
}

static void 
f2p_vec(vector *dest, fvector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

/* Convert (or copy) vector from prevailing to single precision */
static void 
p2f_vec(fvector *dest, vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
d2p_vec(vector *dest, dvector *src){
  memcpy((void *)dest, (void *)src, sizeof(dvector));
}

static void 
p2d_vec(dvector *dest, vector *src){
  memcpy((void *)dest, (void *)src, sizeof(dvector));
}

#endif

/* Conversions from prevailing MILC formats to specified formats */

#define copy_milc_to_F_G(d,s) p2f_mat(d,s);
#define copy_milc_to_D_G(d,s) p2d_mat(d,s);

static void
copy_milc_to_F_F(fmatrix *dest, anti_hermitmat *src){
  matrix t;
  uncompress_anti_hermitian( src, &t );
  p2f_mat( dest, &t );
}

static void
copy_milc_to_D_F(dmatrix *dest, anti_hermitmat *src){
  matrix t;
  uncompress_anti_hermitian( src, &t );
  p2d_mat( dest, &t );
}

#define copy_milc_to_F_V(d,s) p2f_vec(d,s);
#define copy_milc_to_D_V(d,s) p2d_vec(d,s);

/* Conversions from specified formats to prevailing MILC formats */

#define copy_F_G_to_milc(d,s) f2p_mat(d,s);
#define copy_D_G_to_milc(d,s) d2p_mat(d,s);

static void
copy_F_F_to_milc(anti_hermitmat *dest, fmatrix *src){
  matrix t;
  f2p_mat(&t, src);
  make_anti_hermitian( &t, dest );
}

static void
copy_D_F_to_milc(anti_hermitmat *dest, dmatrix *src){
  matrix t;
  d2p_mat(&t, src);
  make_anti_hermitian( &t, dest );
}

#define copy_F_V_to_milc(d,s) f2p_vec(d,s);
#define copy_D_V_to_milc(d,s) d2p_vec(d,s);

void site_coords(int coords[4],site *s){
  coords[0] = s->x;
  coords[1] = s->y;
  coords[2] = s->z;
  coords[3] = s->t;
}

/* Map MILC site links to raw order */

#define make_create_raw4_from_site(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE ** \
create_raw4_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE **raw; \
  MILC_SRCTYPE *tmp; \
  raw = create_raw4_##P##_##T (); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw4_from_site: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = (MILC_SRCTYPE *)F_PT(s, src); \
      copy_milc_to_##P##_##T(raw[dir] + j, &tmp[dir]); \
    } \
  } \
  return raw; \
}

/* Map MILC field links to raw order */

#define make_create_raw4_from_field(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE ** \
create_raw4_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE **raw = NULL; \
  MILC_SRCTYPE *tmp; \
  raw = create_raw4_##P##_##T (); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw4_from_field: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T (coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = src + 4*i; \
      copy_milc_to_##P##_##T(raw[dir] + j, &tmp[dir]); \
    } \
  } \
  return raw; \
}

/* Map MILC site field to raw */

#define make_create_raw_from_site(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE * \
create_raw_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE *raw; \
  MILC_SRCTYPE *tmp; \
  raw = create_raw_##P##_##T(); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw_from_site: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T (coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = (MILC_SRCTYPE *)F_PT(s, src); \
      copy_milc_to_##P##_##T(raw + j, tmp); \
    } \
  } \
  return raw; \
}

/* Map MILC field field to raw */

#define make_create_raw_from_field(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE * \
create_raw_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE *raw; \
  raw = create_raw_##P##_##T(); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw_from_field: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T (coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      copy_milc_to_##P##_##T(raw + j, src + i); \
    } \
  } \
  return raw; \
}

/* Map raw links to MILC site structure */

#define make_unload_raw4_to_site(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw4_##P##_##T##_to_site(field_offset dest, MILC_RAWTYPE *raw[], \
         int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_DSTTYPE *tmp; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw4_to_site: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = (MILC_DSTTYPE *)F_PT(s, dest); \
      copy_##P##_##T##_to_milc(&tmp[dir], raw[dir] + j); \
    } \
  } \
}

#define make_unload_raw4_to_field(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw4_##P##_##T##_to_field(MILC_DSTTYPE *dest, MILC_RAWTYPE *raw[], \
         int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_DSTTYPE *tmp; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw4_to_field: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = dest + 4*i; \
      copy_##P##_##T##_to_milc(&tmp[dir], raw[dir] + j); \
    } \
  } \
}

#define make_unload_raw_to_site(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw_##P##_##T##_to_site(field_offset dest, MILC_RAWTYPE *raw, \
       int milc_parity){ \
  int coords[4]; \
  int i,j; \
  site *s; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw_to_site: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    copy_##P##_##T##_to_milc((MILC_DSTTYPE *)F_PT(s,dest), raw + j); \
  } \
}

#define make_unload_raw_to_field(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw_##P##_##T##_to_field(MILC_DSTTYPE *dest, MILC_RAWTYPE *raw, \
       int milc_parity){ \
  int coords[4]; \
  int i,j; \
  site *s; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw_V_to_field: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    copy_##P##_##T##_to_milc(dest + i, raw + j); \
  } \
}

#define make_load_links_and_mom_site(P, MILCTYPE, MILCFLOAT) \
void \
load_##P##_links_and_mom_site(QOP_##P##3_GaugeField **links, \
   QOP_##P##3_Force **mom, MILCTYPE ***rawlinks, MILCTYPE ***rawmom) \
{ \
  *rawlinks = create_raw4_##P##_G_from_site(F_OFFSET(link), EVENANDODD); \
  if(*rawlinks == NULL)terminate(1); \
  *links = QOP_##P##3_create_G_from_raw((MILCFLOAT **)(*rawlinks),QOP_EVENODD); \
  *rawmom = create_raw4_##P##_F_from_site(F_OFFSET(mom), EVENANDODD); \
  if(*rawmom == NULL)terminate(1); \
  *mom = QOP_##P##3_create_F_from_raw((MILCFLOAT **)(*rawmom),QOP_EVENODD); \
}

#define make_load_links_and_mom_field(P, MILCTYPE, MILCFLOAT) \
void \
load_##P##_links_and_mom_field(QOP_##P##3_GaugeField **links, \
   QOP_##P##3_Force **mom, MILCTYPE ***rawlinks, MILCTYPE ***rawmom, \
   matrix *srclink, anti_hermitmat *srcmom) \
{ \
  *rawlinks = create_raw4_##P##_G_from_field(srclink, EVENANDODD); \
  if(*rawlinks == NULL)terminate(1); \
  *links = QOP_##P##3_create_G_from_raw((MILCFLOAT **)(*rawlinks),QOP_EVENODD); \
  *rawmom = create_raw4_##P##_F_from_field(srcmom, EVENANDODD); \
  if(*rawmom == NULL)terminate(1); \
  *mom = QOP_##P##3_create_F_from_raw((MILCFLOAT **)(*rawmom),QOP_EVENODD); \
}

#define make_load_from_site(P, T, TYPE, MILCTYPE, MILCFLOAT) \
void \
load_##P##_##T##_from_site( TYPE** qop, field_offset src, int parity) \
{ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T##_from_site(src, parity); \
  if(raw == NULL)terminate(1); \
  *qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
					  milc2qop_parity(parity)); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
  return; \
}

/* Map color vector field from MILC field to QOP field */

#define make_load_from_field(P, T, TYPE, MILCTYPE, MILC_SRC_TYPE, MILCFLOAT) \
void \
load_##P##_##T##_from_field( TYPE** qop, MILC_SRC_TYPE *src, int parity) \
{ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T##_from_field(src, parity); \
  if(raw == NULL)terminate(1); \
  *qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
                                        milc2qop_parity(parity)); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
  return; \
}

#define make_unload_links_and_mom_site(P, MILCTYPE, MILCFLOAT) \
void \
unload_##P##_links_and_mom_site(QOP_##P##3_GaugeField **links,  \
   QOP_##P##3_Force **mom, MILCTYPE ***rawlinks, MILCTYPE ***rawmom) \
{ \
  destroy_raw4_##P##_G (*rawlinks);   *rawlinks = NULL; \
  QOP_##P##3_destroy_G (*links);      *links = NULL; \
  QOP_##P##3_extract_F_to_raw((MILCFLOAT **)(*rawmom), *mom, QOP_EVENODD); \
  unload_raw4_##P##_F_to_site(F_OFFSET(mom), *rawmom, EVENANDODD); \
  destroy_raw4_##P##_F (*rawmom);   *rawmom = NULL; \
  QOP_##P##3_destroy_F (*mom);      *mom = NULL; \
}

#define make_unload_links_and_mom_field(P, MILCTYPE, MILCFLOAT) \
void \
unload_##P##_links_and_mom_field(matrix *dstlink, anti_hermitmat *dstmom, \
   QOP_##P##3_GaugeField **links, QOP_##P##3_Force **mom, \
   MILCTYPE ***rawlinks, MILCTYPE ***rawmom) \
{ \
  destroy_raw4_##P##_G (*rawlinks);   *rawlinks = NULL; \
  QOP_##P##3_destroy_G (*links);      *links = NULL; \
  QOP_##P##3_extract_F_to_raw((MILCFLOAT **)(*rawmom), *mom, QOP_EVENODD); \
  unload_raw4_##P##_F_to_field(dstmom, *rawmom, EVENANDODD); \
  destroy_raw4_##P##_F (*rawmom);   *rawmom = NULL; \
  QOP_##P##3_destroy_F (*mom);      *mom = NULL; \
}

/* Map color vector from QOP field to site */

#define make_unload_to_site(P, T, TYPE, MILCTYPE, MILCFLOAT) \
void \
unload_##P##_##T##_to_site( field_offset dest, TYPE *qop, int parity){ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T(); \
  QOP_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qop, milc2qop_parity(parity)); \
  unload_raw_##P##_##T##_to_site(dest, raw, parity); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
}

/* Map color vector from QOP field to MILC field */

#define make_unload_to_field(P, T, TYPE, MILCTYPE, MILC_DSTTYPE, MILCFLOAT) \
void \
unload_##P##_##T##_to_field( MILC_DSTTYPE *dest, TYPE *qop, int parity){ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T(); \
  QOP_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qop, milc2qop_parity(parity)); \
  unload_raw_##P##_##T##_to_field(dest, raw, parity); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
}



/* Storage for raw gauge field */

make_create_raw4(F, G, fmatrix);
make_create_raw4(D, G, dmatrix);

make_destroy_raw4(F, G, fmatrix);
make_destroy_raw4(D, G, dmatrix);

/* Storage for raw gauge momentum */

make_create_raw4(F, F, fmatrix);
make_create_raw4(D, F, dmatrix);

make_destroy_raw4(F, F, fmatrix);
make_destroy_raw4(D, F, dmatrix);

/* Storage for raw su3 vector field */

make_create_raw(F, V, fvector);
make_create_raw(D, V, dvector);

make_destroy_raw(F, V, fvector);
make_destroy_raw(D, V, dvector);

/* Map gauge field from site to raw */

make_create_raw4_from_site(F, G, fmatrix, matrix);
make_create_raw4_from_site(D, G, dmatrix, matrix);

/* Map gauge field from field to raw */

make_create_raw4_from_field(F, G, fmatrix, matrix);
make_create_raw4_from_field(D, G, dmatrix, matrix);

/* Map gauge momentum from site to raw */

make_create_raw4_from_site(F, F, fmatrix, anti_hermitmat);
make_create_raw4_from_site(D, F, dmatrix, anti_hermitmat);

/* Map gauge momentum from field to raw */

make_create_raw4_from_field(F, F, fmatrix, anti_hermitmat);
make_create_raw4_from_field(D, F, dmatrix, anti_hermitmat);

/* Map color vector from site to raw */

make_create_raw_from_site(F, V, fvector, vector);
make_create_raw_from_site(D, V, dvector, vector);

/* Map color vector from field to raw */

make_create_raw_from_field(F, V, fvector, vector);
make_create_raw_from_field(D, V, dvector, vector);

/* Map gauge field from raw to site */

make_unload_raw4_to_site(F, G, matrix, fmatrix);
make_unload_raw4_to_site(D, G, matrix, dmatrix);

/* Map gauge field from raw to field */

make_unload_raw4_to_field(F, G, matrix, fmatrix);
make_unload_raw4_to_field(D, G, matrix, dmatrix);

/* Map gauge momentum from raw to site */

make_unload_raw4_to_site(F, F, anti_hermitmat, fmatrix);
make_unload_raw4_to_site(D, F, anti_hermitmat, dmatrix);

/* Map gauge momentum from raw to field */

make_unload_raw4_to_field(F, F, anti_hermitmat, fmatrix);
make_unload_raw4_to_field(D, F, anti_hermitmat, dmatrix);

/* Map color vector from raw to site */

make_unload_raw_to_site(F, V, vector, fvector);
make_unload_raw_to_site(D, V, vector, dvector);

/* Map color vector from raw to field */

make_unload_raw_to_field(F, V, vector, fvector);
make_unload_raw_to_field(D, V, vector, dvector);

/* Composite mapping */

make_load_links_and_mom_site(F, fmatrix, float);
make_load_links_and_mom_site(D, dmatrix, double);

make_load_links_and_mom_field(F, fmatrix, float);
make_load_links_and_mom_field(D, dmatrix, double);

make_unload_links_and_mom_site(F, fmatrix, float);
make_unload_links_and_mom_site(D, dmatrix, double);

make_unload_links_and_mom_field(F, fmatrix, float);
make_unload_links_and_mom_field(D, dmatrix, double);

make_load_from_site(F, V, QOP_F3_ColorVector, fvector, float);
make_load_from_site(D, V, QOP_D3_ColorVector, dvector, double);

make_load_from_field(F, V, QOP_F3_ColorVector, fvector, vector ,float);
make_load_from_field(D, V, QOP_D3_ColorVector, dvector, vector ,double);

make_unload_to_site(F, V, QOP_F3_ColorVector, fvector, float);
make_unload_to_site(D, V, QOP_D3_ColorVector, dvector, double);

make_unload_to_field(F, V, QOP_F3_ColorVector, fvector, vector, float);
make_unload_to_field(D, V, QOP_D3_ColorVector, dvector, vector, double);

/* map_milc_to_qop.c */


/******************* qop_milc_utilities.c ************************/
/* MILC Version 7 */

/* Functions for mapping QOP-MILC data of a specific precision to flat
   MILC arrays of the prevailing precision specified by the PRECISION
   macro.  Supports the MILC test version of QOP in mixed precision
   calculation */

/* C. DeTar 12/15/2006                                              */

#include "generic_includes.h"
#include <string.h>

/* Convert (or copy) MILC types between specific and prevailing precision */

#if (PRECISION == 1)

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

/********************************************************************/
/* matrix field conversion                                      */
/********************************************************************/

#if ( PRECISION == 1 )

/* No copy necessary if the prevailing precision matches the input array */

matrix *
create_links_from_qop_milc_F(fmatrix *src)
{
  return src;
}

void
destroy_links_from_qop_milc_F(matrix *g){
}

matrix *
create_links_from_qop_milc_D(dmatrix *src)
{
  matrix *g;
  int i, dir;
  site *s;

  g = (matrix *)malloc(sizeof(matrix)*sites_on_node*4);
  if(g == NULL){
    printf("create_links_from_qop_milc_D(%d): No room\n",this_node);
    return NULL;
  }
  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      d2p_mat(g+dir*sites_on_node+i, src+dir*sites_on_node+i);
    }
  }
  return g;
}

void
destroy_links_from_qop_milc_D(matrix *g){
  if(g != NULL) free(g);
}


#else

matrix *
create_links_from_qop_milc_F(fmatrix *src)
{
  matrix *g;
  int i, dir;
  site *s;

  g = (matrix *)malloc(sizeof(matrix)*sites_on_node*4);
  if(g == NULL){
    printf("create_links_from_qop_milc_F(%d): No room\n",this_node);
    return NULL;
  }
  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      f2p_mat(g+dir*sites_on_node+i, src+dir*sites_on_node+i);
    }
  }
  return g;
}

void
destroy_links_from_qop_milc_F(matrix *g){
  if(g != NULL) free(g);
}

matrix *
create_links_from_qop_milc_D(dmatrix *src)
{
  return src;
}

void
destroy_links_from_qop_milc_D(matrix *g){
}


#endif

/********************************************************************/
/* vector field conversion                                      */
/********************************************************************/

#if ( PRECISION == 1 )

vector *
create_latvec_from_qop_milc_F(fvector *src)
{
  return src;
}

void
destroy_latvec_from_qop_milc_F(vector *v){
  return;
}

vector *
create_latvec_from_qop_milc_D(dvector *src)
{
  vector *v;
  int i;
  site *s;

  v = (vector *)malloc(sizeof(vector)*sites_on_node);
  if(v == NULL){
    printf("create_latvec_from_qop_milc_D(%d): No room\n",this_node);
    return NULL;
  }
  FORALLSITES(i,s){
    d2p_vec(v + i, src + i);
  }
  return v;
}

void
destroy_latvec_from_qop_milc_D(vector *v){
  if(v != NULL) free(v);
}


#else

vector *
create_latvec_from_qop_milc_F(fvector *src)
{
  vector *v;
  int i;
  site *s;

  v = (vector *)malloc(sizeof(vector)*sites_on_node);
  if(v == NULL){
    printf("create_latvec_from_qop_milc_F(%d): No room\n",this_node);
    return NULL;
  }
  FORALLSITES(i,s){
    f2p_vec(v + i, src + i);
  }
  return v;
}

void
destroy_latvec_from_qop_milc_F(vector *v){
  if(v != NULL) free(v);
}

vector *
create_latvec_from_qop_milc_D(dvector *src)
{
  return src;
}

void
destroy_latvec_from_qop_milc_D(vector *v){
  return;
}

#endif

/********************************************************************/
/* vector field copy from prevailing MILC precision to specific
   QOP precision */
/********************************************************************/

#if ( PRECISION == 1 )
void
copy_latvec_to_qop_milc_F( fvector *dest, vector *src)
{
  if(dest != src)
    memcpy(dest, src, sites_on_node*sizeof(fvector));
}

void
copy_latvec_to_qop_milc_D( dvector *dest, vector *src)
{
  int i;
  site *s;

  FORALLSITES(i,s){
    p2d_vec(dest + i, src + i);
  }
}

#else

void
copy_latvec_to_qop_milc_F(fvector *dest, vector *src)
{
  int i;
  site *s;

  FORALLSITES(i,s){
    p2f_vec(dest + i, src + i);
  }
}

void
copy_latvec_to_qop_milc_D(dvector *dest, vector *src)
{
  if(dest != src)
    memcpy(dest, src, sites_on_node*sizeof(dvector));
}

#endif

/* qop_milc_utilities */

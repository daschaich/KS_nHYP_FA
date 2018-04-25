// -----------------------------------------------------------------
// Defines and subroutine declarations for SU(3)
#ifndef _SU3_H
#define _SU3_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge group is SU(3) with three-component fundamental fermions
typedef struct { fcomplex e[3][3]; } fmatrix;
typedef struct { fcomplex c[3]; } fvector;
typedef struct {
  fcomplex m01, m02, m12;
  float m00im, m11im, m22im;
  float space;
} fanti_hermitmat;

typedef struct { dcomplex e[3][3]; } dmatrix;
typedef struct { dcomplex c[3]; } dvector;
typedef struct {
  dcomplex m01, m02, m12;
  double m00im, m11im, m22im;
  double space;
} danti_hermitmat;

#if PRECISION == 1
#define matrix          fmatrix
#define vector          fvector
#define anti_hermitmat  fanti_hermitmat
#else
#define matrix          dmatrix
#define vector          dvector
#define anti_hermitmat  danti_hermitmat
#endif

/* For KS spectroscopy */
typedef vector **ks_prop_field;

/* SU(2) */
typedef struct { complex e[2][2]; } su2_matrix;

#define GAMMAFIVE -1    // Some integer which is not a direction

// Flags for selecting M or Mdag
#define PLUS 1
#define MINUS -1

/*
* ROUTINES FOR SU(3) MATRIX OPERATIONS
*
* void mult_nn(matrix *a, matrix *b, matrix *c)
* matrix multiply, no adjoints
* files "m_mat_nn.c"
* void mult_na(matrix *a, matrix *b, matrix *c)
* matrix multiply, second matrix is adjoint
* files "m_mat_na.c"
* void mult_an(matrix *a, matrix *b, matrix *c)
* matrix multiply, first matrix is adjoint
* files "m_mat_an.c"
* complex complextrace(matrix *a, matrix *b)
* (Tr(A_adjoint*B))
* file "complextr.c"
* complex det_su3(matrix *a)
* file "det_su3.c"
* void add_matrix(matrix *a, matrix *b, matrix *c)
* file "addmat.c"
* void sub_matrix(matrix *a, matrix *b, matrix *c)
* file "submat.c"
* void scalar_mult_matrix(matrix *a, Real s, matrix *b)
* file "s_m_mat.c"
* void scalar_mult_add_matrix(matrix *a, matrix *b,
* Real s, matrix *c)
* file "s_m_a_mat.c"
* void scalar_mult_sub_matrix(matrix *a, matrix *b,
* Real s, matrix *c)
* file "s_m_s_mat.c"
* void c_scalar_mult_mat(matrix *src, complex *phase, matrix *dest)
* file "cs_m_mat.c"
* void c_scalar_mult_add_mat(matrix *m1, matrix *m2,
* complex *phase, matrix *m3)
* file "cs_m_a_mat.c"
* void adjoint(matrix *a, matrix *b)
* file "adjoint.c"
* void make_anti_hermitian(matrix *m3,  anti_hermitmat *ah3)
* file "make_ahmat.c"
* void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt)
* (prn_pt passed through to myrand())
* file "rand_ahmat.c"
* void uncompress_anti_hermitian(anti_hermitmat *mat_anti, matrix *mat)
* file "uncmp_ahmat.c"
* void compress_anti_hermitian(matrix *mat, anti_hermitmat *mat_anti)
* file "cmp_ahmat.c"
* void clear_mat(matrix *dest);
*       file clear_mat.c
*          dest <- 0.0
* void mat_copy(matrix *a, matrix *b)
* file "mat_copy.c"
* void dumpmat(matrix *m)
*       file "dumpmat.c"
*
* ROUTINES FOR vector OPERATIONS (3 COMPONENT COMPLEX)
*
* void su3_projector(vector *a, vector *b, matrix *c)
* (outer product of A and B)
* file "su3_proj.c"
  *
* complex su3_dot(vector *a, vector *b)
* file "su3_dot.c"
  *
* Real su3_rdot(vector *a, vector *b)
* file "su3_rdot.c"
  *
* Real magsq_vec(vector *a)
* file "msq_vec.c"
  *
* void vec_copy(vector *a, vector *b)
* file "vec_copy.c"
*
* void mult_mat_vec(matrix *a, vector *b, vector *c)
*  C  <-  A*B
* file "m_matvec.c"
  *
* void mult_mat_vec_sum(matrix *a, vector *b, vector *c)
*  C  <-  C + A*B
* file "m_matvec_s.c"
  *
* void mult_mat_vec_sum_4dir(matrix *a, vector *b0,
* vector *b1, vector *b2, vector *b3, vector *c)
* file "m_mv_s_4dir.c"
* Multiply four vectors by elements of an array of matrices,
* sum results.
  *
* void mult_adj_mat_vec(matrix *a, vector *b, vector *c)
* file "m_amatvec.c"
  *
* void mult_adj_mat_vec_4dir(matrix *a, vector *b, vector *c)
* Same as above, but result vectors need not be in an array.
* file "m_amv_4dir.c"
*
* void add_vector(vector *a, vector *b, vector *c)
* file "addvec.c"
* void sub_four_vecs(vector *a, vector *b1, vector *b2,
*   vector *b3, vector *b4)
* file "sub4vecs.c"
*
* void scalar_mult_vector(vector *a, Real s, vector *c)
* file "s_m_vec.c"
* void scalar_mult_add_vector(vector *a, vector *b, Real s,
* vector *c)
* file "s_m_a_vec.c"
* void scalar_mult_sum_vector(vector *a, vector *b, Real s)
* file "s_m_sum_vec.c"
* void scalar_mult_sub_vector(vector *a, vector *b, Real s,
* vector *c)
* file "s_m_s_vec.c"
* void c_scalar_mult_vec(vector *src, complex *phase, vector *dest)
* file "cs_m_vec.c"
* void c_scalar_mult_add_vec(vector *v1, complex *phase, vector *v2)
* file "cs_m_a_vec.c"
* void c_scalar_mult_sub_vec(vector *v1, complex *phase, vector *v2)
* file "cs_m_s_vec.c"
* void dumpvec(vector *v)
*       file "dumpvec.c"
* void clearvec(vector *v)
*       file "clearvec.c"
*
* ROUTINES MIXING SU(2) and SU(3)
*
* void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link)
*       file "l_su2_hit_n.c"
* void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link)
*       file "r_su2_hit_a.c"
* void dumpsu2(su2_matrix *u)
*       file "dumpsu2.c"
* void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1)
*       file "m_su2_mat_vec_n.c"
* void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);
*       file "m_su2_mat_vec_a.c"
*
* MISCELLANEOUS ROUTINES
*
* Real gaussian_rand_no(double_prn *prn_pt)
* file "gaussrand.c"
* Real z2_rand_no(double_prn *prn_pt);
* file "z2rand.c"
* void byterevn(int32type w[], int n)
* void byterevn64(int32type w[], int n)
*/
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Inline macros

/* We have optional SSE and C inline macros for selected library
   routines */

/* All macros are defined and available for selective inlining.  To
   invoke them selectively, define the macro SSE_INLINE and then
   modify the code by replacing the function call with the appropriate
   macro */

/* Inlining can also be implemented globally without requiring any
   changes in the code.  SSE macros are inlined globally by defining
   SSE_GLOBAL_INLINE.  */

/* SSE macros for selected library functions are available in single
   and double precision */

#if defined SSE_INLINE || defined SSE_GLOBAL_INLINE
#if PRECISION == 2
#include "../sse2/include/inline_sse.h"
#endif
#endif

/* The following definitions cause the macros to be invoked globally
   if SSE_INLINE is defined.  Note that in each stanza
   the same library routines are repeated, but with either a macro
   definition when available or a prototype definition. */

#if defined SSE_GLOBAL_INLINE

// Our available double-precision SSE macros */

#define mult_nn(...) _inline_sse_mult_nn(__VA_ARGS__)
#define mult_na(...) _inline_sse_mult_na(__VA_ARGS__)
#define mult_an(...) _inline_sse_mult_an(__VA_ARGS__)
#define mult_mat_vec(...) _inline_sse_mult_mat_vec(__VA_ARGS__)
#define mult_adj_mat_vec(...) _inline_sse_mult_adj_mat_vec(__VA_ARGS__)
#define mult_adj_mat_vec_4dir(...) _inline_sse_mult_adj_mat_vec_4dir(__VA_ARGS__)
#define mult_mat_vec_sum_4dir(...) _inline_sse_mult_mat_vec_sum_4dir(__VA_ARGS__)
#define su3_projector(...) _inline_sse_su3_projector(__VA_ARGS__)

/********************************************************************/
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Matrix operations
// In file realtr.c
Real realtrace_nn(matrix *a, matrix *b);
Real realtrace(matrix *a, matrix *b);

// In file trace.c
complex trace(matrix *a);

complex complextrace(matrix *a, matrix *b);
complex det_su3(matrix *a);

// In file addmat.c
void sum_matrix(matrix *b, matrix *c);
void add_matrix(matrix *a, matrix *b, matrix *c);

// In file submat.c
void sub_matrix(matrix *a, matrix *b, matrix *c);

void scalar_mult_matrix(matrix *b, Real s, matrix *c);
void scalar_mult_sub_matrix(matrix *a, matrix *b, Real s, matrix *c);

void c_scalar_mult_mat(matrix *b, complex *s, matrix *c);
void scalar_add_diag(matrix *a, Real s);
void c_scalar_add_diag(matrix *a, complex *s);

// In file s_m_a_mat.c
void scalar_mult_sum_matrix(matrix *b, Real s, matrix *c);
void scalar_mult_add_matrix(matrix *a, matrix *b, Real scalar, matrix *c);

// In file cs_m_a_mat.c
void c_scalar_mult_sum_mat(matrix *b, complex *s, matrix *c);
void c_scalar_mult_add_mat(matrix *a, matrix *b, complex *s, matrix *c);

void adjoint(matrix *a, matrix *b);

void make_anti_hermitian(matrix *m3, anti_hermitmat *ah3);
void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt);
void uncompress_anti_hermitian(anti_hermitmat *mat_anti, matrix *mat);
void compress_anti_hermitian(matrix *mat, anti_hermitmat *mat_anti);

void clear_mat(matrix *dest);
void mat_copy(matrix *a, matrix *b);
void dumpmat(matrix *m);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Vector operations
complex su3_dot(vector *a, vector *b);
void vec_copy(vector *a, vector *b);
void dumpvec(vector *v);
void clearvec(vector *v);

void mult_mat_vec_sum(matrix *a, vector *b, vector *c);

// In file addvec.c
void sum_vector(vector *b, vector *c);
void add_vector(vector *a, vector *b, vector *c);

// In file subvec.c
void dif_vector(vector *b, vector *c);
void sub_vector(vector *a, vector *b, vector *c);

void scalar_mult_vector(vector *b, Real s, vector *c);
void scalar_mult_sub_vector(vector *a, vector *b, Real s, vector *c);
void c_scalar_mult_vec(vector *b, complex *s, vector *c);
void c_scalar_mult_sum_vec(vector *b, complex *s, vector *c);
void c_scalar_mult_dif_vec(vector *b, complex *s, vector *c);

void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link);
void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link);
void dumpsu2(su2_matrix *u);
void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1);
void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);

Real gaussian_rand_no(double_prn *prn_pt);

Real z2_rand_no(double_prn *prn_pt);

#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

Real magsq_vec(vector *a);

// In file m_mat_nn.c
#ifndef mult_nn
void mult_nn(matrix *a, matrix *b, matrix *c);
#endif

// In file m_mat_na.c
void mult_na_sum(matrix *a, matrix *b, matrix *c);
#ifndef mult_na
void mult_na(matrix *a, matrix *b, matrix *c);
#endif

// In file m_mat_an.c
void mult_an_sum(matrix *a, matrix *b, matrix *c);
#ifndef mult_an
void mult_an(matrix *a, matrix *b, matrix *c);
#endif

#ifndef mult_mat_vec
void mult_mat_vec(matrix *a, vector *b, vector *c);
#endif

#ifndef mult_adj_mat_vec
void mult_adj_mat_vec(matrix *a, vector *b, vector *c);
#endif

#ifndef mult_adj_mat_vec_4dir
void mult_adj_mat_vec_4dir(matrix *a, vector *b, vector *c);
#endif

#ifndef mult_mat_vec_sum_4dir
void mult_mat_vec_sum_4dir(matrix *a, vector *b0, vector *b1,
                           vector *b2, vector *b3, vector *c);
#endif

// In file s_m_a_vec.c
void scalar_mult_sum_vector(vector *b, Real s, vector *c);
void scalar_mult_add_vector(vector *a, vector *b, Real s, vector *c);

#ifndef su3_projector
void su3_projector(vector *a, vector *b, matrix *c);
#endif

Real su3_rdot(vector *a, vector *b);

void sub_four_vecs(vector *a, vector *b1, vector *b2, vector *b3, vector *b4);

#endif
// -----------------------------------------------------------------

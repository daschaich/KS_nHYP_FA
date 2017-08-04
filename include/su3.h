// -----------------------------------------------------------------
// Defines and subroutine declarations for SU(3)
#ifndef _SU3_H
#define _SU3_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge group is SU(3) with three-component fundamental fermions
typedef struct { fcomplex e[3][3]; } fsu3_matrix;
typedef struct { fcomplex c[3]; } fsu3_vector;
typedef struct {
  fcomplex m01, m02, m12;
  float m00im, m11im, m22im;
  float space;
} fanti_hermitmat;

typedef struct { dcomplex e[3][3]; } dsu3_matrix;
typedef struct { dcomplex c[3]; } dsu3_vector;
typedef struct {
  dcomplex m01, m02, m12;
  double m00im, m11im, m22im;
  double space;
} danti_hermitmat;

#if PRECISION == 1
#define su3_matrix      fsu3_matrix
#define su3_vector      fsu3_vector
#define anti_hermitmat  fanti_hermitmat
#else
#define su3_matrix      dsu3_matrix
#define su3_vector      dsu3_vector
#define anti_hermitmat  danti_hermitmat
#endif

/* For KS spectroscopy */
typedef su3_vector **ks_prop_field;

/* SU(2) */
typedef struct { complex e[2][2]; } su2_matrix;

#define GAMMAFIVE -1    // Some integer which is not a direction

// Flags for selecting M or Mdag
#define PLUS 1
#define MINUS -1

/*
* ROUTINES FOR SU(3) MATRIX OPERATIONS
*
* void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* matrix multiply, no adjoints
* files "m_mat_nn.c"
* void mult_su3_na(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* matrix multiply, second matrix is adjoint
* files "m_mat_na.c"
* void mult_su3_an(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* matrix multiply, first matrix is adjoint
* files "m_mat_an.c"
* complex complextrace_su3(su3_matrix *a, su3_matrix *b)
* (Tr(A_adjoint*B))
* file "complextr.c"
* complex det_su3(su3_matrix *a)
* file "det_su3.c"
* void add_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* file "addmat.c"
* void sub_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* file "submat.c"
* void scalar_mult_su3_matrix(su3_matrix *a, Real s, su3_matrix *b)
* file "s_m_mat.c"
* void scalar_mult_add_su3_matrix(su3_matrix *a, su3_matrix *b,
* Real s, su3_matrix *c)
* file "s_m_a_mat.c"
* void scalar_mult_sub_su3_matrix(su3_matrix *a, su3_matrix *b,
* Real s, su3_matrix *c)
* file "s_m_s_mat.c"
* void c_scalar_mult_su3mat(su3_matrix *src, complex *phase, su3_matrix *dest)
* file "cs_m_mat.c"
* void c_scalar_mult_add_su3mat(su3_matrix *m1, su3_matrix *m2,
* complex *phase, su3_matrix *m3)
* file "cs_m_a_mat.c"
* void su3_adjoint(su3_matrix *a, su3_matrix *b)
* file "su3_adjoint.c"
* void make_anti_hermitian(su3_matrix *m3,  anti_hermitmat *ah3)
* file "make_ahmat.c"
* void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt)
* (prn_pt passed through to myrand())
* file "rand_ahmat.c"
* void uncompress_anti_hermitian(anti_hermitmat *mat_anti, su3_matrix *mat)
* file "uncmp_ahmat.c"
* void compress_anti_hermitian(su3_matrix *mat, anti_hermitmat *mat_anti)
* file "cmp_ahmat.c"
* void clear_su3mat(su3_matrix *dest);
*       file clear_mat.c
*          dest <- 0.0
* void su3mat_copy(su3_matrix *a, su3_matrix *b)
* file "su3mat_copy.c"
* void dumpmat(su3_matrix *m)
*       file "dumpmat.c"
*
* ROUTINES FOR su3_vector OPERATIONS (3 COMPONENT COMPLEX)
*
* void su3_projector(su3_vector *a, su3_vector *b, su3_matrix *c)
* (outer product of A and B)
* file "su3_proj.c"
  *
* complex su3_dot(su3_vector *a, su3_vector *b)
* file "su3_dot.c"
  *
* Real su3_rdot(su3_vector *a, su3_vector *b)
* file "su3_rdot.c"
  *
* Real magsq_su3vec(su3_vector *a)
* file "msq_su3vec.c"
  *
* void su3vec_copy(su3_vector *a, su3_vector *b)
* file "su3vec_copy.c"
*
* void mult_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c)
*  C  <-  A*B
* file "m_matvec.c"
  *
* void mult_su3_mat_vec_sum(su3_matrix *a, su3_vector *b, su3_vector *c)
*  C  <-  C + A*B
* file "m_matvec_s.c"
  *
* void mult_su3_mat_vec_sum_4dir(su3_matrix *a, su3_vector *b0,
* su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c)
* file "m_mv_s_4dir.c"
* Multiply four su3_vectors by elements of an array of su3_matrices,
* sum results.
  *
* void mult_adj_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c)
* file "m_amatvec.c"
  *
* void mult_adj_su3_mat_vec_4dir(su3_matrix *a, su3_vector *b, su3_vector *c)
* Same as above, but result vectors need not be in an array.
* file "m_amv_4dir.c"
*
* void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c)
* file "addvec.c"
* void sub_four_su3_vecs(su3_vector *a, su3_vector *b1, su3_vector *b2,
*   su3_vector *b3, su3_vector *b4)
* file "sub4vecs.c"
*
* void scalar_mult_su3_vector(su3_vector *a, Real s, su3_vector *c)
* file "s_m_vec.c"
* void scalar_mult_add_su3_vector(su3_vector *a, su3_vector *b, Real s,
* su3_vector *c)
* file "s_m_a_vec.c"
* void scalar_mult_sum_su3_vector(su3_vector *a, su3_vector *b, Real s)
* file "s_m_sum_vec.c"
* void scalar_mult_sub_su3_vector(su3_vector *a, su3_vector *b, Real s,
* su3_vector *c)
* file "s_m_s_vec.c"
* void c_scalar_mult_su3vec(su3_vector *src, complex *phase, su3_vector *dest)
* file "cs_m_vec.c"
* void c_scalar_mult_add_su3vec(su3_vector *v1, complex *phase, su3_vector *v2)
* file "cs_m_a_vec.c"
* void c_scalar_mult_sub_su3vec(su3_vector *v1, complex *phase, su3_vector *v2)
* file "cs_m_s_vec.c"
* void dumpvec(su3_vector *v)
*       file "dumpvec.c"
* void clearvec(su3_vector *v)
*       file "clearvec.c"
*
* ROUTINES MIXING SU(2) and SU(3)
*
* void left_su2_hit_n(su2_matrix *u, int p, int q, su3_matrix *link)
*       file "l_su2_hit_n.c"
* void right_su2_hit_a(su2_matrix *u, int p, int q, su3_matrix *link)
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

#define mult_su3_nn(...) _inline_sse_mult_su3_nn(__VA_ARGS__)
#define mult_su3_na(...) _inline_sse_mult_su3_na(__VA_ARGS__)
#define mult_su3_an(...) _inline_sse_mult_su3_an(__VA_ARGS__)
#define mult_su3_mat_vec(...) _inline_sse_mult_su3_mat_vec(__VA_ARGS__)
#define mult_adj_su3_mat_vec(...) _inline_sse_mult_adj_su3_mat_vec(__VA_ARGS__)
#define mult_adj_su3_mat_vec_4dir(...) _inline_sse_mult_adj_su3_mat_vec_4dir(__VA_ARGS__)
#define mult_su3_mat_vec_sum_4dir(...) _inline_sse_mult_su3_mat_vec_sum_4dir(__VA_ARGS__)
#define su3_projector(...) _inline_sse_su3_projector(__VA_ARGS__)

/********************************************************************/
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Matrix operations
// In file realtr.c
Real realtrace_su3_nn(su3_matrix *a, su3_matrix *b);
Real realtrace_su3(su3_matrix *a, su3_matrix *b);

// In file trace_su3.c
complex trace_su3(su3_matrix *a);

complex complextrace_su3(su3_matrix *a, su3_matrix *b);
complex det_su3(su3_matrix *a);

// In file addmat.c
void sum_su3_matrix(su3_matrix *b, su3_matrix *c);
void add_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file submat.c
void sub_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);

void scalar_mult_su3_matrix(su3_matrix *b, Real s, su3_matrix *c);
void scalar_mult_sub_su3_matrix(su3_matrix *a, su3_matrix *b, Real s,
                                su3_matrix *c);
void c_scalar_mult_su3mat(su3_matrix *b, complex *s, su3_matrix *c);
void scalar_add_diag_su3(su3_matrix *a, Real s);
void c_scalar_add_diag_su3(su3_matrix *a, complex *s);

// In file cs_m_a_mat.c
void c_scalar_mult_sum_su3mat(su3_matrix *b, complex *s, su3_matrix *c);
void c_scalar_mult_add_su3mat(su3_matrix *a, su3_matrix *b, complex *s,
                              su3_matrix *c);

void su3_adjoint(su3_matrix *a, su3_matrix *b);

void make_anti_hermitian(su3_matrix *m3, anti_hermitmat *ah3);
void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt);
void uncompress_anti_hermitian(anti_hermitmat *mat_anti, su3_matrix *mat);
void compress_anti_hermitian(su3_matrix *mat, anti_hermitmat *mat_anti);

void clear_su3mat(su3_matrix *dest);
void su3mat_copy(su3_matrix *a, su3_matrix *b);
void dumpmat(su3_matrix *m);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Vector operations
complex su3_dot(su3_vector *a, su3_vector *b);
void su3vec_copy(su3_vector *a, su3_vector *b);
void dumpvec(su3_vector *v);
void clearvec(su3_vector *v);

void mult_su3_mat_vec_sum(su3_matrix *a, su3_vector *b, su3_vector *c);

// In file addvec.c
void sum_su3_vector(su3_vector *b, su3_vector *c);
void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c);

// In file subvec.c
void dif_su3_vector(su3_vector *b, su3_vector *c);
void sub_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c);

void scalar_mult_su3_vector(su3_vector *b, Real s, su3_vector *c);
void scalar_mult_sub_su3_vector(su3_vector *a, su3_vector *b, Real s,
                                su3_vector *c);
void c_scalar_mult_su3vec(su3_vector *b, complex *s, su3_vector *c);
void c_scalar_mult_sum_su3vec(su3_vector *b, complex *s, su3_vector *c);
void c_scalar_mult_dif_su3vec(su3_vector *b, complex *s, su3_vector *c);

void left_su2_hit_n(su2_matrix *u, int p, int q, su3_matrix *link);
void right_su2_hit_a(su2_matrix *u, int p, int q, su3_matrix *link);
void dumpsu2(su2_matrix *u);
void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1);
void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);

Real gaussian_rand_no(double_prn *prn_pt);

Real z2_rand_no(double_prn *prn_pt);

#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

Real magsq_su3vec(su3_vector *a);

// In file m_mat_nn.c
#ifndef mult_su3_nn
void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c);
#endif

// In file m_mat_na.c
void mult_su3_na_sum(su3_matrix *a, su3_matrix *b, su3_matrix *c);
#ifndef mult_su3_na
void mult_su3_na(su3_matrix *a, su3_matrix *b, su3_matrix *c);
#endif

// In file m_mat_an.c
void mult_su3_an_sum(su3_matrix *a, su3_matrix *b, su3_matrix *c);
#ifndef mult_su3_an
void mult_su3_an(su3_matrix *a, su3_matrix *b, su3_matrix *c);
#endif

#ifndef mult_su3_mat_vec
void mult_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);
#endif

#ifndef mult_adj_su3_mat_vec
void mult_adj_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);
#endif

#ifndef mult_adj_su3_mat_vec_4dir
void mult_adj_su3_mat_vec_4dir(su3_matrix *a, su3_vector *b, su3_vector *c);
#endif

#ifndef mult_su3_mat_vec_sum_4dir
void mult_su3_mat_vec_sum_4dir(su3_matrix *a, su3_vector *b0,
  su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c);
#endif

// c <-- c + s * b, in "s_m_a_mat.c"
void scalar_mult_sum_su3_matrix(su3_matrix *b, Real s, su3_matrix *c);

// c <-- a + s * b, in "s_m_a_mat.c"
void scalar_mult_add_su3_matrix(su3_matrix *src1, su3_matrix *src2,
                                Real scalar, su3_matrix *dest);

void scalar_mult_sum_su3_vector(su3_vector *dest, su3_vector *src, Real s);
void scalar_mult_add_su3_vector(su3_vector *src1, su3_vector *src2,
                                Real s, su3_vector *dest);

#ifndef su3_projector
void su3_projector(su3_vector *a, su3_vector *b, su3_matrix *c);
#endif

Real su3_rdot(su3_vector *a, su3_vector *b);

void sub_four_su3_vecs(su3_vector *a, su3_vector *b1, su3_vector *b2,
                       su3_vector *b3, su3_vector *b4);

#endif
// -----------------------------------------------------------------

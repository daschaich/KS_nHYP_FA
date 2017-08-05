// -----------------------------------------------------------------
// Conjugate gradient for staggered fermions
// Traditional algorithm based on v6/generic_ks/d_congrad5.c
// Solves (M^dag.M) psi = chi
// Like v6 code, includes only dslash!

// The dslash code in this version overlaps computation and gathers
// from negative directions, and has an extra lattice loop devoted to
// exclusively to sub_four_vectors.  I stripped the dslash timing

// This version looks at the initial vector every "niter" passes
// "chi" is the source vector
// "psi" is the initial guess and answer
// "r" is the residual vector
// "p" and "mp" are working vectors for the conjugate gradient
// niter = maximum number of iterations
// rsqmin = desired rsq, quit when we reach rsq <= rsqmin * source_norm

// The source is obtained from a random vector with
// average squared magnitude 3 on each site
// Then, on half the sites, we gather and sum the
// eight neighboring random vectors and add 2 * m times the local vector

// It looks like the routine never restarts
// The restart criterion is five times the maximum number of iterations

// parity = EVEN:       do only even sites
// parity = ODD:        do only odd sites
// parity = EVENANDODD: do all sites

// Includes and definitions
#include "eig_includes.h"
#include "../generic_ks/generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND   // For loopend.h
#include "../include/loopend.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Dslash routine
// Set psi on each site to sum of sources parallel transported to site,
// with minus sign for transport from negative directions
void dslash(field_offset chi, field_offset psi, int parity) {
  register int i, dir, otherparity = EVEN;
  register site *s;
  msg_tag *tag[8];
#ifdef INLINE
  register vector *a, *b1, *b2, *b3, *b4;
#endif

  switch(parity) {
    case EVEN:
      otherparity = ODD;
      break;
    case ODD:
      otherparity = EVEN;
      break;
    case EVENANDODD:
      otherparity = EVENANDODD;
      break;
  }

  // Start gathers from positive directions
  for (dir = XUP; dir <= TUP; dir++)
      tag[dir] = start_gather_site(chi, sizeof(vector), dir, parity,
                                   gen_pt[dir]);

  // Multiply by adjoint matrix at other sites
  FORSOMEPARITY(i, s, otherparity) {
    if (i < loopend - FETCH_UP)
      prefetch_4MV4V(&((s + FETCH_UP)->link[XUP]),
                     (vector *)F_PT((s + FETCH_UP), chi),
                     (s+FETCH_UP)->tempvec);

    mult_adj_mat_vec_4dir(s->link, (vector *)F_PT(s, chi),
                              s->tempvec);
  } END_LOOP

  // Start gathers from negative directions
  for (dir = XUP; dir <= TUP; dir++)
    tag[OPP_DIR(dir)] = start_gather_site(F_OFFSET(tempvec[dir]),
                                          sizeof(vector), OPP_DIR(dir),
                                          parity, gen_pt[OPP_DIR(dir)]);

  // Wait gathers from positive directions
  for (dir = XUP; dir <= TUP; dir++)
    wait_gather(tag[dir]);

  // Multiply by matrix and accumulate
  FORSOMEPARITY(i, s, parity) {
    if (i < loopend - FETCH_UP) {
      prefetch_V((vector *)F_PT(s + FETCH_UP, psi));
      prefetch_4MVVVV(&((s + FETCH_UP)->link[XUP]),
                      (vector *)gen_pt[XUP][i + FETCH_UP],
                      (vector *)gen_pt[YUP][i + FETCH_UP],
                      (vector *)gen_pt[ZUP][i + FETCH_UP],
                      (vector *)gen_pt[TUP][i + FETCH_UP]);
    }
    mult_mat_vec_sum_4dir(s->link, (vector *)gen_pt[XUP][i],
                              (vector *)gen_pt[YUP][i],
                              (vector *)gen_pt[ZUP][i],
                              (vector *)gen_pt[TUP][i],
                              (vector *)F_PT(s,psi));
  } END_LOOP

  // Wait gathers from negative directions
  for (dir = XUP; dir <= TUP; dir++)
    wait_gather(tag[OPP_DIR(dir)]);

  // Accumulate (negative)
  FORSOMEPARITY(i, s, parity) {
    if (i < loopend - FETCH_UP)
      prefetch_VVVV((vector *)gen_pt[XDOWN][i + FETCH_UP],
                    (vector *)gen_pt[YDOWN][i + FETCH_UP],
                    (vector *)gen_pt[ZDOWN][i + FETCH_UP],
                    (vector *)gen_pt[TDOWN][i + FETCH_UP] );

#ifndef INLINE
    // Non-inline version
    sub_four_vecs((vector *)F_PT(s,psi),
                      (vector *)(gen_pt[XDOWN][i]),
                      (vector *)(gen_pt[YDOWN][i]),
                      (vector *)(gen_pt[ZDOWN][i]),
                      (vector *)(gen_pt[TDOWN][i]));
#else
    // Inline version
    a  = (vector *)F_PT(s,psi);
    b1 = (vector *)(gen_pt[XDOWN][i]);
    b2 = (vector *)(gen_pt[YDOWN][i]);
    b3 = (vector *)(gen_pt[ZDOWN][i]);
    b4 = (vector *)(gen_pt[TDOWN][i]);

    CSUB(a->c[0], b1->c[0], a->c[0]);
    CSUB(a->c[1], b1->c[1], a->c[1]);
    CSUB(a->c[2], b1->c[2], a->c[2]);

    CSUB(a->c[0], b2->c[0], a->c[0]);
    CSUB(a->c[1], b2->c[1], a->c[1]);
    CSUB(a->c[2], b2->c[2], a->c[2]);

    CSUB(a->c[0], b3->c[0], a->c[0]);
    CSUB(a->c[1], b3->c[1], a->c[1]);
    CSUB(a->c[2], b3->c[2], a->c[2]);

    CSUB(a->c[0], b4->c[0], a->c[0]);
    CSUB(a->c[1], b4->c[1], a->c[1]);
    CSUB(a->c[2], b4->c[2], a->c[2]);
#endif
  } END_LOOP

  // Free buffers
  for (dir = XUP; dir <= TUP; dir++) {
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[OPP_DIR(dir)]);
  }
}
// -----------------------------------------------------------------

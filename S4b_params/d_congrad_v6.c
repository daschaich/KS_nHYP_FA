// -----------------------------------------------------------------
// Conjugate gradient for staggered fermions
// Uses all s->tempvec components for temporary storage
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

// parity = EVEN:       do only even sites
// parity = ODD:        do only odd sites
// parity = EVENANDODD: do all sites
#include "S4b_includes.h"
#include "../generic_ks/generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND   // For loopend.h
#include "../include/loopend.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Stripped final_rsq and global variables from argument list
int ks_congrad(field_offset chi, field_offset psi, Real m, int parity) {
  register int i;
  register site *s;
  int iteration = 0;              // Counter for iterations
  int restarts = 0;               // Counter for restarts
  Real a, b;                      // Alpha and beta
  double rsq;                     // r^2
  double oldrsq, pkp;             // Previous r^2, pkp = p.K.p
  Real msq_x4 = 4.0 * m * m;
  double source_norm;             // Squared magnitude of source vector
  double rsqstop;                 // To be normalized by source norm
  int l_parity = EVEN;            // Local parity, reset for EVENANDODD
  int l_otherparity = ODD;        // The parity we are not doing
  msg_tag *tags1[8], *tags2[8];   // Tags for gathers to parity and opposite
  int special_started = 0;        // 1 if dslash_special has been called
  double nflop = 606;

#ifdef CGTIME
  double dtimec = -dclock();
#endif

  if (parity == EVENANDODD)
    nflop *= 2;

  // If we want both parities, we do even first
  switch(parity) {
    case(EVEN):
      l_parity = EVEN;
      l_otherparity = ODD;
      break;
    case(ODD):
      l_parity = ODD;
      l_otherparity = EVEN;
      break;
    case(EVENANDODD):
      l_parity = EVEN;
      l_otherparity = ODD;
      break;
  }

  // Initialization
start:
#ifdef CG_DEBUG
  node0_printf("CONGRAD: start, parity = %d\n", parity);
  if (restarts > 0)
    node0_printf("CONGRAD: pre-restart rsq/src = %.8g\n", rsq / source_norm);
#endif
  // mp <- -(M^dag.M) psi
  // (Negative sign leads to adds instead of subtracts,
  //  compared to clover code)
  // r <- chi + mp
  // p <- chi + mp
  // rsq = |r|^2
  // source_norm = |chi|^2

  if (special_started == 1) {    // Clean up gathers
    cleanup_gathers(tags1, tags2);
    special_started = 0;
  }
  rsq = 0;
  source_norm = 0;
  dslash(psi, F_OFFSET(mp), l_otherparity);
  dslash(F_OFFSET(mp), F_OFFSET(mp), l_parity);

  // mp  <- mp - msq_x4 * psi
  FORSOMEPARITY(i, s, l_parity) {
    scalar_mult_add_vector(&(s->mp), (vector *)F_PT(s, psi),
                               -msq_x4, &(s->mp));

    // mp contains -(M^dag M) psi
    add_vector((vector *)F_PT(s, chi), &(s->mp), &(s->r));

    s->p = s->r;
    rsq += (double)magsq_vec(&(s->r));
    source_norm += (double)magsq_vec((vector *)F_PT(s, chi));
  } END_LOOP
  g_doublesum(&source_norm);
  g_doublesum(&rsq);
#ifdef CG_DEBUG
  node0_printf("CONGRAD: source_norm = %.8g\n", source_norm);
  node0_printf("CONGRAD: (re)start %d rsq/src = %.8g\n",
               restarts, rsq / source_norm);
#endif

  rsqstop = rsqmin * source_norm;
  iteration++;    // Count number of multiplications by (M^dag M)
  total_iters++;

  if (rsq <= rsqstop) {
    // If we want both parities, reset l_parity to ODD and repeat
    if (parity == EVENANDODD) {
      l_parity = ODD;
      l_otherparity = EVEN;
      parity = EVEN;  // So we won't loop endlessly
      iteration = 0;
      goto start;
    }

    // Otherwise, we're done
    if (special_started == 1) {   // Clean up gathers
      cleanup_gathers(tags1, tags2);
      special_started = 0;
    }

#ifdef CGTIME
    dtimec += dclock();
    node0_printf("CONGRAD: time = %.4g iters = %d mflops = %.4g\n",
                 dtimec, iteration, (double)(nflop * volume * iteration
                                     / (1.0e6 * dtimec * numnodes())));
    fflush(stdout);
#endif
    return iteration;
  }

  // Main loop -- do until convergence or time to restart
  // Note negative signs in mp, pkp and a cancel out with adds
  // oldrsq <- rsq
  // mp <- -M^dag.M.p
  // pkp <- p.mp
  // a <- -rsq / pkp
  // psi <- psi + a * p
  // r <- r + a * mp
  // rsq <- |r|^2
  // b <- rsq / oldrsq
  // p <- r + b * p
  do { // while (iteration % niter != 0);
    oldrsq = rsq;
    if (special_started == 0) {
      dslash_special(F_OFFSET(p), F_OFFSET(mp), l_otherparity, tags2, 1);
      dslash_special(F_OFFSET(mp), F_OFFSET(mp), l_parity, tags1, 1);
      special_started = 1;
    }
    else {
      dslash_special(F_OFFSET(p), F_OFFSET(mp), l_otherparity, tags2, 0);
      dslash_special(F_OFFSET(mp), F_OFFSET(mp), l_parity, tags1, 0);
    }

    // Finish computation of M^dag.M.p and p.M^dag.M.p
    // mp <- mp - msq_x4 * p
    // pkp <- p.mp
    pkp = 0;
    FORSOMEPARITY(i, s, l_parity) {
      if (i < loopend - FETCH_UP)
        prefetch_VV(&((s + FETCH_UP)->mp), &((s + FETCH_UP)->p));

      scalar_mult_add_vector(&(s->mp), &(s->p), -msq_x4, &(s->mp));
      pkp += (double)su3_rdot(&(s->p), &(s->mp));
    } END_LOOP
    g_doublesum(&pkp);
    iteration++;
    total_iters++;

    // Note negative sign; this is -a
    a = (Real)(-rsq/pkp);

    // psi <- psi - a * p
    // r <- r - a * mp
    rsq = 0;
    FORSOMEPARITY(i, s, l_parity) {
      if (i < loopend - FETCH_UP)
        prefetch_VVVV((vector *)F_PT(s + FETCH_UP, psi),
                      &((s+FETCH_UP)->p), &((s+FETCH_UP)->r),
                      &((s+FETCH_UP)->mp));

      scalar_mult_add_vector((vector *)F_PT(s, psi),
                                 &(s->p), a, (vector *)F_PT(s, psi));
      scalar_mult_add_vector(&(s->r), &(s->mp), a, &(s->r));
      rsq += (double)magsq_vec(&(s->r));
    } END_LOOP
    g_doublesum(&rsq);

#ifdef CG_DEBUG
    node0_printf("CONGRAD: iter %d rsq/src = %.8g, pkp = %.8g\n",
                 iteration, rsq / source_norm, pkp);
#endif

    if (rsq <= rsqstop) {
      // If we want both parities, reset l_parity to ODD and repeat
      if (parity == EVENANDODD) {
        l_parity = ODD;
        l_otherparity = EVEN;
        parity = EVEN;  // So we won't loop endlessly
        iteration = 0;
        goto start;
      }

      // Otherwise, we're done
      if (special_started == 1) {   // Clean up gathers
        cleanup_gathers(tags1, tags2);
        special_started = 0;
      }

#ifdef CGTIME
      dtimec += dclock();
      node0_printf("CONGRAD: time = %.4g iters = %d mflops = %.4g\n",
                   dtimec, iteration, (double)(nflop * volume * iteration
                                       / (1.0e6 * dtimec * numnodes())));
      fflush(stdout);
#endif
      return iteration;
    }

    b = (Real)rsq / oldrsq;
    // p <- r + b * p
    scalar_mult_add_latvec(F_OFFSET(r), F_OFFSET(p), b,
                           F_OFFSET(p), l_parity);
  } while (iteration % niter != 0);

  if (iteration < nrestart * niter) {
#ifdef CG_DEBUG
    node0_printf("CONGRAD: restarting after %d iterations\n", iteration);
#endif
    restarts++;
    goto start;
  }

  // If we want both parities, reset l_parity to ODD and repeat
  if (parity == EVENANDODD) {
    l_parity = ODD;
    l_otherparity = EVEN;
    parity = EVEN;          // So we won't loop endlessly */
    iteration = 0;
    goto start;
  }

  if (special_started == 1) {   // Clean up gathers
    cleanup_gathers(tags1, tags2);
    special_started = 0;
  }

  node0_printf("CONGRAD: did not converge after %d iterations, ",
               iteration);
  node0_printf("rsq/src = %.8g wanted %.8g\n",
               rsq / source_norm, rsqmin);
  fflush(stdout);
  return iteration;
}
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
                    (vector *)gen_pt[TDOWN][i + FETCH_UP]);

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



// -----------------------------------------------------------------
// Special dslash for use by CG
// Uses restart_gather() when possible
// Last argument is an array of message tags,
// to be set if this is the first use, otherwise reused
// If start = 1, use start_gather, otherwise use restart_gather
// The calling program must clean up the gathers!
void dslash_special(field_offset chi, field_offset psi, int parity,
                    msg_tag **tag, int start) {

  register int i, dir, otherparity = EVEN;
  register site *s;
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
  for (dir = XUP; dir <= TUP; dir++) {
    if (start == 1)
      tag[dir] = start_gather_site(chi, sizeof(vector), dir,
                                   parity, gen_pt[dir]);
    else
      restart_gather_site(chi, sizeof(vector), dir,
                          parity, gen_pt[dir], tag[dir]);
  }

  // Multiply by adjoint matrix at other sites
  FORSOMEPARITY(i, s, otherparity) {
    if (i < loopend-FETCH_UP)
      prefetch_4MV4V(&((s + FETCH_UP)->link[XUP]),
                     (vector *)F_PT((s + FETCH_UP), chi),
                     (s + FETCH_UP)->tempvec);

    mult_adj_mat_vec_4dir(s->link, (vector *)F_PT(s, chi),
                              s->tempvec);
  } END_LOOP

  // Start gathers from negative directions
  for (dir = XUP; dir <= TUP; dir++) {
    if (start == 1)
      tag[OPP_DIR(dir)] = start_gather_site(F_OFFSET(tempvec[dir]),
                                            sizeof(vector), OPP_DIR(dir),
                                            parity, gen_pt[OPP_DIR(dir)]);
    else
      restart_gather_site(F_OFFSET(tempvec[dir]), sizeof(vector),
                          OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)],
                          tag[OPP_DIR(dir)]);
    }

  // Wait gathers from positive directions
  for (dir = XUP; dir <= TUP; dir++)
    wait_gather(tag[dir]);

  // Multiply by matrix and accumulate
  FORSOMEPARITY(i, s, parity) {
    if (i < loopend - FETCH_UP)
      prefetch_4MVVVV(&((s + FETCH_UP)->link[XUP]),
                      (vector *)gen_pt[XUP][i + FETCH_UP],
                      (vector *)gen_pt[YUP][i + FETCH_UP],
                      (vector *)gen_pt[ZUP][i + FETCH_UP],
                      (vector *)gen_pt[TUP][i + FETCH_UP]);

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
    if (i < loopend-FETCH_UP)
      prefetch_VVVV((vector *)gen_pt[XDOWN][i + FETCH_UP],
                    (vector *)gen_pt[YDOWN][i + FETCH_UP],
                    (vector *)gen_pt[ZDOWN][i + FETCH_UP],
                    (vector *)gen_pt[TDOWN][i + FETCH_UP]);

#ifndef INLINE
    // Non-inline version
    sub_four_vecs((vector *)F_PT(s, psi),
                      (vector *)(gen_pt[XDOWN][i]),
                      (vector *)(gen_pt[YDOWN][i]),
                      (vector *)(gen_pt[ZDOWN][i]),
                      (vector *)(gen_pt[TDOWN][i]));
#else
    // Inline version
    a  = (vector *)F_PT(s, psi);
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
}
// -----------------------------------------------------------------

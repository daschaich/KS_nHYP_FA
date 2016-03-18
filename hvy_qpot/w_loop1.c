// -----------------------------------------------------------------
// Compute time-like, planar Wilson loops in axial gauge
// with non-trivial timelike links copied to all timeslices

// Includes and definitions
#include "hvy_qpot_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// tot_smear will label the output; currently we only use 1
void w_loop1(int tot_smear) {
  register int i, dir, r, t;
  register site *s;
  int nth, nxh;
  double *wils_loop1, td;
  su3_matrix tmat1, tmat2;
  msg_tag *mtag[4];

  if (nx != ny || nx != ny) {
    node0_printf("w_loop1 gives wrong results for nx != ny != nz");
    return;
  }

  nth = nt / 2;
  nxh = nx / 2;
  wils_loop1 = (double *)malloc(nth * nxh * sizeof(double));
  for (t = 0; t < nth; t++) {
    for (r = 0; r < nxh; r++)
      wils_loop1[r + nxh * t] = 0;
  }

  for (dir = XUP; dir <= ZUP; dir++) {
    FORALLSITES(i, s) {
      su3mat_copy(&(s->link[dir]), &(s->s_link));
      su3mat_copy(&(s->link[TUP]), &(s->t_link_f));
    }

    // Start gather of forward temporal links
    mtag[1] = NULL;    // Avoid compiler warning
    mtag[0] = start_gather_site(F_OFFSET(t_link_f), sizeof(su3_matrix),
                                dir, EVENANDODD, gen_pt[0]);

    // Recursively construct the spatial segments and compute
    // the Wilson loops with that segment
    for (r = 0; r < nxh; r++) {
      if (r > 0) {
        wait_gather(mtag[1]);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[1][i]), &(s->staple));

        FORALLSITES(i, s)
          mult_su3_nn(&(s->link[dir]), &(s->staple), &(s->s_link));
      }

      FORALLSITES(i, s)
        su3mat_copy(&(s->s_link), &(s->s_link_f));

      // Start gather of forward spatial segments
      mtag[TUP] = start_gather_site(F_OFFSET(s_link_f), sizeof(su3_matrix),
                                    TUP, EVENANDODD, gen_pt[TUP]);

      // Concurrently gather spatial links for next r, if still needed
      if (r == 0)
        mtag[1] = start_gather_site(F_OFFSET(s_link), sizeof(su3_matrix),
                                    dir, EVENANDODD, gen_pt[1]);
      else if (r < nxh - 1)
        restart_gather_site(F_OFFSET(s_link), sizeof(su3_matrix),
                            dir, EVENANDODD, gen_pt[1], mtag[1]);
      else
        cleanup_gather(mtag[1]);

      // Collect forward temporal links
      wait_gather(mtag[0]);
      FORALLSITES(i, s)
        su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->staple));

      FORALLSITES(i, s)
        su3mat_copy(&(s->staple), &(s->t_link_f));

      // Recursively compute the Wilson loops of different time extent
      for (t = 0; t < nth; t++) {
        // Collect forward spatial segments
        wait_gather(mtag[TUP]);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[TUP][i]), &(s->staple));

        FORALLSITES(i, s)
          su3mat_copy(&(s->staple), &(s->s_link_f));

        // Start gather for next t, if still needed.
        if (t < nth - 1)
          restart_gather_site(F_OFFSET(s_link_f), sizeof(su3_matrix),
                            TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
        else
          cleanup_gather(mtag[TUP]);

        // Finally, compute the Wilson loops
        FORALLSITES(i, s) {
          if ((s->t) + t + 1 >= nt) {
            mult_su3_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
            mult_su3_na(&tmat1, &(s->t_link_f), &tmat2);
            wils_loop1[r + nxh * t] += (double)realtrace_su3(&tmat2,
                                                             &(s->s_link));
          }
          else {
            wils_loop1[r + nxh * t] += (double)realtrace_su3(&(s->s_link_f),
                                                             &(s->s_link));
          }
        }
      } // End loop over t

      // Start gather of forward time-like links for next r
      if (r < (nxh - 1))
        restart_gather_site(F_OFFSET(t_link_f), sizeof(su3_matrix),
                            dir, EVENANDODD, gen_pt[0], mtag[0]);
      else
        cleanup_gather(mtag[0]);
    } // End loop over r
  } // End loop over dir

  // Normalize and print the Wilson loops
  for (t = 0; t < nth; t++) {
    for (r = 0; r < nxh; r++) {
      td = wils_loop1[r + nxh * t];
      g_doublesum(&td);
      td /= (double)(9 * volume);
      node0_printf("WILS_LOOP1_%d %d %d %.8g\n", tot_smear, r, t, td);
    }
  }
  free(wils_loop1);
}
// -----------------------------------------------------------------

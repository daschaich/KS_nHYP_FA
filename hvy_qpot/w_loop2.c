// -----------------------------------------------------------------
// Compute temporal, off-axis Wilson loops in axial gauge
// with non-trivial timelike links copied to all timeslices

// Includes and definitions
#include "hvy_qpot_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// tot_smear will label the output; currently we only use 1
void w_loop2(int tot_smear) {
  register int i,dir1,dir2,dir3,r,t,r_off;
  register site *s;
  int nth,nxh,nrmax;
  int disp[4];      // Displacement vector for general gather
  double *wils_loop2, dt;
  matrix tmat1, tmat2;
  msg_tag *mtag[8], *gmtag;

  if (nx != ny || nx != nz) {
    node0_printf("w_loop2 gives wrong results for nx != ny != nz");
    return;
  }

  nth = nt / 2;
  nxh = nx / 2;
  nrmax = 2 * nxh + nxh / 2;
  wils_loop2 = (double *)malloc(nth * nrmax * sizeof(double));
  for (t = 0; t < nth; t++) {
    for (r = 0; r < nrmax; r++) {
      wils_loop2[r + nrmax * t] = 0;
    }
  }

  for (dir1 = XUP; dir1 <= YUP; dir1++) {
    for (dir2 = dir1 + 1; dir2 <= ZUP; dir2++) {
      // First "sqrt(2)" loops in the (dir1, dir2) plane
      // Construct the "diagonal" link in the (dir1, dir2) direction
      mtag[dir1] = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                     dir1, EVENANDODD, gen_pt[dir1]);
      mtag[dir2] = start_gather_site(F_OFFSET(link[dir1]), sizeof(matrix),
                                     dir2, EVENANDODD, gen_pt[dir2]);

      FORALLSITES(i,s)
        mat_copy(&(s->link[TUP]), &(s->t_link_f));

      // Concurrently start gather of temporal links across the diagonal
      for (i = XUP; i <= TUP; i++)
        disp[i] = 0;

      disp[dir1] = 1;
      disp[dir2] = 1;
      gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

      wait_gather(mtag[dir1]);
      wait_gather(mtag[dir2]);
      dt = 0.5;
      FORALLSITES(i, s) {
        mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[dir1][i]),
                    &(s->diag));
        mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[dir2][i]),
                    &tmat1);
        add_matrix(&(s->diag), &tmat1, &(s->diag));
        scalar_mult_matrix(&(s->diag), dt, &(s->diag));
      }
      cleanup_gather(mtag[dir1]);
      cleanup_gather(mtag[dir2]);

      // Recursively construct the spatial segments
      // Compute the Wilson loops with that segment
      for (r = 0; r < nxh; r++) {
        if (r == 0) {
          FORALLSITES(i, s)
            mat_copy(&(s->diag), &(s->s_link));
        }
        else {
          wait_general_gather(gmtag);
          FORALLSITES(i, s)
            mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

          cleanup_general_gather(gmtag);

          // Concurrently gather temporal links across the diagonal
          gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                            sizeof(matrix), disp,
                                            EVENANDODD, gen_pt[4]);

          FORALLSITES(i, s)
            mult_nn(&(s->diag), &(s->staple), &(s->s_link));
        }
        FORALLSITES(i, s)
          mat_copy(&(s->s_link), &(s->s_link_f));

        // Start gather of forward spatial segments
        mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                      sizeof(matrix), TUP,
                                      EVENANDODD, gen_pt[TUP]);

        // Collect forward temporal links
        wait_general_gather(gmtag);
        FORALLSITES(i, s)
          mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

        FORALLSITES(i, s)
          mat_copy(&(s->staple), &(s->t_link_f));

        cleanup_general_gather(gmtag);

        // Concurrently gather spatial links across the diagonal for next r
        if (r < nxh - 1)
          gmtag = start_general_gather_site(F_OFFSET(s_link),
                                            sizeof(matrix), disp,
                                            EVENANDODD, gen_pt[4]);

        // Recursively compute the Wilson loops of different time extent
        for (t = 0; t < nth; t++) {
          // Collect forward spatial segments
          wait_gather(mtag[TUP]);
          FORALLSITES(i, s)
            mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

          FORALLSITES(i, s)
            mat_copy(&(s->staple), &(s->s_link_f));

          // Start gather for next t, if still needed
          if (t < nth - 1)
            restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                                TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
          else
            cleanup_gather(mtag[TUP]);

          // Finally, compute the Wilson loops
          FORALLSITES(i, s) {
            if ((s->t) + t + 1 >= nt) {
              mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
              mult_na(&tmat1, &(s->t_link_f), &tmat2);
              wils_loop2[r + nrmax * t]
                += (double)realtrace(&tmat2, &(s->s_link));
            }
            else
              wils_loop2[r + nrmax * t]
                += (double)realtrace(&(s->s_link_f), &(s->s_link));
          }
        } // End loop over t
      } // End loop over r

      // Next, "sqrt(2)" loops in the (dir1, -dir2) plane
      // Construct the "diagonal" link in the (dir1, -dir2) direction
      // Gather dir2 link from across the diagonal
      for (i = XUP; i <= TUP; i++)
        disp[i] = 0;

      disp[dir1] = 1;
      disp[dir2] = -1;
      gmtag = start_general_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                        disp, EVENANDODD, gen_pt[4]);

      // Multiply one corner and then gather it
      FORALLSITES(i, s)
        mult_an(&(s->link[dir2]), &(s->link[dir1]), &(s->staple));

      mtag[dir2] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                                     OPP_DIR(dir2), EVENANDODD, gen_pt[dir2]);

      FORALLSITES(i, s)
        mat_copy(&(s->link[TUP]), &(s->t_link_f));

      // Make second corner
      wait_general_gather(gmtag);
      FORALLSITES(i, s)
        mult_na(&(s->link[dir1]), (matrix *)(gen_pt[4][i]),
                    &(s->diag));

      cleanup_general_gather(gmtag);

      // Collect first corner and add to second
      wait_gather(mtag[dir2]);
      dt = 0.5;
      FORALLSITES(i, s) {
        add_matrix(&(s->diag), (matrix *)(gen_pt[dir2][i]),
                       &(s->diag));
        scalar_mult_matrix(&(s->diag), dt, &(s->diag));
      }
      cleanup_gather(mtag[dir2]);

      // Start gather of temporal links across the diagonal
      gmtag = start_general_gather_site(F_OFFSET(t_link_f), sizeof(matrix),
                                        disp, EVENANDODD, gen_pt[4]);

      // Recursively construct the spatial segments
      // Compute the Wilson loops with that segment
      for (r = 0; r < nxh; r++) {
        if (r == 0) {
          FORALLSITES(i, s)
            mat_copy(&(s->diag), &(s->s_link));
        }
        else {
          wait_general_gather(gmtag);
          FORALLSITES(i, s)
            mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

          cleanup_general_gather(gmtag);

          // Concurrently gather temporal links across the diagonal
          gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                            sizeof(matrix), disp,
                                            EVENANDODD, gen_pt[4]);

          FORALLSITES(i, s)
            mult_nn(&(s->diag), &(s->staple), &(s->s_link));
        }
        FORALLSITES(i, s)
          mat_copy(&(s->s_link), &(s->s_link_f));

        // Start gather of forward spatial segments
        mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                      sizeof(matrix), TUP,
                                      EVENANDODD, gen_pt[TUP]);

        // Collect forward temporal links
        wait_general_gather(gmtag);
        FORALLSITES(i, s)
          mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

        FORALLSITES(i, s)
          mat_copy(&(s->staple), &(s->t_link_f));

        cleanup_general_gather(gmtag);

        // Concurrently gather spatial links across the diagonal for next r
        if (r < nxh - 1)
          gmtag = start_general_gather_site(F_OFFSET(s_link),
                                            sizeof(matrix), disp,
                                            EVENANDODD, gen_pt[4]);

        // Recursively compute the Wilson loops of different time extent
        for (t = 0; t < nth; t++) {
          // Collect forward spatial segments
          wait_gather(mtag[TUP]);
          FORALLSITES(i, s)
            mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

          FORALLSITES(i, s)
            mat_copy(&(s->staple), &(s->s_link_f));

          // Start gather for next t, if still needed
          if (t < nth - 1)
            restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                                TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
          else
            cleanup_gather(mtag[TUP]);

          // Finally, compute the Wilson loops
          FORALLSITES(i, s) {
            if ((s->t) + t + 1 >= nt) {
              mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
              mult_na(&tmat1, &(s->t_link_f), &tmat2);
              wils_loop2[r + nrmax * t] += (double)realtrace(&tmat2,
                                                                 &(s->s_link));
            }
            else
              wils_loop2[r + nrmax * t] += (double)realtrace(&(s->s_link_f),
                                                                 &(s->s_link));
          }
        } // End loop over t
      } // End loop over r
    } // End loop over dir2 > dir1
  } // End loop over dir1

  // Offset for next bunch of Wilson loops
  r_off = nxh;

  for (dir1 = XUP; dir1 <= ZUP; dir1++) {
    for (dir2 = XUP; dir2 <= ZUP; dir2++) {
      if (dir1 != dir2) {
        // First, "sqrt(5)" loops in (dir1, dir2) plane
        // Construct the "diagonal" link in the (2 * dir1, dir2) direction
        mtag[dir1] = start_gather_site(F_OFFSET(link[dir1]),
                                       sizeof(matrix), dir1,
                                       EVENANDODD, gen_pt[dir1]);

        // Start gather of dir2-link from 2 * dir1
        for (i = XUP; i <= TUP; i++)
          disp[i] = 0;

        disp[dir1] = 2;
        gmtag = start_general_gather_site(F_OFFSET(link[dir2]),
                                          sizeof(matrix), disp,
                                          EVENANDODD, gen_pt[4]);

        FORALLSITES(i, s)
          mat_copy(&(s->link[TUP]), &(s->t_link_f));

        // Make double links in dir1 direction
        wait_gather(mtag[dir1]);
        FORALLSITES(i, s)
          mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[dir1][i]),
                      &(s->s_link));

        cleanup_gather(mtag[dir1]);

        // Gather the double links from dir2 direction
        mtag[dir2] = start_gather_site(F_OFFSET(s_link), sizeof(matrix),
                                       dir2, EVENANDODD, gen_pt[dir2]);

        // Make first corner
        wait_general_gather(gmtag);
        FORALLSITES(i, s)
          mult_nn(&(s->s_link), (matrix *)(gen_pt[4][i]),
                      &(s->diag));

        cleanup_general_gather(gmtag);

        // Make second corner and add to first
        wait_gather(mtag[dir2]);
        dt = 0.5;
        FORALLSITES(i, s) {
          mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[dir2][i]),
                      &tmat1);
          add_matrix(&(s->diag), &tmat1, &(s->diag));
          scalar_mult_matrix(&(s->diag), dt, &(s->diag));
        }
        cleanup_gather(mtag[dir2]);

        // Start gather of temporal links across the diagonal
        for (i = XUP; i <= TUP; i++)
          disp[i] = 0;

        disp[dir1] = 2;
        disp[dir2] = 1;
        gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                          sizeof(matrix), disp,
                                          EVENANDODD, gen_pt[4]);

        // Recursively construct the spatial segments
        // Compute the Wilson loops with that segment
        for (r = 0; r < nxh / 2; r++) {
          if (r == 0) {
            FORALLSITES(i, s)
              mat_copy(&(s->diag), &(s->s_link));
          }
          else {
            wait_general_gather(gmtag);
            FORALLSITES(i, s)
              mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

            cleanup_general_gather(gmtag);

            // Concurrently gather temporal links across the diagonal
            gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                              sizeof(matrix), disp,
                                              EVENANDODD, gen_pt[4]);

            FORALLSITES(i, s)
              mult_nn(&(s->diag), &(s->staple), &(s->s_link));
          }
          FORALLSITES(i, s)
            mat_copy(&(s->s_link), &(s->s_link_f));

          // Start gather of forward spatial segments
          mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                        sizeof(matrix), TUP,
                                        EVENANDODD, gen_pt[TUP]);

          // Collect forward temporal links
          wait_general_gather(gmtag);
          FORALLSITES(i, s)
            mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

          FORALLSITES(i, s)
            mat_copy(&(s->staple), &(s->t_link_f));

          cleanup_general_gather(gmtag);

          // Concurrently gather spatial links across the diagonal for next r
          if (r < nxh / 2 - 1)
            gmtag = start_general_gather_site(F_OFFSET(s_link),
                                              sizeof(matrix), disp,
                                              EVENANDODD, gen_pt[4]);

          // Recursively compute the Wilson loops of different time extent
          for (t = 0; t < nth; t++) {
            // Collect forward spatial segments
            wait_gather(mtag[TUP]);
            FORALLSITES(i, s)
              mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

            FORALLSITES(i, s)
              mat_copy(&(s->staple), &(s->s_link_f));

            // Start gather for next t, if still needed
            if (t < nth - 1)
              restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                                  TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
            else
              cleanup_gather(mtag[TUP]);

            // Finally, compute the Wilson loops
            FORALLSITES(i, s) {
              if ((s->t) + t + 1 >= nt) {
                mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
                mult_na(&tmat1, &(s->t_link_f), &tmat2);
                wils_loop2[r + r_off + nrmax * t]
                  += (double)realtrace(&tmat2, &(s->s_link));
              }
              else
                wils_loop2[r + r_off + nrmax * t]
                  += (double)realtrace(&(s->s_link_f), &(s->s_link));
            }
          } // End loop over t
        } // End loop over r

        // Next, "sqrt(5)" loops in the (dir1, -dir2) plane
        // Construct the "diagonal" link in the (2 * dir1, -dir2) direction
        mtag[dir1] = start_gather_site(F_OFFSET(link[dir1]), sizeof(matrix),
                                       dir1, EVENANDODD, gen_pt[dir1]);

        // Gather dir2-link from across the diagonal
        for (i = XUP; i <= TUP; i++)
          disp[i] = 0;

        disp[dir1] = 2;
        disp[dir2] = -1;
        gmtag = start_general_gather_site(F_OFFSET(link[dir2]),
                                          sizeof(matrix), disp,
                                          EVENANDODD, gen_pt[4]);

        FORALLSITES(i, s)
          mat_copy(&(s->link[TUP]), &(s->t_link_f));

        // Make double links in dir1 direction
        wait_gather(mtag[dir1]);
        FORALLSITES(i, s)
          mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[dir1][i]),
                      &(s->s_link));

        cleanup_gather(mtag[dir1]);

        // Make one corner and then gather it
        FORALLSITES(i, s)
        mult_an(&(s->link[dir2]), &(s->s_link), &(s->staple));

        mtag[dir2] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                                       OPP_DIR(dir2), EVENANDODD, gen_pt[dir2]);

        // Make second corner
        wait_general_gather(gmtag);
        FORALLSITES(i, s)
          mult_na(&(s->s_link), (matrix *)(gen_pt[4][i]),
                      &(s->diag));

        cleanup_general_gather(gmtag);

        // Collect first corner and add to second
        wait_gather(mtag[dir2]);
        dt = 0.5;
        FORALLSITES(i, s) {
          add_matrix(&(s->diag), (matrix *)(gen_pt[dir2][i]),
                         &(s->diag));
          scalar_mult_matrix(&(s->diag), dt, &(s->diag));
        }
        cleanup_gather(mtag[dir2]);

        // Start gather of temporal links across the diagonal
        gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                          sizeof(matrix), disp,
                                          EVENANDODD, gen_pt[4]);

        // Recursively construct the spatial segments
        // Compute the Wilson loops with that segment
        for (r = 0; r < nxh / 2; r++) {
          if (r == 0) {
            FORALLSITES(i, s)
              mat_copy(&(s->diag), &(s->s_link));
          }
          else {
            wait_general_gather(gmtag);
            FORALLSITES(i, s)
              mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

            cleanup_general_gather(gmtag);

            // Concurrently gather temporal links across the diagonal
            gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                              sizeof(matrix), disp,
                                              EVENANDODD, gen_pt[4]);

            FORALLSITES(i, s)
              mult_nn(&(s->diag), &(s->staple), &(s->s_link));
          }
          FORALLSITES(i, s)
            mat_copy(&(s->s_link), &(s->s_link_f));

          // Start gather of forward spatial segments
          mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                        sizeof(matrix), TUP,
                                        EVENANDODD, gen_pt[TUP]);

          // Collect forward temporal links
          wait_general_gather(gmtag);
          FORALLSITES(i, s)
            mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

          FORALLSITES(i, s)
            mat_copy(&(s->staple), &(s->t_link_f));

          cleanup_general_gather(gmtag);

          // Concurrently gather spatial links across the diagonal for next r
          if (r < nxh / 2 - 1)
            gmtag = start_general_gather_site(F_OFFSET(s_link),
                                              sizeof(matrix), disp,
                                              EVENANDODD, gen_pt[4]);

          // Recursively compute the Wilson loops of different time extent
          for (t = 0; t < nth; t++) {
            // Collect forward spatial segments
            wait_gather(mtag[TUP]);
            FORALLSITES(i, s)
              mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

            FORALLSITES(i, s)
              mat_copy(&(s->staple), &(s->s_link_f));

            // Start gather for next t, if still needed
            if (t < nth - 1)
              restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                                  TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
            else
              cleanup_gather(mtag[TUP]);

            // Finally, compute the Wilson loops
            FORALLSITES(i, s) {
              if ((s->t) + t + 1 >= nt) {
                mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
                mult_na(&tmat1, &(s->t_link_f), &tmat2);
                wils_loop2[r + r_off + nrmax * t]
                  += (double)realtrace(&tmat2, &(s->s_link));
              }
              else
                wils_loop2[r + r_off + nrmax * t]
                  += (double)realtrace(&(s->s_link_f), &(s->s_link));
            }
          } // End loop over t
        } // End loop over r
      }
    } // End loop over dir2 != dir1
  } // End loop over dir1

  // Offset for next bunch of Wilson loops
  r_off = nxh + nxh / 2;
  dir1 = XUP;
  dir2 = YUP;
  dir3 = ZUP;

  // First, "sqrt(3)" loops in (x, y, z) space
  // Construct the "body diagonal" link in the (x, y, z) direction
  // Gather for first "plaquette"
  mtag[0] = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                              dir1, EVENANDODD, gen_pt[0]);
  mtag[1] = start_gather_site(F_OFFSET(link[dir1]), sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[1]);

  // Gather for second "plaquette"
  mtag[7] = start_gather_site(F_OFFSET(link[dir3]), sizeof(matrix),
                              dir1, EVENANDODD, gen_pt[7]);
  mtag[5] = start_gather_site(F_OFFSET(link[dir1]), sizeof(matrix),
                              dir3, EVENANDODD, gen_pt[5]);

  FORALLSITES(i, s)
    mat_copy(&(s->link[TUP]), &(s->t_link_f));

  // Make diagonal link in (x, y) direction and gather it
  wait_gather(mtag[0]);
  wait_gather(mtag[1]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[0][i]),
                &(s->s_link));
    mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[1][i]), &tmat1);
    add_matrix(&(s->s_link), &tmat1, &(s->s_link));
  }
  cleanup_gather(mtag[0]);
  cleanup_gather(mtag[1]);
  mtag[2] = start_gather_site(F_OFFSET(s_link), sizeof(matrix),
                              dir3, EVENANDODD, gen_pt[2]);

  // Gather for third "plaquette"
  mtag[1] = start_gather_site(F_OFFSET(link[dir3]), sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[1]);
  mtag[3] = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                              dir3, EVENANDODD, gen_pt[3]);

  // Make diagonal link in (x, z) direction and gather it
  wait_gather(mtag[7]);
  wait_gather(mtag[5]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[7][i]),
                &(s->s_link_f));
    mult_nn(&(s->link[dir3]), (matrix *)(gen_pt[5][i]), &tmat1);
    add_matrix(&(s->s_link_f), &tmat1, &(s->s_link_f));
  }
  cleanup_gather(mtag[7]);
  cleanup_gather(mtag[5]);
  mtag[6] = start_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[6]);

  // Make first body diagonal
  wait_gather(mtag[2]);
  FORALLSITES(i, s)
    mult_nn(&(s->link[dir3]), (matrix *)(gen_pt[2][i]), &(s->diag));

  cleanup_gather(mtag[2]);

  // Make diagonal link in (y, z) direction and gather it
  wait_gather(mtag[1]);
  wait_gather(mtag[3]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[1][i]),
                &(s->s_link));
    mult_nn(&(s->link[dir3]), (matrix *)(gen_pt[3][i]), &tmat1);
    add_matrix(&(s->s_link), &tmat1, &(s->s_link));
  }
  cleanup_gather(mtag[1]);
  cleanup_gather(mtag[3]);
  mtag[0] = start_gather_site(F_OFFSET(s_link), sizeof(matrix),
                              dir1, EVENANDODD, gen_pt[0]);

  // Make second body diagonal and add to first
  wait_gather(mtag[6]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[6][i]), &tmat1);
    add_matrix(&(s->diag), &tmat1, &(s->diag));
  }
  cleanup_gather(mtag[6]);

  // Make third body diagonal and add
  wait_gather(mtag[0]);
  dt = 1.0 / 6.0;
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[0][i]), &tmat1);
    add_matrix(&(s->diag), &tmat1, &(s->diag));
    scalar_mult_matrix(&(s->diag), dt, &(s->diag));
  }
  cleanup_gather(mtag[0]);

  // Start gather of temporal links across the body diagonal
  disp[dir1] = 1;
  disp[dir2] = 1;
  disp[dir3] = 1;
  disp[TUP] = 0;
  gmtag = start_general_gather_site(F_OFFSET(t_link_f), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);

  // Recursively construct the spatial segments
  // Compute the Wilson loops with that segment
  for (r = 0; r < nxh; r++) {
    if (r == 0) {
      FORALLSITES(i, s)
        mat_copy(&(s->diag), &(s->s_link));
    }
    else {
      wait_general_gather(gmtag);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

      cleanup_general_gather(gmtag);

      // Concurrently gather temporal links across the diagonal
      gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

      FORALLSITES(i, s)
        mult_nn(&(s->diag), &(s->staple), &(s->s_link));
    }
    FORALLSITES(i, s)
      mat_copy(&(s->s_link), &(s->s_link_f));

    // Start gather of forward spatial segments
    mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                  sizeof(matrix), TUP,
                                  EVENANDODD, gen_pt[TUP]);

    // Collect forward temporal links
    wait_general_gather(gmtag);
    FORALLSITES(i, s)
      mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

    FORALLSITES(i, s)
      mat_copy(&(s->staple), &(s->t_link_f));

    cleanup_general_gather(gmtag);

    // Concurrently gather spatial links across the diagonal for next r
    if (r < nxh - 1)
      gmtag = start_general_gather_site(F_OFFSET(s_link),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

    // Recursively compute the Wilson loops of different time extent
    for (t = 0; t < nth; t++) {
      // Collect forward spatial segments
      wait_gather(mtag[TUP]);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

      FORALLSITES(i, s)
        mat_copy(&(s->staple), &(s->s_link_f));

      // Start gather for next t, if still needed
      if (t < nth - 1)
        restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                            TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
      else
        cleanup_gather(mtag[TUP]);

      // Finally, compute the Wilson loops
      FORALLSITES(i, s) {
        if ((s->t) + t + 1 >= nt) {
          mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
          mult_na(&tmat1, &(s->t_link_f), &tmat2);
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&tmat2, &(s->s_link));
        }
        else
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&(s->s_link_f), &(s->s_link));
      }
    } // End loop over t
  } // End loop over r

  // Next, "sqrt(3)" loops in (x, y, -z) space
  // Construct the "body diagonal" link in the (x, y, -z) direction
  // Gather for first "plaquette"
  mtag[0] = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                              dir1, EVENANDODD, gen_pt[0]);
  mtag[1] = start_gather_site(F_OFFSET(link[dir1]), sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[1]);

  // Start one corner for second "plaquette" and gather
  FORALLSITES(i, s)
      mult_an(&(s->link[dir3]), &(s->link[dir1]), &(s->staple));

  for (i = XUP; i <= TUP; i++)
    disp[i] = 0;

  disp[dir1] = 1;
  disp[dir3] = -1;
  gmtag = start_general_gather_site(F_OFFSET(link[dir3]), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);
  mtag[5] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                              OPP_DIR(dir3), EVENANDODD, gen_pt[5]);

  // Make diagonal link in (x, y) direction
  // Multiply to get first body diagonal and gather it
  wait_gather(mtag[0]);
  wait_gather(mtag[1]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[0][i]), &tmat1);
    mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[1][i]), &tmat2);
    add_matrix(&tmat1, &tmat2, &tmat1);
    mult_an(&(s->link[dir3]), &tmat1, &(s->s_link));
  }
  cleanup_gather(mtag[0]);
  cleanup_gather(mtag[1]);
  mtag[2] = start_gather_site(F_OFFSET(s_link), sizeof(matrix),
                              OPP_DIR(dir3), EVENANDODD, gen_pt[2]);

  // Make diagonal link in (x, -z) direction and gather it
  wait_general_gather(gmtag);
  wait_gather(mtag[5]);
  FORALLSITES(i, s) {
    mult_na(&(s->link[dir1]), (matrix *)(gen_pt[4][i]),
                &(s->s_link_f));
    add_matrix(&(s->s_link_f), (matrix *)(gen_pt[5][i]),
                   &(s->s_link_f));
  }
  cleanup_general_gather(gmtag);
  cleanup_gather(mtag[5]);
  mtag[6] = start_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[6]);

  // Start one corner for third "plaquette" and gather
  FORALLSITES(i, s)
    mult_an(&(s->link[dir3]), &(s->link[dir2]), &(s->diag));

  for (i = XUP; i <= TUP; i++)
    disp[i] = 0;

  disp[dir2] = 1;
  disp[dir3] = -1;
  gmtag = start_general_gather_site(F_OFFSET(link[dir3]), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);
  mtag[3] = start_gather_site(F_OFFSET(diag), sizeof(matrix),
                              OPP_DIR(dir3), EVENANDODD, gen_pt[3]);

  FORALLSITES(i, s)
    mat_copy(&(s->link[TUP]), &(s->t_link_f));

  // Make diagonal link in (y, -z) direction and gather it
  wait_general_gather(gmtag);
  wait_gather(mtag[3]);
  FORALLSITES(i, s) {
    mult_na(&(s->link[dir2]), (matrix *)(gen_pt[4][i]),
                &(s->staple));
    add_matrix(&(s->staple), (matrix *)(gen_pt[3][i]),
                   &(s->staple));
  }
  cleanup_general_gather(gmtag);
  cleanup_gather(mtag[3]);
  mtag[0] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                              dir1, EVENANDODD, gen_pt[0]);

  // Make second body diagonal and add to that gathered first
  wait_gather(mtag[2]);
  wait_gather(mtag[6]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[6][i]), &(s->diag));
    add_matrix(&(s->diag), (matrix *)(gen_pt[2][i]), &(s->diag));
  }
  cleanup_gather(mtag[2]);
  cleanup_gather(mtag[6]);

  // Make third body diagonal and add
  wait_gather(mtag[0]);
  dt = 1.0 / 6.0;
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[0][i]), &tmat1);
    add_matrix(&(s->diag), &tmat1, &(s->diag));
    scalar_mult_matrix(&(s->diag), dt, &(s->diag));
  }
  cleanup_gather(mtag[0]);

  // Start gather of temporal links across the body diagonal
  disp[dir1] = 1;
  disp[dir2] = 1;
  disp[dir3] = -1;
  disp[TUP] = 0;
  gmtag = start_general_gather_site(F_OFFSET(t_link_f), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);

  // Recursively construct the spatial segments
  // Compute the Wilson loops with that segment
  for (r = 0; r < nxh; r++) {
    if (r == 0) {
      FORALLSITES(i, s)
        mat_copy(&(s->diag), &(s->s_link));
    }
    else {
      wait_general_gather(gmtag);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

      cleanup_general_gather(gmtag);

      // Concurrently gather temporal links across the diagonal
      gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

      FORALLSITES(i, s)
        mult_nn(&(s->diag), &(s->staple), &(s->s_link));
    }
    FORALLSITES(i, s)
      mat_copy(&(s->s_link), &(s->s_link_f));

    // Start gather of forward spatial segments
    mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                  sizeof(matrix), TUP,
                                  EVENANDODD, gen_pt[TUP]);

    // Collect forward temporal links
    wait_general_gather(gmtag);
    FORALLSITES(i, s)
      mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

    FORALLSITES(i, s)
      mat_copy(&(s->staple), &(s->t_link_f));

    cleanup_general_gather(gmtag);

    // Concurrently gather spatial links across the diagonal for next r
    if (r < nxh - 1)
      gmtag = start_general_gather_site(F_OFFSET(s_link),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

    // Recursively compute the Wilson loops of different time extent
    for (t = 0; t < nth; t++) {
      // Collect forward spatial segments
      wait_gather(mtag[TUP]);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

      FORALLSITES(i, s)
        mat_copy(&(s->staple), &(s->s_link_f));

      // Start gather for next t, if still needed
      if (t < nth - 1)
        restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                            TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
      else
        cleanup_gather(mtag[TUP]);

      // Finally, compute the Wilson loops
      FORALLSITES(i, s) {
        if ((s->t) + t + 1 >= nt) {
          mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
          mult_na(&tmat1, &(s->t_link_f), &tmat2);
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&tmat2, &(s->s_link));
        }
        else
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&(s->s_link_f), &(s->s_link));
      }
    } // End loop over t
  } // End loop over r

  // Next, "sqrt(3)" loops in (x, -y, z) space
  // Construct the "body diagonal" link in the (x, -y, z) direction
  // Gather for first "plaquette"
  mtag[0] = start_gather_site(F_OFFSET(link[dir3]), sizeof(matrix),
                              dir1, EVENANDODD, gen_pt[0]);
  mtag[2] = start_gather_site(F_OFFSET(link[dir1]), sizeof(matrix),
                              dir3, EVENANDODD, gen_pt[2]);

  // Start one corner for second "plaquette" and gather
  FORALLSITES(i, s)
    mult_an(&(s->link[dir2]), &(s->link[dir1]), &(s->staple));

  for (i = XUP; i <= TUP; i++)
    disp[i]=0;

  disp[dir1] = 1;
  disp[dir2] = -1;
  gmtag = start_general_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);
  mtag[6] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[6]);

  // Make diagonal link in (x, z) direction
  // Multiply to get first body diagonal and gather it
  wait_gather(mtag[0]);
  wait_gather(mtag[2]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[0][i]), &tmat1);
    mult_nn(&(s->link[dir3]), (matrix *)(gen_pt[2][i]), &tmat2);
    add_matrix(&tmat1, &tmat2, &tmat1);
    mult_an(&(s->link[dir2]), &tmat1, &(s->s_link));
  }
  cleanup_gather(mtag[0]);
  cleanup_gather(mtag[2]);
  mtag[1] = start_gather_site(F_OFFSET(s_link), sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

  // Make diagonal link in (x, -y) direction and gather it
  wait_general_gather(gmtag);
  wait_gather(mtag[6]);
  FORALLSITES(i, s) {
    mult_na(&(s->link[dir1]), (matrix *)(gen_pt[4][i]),
                &(s->s_link_f));
    add_matrix(&(s->s_link_f), (matrix *)(gen_pt[6][i]),
                   &(s->s_link_f));
  }
  cleanup_general_gather(gmtag);
  cleanup_gather(mtag[6]);
  mtag[5] = start_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                              dir3, EVENANDODD, gen_pt[5]);

  // Start one corner for third "plaquette" and gather
  FORALLSITES(i, s)
    mult_an(&(s->link[dir2]), &(s->link[dir3]), &(s->diag));

  for (i = XUP; i <= TUP; i++)
    disp[i] = 0;

  disp[dir2] = -1;
  disp[dir3] = 1;
  gmtag = start_general_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);
  mtag[3] = start_gather_site(F_OFFSET(diag), sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[3]);

  FORALLSITES(i, s)
    mat_copy(&(s->link[TUP]), &(s->t_link_f));

  // Make diagonal link in (-y, z) direction and gather it
  wait_general_gather(gmtag);
  wait_gather(mtag[3]);
  FORALLSITES(i, s) {
    mult_na(&(s->link[dir3]), (matrix *)(gen_pt[4][i]),
                &(s->staple));
    add_matrix(&(s->staple), (matrix *)(gen_pt[3][i]),
                   &(s->staple));
  }
  cleanup_general_gather(gmtag);
  cleanup_gather(mtag[3]);
  mtag[0] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                              dir1, EVENANDODD, gen_pt[0]);

  // Make second body diagonal and add to that gathered first
  wait_gather(mtag[1]);
  wait_gather(mtag[5]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir3]), (matrix *)(gen_pt[5][i]), &(s->diag));
    add_matrix(&(s->diag), (matrix *)(gen_pt[1][i]), &(s->diag));
  }
  cleanup_gather(mtag[1]);
  cleanup_gather(mtag[5]);

  // Make third body diagonal and add
  wait_gather(mtag[0]);
  dt = 1.0 / 6.0;
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir1]), (matrix *)(gen_pt[0][i]), &tmat1);
    add_matrix(&(s->diag), &tmat1, &(s->diag));
    scalar_mult_matrix(&(s->diag), dt, &(s->diag));
  }
  cleanup_gather(mtag[0]);

  // Start gather of temporal links across the body diagonal
  disp[dir1] = 1;
  disp[dir2] = -1;
  disp[dir3] = 1;
  disp[TUP] = 0;
  gmtag = start_general_gather_site(F_OFFSET(t_link_f), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);

  // Recursively construct the spatial segments
  // Compute the Wilson loops with that segment
  for (r = 0; r < nxh; r++) {
    if (r == 0) {
      FORALLSITES(i, s)
        mat_copy(&(s->diag), &(s->s_link));
    }
    else {
      wait_general_gather(gmtag);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

      cleanup_general_gather(gmtag);

      // Concurrently gather temporal links across the diagonal
      gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

      FORALLSITES(i, s)
        mult_nn(&(s->diag), &(s->staple), &(s->s_link));
    }
    FORALLSITES(i, s)
      mat_copy(&(s->s_link), &(s->s_link_f));

    // Start gather of forward spatial segments
    mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                  sizeof(matrix), TUP,
                                  EVENANDODD, gen_pt[TUP]);

    // Collect forward temporal links
    wait_general_gather(gmtag);
    FORALLSITES(i, s)
      mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

    FORALLSITES(i, s)
      mat_copy(&(s->staple), &(s->t_link_f));

    cleanup_general_gather(gmtag);

    // Concurrently gather spatial links across the diagonal for next r
    if (r < nxh - 1)
      gmtag = start_general_gather_site(F_OFFSET(s_link),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

    // Recursively compute the Wilson loops of different time extent
    for (t = 0; t < nth; t++) {
      // Collect forward spatial segments
      wait_gather(mtag[TUP]);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

      FORALLSITES(i, s)
        mat_copy(&(s->staple), &(s->s_link_f));

      // Start gather for next t, if still needed
      if (t < nth - 1)
        restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                            TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
      else
        cleanup_gather(mtag[TUP]);

      // Finally, compute the Wilson loops
      FORALLSITES(i, s) {
        if ((s->t) + t + 1 >= nt) {
          mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
          mult_na(&tmat1, &(s->t_link_f), &tmat2);
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&tmat2, &(s->s_link));
        }
        else
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&(s->s_link_f), &(s->s_link));
      }
    } // End loop over t
  } // End loop over r

  // Next, "sqrt(3)" loops in (x, -y, -z) space
  // Construct the "body diagonal" link in the (x, -y, -z) direction
  // Gather for first "plaquette"
  mtag[1] = start_gather_site(F_OFFSET(link[dir3]), sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[1]);
  mtag[2] = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                              dir3, EVENANDODD, gen_pt[2]);

  // Start one corner for second "plaquette" and gather
  FORALLSITES(i, s)
    mult_an(&(s->link[dir2]), &(s->link[dir1]), &(s->staple));

  for (i = XUP; i <= TUP; i++)
    disp[i] = 0;

  disp[dir1] = 1;
  disp[dir2] = -1;
  gmtag = start_general_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);
  mtag[6] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[6]);

  // Make diagonal link in (y, z) direction
  wait_gather(mtag[1]);
  wait_gather(mtag[2]);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[dir2]), (matrix *)(gen_pt[1][i]),
                &(s->s_link));
    mult_nn(&(s->link[dir3]), (matrix *)(gen_pt[2][i]), &tmat1);
    add_matrix(&(s->s_link), &tmat1, &(s->s_link));
  }
  cleanup_gather(mtag[1]);
  cleanup_gather(mtag[2]);

  // Make diagonal link in (x, -y) direction
  // Multiply to get second body diagonal and gather it
  wait_general_gather(gmtag);
  wait_gather(mtag[6]);
  FORALLSITES(i, s) {
    mult_na(&(s->link[dir1]), (matrix *)(gen_pt[4][i]), &tmat1);
    add_matrix(&tmat1, (matrix *)(gen_pt[6][i]), &tmat1);
    mult_an(&(s->link[dir3]), &tmat1, &(s->s_link_f));
  }
  cleanup_general_gather(gmtag);
  cleanup_gather(mtag[6]);
  mtag[2] = start_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                              OPP_DIR(dir3), EVENANDODD, gen_pt[2]);

  // Start one corner for third "plaquette" and gather
  FORALLSITES(i, s)
      mult_an(&(s->link[dir3]), &(s->link[dir1]), &(s->diag));

  for (i = XUP; i <= TUP; i++)
    disp[i] = 0;

  disp[dir1] = 1;
  disp[dir3] = -1;
  gmtag = start_general_gather_site(F_OFFSET(link[dir3]), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);
  mtag[3] = start_gather_site(F_OFFSET(diag), sizeof(matrix),
                              OPP_DIR(dir3), EVENANDODD, gen_pt[3]);

  FORALLSITES(i, s)
      mat_copy(&(s->link[TUP]), &(s->t_link_f));

  // Make diagonal link in (x, -z) direction
  // Multiply to get third body diagonal and gather it
  wait_general_gather(gmtag);
  wait_gather(mtag[3]);
  FORALLSITES(i, s) {
    mult_na(&(s->link[dir1]), (matrix *)(gen_pt[4][i]), &tmat1);
    add_matrix(&tmat1, (matrix *)(gen_pt[3][i]), &tmat1);
    mult_an(&(s->link[dir2]), &tmat1, &(s->staple));
  }
  cleanup_general_gather(gmtag);
  cleanup_gather(mtag[3]);
  mtag[1] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

  // Finally gather first "plaquette"
  disp[dir1] = 1;
  disp[dir2] = -1;
  disp[dir3] = -1;
  disp[TUP] = 0;
  gmtag = start_general_gather_site(F_OFFSET(s_link), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);

  // Collect second and third body diagonal and add them
  wait_gather(mtag[1]);
  wait_gather(mtag[2]);
  FORALLSITES(i, s)
    add_matrix((matrix *)(gen_pt[1][i]),
                   (matrix *)(gen_pt[2][i]), &(s->diag));

  cleanup_gather(mtag[1]);
  cleanup_gather(mtag[2]);

  // Make first body diagonal and add
  wait_general_gather(gmtag);
  dt = 1.0 / 6.0;
  FORALLSITES(i, s) {
    mult_na(&(s->link[dir1]), (matrix *)(gen_pt[4][i]), &tmat1);
    add_matrix(&(s->diag), &tmat1, &(s->diag));
    scalar_mult_matrix(&(s->diag), dt, &(s->diag));
  }
  cleanup_general_gather(gmtag);

  // Start gather of temporal links across the body diagonal
  gmtag = start_general_gather_site(F_OFFSET(t_link_f), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);

  // Recursively construct the spatial segments
  // Compute the Wilson loops with that segment
  for (r = 0; r < nxh; r++) {
    if (r == 0) {
      FORALLSITES(i, s)
        mat_copy(&(s->diag), &(s->s_link));
    }
    else {
      wait_general_gather(gmtag);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

      cleanup_general_gather(gmtag);

      // Concurrently gather temporal links across the diagonal
      gmtag = start_general_gather_site(F_OFFSET(t_link_f),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

      FORALLSITES(i, s)
        mult_nn(&(s->diag), &(s->staple), &(s->s_link));
    }
    FORALLSITES(i, s)
      mat_copy(&(s->s_link), &(s->s_link_f));

    // Start gather of forward spatial segments
    mtag[TUP] = start_gather_site(F_OFFSET(s_link_f),
                                  sizeof(matrix), TUP,
                                  EVENANDODD, gen_pt[TUP]);

    // Collect forward temporal links
    wait_general_gather(gmtag);
    FORALLSITES(i, s)
      mat_copy((matrix *)(gen_pt[4][i]), &(s->staple));

    FORALLSITES(i, s)
      mat_copy(&(s->staple), &(s->t_link_f));

    cleanup_general_gather(gmtag);

    // Concurrently gather spatial links across the diagonal for next r
    if (r < nxh - 1)
      gmtag = start_general_gather_site(F_OFFSET(s_link),
                                        sizeof(matrix), disp,
                                        EVENANDODD, gen_pt[4]);

    // Recursively compute the Wilson loops of different time extent
    for (t = 0; t < nth; t++) {
      // Collect forward spatial segments
      wait_gather(mtag[TUP]);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));

      FORALLSITES(i, s)
        mat_copy(&(s->staple), &(s->s_link_f));

      // Start gather for next t, if still needed
      if (t < nth - 1)
        restart_gather_site(F_OFFSET(s_link_f), sizeof(matrix),
                            TUP, EVENANDODD, gen_pt[TUP], mtag[TUP]);
      else
        cleanup_gather(mtag[TUP]);

      // Finally, compute the Wilson loops
      FORALLSITES(i, s) {
        if ((s->t) + t + 1 >= nt) {
          mult_nn(&(s->link[TUP]), &(s->s_link_f), &tmat1);
          mult_na(&tmat1, &(s->t_link_f), &tmat2);
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&tmat2, &(s->s_link));
        }
        else
          wils_loop2[r + r_off + nrmax * t]
            += (double)realtrace(&(s->s_link_f), &(s->s_link));
      }
    } // End loop over t
  } // End loop over r
  // End of "sqrt(3)" Wilson loops with dir1 = x, dir2 = y, dir3 = z

  // Normalize and print the Wilson loops
  for (t = 0; t < nth; t++) {
    for (r = 0; r < nxh; r++) {
      dt = wils_loop2[r + nrmax * t];
      g_doublesum(&dt);
      dt /= (double)(18 * volume);
      node0_printf("WILS_LOOP2_%d %d %d %.8g\n", tot_smear, r, t, dt);
    }

    r_off = nxh;
    for (r = 0; r < nxh / 2; r++) {
      dt = wils_loop2[r + r_off + nrmax * t];
      g_doublesum(&dt);
      dt /= (double)(36 * volume);
      node0_printf("WILS_LOOP2_%d %d %d %.8g\n", tot_smear, r + r_off, t, dt);
    }

    r_off = nxh + nxh / 2;
    for (r = 0; r < nxh; r++) {
      dt = wils_loop2[r + r_off+nrmax*t];
      g_doublesum(&dt);
      dt /= (double)(12 * volume);
      node0_printf("WILS_LOOP2_%d %d %d %.8g\n", tot_smear, r + r_off, t, dt);
    }
  }
  free(wils_loop2);
}
// -----------------------------------------------------------------

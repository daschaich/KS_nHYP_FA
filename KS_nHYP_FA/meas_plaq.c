// -----------------------------------------------------------------
// Measure the average plaquettes (space--space and space--time)
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void meas_plaq() {
  register int i, dir1, dir2;
  register site *s;
  register su3_matrix *m1, *m4;
  double plaq_ss = 0, plaq_e[5], plaq_o[5];    // st, x, y, z, a
  su3_matrix mtmp, *tmat;         // Scratch space
  msg_tag *mtag0,*mtag1;

  tmat = (su3_matrix *)malloc(sizeof(su3_matrix) * sites_on_node);
  if (tmat == NULL) {
    node0_printf("ERROR: can't malloc tmat in plaq_diff()\n");
    fflush(stdout);
    terminate(1);
  }

  for (i = 0; i < 5; i++) {
    plaq_e[i] = 0;
    plaq_o[i] = 0;
  }

  for (dir1 = YUP; dir1 <= TUP; dir1++) {
    for (dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(su3_matrix),
                                dir1, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir1]), sizeof(su3_matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->link[dir1]);
        m4 = &(s->link[dir2]);
        mult_su3_an(m4, m1, &tmat[i]);
      }
      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mult_su3_nn(&tmat[i], (su3_matrix *)(gen_pt[0][i]), &mtmp);

        if (dir1 == TUP) {
          if ((s->t) % 2 == 0)
            plaq_e[0] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
          else
            plaq_o[0] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
        }
        else    // Check
          plaq_ss += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                           &mtmp);

        if (dir1 == XUP || dir2 == XUP) {
          if ((s->x) % 2 == 0)
            plaq_e[1] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
          else
            plaq_o[1] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
        }

        if (dir1 == YUP || dir2 == YUP) {
          if ((s->y) % 2 == 0)
            plaq_e[2] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
          else
            plaq_o[2] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
        }

        if (dir1 == ZUP || dir2 == ZUP) {
          if ((s->z) % 2 == 0)
            plaq_e[3] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
          else
            plaq_o[3] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                               &mtmp);
        }

        // "a" is not really even/odd, but leave the labels for consistency
        if ((s->t) % 2 == 1 && (s->x) % 2 == 1
                            && (s->y) % 2 == 1 && (s->z) % 2 == 1)
          plaq_e[4] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                             &mtmp);
        if ((s->t) % 2 == 0 && (s->x) % 2 == 0
                            && (s->y) % 2 == 0 && (s->z) % 2 == 0)
          plaq_o[4] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                             &mtmp);
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    } // End loop over dir2
  } // End loop over dir1

  // Accumulate, normalize and print results
  g_doublesum(&plaq_ss);
  for (i = 0; i < 5; i++) {
    g_doublesum(&plaq_e[i]);
    g_doublesum(&plaq_o[i]);
  }
  node0_printf("Check plaquettes: %.8g %.8g\n",
               plaq_ss / ((double)(3 * volume)),
               (plaq_e[0] + plaq_o[0]) / ((double)(3 * volume)));
  for (i = 0; i < 5; i++) {
  node0_printf("StaggPlaq ");
    switch(i) {
      case 0: node0_printf("t "); break;
      case 1: node0_printf("x "); break;
      case 2: node0_printf("y "); break;
      case 3: node0_printf("z "); break;
      case 4: node0_printf("a "); break;
    }
    node0_printf("%.8g %.8g\n",
                 plaq_e[i] / ((double)(3 * volume)),
                 plaq_o[i] / ((double)(3 * volume)));
  }
  free(tmat);
}
// -----------------------------------------------------------------

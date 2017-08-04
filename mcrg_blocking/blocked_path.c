// -----------------------------------------------------------------
// An arbitrary path walker using ordinary gathers
#include "mcrg_includes.h"

// Put result in tempmat2
// Use tempmat for temporary storage
void blocked_path(int block, int *dir, int *sign, int length) {
  register int i;
  register site *s;
  int j, k, bl;
  msg_tag *mtag;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  if (sign[0] > 0) {
    FORALLSITES(i, s)
      su3mat_copy(&(s->link[dir[0]]), &(tempmat2[i]));

    // Shift block times
    for (k = 0; k < bl; k++) {
      mtag = start_gather_field(tempmat2, sizeof(su3_matrix),
                                OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat[i]));

      cleanup_gather(mtag);
      mtag = start_gather_field(tempmat, sizeof(su3_matrix),
                                OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat2[i]));

      cleanup_gather(mtag);
    }
  }

  if (sign[0] < 0) {
    FORALLSITES(i, s)
      su3_adjoint(&(s->link[dir[0]]), &(tempmat2[i]));
  }

  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {
      FORALLSITES(i, s)
        mult_su3_nn(&(tempmat2[i]), &(s->link[dir[j]]), &(tempmat[i]));

      // Shift block times
      for (k = 0; k < bl; k++) {
        mtag = start_gather_field(tempmat, sizeof(su3_matrix),
                                  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat2[i]));

        cleanup_gather(mtag);
        mtag = start_gather_field(tempmat2, sizeof(su3_matrix),
                                  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat[i]));

        cleanup_gather(mtag);
      }
      FORALLSITES(i, s)
        su3mat_copy(&(tempmat[i]), &(tempmat2[i]));
    }

    if (sign[j] < 0) {
      // Shift block times
      for (k = 0; k < bl; k++) {
        mtag = start_gather_field(tempmat2, sizeof(su3_matrix),
                                  dir[j], EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat[i]));

        cleanup_gather(mtag);
        mtag = start_gather_field(tempmat, sizeof(su3_matrix),
                                  dir[j], EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(tempmat2[i]));

        cleanup_gather(mtag);
      }
      FORALLSITES(i, s)
        su3mat_copy(&(tempmat2[i]), &(tempmat[i]));

      FORALLSITES(i, s)
        mult_su3_na(&(tempmat[i]), &(s->link[dir[j]]), &(tempmat2[i]));
    }
  }
}
// -----------------------------------------------------------------

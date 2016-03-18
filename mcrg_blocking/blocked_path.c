// -----------------------------------------------------------------
// An arbitrary path walker using ordinary gathers
#include "mcrg_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Result is reported in resmat
void blocked_path(int block, int *dir, int *sign,
                  int length, su3_matrix *resmat) {

  register int i;
  register site *s;
  int j, k, bl;
  msg_tag *mtag;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  if (sign[0] > 0) {
    FORALLSITES(i, s)
      su3mat_copy(&(s->link[dir[0]]), &(s->tempmat2));

    // Shift block times
    for (k = 0; k < bl; k++) {
      mtag = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix),
                               OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->tempmat1));

      cleanup_gather(mtag);
      mtag = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix),
                               OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->tempmat2));

      cleanup_gather(mtag);
    }
  }

  if (sign[0] < 0) {
    FORALLSITES(i, s)
      su3_adjoint(&(s->link[dir[0]]), &(s->tempmat2));
  }

  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {
      FORALLSITES(i, s)
        mult_su3_nn(&(s->tempmat2), &(s->link[dir[j]]), &(s->tempmat1));

      // Shift block times
      for (k = 0; k < bl; k++) {
        mtag = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix),
                                 OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->tempmat2));

        cleanup_gather(mtag);
        mtag = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix),
                                 OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->tempmat1));

        cleanup_gather(mtag);
      }
      FORALLSITES(i, s)
        su3mat_copy(&(s->tempmat1), &(s->tempmat2));
    }

    if (sign[j] < 0) {
      // Shift block times
      for (k = 0; k < bl; k++) {
        mtag = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix),
                                 dir[j], EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->tempmat1));

        cleanup_gather(mtag);
        mtag = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix),
                                 dir[j], EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          su3mat_copy((su3_matrix *)(gen_pt[0][i]), &(s->tempmat2));

        cleanup_gather(mtag);
      }
      FORALLSITES(i, s)
        su3mat_copy(&(s->tempmat2), &(s->tempmat1));

      FORALLSITES(i, s)
        mult_su3_na(&(s->tempmat1), &(s->link[dir[j]]), &(s->tempmat2));
    }
  }

  // Copy result into matrix to be returned
  FORALLSITES(i, s)
    su3mat_copy(&(s->tempmat2), &resmat[i]);
}
// -----------------------------------------------------------------

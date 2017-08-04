// -----------------------------------------------------------------
// An arbitrary path walker using ordinary gathers
#include "wflow_includes.h"

// Put result in tempmat2
// Use tempmat for temporary storage
void path(int *dir, int *sign, int length) {
  register int i;
  register site *s;
  int j;
  msg_tag *mtag;

  if (sign[0] > 0)  {
    mtag = start_gather_site(F_OFFSET(link[dir[0]]), sizeof(matrix),
                             OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
    wait_gather(mtag);
    FORALLSITES(i, s)
      mat_copy((matrix *)(gen_pt[0][i]), &(tempmat2[i]));

    cleanup_gather(mtag);
  }

  if (sign[0] < 0) {
    FORALLSITES(i, s)
      su3_adjoint(&(s->link[dir[0]]), &(tempmat2[i]));
  }

  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {
      FORALLSITES(i, s)
        mult_su3_nn(&(tempmat2[i]), &(s->link[dir[j]]), &(tempmat[i]));

      mtag = start_gather_field(tempmat, sizeof(matrix),
                                OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[0][i]), &(tempmat2[i]));

      cleanup_gather(mtag);
    }

    if (sign[j] < 0) {
      mtag = start_gather_field(tempmat2, sizeof(matrix),
                                dir[j], EVENANDODD, gen_pt[1]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        mult_su3_na((matrix *)(gen_pt[1][i]),
                    &(s->link[dir[j]]), &(tempmat[i]));

      FORALLSITES(i, s)
        mat_copy(&(tempmat[i]), &(tempmat2[i]));

      cleanup_gather(mtag);
    }
  }
}
// -----------------------------------------------------------------

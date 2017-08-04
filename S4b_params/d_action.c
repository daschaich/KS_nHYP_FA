// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the conjugate gradient should already
// have been run on both even and odd sites, depending on Nf
// The vectors psi1 and psi2 contain (M^dag.M)^(-1) chi

#include "S4b_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Adds adjoint plaquette term
void d_plaquette_a(double *ss_plaq, double *st_plaq) {
  register int i, dir1, dir2;
  register site *s;
  register matrix *m1, *m4;
  matrix mtmp, *mat;   // mat is scratch space
  double ss_sum, st_sum;
  complex ctt;
  msg_tag *mtag0, *mtag1;

  ss_sum = 0;
  st_sum = 0;
  mat = (matrix *)malloc(sizeof(matrix) * sites_on_node);
  if (mat == NULL) {
    printf("d_plaquette_a: can't malloc mat\n");
    fflush(stdout);
    terminate(1);
  }

  for(dir1 = YUP; dir1 <= TUP; dir1++) {
    for(dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                dir1, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir1]), sizeof(matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->link[dir1]);
        m4 = &(s->link[dir2]);
        mult_su3_an(m4, m1, &mat[i]);
      }
      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mult_su3_nn(&mat[i], (matrix *)(gen_pt[0][i]), &mtmp);

        if (dir1 == TUP) {
          ctt = complextrace_su3((matrix *)(gen_pt[1][i]), &mtmp);
          st_sum += (ctt.real / 3)
                  + beta_a * (ctt.real * ctt.real + ctt.imag * ctt.imag) / 9;
        }
        else {
          ctt = complextrace_su3((matrix *)(gen_pt[1][i]), &mtmp);
          ss_sum += (ctt.real / 3)
                  + beta_a * (ctt.real * ctt.real + ctt.imag * ctt.imag) / 9;
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  *ss_plaq = ss_sum / ((double)(nx * ny * nz * nt));
  *st_plaq = st_sum / ((double)(nx * ny * nz * nt));
  free(mat);
}
// -----------------------------------------------------------------

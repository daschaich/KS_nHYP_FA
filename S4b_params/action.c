// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the conjugate gradient should already
// have been run on both even and odd sites, depending on Nf
// The vectors psi1 and psi2 contain (M^dag.M)^(-1) chi
#include "S4b_includes.h"

// Add adjoint plaquette term
// Use tempmat for temporary storage
void plaquette_a(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  register matrix *m1, *m4;
  double ss_sum = 0.0, st_sum = 0.0;
  complex tc;
  msg_tag *mtag, *mtag2;
  matrix tmat;

  // We can exploit a symmetry under dir<-->dir2
  for (dir = YUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
      mtag = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                dir, EVENANDODD, gen_pt[0]);
      mtag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      // tempmat = Udag_b(x) U_a(x)
      FORALLSITES(i, s) {
        m1 = &(s->link[dir]);
        m4 = &(s->link[dir2]);
        mult_an(m4, m1, &(tempmat[i]));
      }
      wait_gather(mtag);
      wait_gather(mtag2);

      // Compute tc = tr[Udag_a(x+b) Udag_b(x) U_a(x) U_b(x+a)] / 3
      // and combine tc.real + beta_A * |tc|^2
      FORALLSITES(i, s) {
        m1 = (matrix *)(gen_pt[0][i]);
        m4 = (matrix *)(gen_pt[1][i]);
        mult_nn(&(tempmat[i]), m1, &tmat);
        tc = complextrace(m4, &tmat);

        if (dir == TUP)
          st_sum += tc.real / 3.0 + beta_a * cabs_sq(&tc) / 9.0;
        else
          ss_sum += tc.real / 3.0 + beta_a * cabs_sq(&tc) / 9.0;
      }
      cleanup_gather(mtag);
      cleanup_gather(mtag2);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  *ss_plaq = ss_sum / ((double)volume);
  *st_plaq = st_sum / ((double)volume);
}
// -----------------------------------------------------------------

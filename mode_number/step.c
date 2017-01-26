// -----------------------------------------------------------------
// Step function implemented through Clenshaw algorithm
#include "mode_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest = src - 8M^2 (DD^dag + 4M^2)^{-1} src
// Hard-code EVEN parity in inverse
void X(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  Real MSq_x2 = -8 * M * M / starSq;

  // psi = (DD^dag + 4M^2)^{-1} src
  clear_latvec(F_OFFSET(psi), EVENANDODD);    // Zero initial guess
#ifdef CG_DEBUG
  int iters = ks_congrad(src, F_OFFSET(psi), M / star, EVEN);
  node0_printf("%d iters in congrad with M=%.4g\n", iters, M / star);
#else
  ks_congrad(src, F_OFFSET(psi), M / star, EVEN);
#endif

  // dest = src - 8M^2 psi
  FOREVENSITES(i, s)
    scalar_mult_add_su3_vector((su3_vector *)F_PT(s, src), &(s->psi),
                               MSq_x2, (su3_vector *)F_PT(s, dest));
}

// dest = (2X^2 - 1 - epsilon) src / (1 - epsilon)
// Hard-code EVEN parity
void Z(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  double toAdd = -1.0 - epsilon;
  double norm = 1.0 / (1.0 - epsilon);

  X(src, F_OFFSET(R1));
  X(F_OFFSET(R1), F_OFFSET(Xsrc));

  FOREVENSITES(i, s) {
    scalar_mult_su3_vector(&(s->Xsrc), 2.0, &(s->Xsrc));
    scalar_mult_add_su3_vector(&(s->Xsrc), (su3_vector *)F_PT(s, src),
                               toAdd, (su3_vector *)F_PT(s, dest));
    scalar_mult_su3_vector((su3_vector *)F_PT(s, dest), norm,
                           (su3_vector *)F_PT(s, dest));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clenshaw algorithm:
// P(x)R(0) = \sum_i^n c[i]T[i]R(0) = (b[0] - xb[1])R(0),
// where b[i] = c[i] + 2zb[i + 1] - b[i + 2], b[n] = b[n + 1] = 0
// So want to compute b[0] - xb[1];
// Hard-code EVEN parity
void clenshaw(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  int j;

  for (j = Norder; j >= 0; j--) {
    // Construct bj.src = (cj + 2Zbjp1 - bjp2).src
    // Start with bj.src = cj.src
    FOREVENSITES(i, s)
      scalar_mult_su3_vector((su3_vector *)F_PT(s, src), coeffs[j], &(s->bj));

    // Now subtract bjp2.src calculated in previous iterations
    if (j < Norder - 1) {
      FOREVENSITES(i, s)
        sub_su3_vector(&(s->bj), &(s->bjp2), &(s->bj));
    }

    // Finally we need 2Z(bjp1.src)
    // Based on bjp1.src calculated in previous iterations
    if (j < Norder) {
      Z(F_OFFSET(bjp1), F_OFFSET(Zbjp1));
      FOREVENSITES(i, s)
        scalar_mult_add_su3_vector(&(s->bj), &(s->Zbjp1), 2.0, &(s->bj));
    }

    // Now move bjp1-->bjp2 and bj-->bjp1 for next iteration
    if (j > 0) {
      copy_latvec(F_OFFSET(bjp1), F_OFFSET(bjp2), EVENANDODD);
      copy_latvec(F_OFFSET(bj), F_OFFSET(bjp1), EVENANDODD);
    }
  }

  // We now have bj = b[0].src and Zbjp1 = Z(b[1].src)
  // We want to return (b[0].src - Z(b[1].src))
  FOREVENSITES(i, s)
    sub_su3_vector(&(s->bj), &(s->Zbjp1), (su3_vector *)F_PT(s, dest));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Step function approximated by h(x) = [1 - xp(x)^2] / 2
// Hard-code EVEN parity
void step(field_offset src, field_offset dest) {
  register int i;
  register site *s;

  // dest = P(X^2) src temporarily
  clenshaw(src, dest);

  // dest = (src - X P(X^2) src) / 2
  X(dest, F_OFFSET(Xsrc));
  FOREVENSITES(i, s) {
    sub_su3_vector((su3_vector *)F_PT(s, src), &(s->Xsrc),
                   (su3_vector *)F_PT(s, dest));
    scalar_mult_su3_vector((su3_vector *)F_PT(s, dest), 0.5,
                           (su3_vector *)F_PT(s, dest));
  }
}
// -----------------------------------------------------------------

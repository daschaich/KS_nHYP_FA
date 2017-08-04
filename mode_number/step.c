// -----------------------------------------------------------------
// Step function implemented through Clenshaw algorithm
#include "mode_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest = src - 8M^2 (DD^dag + 4M^2)^(-1) src
// Hard-code EVEN parity in inverse
void X(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  Real MSq_x2 = -8.0 * M * M / starSq;

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
    scalar_mult_add_vector((vector *)F_PT(s, src), &(s->psi),
                               MSq_x2, (vector *)F_PT(s, dest));
}

// dest = (2X^2 - 1 - epsilon) src / (1 - epsilon)
// Hard-code EVEN parity
void Z(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  Real scale = 2.0 / (1.0 - epsilon);
  Real toAdd = (-1.0 - epsilon) / (1.0 - epsilon);

  // This is more compact, but the subtraction in X(src)
  // seems to leave it relatively poorly conditioned,
  // increasing CG iterations by almost 10% even for a small 4nt4 test
  // (compared to inverting twice on src itself)
  X(src, F_OFFSET(Xsrc));
  X(F_OFFSET(Xsrc), dest);

  FOREVENSITES(i, s) {
    scalar_mult_vector((vector *)F_PT(s, dest), scale,
                       (vector *)F_PT(s, dest));
    scalar_mult_sum_vector((vector *)F_PT(s, src), toAdd,
                           (vector *)F_PT(s, dest));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clenshaw algorithm:
// P(X) src = sum_i^n c[i] T[i] src = (b[0] - X b[1]) src,
// where b[i] = c[i] + 2zb[i + 1] - b[i + 2], b[n] = b[n + 1] = 0
// Hard-code EVEN parity
void clenshaw(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  int j;

  for (j = Norder; j >= 0; j--) {
    // Construct bj.src = (cj + 2Zbjp1 - bjp2).src
    // Start with bj.src = cj.src
    FOREVENSITES(i, s)
      scalar_mult_vector((vector *)F_PT(s, src), coeffs[j], &(s->bj));

    // Now subtract bjp2.src calculated in previous iterations
    if (j < Norder - 1) {
      FOREVENSITES(i, s)
        sub_vector(&(s->bj), &(s->bjp2), &(s->bj));
    }

    // Finally we need 2Z(bjp1.src)
    // Based on bjp1.src calculated in previous iterations
    if (j < Norder) {
      Z(F_OFFSET(bjp1), F_OFFSET(Zbjp1));
      FOREVENSITES(i, s)
        scalar_mult_add_vector(&(s->bj), &(s->Zbjp1), 2.0, &(s->bj));
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
    sub_vector(&(s->bj), &(s->Zbjp1), (vector *)F_PT(s, dest));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Step function approximated by h(X) = [1 - X p(X)^2] / 2
// Hard-code EVEN parity
void step(field_offset src, field_offset dest) {
  register int i;
  register site *s;

  // dest = P(X^2) src temporarily
  clenshaw(src, dest);

  // dest = (src - X P(X^2) src) / 2
  X(dest, F_OFFSET(Xsrc));
  FOREVENSITES(i, s) {
    sub_vector((vector *)F_PT(s, src), &(s->Xsrc),
                   (vector *)F_PT(s, dest));
    scalar_mult_vector((vector *)F_PT(s, dest), 0.5,
                           (vector *)F_PT(s, dest));
  }
}
// -----------------------------------------------------------------

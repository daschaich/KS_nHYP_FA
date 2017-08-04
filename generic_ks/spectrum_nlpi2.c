// -----------------------------------------------------------------
// Spectrum for lots of Kogut--Susskind pion (and some rho) operators
// Use s->chi for source, s->psi for temporary in mat_invert_uml
// Use s->g_rand, s->prop and s->aprop for temporary storage
// s->tempvec is used by ks_congrad through mat_invert_uml
// tempvec[0] and tempvec[1] are used by some of the operators
// tempvec[0] is also used by zeta_shift, ttt by sym_shift

// Call with rephase(ON) to absorb (anti-periodic) boundary conditions
// into the link matrices

// Coulomb gauge should be fixed before calling this routine

// Build shorthand for our twelve propagators summarized below
// We only consider 'diagonal' propagators with the same source and sink
// Name        summary    KS trans          partner  phase
// pion5       local 0-+  (flavor) gamma_5  0+-      1
// pion05      local 0-+   gamma_0 gamma_5  0++    (-1)^(x + y + z + t)
// pioni5      1link 0-+   gamma_i gamma_5  0+-
// pionij      1link 0-+   gamma_i gamma_j  0++
// pioni       2link 0-+   gamma_i          0++
// pioni0      2link 0-+   gamma_i gamma_0  0+-
// pions       3link 0-+   singlet          0++
// pion0       3link 0-+   gamma_0          0+-
// rhoi (VT)   local 1--   gamma_i          1+-    (-1)^(dir)
// rhoi0 (PV)  local 1--   gamma_i gamma_0  1++    (-1)^(x + y + z + dir)
// rhos        1link 1--   singlet          1+-
// rho0        1link 1--   gamma_0          1++
enum prop_name {
  prop_pion5,
  prop_pion05,
  prop_pioni5,
  prop_pionij,
  prop_pioni,
  prop_pioni0,
  prop_pions,
  prop_pion0,
  prop_rhoi,
  prop_rhoi0,
  prop_rhos,
  prop_rho0,
  nprops
};

#include "generic_ks_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Twelve meson operators
// All contain an extra gamma_5 = (-1)^(x + y + z + t)
// from computing both fermion and anti-fermion props from the same source

// The usual (Goldstone) local pion operator (gamma_5)
// gamma_5 times the gamma_5 explained above leaves us with nothing to do
void mult_pion5(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  FORALLSITES(i, s)
    *(vector *)F_PT(s, dest) = *(vector *)F_PT(s, src);
}

// The second (pi_05) local pion operator (gamma_0 gamma_5)
// gamma_0 = 1 times the extra gamma_5 = gamma_5 = (-1)^(x + y + z + t)
void mult_pion05(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  FORALLSITES(i, s) {
    if (s->parity == EVEN)
      *(vector *)F_PT(s, dest) = *(vector *)F_PT(s, src);
    else {
      scalar_mult_vector((vector *)F_PT(s, src), -1.0,
                             (vector *)F_PT(s, dest));
    }
  }
}

// The pi_i5 one-link pion operator (gamma_5 x gamma_i gamma_5)
void mult_pioni5(int fdir, field_offset src, field_offset dest) {
  register int i;
  register site *s;
  Real sign;

  // Apply the symmetric shift operator
  sym_shift(fdir, src, dest);
  FORALLSITES(i, s) {
    sign = 1.0;
    if (((((short *)&(s->x))[fdir]) & 0x1) == 1)
      sign = -1.0;
    if (s->parity == ODD)
      sign = -sign;     // The extra gamma_5
    scalar_mult_vector((vector *)F_PT(s, dest), sign,
                           (vector *)F_PT(s, dest));
  }
}

// The pi_ij one-link pion operator (gamma_0 gamma_5 x gamma_i gamma_j)
void mult_pionij(int fdir, field_offset src, field_offset dest) {
  register int i;
  register site *s;
  Real sign;

  // Apply the symmetric shift operator
  sym_shift(fdir, src, dest);
  FORALLSITES(i, s) {
    sign = 1.0;
    if (((((short *)&(s->x))[fdir]) & 0x1) == 1)
      sign = -1.0;
    // The extra gamma_5 leaves us with nothing more to do
    scalar_mult_vector((vector *)F_PT(s, dest), sign,
                           (vector *)F_PT(s, dest));
  }
}

// The pi_i two-link pion operator (gamma_i)
// epsilon_sign normalized to average the (1, 2) and (2, 1) paths
void mult_pioni(int fdir, field_offset src, field_offset dest) {
  register int i, d1, d2;
  register site *s;
  short dir[2];
  Real sign, epsilon_sign = 0.5;

  switch(fdir) {
    case 0:  dir[0] = 1; dir[1] = 2; break;
    case 1:  dir[0] = 2; dir[1] = 0; break;
    case 2:  dir[0] = 0; dir[1] = 1; break;
    default: node0_printf("ERROR! invalid direction %i\n", fdir);
  }

  FORALLSITES(i, s)
    clearvec((vector *)F_PT(s, dest));

  for (d1 = 0, d2 = 1; d1 < 2; d1++, d2--) {
//    printf("In mult_pioni: (d1, d2): (%i,%i)\n", dir[d1], dir[d2]);

    // Apply the symmetric shift operator in dir2
    sym_shift(dir[d2], src, F_OFFSET(tempvec[0]));

    // Multiply by zeta_dir2
    // Because the phases are on,
    // multiply by zeta * eta = epsilon * (-1)^coord[dir2]
    FORALLSITES(i, s) {
      if (s->parity == EVEN)
        sign = 1.0;
      else
        sign = -1.0;
      if (((((short *)&(s->x))[d2]) & 0x1) == 1)
        sign = -sign;
      scalar_mult_vector(&(s->tempvec[0]), sign, &(s->tempvec[0]));
    }

    // Apply the symmetric shift operator in dir1
    sym_shift(dir[d1], F_OFFSET(tempvec[0]), F_OFFSET(tempvec[1]));

    // Multiply by zeta_dir1
    // Because the phases are on,
    // multiply by zeta * eta * epsilon = (-1)^coord[dir1]
    // The extra epsilon comes from the anti-fermion propagaton
    FORALLSITES(i, s) {
      sign = 1.0;
      if (((((short *)&(s->x))[d1]) & 0x1) == 1)
        sign = -1.0;
      sign *= epsilon_sign;
      scalar_mult_sum_vector(&(s->tempvec[1]), sign, (vector *)F_PT(s, dest));
    }
    epsilon_sign = -epsilon_sign;
  }
}

// The pi_i0 two-link pion operator (gamma_i gamma_0)
// We can call mult_pioni and hit it with zeta_4 * eta_4 = (-1)^(x + y + z)
void mult_pioni0(int fdir, field_offset src, field_offset dest) {
  register int i;
  register site *s;

  mult_pioni(fdir, src, dest);
  FORALLSITES(i, s) {
    if ((s->x + s->y + s->z) % 2 == 1) {
      scalar_mult_vector((vector *)F_PT(s, dest), -1.0,
                         (vector *)F_PT(s, dest));
    }
  }
}

// The singlet three-link pion operator
void mult_pions(field_offset src, field_offset dest) {
  register int i;
  register site *s;
  int c;
  Real norm = 1.0 / 6.0;    // Account for permutation multiplicity

  // Collect all permutations, with the appropriate sign
  struct {
    int d[3];
    Real sign;
  } p[6] = {{{0, 1, 2},  norm}, {{1, 2, 0},  norm}, {{2, 0, 1},  norm},
            {{0, 2, 1}, -norm}, {{1, 0, 2}, -norm}, {{2, 1, 0}, -norm}};

  FORALLSITES(i, s)
    clearvec((vector *)F_PT(s, dest));

  for (c = 0; c < 6; c++) {
    zeta_shift(3, p[c].d, src, F_OFFSET(tempvec[1]));
    FORALLSITES(i, s) {
      scalar_mult_sum_vector(&(s->tempvec[1]), p[c].sign,
                             (vector *)F_PT(s, dest));
    }
  }

  // Multiply by epsilon for the anti-quark
  FORODDSITES(i, s)
    scalar_mult_vector((vector *)F_PT(s, dest), -1.0, (vector *)F_PT(s, dest));
}

// The three-link gamma_0 pion operator
// We can call mult_pions and hit it with zeta_4 * eta_4 = (-1)^(x + y + z)
void mult_pion0(field_offset src, field_offset dest) {
  register int i;
  register site *s;

  mult_pions(src, dest);
  FORALLSITES(i, s) {
    if ((s->x + s->y + s->z) % 2 == 1) {
      scalar_mult_vector((vector *)F_PT(s, dest), -1.0,
                         (vector *)F_PT(s, dest));
    }
  }
}

// The VT local rho operator (gamma_pdir) plus the extra gamma_5
// End up with phase of just (-1)^(pdir)
void mult_rhoi(int pdir, field_offset src, field_offset dest) {
  register int i;
  register site *s;
  FORALLSITES(i, s) {
    if (((((short *)&(s->x))[pdir]) & 0x1) == 0)
      *(vector *)F_PT(s, dest) = *(vector *)F_PT(s, src);
    else {
      scalar_mult_vector((vector *)F_PT(s, src), -1.0,
                         (vector *)F_PT(s, dest));
    }
  }
}

// The PV local rho operator (gamma_0 gamma_pdir) plus the extra gamma_5
// End up with phase of (-1)^(x + y + z + t + pdir)
void mult_rhoi0(int pdir, field_offset src, field_offset dest) {
  register int i;
  register site *s;
  FORALLSITES(i, s) {
    if (((((short *)&(s->x))[pdir] + s->x + s->y + s->z) & 0x1) == 0)
      *(vector *)F_PT(s, dest) = *(vector *)F_PT(s, src);
    else {
      scalar_mult_vector((vector *)F_PT(s, src), -1.0,
                         (vector *)F_PT(s, dest));
    }
  }
}

// The singlet one-link rho operator
// Because the phases are on, only need to apply the extra gamma_5
void mult_rhos(int fdir, field_offset src, field_offset dest) {
  register int i;
  register site *s;

  // Apply the symmetric shift operator
  sym_shift(fdir, src, dest);
  FORALLSITES(i, s) {
    if (s->parity == ODD) {   // The extra gamma_5
      scalar_mult_vector((vector *)F_PT(s, dest), -1.0,
                         (vector *)F_PT(s, dest));
    }
  }
}

// The one-link gamma_0 rho operator
// Because the phases are on,
// only need to apply gamma_5 times (-1)^(x + y + z) = (-1)^t
void mult_rho0(int fdir, field_offset src, field_offset dest) {
  register int i;
  register site *s;

  // Apply the symmetric shift operator
  sym_shift(fdir, src, dest);
  FORALLSITES(i, s) {
    if ((s->t) % 2 == 1) {
      scalar_mult_vector((vector *)F_PT(s, dest), -1.0,
                         (vector *)F_PT(s, dest));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// fmass and amass are masses for the fermion and the anti-fermion
// n_src, src_start and src_inc need to be set by each application's setup.c
// Return the CG iteration count
int spectrum_nlpi2(Real fmass, Real amass) {
  register int i, j, t_src;
  register complex tc = cmplx(0.0, 0.0);
  register site* s;
  int icol, isrc, cgn = 0;
  Real norm = 1.0 / (nx * ny * nz * n_src);
  complex **props = malloc(nprops * sizeof(**props));

  // Set up and clear arrays to accumulate propagators
  for (i = 0; i < nprops; i++) {
    props[i] = malloc(nt * sizeof(complex));
    for (j = 0; j < nt; j++)
      props[i][j] = tc;
  }

  // Loop over sources
  for (isrc = 0, t_src = src_start;
       t_src < nt && isrc < n_src;
       t_src += src_inc, isrc++) {

    node0_printf("spectrum_nlpi2: t_src = %d\n", t_src);

    // Set up wall source in chi
    for (icol = 0; icol < 3; icol++) {
      FORALLSITES(i, s) {
        clearvec(&(s->chi));
        if (s->t == t_src)
          s->chi.c[icol].real = 1.0;
      }

      // Compute prop <-- M^(-1) chi
      cgn += mat_invert_uml(F_OFFSET(chi), F_OFFSET(prop),
                            F_OFFSET(psi), fmass);

      // Make source g_rand for aprop
      // by summing over desired operators times quark source
      // First source couples to pion5, pioni5, pioni, pions, rhoi, rhos
      mult_pion5(F_OFFSET(chi), F_OFFSET(g_rand));
      mult_pioni(ZUP, F_OFFSET(chi), F_OFFSET(aprop));
      FORALLSITES(i, s)
        sum_vector(&(s->aprop), &(s->g_rand));

      mult_rhoi(ZUP, F_OFFSET(chi), F_OFFSET(aprop));
      FORALLSITES(i, s)
        sum_vector(&(s->aprop), &(s->g_rand));

      // Compute aprop <-- M^(-1) g_rand
      cgn += mat_invert_uml(F_OFFSET(g_rand), F_OFFSET(aprop),
                            F_OFFSET(psi), amass);

      // Loop over all desired sink operators
      // Tie propagators together at sink end to project out desired mesons
      // Add to propagator arrays at distance from source
      mult_pion5(F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pion5][(s->t + nt - t_src) % nt], tc);
      }

      mult_pioni5(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pioni5][(s->t + nt - t_src) % nt], tc);
      }

      mult_pioni(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pioni][(s->t + nt - t_src) % nt], tc);
      }

      mult_pions(F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pions][(s->t + nt - t_src) % nt], tc);
      }

      // Two rho correlators
      mult_rhoi(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_rhoi][(s->t + nt - t_src) % nt], tc);
      }

      mult_rhos(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_rhos][(s->t + nt - t_src) % nt], tc);
      }

      // Make source g_rand for aprop
      // by summing over desired operators times quark source
      // Second source couples to pion05, pionij, pioni0, pion0, rhoi0, rho0
      mult_pion05(F_OFFSET(chi), F_OFFSET(g_rand));
      mult_pioni0(ZUP, F_OFFSET(chi), F_OFFSET(aprop));
      FORALLSITES(i, s)
        sum_vector(&(s->aprop), &(s->g_rand));

      mult_rhoi0(ZUP, F_OFFSET(chi), F_OFFSET(aprop));
      FORALLSITES(i, s)
        sum_vector(&(s->aprop), &(s->g_rand));

      // Compute aprop <-- M^(-1) g_rand
      cgn += mat_invert_uml(F_OFFSET(g_rand), F_OFFSET(aprop),
                            F_OFFSET(psi), amass);

      mult_pion05(F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pion05][(s->t + nt - t_src) % nt], tc);
      }

      mult_pionij(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pionij][(s->t + nt - t_src) % nt], tc);
      }

      mult_pioni0(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pioni0][(s->t + nt - t_src) % nt], tc);
      }

      mult_pion0(F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_pion0][(s->t + nt - t_src) % nt], tc);
      }

      mult_rhoi0(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_rhoi0][(s->t + nt - t_src) % nt], tc);
      }

      mult_rho0(ZUP, F_OFFSET(prop), F_OFFSET(g_rand));
      FORALLSITES(i, s) {
        tc = su3_dot(&(s->aprop), &(s->g_rand));
        CSUM(props[prop_rho0][(s->t + nt - t_src) % nt], tc);
      }
    }
  }

  // Sum propagator arrays over nodes and normalize
  for (i = 0; i < nprops; i++) {
    g_veccomplexsum(props[i], nt);
    for (j = 0; j < nt; j++)
      CMULREAL(props[i][j], norm, props[i][j]);
  }

  // Dump the (complex) propagators
  if (this_node == 0) {
    // First source couples to pion5, pioni5, pioni, pions, rhoi, rhos
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n", fmass, amass);
    printf("SOURCE: FUNNYWALL1\n");
    printf("SINKS: PION_5 PION_i5 PION_i PION_s RHO_i RHO_s\n");
    for (j = 0; j < nt; j++) {
      printf("%d %e %e %e %e %e %e %e %e %e %e %e %e\n",j,
             props[prop_pion5][j].real, props[prop_pion5][j].imag,
             props[prop_pioni5][j].real, props[prop_pioni5][j].imag,
             props[prop_pioni][j].real, props[prop_pioni][j].imag,
             props[prop_pions][j].real, props[prop_pions][j].imag,
             props[prop_rhoi][j].real, props[prop_rhoi][j].imag,
             props[prop_rhos][j].real, props[prop_rhos][j].imag);
    }
    printf("ENDPROP\n");

    // Second source couples to pion05, pionij, pioni0, pion0, rhoi0, rho0
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n", fmass, amass);
    printf("SOURCE: FUNNYWALL2\n");
    printf("SINKS: PION_05 PION_ij PION_i0 PION_0 RHO_i0 RHO_0\n");
    for (j = 0; j < nt; j++) {
      printf("%d %e %e %e %e %e %e %e %e %e %e %e %e\n",j,
             props[prop_pion05][j].real, props[prop_pion05][j].imag,
             props[prop_pionij][j].real, props[prop_pionij][j].imag,
             props[prop_pioni0][j].real, props[prop_pioni0][j].imag,
             props[prop_pion0][j].real, props[prop_pion0][j].imag,
             props[prop_rhoi0][j].real, props[prop_rhoi0][j].imag,
             props[prop_rho0][j].real, props[prop_rho0][j].imag);
    }
    printf("ENDPROP\n");
    fflush(stdout);
  }

  // Free temporaries
  for (i = 0; i < nprops; i++)
    free(props[i]);
  free(props);
  return cgn;
}
// -----------------------------------------------------------------

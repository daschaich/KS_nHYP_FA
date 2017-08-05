// -----------------------------------------------------------------
// Measure fermionic observables for general "even plus odd" quark actions:
//   pbp separately on even and odd sites
//   fermion action
// ferm_links_t, DM_DU0 and CHEM_POT stuff stripped out

// Stuck CG helpers in this directory
#include "generic_ks_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Modified to return total number of iterations
int f_meas_imp(field_offset chi_off, field_offset psi_off, Real mass) {
  register int i;
  register site *s;
  int jpbp_reps, tot_iters = 0, miters, npbp_reps = 1;
  Real r_pbp_even, i_pbp_even, r_pbp_odd, i_pbp_odd, r_ferm_action;
  double rfaction;
  complex tc;
  double_complex pbp_e, pbp_o;

#ifdef NPBP_REPS
  double pbp_pbp;
  npbp_reps = NPBP_REPS;    // Number of stochastic estimations
#endif

  for (jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++) {
    rfaction = 0;
    pbp_e = dcmplx(0, 0);
    pbp_o = dcmplx(0, 0);

    // Make random source and do inversion
    // Generate (one-mass) g_rand; chi_off = M g_rand
    grsource_imp(chi_off, mass, EVENANDODD);
    // chi_off = M g_rand (still)
    // psi_off = M^{-1} g_rand
    clear_latvec(psi_off, EVENANDODD);
    miters = mat_invert_uml(F_OFFSET(g_rand), psi_off, chi_off, mass);
    tot_iters += miters;

    // Fermion action = chi.psi
    // pbp on even sites = g_rand.psi
    FOREVENSITES(i, s) {
      tc = su3_dot((vector *)F_PT(s, chi_off),
                   (vector *)F_PT(s, psi_off));
      rfaction += tc.real;
      tc = su3_dot(&(s->g_rand), (vector *)F_PT(s, psi_off));
      CSUM(pbp_e, tc);
    }

    // pbp on odd sites
    FORODDSITES(i, s) {
      tc = su3_dot(&(s->g_rand), (vector *)F_PT(s, psi_off));
      CSUM(pbp_o, tc);
    }
    g_dcomplexsum(&pbp_o);
    g_dcomplexsum(&pbp_e);
    g_doublesum(&rfaction);

    r_pbp_odd = pbp_o.real * (2 / (double)volume);
    i_pbp_odd = pbp_o.imag * (2 / (double)volume);
    r_pbp_even = pbp_e.real * (2 / (double)volume);
    i_pbp_even = pbp_e.imag * (2 / (double)volume);
    r_ferm_action = rfaction * (1 / (double)volume);
    node0_printf("PBP: mass %.4g, %.8g %.8g %.8g %.8g ( %d of %d ) %d\n",
                 mass, r_pbp_even, r_pbp_odd, i_pbp_even, i_pbp_odd,
                 jpbp_reps + 1, npbp_reps, miters);
    node0_printf("FACTION: mass = %.4g, %.8g ( %d of %d )\n",
                 mass, r_ferm_action, jpbp_reps+1, npbp_reps);

#ifdef NPBP_REPS
    pbp_pbp = 0;
    FORALLSITES(i, s)
      vec_copy((vector *)F_PT(s, psi_off), &(s->M_inv));

    clear_latvec(psi_off, EVENANDODD);
    mat_invert_uml(F_OFFSET(M_inv), psi_off, chi_off, mass);
    FORALLSITES(i, s) {
      tc = su3_dot(&(s->g_rand), (vector *)F_PT(s, psi_off));
      pbp_pbp += tc.real;
    }
    g_doublesum(&pbp_pbp);
    pbp_pbp /= (double)volume;
    node0_printf("TR_MM_INV: mass %.4g, %.8g ( %d of %d )\n",
                 mass, pbp_pbp, jpbp_reps + 1, npbp_reps);
#endif
  }
  return tot_iters;
}
// -----------------------------------------------------------------

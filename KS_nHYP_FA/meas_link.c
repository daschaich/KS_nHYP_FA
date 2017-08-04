// -----------------------------------------------------------------
// Measure fermion observables for general "even plus odd" quark actions:
//   pbp separately on even and odd sites
//   fermion action
//   links in all directions separately on even and odd sites
#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return total number of iterations
int meas_link(field_offset chi_off, field_offset psi_off, Real mass) {
  register int i;
  register site *s;
  int jpbp_reps, tot_iters = 0, miters;
  double r_pbp_even, i_pbp_even, r_pbp_odd, i_pbp_odd, faction, pbp_pbp;
  complex cc;
  double_complex pbp_e[6], pbp_o[6];  // PBP, 1, t, x, y, z
  double_complex ave_e[6], ave_o[6];
  msg_tag *tag0, *tag1, *tag2, *tag3;

  for (i = 0; i < 6; i++) {
    ave_e[i] = dcmplx(0, 0);
    ave_o[i] = dcmplx(0, 0);
  }

  for (jpbp_reps = 0; jpbp_reps < npbp; jpbp_reps++) {
    faction = 0;
    for (i = 0; i < 6; i++) {
      pbp_e[i] = dcmplx(0, 0);
      pbp_o[i] = dcmplx(0, 0);
    }

    // Make random source and do inversion
    // Generate (one-mass) g_rand; chi = Mdag g_rand; psi = M^{-1} g_rand
    grsource_imp(chi_off, mass, EVENANDODD);
    clear_latvec(psi_off, EVENANDODD);
//    miters = mat_invert_uml(F_OFFSET(g_rand), psi_off, chi_off, mass);
    // Adapted from mat_invert.c "preconditioning": chi <- M^dag * g_rand
    dslash(F_OFFSET(g_rand), F_OFFSET(ttt[0][0]), EVENANDODD);
    scalar_mult_add_latvec(F_OFFSET(ttt[0][0]), F_OFFSET(g_rand),
                           -2.0 * mass, chi_off, EVENANDODD);
    scalar_mult_latvec(chi_off, -1.0, chi_off, EVENANDODD);
    miters = ks_congrad(chi_off, psi_off, mass, EVENANDODD);
    tot_iters += miters;

    // Fermion action = chi.psi
    // pbp on even sites = g_rand.psi
    FOREVENSITES(i, s) {
      cc = su3_dot((vector *)F_PT(s, chi_off),
                   (vector *)F_PT(s, psi_off));
      faction += cc.real;
      cc = su3_dot(&(s->g_rand), (vector *)F_PT(s, psi_off));
      CSUM(pbp_e[0], cc);
    }

    // pbp on odd sites
    FORODDSITES(i, s) {
      cc = su3_dot(&(s->g_rand), (vector *)F_PT(s, psi_off));
      CSUM(pbp_o[0], cc);
    }

    // Now calculate the link differences
    tag0 = start_gather_site(psi_off, sizeof(vector), TUP,
                             EVENANDODD, gen_pt[0]);
    tag1 = start_gather_site(psi_off, sizeof(vector), XUP,
                             EVENANDODD, gen_pt[1]);
    tag2 = start_gather_site(psi_off, sizeof(vector), YUP,
                             EVENANDODD, gen_pt[2]);
    tag3 = start_gather_site(psi_off, sizeof(vector), ZUP,
                             EVENANDODD, gen_pt[3]);

    wait_gather(tag0);
    wait_gather(tag1);
    wait_gather(tag2);
    wait_gather(tag3);

    FORALLSITES(i, s) {
      mult_mat_vec(&(s->link[TUP]),
                       (vector *)gen_pt[0][i], &(s->tempvec[0]));
      mult_mat_vec(&(s->link[XUP]),
                       (vector *)gen_pt[1][i], &(s->tempvec[1]));
      mult_mat_vec(&(s->link[YUP]),
                       (vector *)gen_pt[2][i], &(s->tempvec[2]));
      mult_mat_vec(&(s->link[ZUP]),
                       (vector *)gen_pt[3][i], &(s->tempvec[3]));
    }

    FORALLSITES(i, s) {
      cc = su3_dot(&(s->g_rand), &(s->tempvec[0]));
      if ((s->x) % 2 == 0)
        CSUM(pbp_e[1], cc);
      if ((s->x) % 2 == 1)
        CSUM(pbp_o[1], cc);
      if ((s->t) % 2 == 0)
        CSUM(pbp_e[2], cc);
      if ((s->t) % 2 == 1)
        CSUM(pbp_o[2], cc);

      cc = su3_dot(&(s->g_rand), &(s->tempvec[1]));
      if ((s->x) % 2 == 0)
        CSUM(pbp_e[3], cc);
      if ((s->x) % 2 == 1)
        CSUM(pbp_o[3], cc);

      cc = su3_dot(&(s->g_rand), &(s->tempvec[2]));
      if ((s->y) % 2 == 0)
        CSUM(pbp_e[4], cc);
      if ((s->y) % 2 == 1)
        CSUM(pbp_o[4], cc);

      cc = su3_dot(&(s->g_rand), &(s->tempvec[3]));
      if ((s->z) % 2 == 0)
        CSUM(pbp_e[5], cc);
      if ((s->z) % 2 == 1)
        CSUM(pbp_o[5], cc);
    }

    // Accumulate, normalize and print
    for (i = 0; i < 6; i++) {
      g_dcomplexsum(&pbp_e[i]);
      g_dcomplexsum(&pbp_o[i]);
      ave_e[i].real += pbp_e[i].real;
      ave_e[i].imag += pbp_e[i].imag;
      ave_o[i].real += pbp_o[i].real;
      ave_o[i].imag += pbp_o[i].imag;
      switch(i) {
        case 0: node0_printf("PBP: "); break;
        case 1: node0_printf("PBP1: "); break;
        case 2: node0_printf("PBPt: "); break;
        case 3: node0_printf("PBPx: "); break;
        case 4: node0_printf("PBPy: "); break;
        case 5: node0_printf("PBPz: "); break;
      }
      r_pbp_even = pbp_e[i].real * (2 / (double)volume);
      i_pbp_even = pbp_e[i].imag * (2 / (double)volume);
      r_pbp_odd = pbp_o[i].real * (2 / (double)volume);
      i_pbp_odd = pbp_o[i].imag * (2 / (double)volume);
      node0_printf("mass %.4g, %.8g %.8g %.8g %.8g ( %d of %d ) %d\n",
                   mass, r_pbp_even, r_pbp_odd, i_pbp_even, i_pbp_odd,
                   jpbp_reps + 1, npbp, miters);
    }
    g_doublesum(&faction);
    node0_printf("FACTION: mass = %.4g, %.8g ( %d of %d )\n",
                 mass, faction / (double)volume, jpbp_reps + 1, npbp);

    if (npbp > 1) {
      pbp_pbp = 0;
      FORALLSITES(i, s)
        su3vec_copy((vector *)F_PT(s, psi_off), &(s->M_inv));

      clear_latvec(psi_off, EVENANDODD);
//      mat_invert_uml(F_OFFSET(M_inv), psi_off, chi_off, mass);
      // Adapted from mat_invert.c "preconditioning": temp <- M^dag * M_inv
      dslash(F_OFFSET(M_inv), F_OFFSET(ttt[0][0]), EVENANDODD);
      scalar_mult_add_latvec(F_OFFSET(ttt[0][0]), F_OFFSET(M_inv),
                             -2.0 * mass, chi_off, EVENANDODD);
      scalar_mult_latvec(chi_off, -1.0, chi_off, EVENANDODD);
      ks_congrad(chi_off, psi_off, mass, EVENANDODD);
      FORALLSITES(i, s) {
        cc = su3_dot(&(s->g_rand), (vector *)F_PT(s, psi_off));
        pbp_pbp += cc.real;
      }
      g_doublesum(&pbp_pbp);
      pbp_pbp /= (double)volume;
      node0_printf("TR_MM_INV: mass %.4g, %.8g ( %d of %d )\n",
                   mass, pbp_pbp, jpbp_reps + 1, npbp);
    }
  }

  // Print out averages, summed over the entire (half-)volume
  for (i = 0; i < 6; i++) {
    switch(i) {
      case 0: node0_printf("pbp: "); break;
      case 1: node0_printf("pbp1: "); break;
      case 2: node0_printf("pbpt: "); break;
      case 3: node0_printf("pbpx: "); break;
      case 4: node0_printf("pbpy: "); break;
      case 5: node0_printf("pbpz: "); break;
    }
    r_pbp_even = ave_e[i].real / npbp;
    i_pbp_even = ave_e[i].imag / npbp;
    r_pbp_odd = ave_o[i].real / npbp;
    i_pbp_odd = ave_o[i].imag / npbp;
    node0_printf("mass %.4g, %.8g %.8g %.8g %.8g ( ave over %d )\n",
                 mass, r_pbp_even, r_pbp_odd, i_pbp_even, i_pbp_odd, npbp);
  }
  return tot_iters;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Spectrum for Kogut--Susskind pointlike hadrons
// Can do both "fat plus Naik" and "even plus odd" spectrum by defining
// "dslash_site" to be the appropriate "dslash_fn_site" or "dslash_eo_site"

// Does arbitrary number of wall sources
// Does NOT fix the gauge -- do that before calling

// Source information:
//   src_start >= 0  with src_start + n_src * src_inc < nt
//     invokes corner wall src
//   src_start >= nt with src_start + n_src * src_inc < 2 * nt
//     invokes point src (zero momentum only)
//   Caution: any other choice mixes point and corner wall sources

// Output information:
//   rho is VT (rho-b_1), rho2 is PV (rho-a_1)
//   also pi and pi2 -> pi_ps_prop, pi_sc_prop

// Stuck CG wrappers in this directory
#include "../ks_spectrum/spectrum_includes.h"
#include "generic_ks_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// vmass is valence mass, temp1 is source for inversion, temp2 is result
// n_src, src_start and src_inc need to be set by each application's setup.c
// Return the CG iteration count
int spectrum2(Real vmass, field_offset temp1, field_offset temp2) {
  register int i, x, y, z, t, icol, cgn = 0;
  register int t_src, t_off;
  register complex tc;
  int isrc, source_type = 0;    // 1 corner wall, 2 point, 3 mixed
  Real vmass_x2 = 2.0 * vmass, one_ov_N = 1.0 / (Real)n_src;
  Real pt_norm = -0.125 * nx * ny * nz;
  double *pi_ps_prop = malloc(nt * sizeof(*pi_ps_prop));
  double *pi_sc_prop = malloc(nt * sizeof(*pi_sc_prop));
  double *rho_pv_prop = malloc(nt * sizeof(*rho_pv_prop));
  double *rho_vt_prop = malloc(nt * sizeof(*rho_vt_prop));
  double *barprop = malloc(nt * sizeof(*barprop));
  char *source_string[4] = {"GOOFED", "CORNER", "POINT", "MIXED"};

  for (t = 0; t < nt; t++) {
    pi_ps_prop[t] = 0.0;
    pi_sc_prop[t] = 0.0;
    rho_pv_prop[t] = 0.0;
    rho_vt_prop[t] = 0.0;
    barprop[t] = 0.0;
  }

  for (t_src = src_start, isrc = 0;
       t_src < 2 * nt && isrc < n_src;
       ++isrc, t_src += src_inc) {

    node0_printf("spectrum(): source time = %d\n", t_src);
    for (icol = 0; icol < 3; icol++) {
      // Set up either wall or point source
      clear_latvec(temp1, EVENANDODD);
      if (t_src < nt) {       // Wall source
        for (x = 0; x < nx; x += 2) {
          for (y = 0; y < ny; y += 2) {
            for (z = 0; z < nz; z += 2) {
              if (node_number(x, y, z, t_src) != mynode())
                continue;
              i = node_index(x, y, z, t_src);
              ((vector *)(F_PT(&lattice[i], temp1)))->c[icol].real = -1.0;
            }
          }
        }
        source_type |= 1;
      }
      else {                  // Point source at origin
        if (node_number(0, 0, 0, t_src % nt) == mynode()) {
          i = node_index(0, 0, 0, t_src % nt);
          ((vector *)(F_PT(&lattice[i], temp1)))->c[icol].real = pt_norm;
        }
        source_type |= 2;
      }

      // Invert with source in temp1 and result in temp2
      // Then multiply by -M^dag to define propmat[icol]
      clear_latvec(temp2, EVENANDODD);
      if (t_src % 2 == 0) {
        cgn += ks_congrad(temp1, temp2, vmass, EVEN);
        dslash(temp2, F_OFFSET(propmat[icol]), ODD);
        scalar_mult_latvec(temp2, -vmass_x2, F_OFFSET(propmat[icol]), EVEN);
      }
      else {
        cgn += ks_congrad(temp1, temp2, vmass, ODD);
        dslash(temp2, F_OFFSET(propmat[icol]), EVEN);
        scalar_mult_latvec(temp2, -vmass_x2, F_OFFSET(propmat[icol]), ODD);
      }
    }

    // Measure the meson propagator
    for (t = 0; t < nt; t++) {
      // Define the time value offset t from t_src
      t_off = (t + t_src) % nt;

      for (x = 0; x < nx; x++) {
        for (y = 0; y < ny; y++) {
          for (z = 0; z < nz; z++) {
            if (node_number(x, y, z, t_off) != mynode())
              continue;
            for (icol = 0; icol < 3; icol++) {
              i = node_index(x, y, z, t_off);
              tc = su3_dot(&lattice[i].propmat[icol],
                           &lattice[i].propmat[icol]);

              pi_ps_prop[t] += tc.real;
              if ((x + y) % 2 == 0) rho_pv_prop[t] += tc.real;
              else                  rho_pv_prop[t] -= tc.real;
              if ((y + z) % 2 == 0) rho_pv_prop[t] += tc.real;
              else                  rho_pv_prop[t] -= tc.real;
              if ((z + x) % 2 == 0) rho_pv_prop[t] += tc.real;
              else                  rho_pv_prop[t] -= tc.real;

              if (x % 2 == 0) rho_vt_prop[t] += tc.real;
              else            rho_vt_prop[t] -= tc.real;
              if (y % 2 == 0) rho_vt_prop[t] += tc.real;
              else            rho_vt_prop[t] -= tc.real;
              if (z % 2 == 0) rho_vt_prop[t] += tc.real;
              else            rho_vt_prop[t] -= tc.real;

              if ((x + y + z) % 2 == 0) pi_sc_prop[t] += tc.real;
              else                      pi_sc_prop[t] -= tc.real;
            }
          }
        }
      }
    }

    // Measure the baryon propagator
    for (t = 0; t < nt; t++) {
      // Define the time value offset t from t_src
      t_off = (t + t_src) % nt;

      for (x = 0; x < nx; x += 2) {
        for (y = 0; y < ny; y += 2) {
          for (z = 0; z < nz; z += 2) {
            if (node_number(x, y, z, t_off) != mynode())
              continue;
            i = node_index(x, y, z, t_off);

            // typecast trick: propmat is three vectors
            tc = det_su3((matrix *)lattice[i].propmat);

            // Must get sign right
            // This looks to see if we have wrapped around the lattice
            // "t" is the distance from the source to the measurement,
            // so we are trying to find out if (t_src + t) >= nt
            // The t_src / nt subtraction is needed in case t_src = nt
            // Otherwise anti-periodic BCs change sign
            // (Sink point should really be in a copy of the lattice)
            if ((((t + t_src) / nt - t_src / nt) % 2) == 0)
              barprop[t] += tc.real;
            else
              barprop[t] -= tc.real;
          }
        }
      }
    }
  }

  // Dump the propagators, including vanishing imaginary components
  g_vecdoublesum(pi_ps_prop, nt);
  g_vecdoublesum(rho_pv_prop, nt);
  g_vecdoublesum(rho_vt_prop, nt);
  g_vecdoublesum(pi_sc_prop, nt);
  g_vecdoublesum(barprop, nt);
  for (t = 0; t < nt; t++) {
    pi_ps_prop[t] *= one_ov_N;
    pi_sc_prop[t] *= one_ov_N;
    rho_pv_prop[t] *= one_ov_N;
    rho_vt_prop[t] *= one_ov_N;
    barprop[t] *= one_ov_N;
  }
  if (this_node == 0) {
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n", vmass, vmass);
    printf("SOURCE: %s\n", source_string[source_type]);
    printf("SINKS: PION_PS PION_SC RHO_VT RHO_PV\n");
    for (t = 0; t < nt; t++) {
      printf("%d %e 0.0 %e 0.0 %e 0.0 %e 0.0\n", t,
             pi_ps_prop[t], pi_sc_prop[t], rho_vt_prop[t], rho_pv_prop[t]);
    }
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %e   %e  %e\n", vmass, vmass, vmass);
    printf("SOURCE: %s\n", source_string[source_type]);
    printf("SINKS: NUCLEON\n");
    for (t = 0; t < nt; t++)
      printf("%d %e 0.0\n", t, barprop[t]);
    printf("ENDPROP\n");

    fflush(stdout);
  }
  free(pi_ps_prop);
  free(pi_sc_prop);
  free(rho_pv_prop);
  free(rho_vt_prop);
  free(barprop);
  return cgn;
}
// -----------------------------------------------------------------

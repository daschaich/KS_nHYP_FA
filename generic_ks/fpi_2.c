// -----------------------------------------------------------------
// Wall and point source/sink propagators for f_pi, f_K, ...
// Use unit random vectors for point sources
// Use s->chi for source, allocate temporary propagators

// Call with rephase(ON) to absorb (anti-periodic) boundary conditions
// into the link matrices

// Coulomb gauge should be fixed before calling this routine
// Results are NOT gauge invariant for any given source
// Should get gauge-invariant results by averaging over sources

// Sources normalized by a factor of 2 (twice the usual convention)
// since ks_congrad inverts matrix with 2m on diagonal,

#include "generic_ks_includes.h"
#include <string.h>

// Number of vectors used to construct each point (RANDOM_WALL) source
#define NSRCVECS 3

// Build shorthand for our four propagators, and count them with NPROP
enum prop_name {
  PP,       // Point source point sink
  PW,       // Point source wall sink
  WP,       // Wall source point sink
  WW,       // Wall source wall sink
  NPROP     // Number of propagators
};
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrappers for dslash and CG
// Just call dslash on src and copy result into dest
// Use psi in site for temporary storage
void dslash_wrapper(vector *src, vector *dest, int parity) {
  register int i;
  register site *s;

  FORALLSITES(i, s)
    vec_copy(&(src[i]), &(s->chi));
  dslash(F_OFFSET(chi), F_OFFSET(psi), parity);
  FORALLSITES(i, s)
    vec_copy(&(s->psi), &(dest[i]));
}

// Just call CG on chi and copy result into props[0]
// Use psi in site for temporary storage
int CG_wrapper(field_offset chi, vector **props, Real m, int parity) {
  register int i, iter;
  register site *s;

  iter = ks_congrad(chi, F_OFFSET(psi), m, parity);
  FORALLSITES(i, s)
    vec_copy(&(s->psi), &(props[0][i]));
  return iter;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Single mass hard-coded at the moment
// n_src, src_start and src_inc need to be set by each application's setup.c
// In general only propagators with m1 <= m2 are printed
// Return number of CG iterations
int fpi_2(Real *masses) {
  register int i, j, m1, m2, t_src;
  register site *s;
  register complex tc = cmplx(0.0, 0.0);
  int cgn = 0, icol, isrc, nmasses = 1, offset;
  Real tr;
  Real pt_norm = 1.0 / ((Real)NSRCVECS * n_src * nx * ny * nz * nx * ny * nz);
  Real wall_norm = 1.0 / (3.0 * n_src * nx * ny * nz * nx * ny * nz);
  complex **props = malloc(NPROP * nmasses * nmasses * sizeof(**props));
  vector **fprops = malloc(nmasses * sizeof(**fprops));
  vector *temp_prop = malloc(sites_on_node * sizeof(*temp_prop));
  vector *wall_sink_m1 = malloc(nt * sizeof(*wall_sink_m1));
  vector *wall_sink_m2 = malloc(nt * sizeof(*wall_sink_m2));

  // Set up and clear arrays to accumulate propagators
  for (i = 0; i < NPROP * nmasses * nmasses; i++) {
    props[i] = malloc(nt * sizeof(complex));
    for (j = 0; j < nt; j++)
      props[i][j] = tc;
  }
  for (i = 0; i < nmasses; i++) {
    fprops[i] = malloc(sites_on_node * sizeof(vector));
    memset((void *)fprops[i], '\0', sites_on_node * sizeof(vector));
  }

  // Loop over sources
  for (isrc = 0, t_src = src_start;
       t_src < nt && isrc < n_src;
       t_src += src_inc, isrc++) {

    node0_printf("fpi_2: t_src = %d, num_sourcevecs = %d\n",
                 t_src, NSRCVECS);

    // Point (RANDOM_WALL) source
    for (j = 0; j < NSRCVECS; j++) {
      FORALLSITES(i, s) {
        clearvec(&(s->chi));
        if (s->t == t_src) {
          for (icol = 0; icol < 3; icol++) {
            s->chi.c[icol].real = gaussian_rand_no(&(s->site_prn));
            s->chi.c[icol].imag = gaussian_rand_no(&(s->site_prn));
          }
          // Factor of 2 because we invert 2m + 2 Dslash
          tr = 2.0 / sqrt(magsq_vec(&(s->chi)));
          scalar_mult_vector(&(s->chi), tr, &(s->chi));
        }
      }

      // Compute M^-1 * chi
      cgn += CG_wrapper(F_OFFSET(chi), fprops, masses[0], EVEN);
      cgn += CG_wrapper(F_OFFSET(chi), fprops, masses[0], ODD);
      // Multiply by M^dag
      dslash_wrapper(fprops[0], temp_prop, EVENANDODD);
      FORALLSITES(i, s) {
        scalar_mult_vector(&(fprops[0][i]), 2 * masses[0],
                               &(fprops[0][i]));
        dif_vector(&(temp_prop[i]), &(fprops[0][i]));
      }

      /* 0-+ (kaon) propagators */
      for (m1 = 0; m1 < nmasses; m1++) {
        for (m2 = m1; m2 < nmasses; m2++) {
          offset = NPROP * (nmasses * m1 + m2);
          for (i = 0; i < nt; i++) {
            clearvec(&(wall_sink_m1[i]));
            clearvec(&(wall_sink_m2[i]));
          }
          FORALLSITES(i, s) {
            sum_vector(&(fprops[m1][i]), &(wall_sink_m1[s->t]));
            sum_vector(&(fprops[m2][i]), &(wall_sink_m2[s->t]));
            tc = su3_dot(&(fprops[m1][i]), &(fprops[m2][i]));
            CSUM(props[offset + PP][(s->t + nt - t_src) % nt], tc);
          }
          g_veccomplexsum((complex *)wall_sink_m1, 3 * nt);
          g_veccomplexsum((complex *)wall_sink_m2, 3 * nt);
          for (i = 0; i < nt; i++) {
            tc = su3_dot(&(wall_sink_m1[i]), &(wall_sink_m2[i]));
            CSUM(props[offset + PW][(i + nt - t_src) % nt], tc);
          }
        }
      }
    }

    // Wall source
    for (icol = 0; icol < 3; icol++) {
      FORALLSITES(i, s) {
        clearvec(&(s->chi));
        if (s->t==t_src)   // Factor of 2 from inverting 2(m + Dslash)
          s->chi.c[icol] = cmplx(0.0, 2.0);
      }

      // Compute M^-1 * chi
      cgn += CG_wrapper(F_OFFSET(chi), fprops, masses[0], EVEN);
      cgn += CG_wrapper(F_OFFSET(chi), fprops, masses[0], ODD);
      // Multiply by M^dag
      for (j = 0; j < nmasses; j++) {
        dslash_wrapper(fprops[j], temp_prop, EVENANDODD);
        FORALLSITES(i,s) {
          scalar_mult_vector(&(fprops[j][i]), 2.0 * masses[j],
                                 &(fprops[j][i]));
          dif_vector(&(temp_prop[i]), &(fprops[j][i]));
        }
      }

      /* 0-+ (kaon) propagators */
      for (m1 = 0; m1 < nmasses; m1++) {
        for (m2 = m1; m2 < nmasses; m2++) {
          offset = NPROP * (nmasses * m1 + m2);
          for (i = 0; i < nt; i++) {
            clearvec(&(wall_sink_m1[i]));
            clearvec(&(wall_sink_m2[i]));
          }
          FORALLSITES(i,s) {
            sum_vector(&(fprops[m1][i]), &(wall_sink_m1[s->t]));
            sum_vector(&(fprops[m2][i]), &(wall_sink_m2[s->t]));
            tc = su3_dot(&(fprops[m1][i]), &(fprops[m2][i]));
            CSUM(props[offset + WP][(s->t + nt - t_src) % nt], tc);
          }
          g_veccomplexsum((complex *)wall_sink_m1, 3*nt);
          g_veccomplexsum((complex *)wall_sink_m2, 3*nt);
          for (i = 0; i < nt; i++) {
            tc = su3_dot(&(wall_sink_m1[i]), &(wall_sink_m2[i]));
            CSUM(props[offset + WW][(i + nt - t_src) % nt], tc);
          }
        }
      }
    }
  }

  // Sum propagator arrays over nodes and normalize
  for (i = 0; i < NPROP * nmasses * nmasses; i += NPROP) {
    g_veccomplexsum(props[i + PP], nt);
    g_veccomplexsum(props[i + WP], nt);
  }
  for (i = 0; i < nmasses * nmasses; i++) {
    offset = NPROP * i;
    for (j = 0; j < nt; j++) {
      CMULREAL(props[offset + PP][j], pt_norm, props[offset + PP][j]);
      CMULREAL(props[offset + PW][j], pt_norm, props[offset + PW][j]);
      CMULREAL(props[offset + WP][j], wall_norm, props[offset + WP][j]);
      CMULREAL(props[offset + WW][j], wall_norm, props[offset + WW][j]);
    }
  }

  // Dump the propagators
  // Seem to be purely real, possibly due to single mass being used
  if (this_node == 0) {
    for (m1 = 0; m1 < nmasses; m1++) {
      for (m2 = m1; m2 < nmasses; m2++) {
        offset = NPROP * (nmasses * m1 + m2);
        printf("STARTPROP\n");
        printf("MASSES:  %.5e   %.5e\n", masses[m1], masses[m2]);
        printf("SOURCE: RANDOM_WALL\n");
        printf("SINKS: POINT_KAON_5 WALL_KAON_5\n");
        for (j = 0; j < nt; j++) {
          printf("%d %e %g %e %g\n",j,
                 props[offset + PP][j].real, props[offset + PP][j].imag,
                 props[offset + PW][j].real, props[offset + PW][j].imag);
        }
        printf("ENDPROP\n");

        printf("STARTPROP\n");
        printf("MASSES:  %.5e   %.5e\n", masses[m1], masses[m2]);
        printf("SOURCE: FULL_WALL\n");
        printf("SINKS: POINT_KAON_5 WALL_KAON_5\n");
        for (j = 0; j < nt; j++) {
          printf("%d %e %g %e %g\n",j,
                 props[offset + WP][j].real, props[offset + WP][j].imag,
                 props[offset + WW][j].real, props[offset + WW][j].imag);
        }
        printf("ENDPROP\n");
      }
    }
    fflush(stdout);
  }

  // Free temporaries
  for (i = 0; i < NPROP * nmasses * nmasses; i++)
    free(props[i]);
  free(props);
  for (i = 0; i < nmasses; i++)
    free(fprops[i]);
  free(fprops);
  free(temp_prop);
  free(wall_sink_m1);
  free(wall_sink_m2);
  return cgn;
}
// -----------------------------------------------------------------

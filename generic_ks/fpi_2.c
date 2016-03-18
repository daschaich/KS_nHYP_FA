// -----------------------------------------------------------------
// Wall and point source/sink propagators for f_pi, f_K, ...
// Use unit random vectors for point sources
//
// Call with rephase(ON) to absorb (antiperiodic) boundary conditions
// into the link matrices
//
// Coulomb gauge should be fixed before calling this routine
// Results are NOT gauge invariant for any given source
// Should get gauge-invariant results by averaging over sources

// Sources normalized by a factor of 2 (twice the usual convention)
// since ks_congrad inverts matrix with 2m on diagonal,

// Only propagators with m1 <= m2 are printed
#define NSOURCEVECS 3
enum prop_name {
    POINTPOINT,   // Point source point sink
    POINTWALL,    // Point source wall sink
    WALLPOINT,    // Wall source point sink
    WALLWALL,     // Wall source wall sink
    NPROPS        // Number of propagators
};

// Stuck CG wrappers in this directory
#include "../ks_spectrum/spectrum_includes.h"
#include "generic_ks_includes.h"
#include <string.h>
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Stripped ferm_links_t
// Only works for a single mass at the moment
// Return number of CG iterations
int fpi_2(Real *masses) {
  register int i, j, m1, m2;
  register site* s;
  register complex cc, czero;
  register int t_source;
  int cgn = 0, sourcevec, color;
  int src_count; /* number of source time slices used */
  int nmasses = 1;
  Real x;
  complex **props;  /* arrays of propagators */
  su3_vector **quark_props, *temp_prop;
  su3_vector *wall_sink_m1, *wall_sink_m2;

  czero.real = 0.0;
  czero.imag = 0.0;

  // Allocate space for quark propagators
  quark_props = (su3_vector **)malloc(nmasses*sizeof(su3_vector *));
  for (i = 0; i < nmasses; i++) {
    quark_props[i] = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
    memset((void *)quark_props[i], '\0', sites_on_node * sizeof(su3_vector));
  }
  temp_prop = (su3_vector *)malloc(sites_on_node * sizeof(su3_vector));
  wall_sink_m1 = (su3_vector *)malloc(nt * sizeof(su3_vector));
  wall_sink_m2 = (su3_vector *)malloc(nt * sizeof(su3_vector));

  /* allocate space for meson propagators (NPROPS numbers per timeslice) */
  /* for meson with masses number m1 and m2, use props[NPROPS*(nmasses*m1+m2)+TYPE] */
  props = (complex **)malloc(NPROPS*nmasses*nmasses*sizeof(complex *) );
  props[0] = (complex *)malloc(NPROPS*nmasses*nmasses*nt*sizeof(complex) );
  for (i = 1; i < NPROPS * nmasses * nmasses; i++)
    props[i] = props[i - 1] + nt;

  // Zero propagators
  for (i = 0; i < NPROPS * nmasses * nmasses; i++) {
    for (j = 0; j < nt; j++) {
      props[i][j] = czero;
    }
  }

  // Loop over sources
  for (src_count = 0, t_source = src_start;
       t_source < nt && src_count < n_src;
       t_source += src_inc, src_count++) {

    node0_printf("fpi_2: t_src = %d, num_sourcevecs = %d\n",
                 t_source, NSOURCEVECS);

    // Point (RANDOM_WALL) source
    for (sourcevec = 0; sourcevec < NSOURCEVECS; sourcevec++) {
      FORALLSITES(i, s) {
        clearvec(&(s->quark_source));
        if (s->t == t_source) {
          for (color = 0; color < 3; color++) {
            s->quark_source.c[color].real = gaussian_rand_no(&(s->site_prn));
            s->quark_source.c[color].imag = gaussian_rand_no(&(s->site_prn));
          }
          // Factor of 2 because we invert 2m + 2 Dslash
          x = 2.0 / sqrt(magsq_su3vec(&(s->quark_source)));
          scalar_mult_su3_vector(&(s->quark_source), x, &(s->quark_source));
        }
      }

      // Compute M^-1 * quark_source
      cgn += CG_wrapper(F_OFFSET(quark_source), quark_props, masses[0], EVEN);
      cgn += CG_wrapper(F_OFFSET(quark_source), quark_props, masses[0], ODD);
      // Multiply by M^dag
      dslash_wrapper(quark_props[0], temp_prop, EVENANDODD);
      FORALLSITES(i, s) {
        scalar_mult_su3_vector(&(quark_props[0][i]), 2 * masses[0],
                               &(quark_props[0][i]));
        scalar_mult_add_su3_vector(&(quark_props[0][i]), &(temp_prop[i]),
                                   -1.0, &(quark_props[0][i]));
      }

      /* 0-+ (kaon) propagators */
      for (m1=0;m1<nmasses;m1++)for (m2=m1;m2<nmasses;m2++) {
        for (i = 0; i < nt; i++) {
          clearvec(&(wall_sink_m1[i]));
          clearvec(&(wall_sink_m2[i]));
        }
        FORALLSITES(i, s) {
          add_su3_vector(&(wall_sink_m1[s->t]), &(quark_props[m1][i]),
                         &(wall_sink_m1[s->t]));
          add_su3_vector(&(wall_sink_m2[s->t]), &(quark_props[m2][i]),
                         &(wall_sink_m2[s->t]));
          cc = su3_dot(&(quark_props[m1][i]), &(quark_props[m2][i]));
          CSUM(props[NPROPS*(nmasses*m1+m2)+POINTPOINT][(s->t+nt-t_source)%nt], cc);
        }
        g_veccomplexsum((complex *)wall_sink_m1, 3 * nt);
        g_veccomplexsum((complex *)wall_sink_m2, 3 * nt);
        for (i = 0; i < nt; i++) {
          cc = su3_dot(&(wall_sink_m1[i]), &(wall_sink_m2[i]));
          CSUM(props[NPROPS*(nmasses*m1+m2)+POINTWALL][(i+nt-t_source)%nt], cc);
        }
      } /* m1,m2=masses */
    } /* sourcevec */

    // Wall source
    for (color = 0; color < 3; color++) {  /* loop over colors instead of random vectors */
      FORALLSITES(i, s) {
        clearvec(&(s->quark_source));
        if (s->t==t_source) {
          s->quark_source.c[color].real  = 0.0;
          s->quark_source.c[color].imag  = 2.0;
          /* "2" because we invert 2m + 2 Dslash */
        }
      }

      // Compute M^-1 * quark_source
      cgn += CG_wrapper(F_OFFSET(quark_source), quark_props, masses[0], EVEN);
      cgn += CG_wrapper(F_OFFSET(quark_source), quark_props, masses[0], ODD);
      // Multiply by M^dag
      for (j=0;j<nmasses;j++) {
        dslash_wrapper(quark_props[j], temp_prop, EVENANDODD);
        FORALLSITES(i,s) {
          scalar_mult_su3_vector( &(quark_props[j][i]), 2.0*masses[j], &(quark_props[j][i]) );
          scalar_mult_add_su3_vector( &(quark_props[j][i]), &(temp_prop[i]), -1.0,
              &(quark_props[j][i]) );
        }
      } /* j=masses */

      /* 0-+ (kaon) propagators */
      for (m1=0;m1<nmasses;m1++)for (m2=m1;m2<nmasses;m2++) {
        for (i = 0; i < nt; i++) {
          clearvec( &(wall_sink_m1[i]) );
          clearvec( &(wall_sink_m2[i]) );
        }
        FORALLSITES(i,s) {
          add_su3_vector( &(wall_sink_m1[s->t]),&(quark_props[m1][i]), &(wall_sink_m1[s->t]) );
          add_su3_vector( &(wall_sink_m2[s->t]),&(quark_props[m2][i]), &(wall_sink_m2[s->t]) );
          cc = su3_dot( &(quark_props[m1][i]), &(quark_props[m2][i]) );
          CSUM( props[NPROPS*(nmasses*m1+m2)+WALLPOINT][(s->t+nt-t_source)%nt], cc );
        }
        g_veccomplexsum( (complex *)wall_sink_m1, 3*nt);
        g_veccomplexsum( (complex *)wall_sink_m2, 3*nt);
        for (i = 0; i < nt; i++) {
          cc = su3_dot( &(wall_sink_m1[i]), &(wall_sink_m2[i]) );
          CSUM( props[NPROPS*(nmasses*m1+m2)+WALLWALL][(i+nt-t_source)%nt], cc );
        }
      } /* m1,m2=masses */

    } /* color */
  } /* end loop on t_source */

  /* Sum propagator arrays over nodes */
  /* print out propagators */
  for (i=0;i<NPROPS*nmasses*nmasses;i+=NPROPS) {
      g_veccomplexsum( props[i+POINTPOINT], nt);
      g_veccomplexsum( props[i+WALLPOINT], nt);
  }
  for (i=0;i<nmasses*nmasses;i++) {
    for (j=0;j<nt;j++) {
      CDIVREAL(props[NPROPS*i+POINTPOINT][j],(Real)NSOURCEVECS*n_src*nx*ny*nz*nx*ny*nz,
          props[NPROPS*i+POINTPOINT][j]);
      CDIVREAL(props[NPROPS*i+POINTWALL ][j],(Real)NSOURCEVECS*n_src*nx*ny*nz*nx*ny*nz,
          props[NPROPS*i+POINTWALL ][j]);
      CDIVREAL(props[NPROPS*i+WALLPOINT ][j],3.0*n_src*nx*ny*nz*nx*ny*nz,
          props[NPROPS*i+WALLPOINT ][j]);
      CDIVREAL(props[NPROPS*i+WALLWALL  ][j],3.0*n_src*nx*ny*nz*nx*ny*nz,
          props[NPROPS*i+WALLWALL  ][j]);
    }
  }
  if (this_node==0) {
    for (m1=0;m1<nmasses;m1++)for (m2=m1;m2<nmasses;m2++) {
        printf("STARTPROP\n");
        printf("MASSES:  %.5e   %.5e\n",masses[m1],masses[m2]);
        printf("SOURCE: RANDOM_WALL\n");
        printf("SINKS: POINT_KAON_5 WALL_KAON_5\n");
        for (j=0;j<nt;j++) {
      printf("%d %e %e %e %e\n",j,
        props[NPROPS*(nmasses*m1+m2)+POINTPOINT][j].real,
        props[NPROPS*(nmasses*m1+m2)+POINTPOINT][j].imag,
        props[NPROPS*(nmasses*m1+m2)+POINTWALL][j].real,
        props[NPROPS*(nmasses*m1+m2)+POINTWALL][j].imag );
        }
        printf("ENDPROP\n");

        printf("STARTPROP\n");
        printf("MASSES:  %.5e   %.5e\n",masses[m1],masses[m2]);
        printf("SOURCE: FULL_WALL\n");
        printf("SINKS: POINT_KAON_5 WALL_KAON_5\n");
        for (j=0;j<nt;j++) {
      printf("%d %e %e %e %e\n",j,
        props[NPROPS*(nmasses*m1+m2)+WALLPOINT][j].real,
        props[NPROPS*(nmasses*m1+m2)+WALLPOINT][j].imag,
        props[NPROPS*(nmasses*m1+m2)+WALLWALL][j].real,
        props[NPROPS*(nmasses*m1+m2)+WALLWALL][j].imag );
        }
        printf("ENDPROP\n");
    } /* masses */
    fflush(stdout);
  }

  // Free arrays
  free(props[0]);
  free(props);
  for (i = 0; i < nmasses; i++)
    free(quark_props[i]);
  free(quark_props);
  free(temp_prop);
  free(wall_sink_m1);
  free(wall_sink_m2);
  return cgn;
}
// -----------------------------------------------------------------

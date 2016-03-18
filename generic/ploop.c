// -----------------------------------------------------------------
// Evaluate Polyakov loops in arbitrary (even-length) direction
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
complex ploop(int dir) {
  register int i, j;
  register site *s;
  msg_tag *tag;
  complex sum, plp;
  int d[4], N = nt;

  sum = cmplx(0, 0);
  for (i = XUP; i <= TUP; i++)
    d[i] = 0;

  switch(dir) {
    case XUP: N = nx; break;
    case YUP: N = ny; break;
    case ZUP: N = nz; break;
    case TUP: N = nt; break;
    default:
      node0_printf("ERROR: UNRECOGNIZED DIRECTION IN BLOCKED_PLOOP()\n");
      terminate(1);
  }

  // First multiply the link on every even site by the next link
  // Compute the loop "at" the even sites in the first two time slices
  tag = start_gather_site(F_OFFSET(link[dir]), sizeof(su3_matrix),
                          dir, EVEN, gen_pt[0]);

  wait_gather(tag);
  FOREVENSITES(i, s)
    mult_su3_nn(&(s->link[dir]), (su3_matrix *)gen_pt[0][i],
                (su3_matrix *)&(s->tempmat1));

  cleanup_gather(tag);

  for (j = 2; j < N; j += 2) {
    d[dir] = j;     // Distance from which to gather
    tag = start_general_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix),
                                    d, EVEN, gen_pt[0]);
    wait_general_gather(tag);
    FOREVENSITES(i, s) {
      // Overwrite tempmat1 on the first two slices
      // Leave other links undisturbed so we can still gather them
      switch(dir) {
        case XUP: if (s->x > 1) continue; break;
        case YUP: if (s->y > 1) continue; break;
        case ZUP: if (s->z > 1) continue; break;
        case TUP: if (s->t > 1) continue; break;
      }
      mult_su3_nn((su3_matrix *)&(s->tempmat1), (su3_matrix *)gen_pt[0][i],
                  (su3_matrix *)&(s->tempmat2));
      su3mat_copy((su3_matrix *)&(s->tempmat2),
                  (su3_matrix *)&(s->tempmat1));
    }
    cleanup_general_gather(tag);
  }
  FOREVENSITES(i, s) {
    switch(dir) {
      case XUP: if (s->x > 1) continue; break;
      case YUP: if (s->y > 1) continue; break;
      case ZUP: if (s->z > 1) continue; break;
      case TUP: if (s->t > 1) continue; break;
    }
    plp = trace_su3((su3_matrix *)&(s->tempmat1));
    CSUM(sum, plp);   // Running complex sum
  }

  // Average all the loops we just calculated
  g_complexsum(&sum);
  plp.real = sum.real * N / ((double)volume);
  plp.imag = sum.imag * N / ((double)volume);
  return(plp);
}
// -----------------------------------------------------------------

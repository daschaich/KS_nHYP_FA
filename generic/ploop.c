// -----------------------------------------------------------------
// Evaluate Polyakov loops in arbitrary (even-length) direction
// Use tempmat and tempmat2 for temporary storage
#include "generic_includes.h"

complex ploop(int dir) {
  register int i, j;
  register site *s;
  int d[4], N = nt;
  msg_tag *tag;
  complex sum = cmplx(0.0, 0.0), plp;

  FORALLUPDIR(i)
    d[i] = 0;

  switch(dir) {
    case XUP: N = nx; break;
    case YUP: N = ny; break;
    case ZUP: N = nz; break;
    case TUP: N = nt; break;
    default:
      node0_printf("ERROR: Unrecognized direction in blocked_ploop\n");
      terminate(1);
  }

  // First multiply the link on every even site by the next link
  // Compute the loop "at" the even sites in the first two time slices
  tag = start_gather_site(F_OFFSET(link[dir]), sizeof(su3_matrix),
                          dir, EVEN, gen_pt[0]);

  wait_gather(tag);
  FOREVENSITES(i, s)
    mult_su3_nn(&(s->link[dir]), (su3_matrix *)gen_pt[0][i], &(tempmat[i]));

  cleanup_gather(tag);

  for (j = 2; j < N; j += 2) {
    d[dir] = j;     // Distance from which to gather
    tag = start_general_gather_field(tempmat, sizeof(su3_matrix),
                                     d, EVEN, gen_pt[0]);
    wait_general_gather(tag);
    FOREVENSITES(i, s) {
      // Overwrite tempmat on the first two slices
      // Leave other links undisturbed so we can still gather them
      switch(dir) {
        case XUP: if (s->x > 1) continue; break;
        case YUP: if (s->y > 1) continue; break;
        case ZUP: if (s->z > 1) continue; break;
        case TUP: if (s->t > 1) continue; break;
      }
      mult_su3_nn(&(tempmat[i]), (su3_matrix *)gen_pt[0][i], &(tempmat2[i]));
      su3mat_copy(&(tempmat2[i]), &(tempmat[i]));
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
    plp = trace_su3(&(tempmat[i]));
    CSUM(sum, plp);   // Running complex sum
  }

  // Average all the loops we just calculated
  g_complexsum(&sum);
  plp.real = sum.real * N / ((double)volume);
  plp.imag = sum.imag * N / ((double)volume);
  return plp;
}
// -----------------------------------------------------------------

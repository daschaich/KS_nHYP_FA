// -----------------------------------------------------------------
// Evaluate Polyakov loops in arbitrary (even-length) direction
// Use tempmat and tempmat2 for temporary storage
// General gathers can be replaced, but probably don't matter
#include "generic_includes.h"

// Compute the Polyakov loop at the even sites in the first two time-slices
complex ploop(int dir) {
  register int i, j;
  register site *s;
  int d[4] = {0, 0, 0, 0}, Ndir = nt;
  complex plp = cmplx(0.0, 0.0);
  msg_tag *tag;

  switch(dir) {
    case XUP: Ndir = nx; break;
    case YUP: Ndir = ny; break;
    case ZUP: Ndir = nz; break;
    case TUP: Ndir = nt; break;
    default:
      node0_printf("ERROR: Unrecognized direction in blocked_ploop\n");
      terminate(1);
  }

  // First multiply the link on every even site by the next link
  tag = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                          dir, EVEN, gen_pt[0]);
  wait_gather(tag);
  FOREVENSITES(i, s)
    mult_nn(&(s->link[dir]), (matrix *)gen_pt[0][i], &(tempmat[i]));
  cleanup_gather(tag);

  for (j = 2; j < Ndir; j += 2) {
    d[dir] = j;     // Distance from which to gather
    tag = start_general_gather_field(tempmat, sizeof(matrix),
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
      mult_nn(&(tempmat[i]), (matrix *)gen_pt[0][i], &(tempmat2[i]));
      mat_copy(&(tempmat2[i]), &(tempmat[i]));
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
    trace_sum(&(tempmat[i]), &plp);
  }
  g_complexsum(&plp);
  plp.real *= Ndir * one_ov_vol;
  plp.imag *= Ndir * one_ov_vol;
  return plp;
}
// -----------------------------------------------------------------

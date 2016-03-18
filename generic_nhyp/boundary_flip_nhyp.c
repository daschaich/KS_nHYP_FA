// -----------------------------------------------------------------
#include "nhyp_includes.h"
static int current_boundary = PLUS;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void boundary_flip_nhyp(int sign) {
  register int i, j, k;
  register site *s;
  int i1;
  if (sign == current_boundary) {
    node0_printf("WARNING: you lost track of the boundary conditions!\n");
    return;
  }

  FORALLSITES(i, s) {
    if (s->t != nt - 1)
      continue; // Hit only last time slice
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++){
        s->link[TUP].e[j][k].real *= -1.0;
        s->link[TUP].e[j][k].imag *= -1.0;
        gauge_field_thin[TUP][i].e[j][k].real *= -1.0;
        gauge_field_thin[TUP][i].e[j][k].imag *= -1.0;
        gauge_field[TUP][i].e[j][k].real *= -1.0;
        gauge_field[TUP][i].e[j][k].imag *= -1.0;
      }
    }
  }
  FORALLSITES(i, s) {
    if (s->t != nt - 1)
      continue; // Hit only last time slice
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        for (i1 = 0; i1 < 4; i1++) {
          if (i1 != TUP) {
            hyplink1[i1][TUP][i].e[j][k].real *= -1.0;
            hyplink1[i1][TUP][i].e[j][k].imag *= -1.0;
            hyplink2[i1][TUP][i].e[j][k].real *= -1.0;
            hyplink2[i1][TUP][i].e[j][k].imag *= -1.0;
            Staple1[i1][TUP][i].e[j][k].real *= -1.0;
            Staple1[i1][TUP][i].e[j][k].imag *= -1.0;
            Staple2[i1][TUP][i].e[j][k].real *= -1.0;
            Staple2[i1][TUP][i].e[j][k].imag *= -1.0;
          }
        }
        Staple3[TUP][i].e[j][k].real *= -1.0;
        Staple3[TUP][i].e[j][k].imag *= -1.0;
      }
    }
  }
  current_boundary = sign;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Set up and apply Kogut--Susskind phase factors
// Default is antiperiodic boundary conditions in time,
// periodic boundary conditions in space
#include "generic_ks_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up Kogut--Susskind phase factors, sum(nu < mu) {-1^i[nu]}
// combined with boundary conditions
void phaseset() {
  register site *s;
  register int i;

#ifdef PERIODICBC
  node0_printf("Using periodic boundary conditions in time\n");
#endif
#ifdef ANTI_PBC
  node0_printf("Using antiperiodic boundary conditions in all directions\n");
#endif

  FORALLSITES(i, s){
    s->phase[TUP] = 1.0;
    if ((s->t) % 2 == 1)
      s->phase[XUP] = -1.0;
    else
      s->phase[XUP] = 1.0;

    if ((s->x) % 2 == 1)
      s->phase[YUP] = -s->phase[XUP];
    else
      s->phase[YUP] = s->phase[XUP];
    if ((s->y) % 2 == 1)
      s->phase[ZUP] = -s->phase[YUP];
    else
      s->phase[ZUP] = s->phase[YUP];

#ifndef PERIODICBC
    // Antiperiodic boundary conditions in time
    // All t phases for t = nt - 1 time slice get extra minus sign
    if (s->t == nt - 1)
      s->phase[TUP] = -s->phase[TUP];
#endif

#ifdef ANTI_PBC
#ifdef PERIODICBC
    node0_printf("DUMMY: can't have periodic and antiperiodic BCs\n");
    terminate(1);
#endif
    // Antiperiodic boundary conditions in space
    // Time already done above
    if (s->x == nx - 1)
      s->phase[XUP] = -s->phase[XUP];
    if (s->y == ny - 1)
      s->phase[YUP] = -s->phase[YUP];
    if (s->z == nz - 1)
      s->phase[ZUP] = -s->phase[ZUP];
#endif
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Apply or remove Kogut--Susskind phase factors and boundary conditions
void rephase(int flag) {
  register int i, j, k, dir;
  register site *s;
  // Check to make sure we are going in expected direction
  if (!((flag == ON && phases_in == OFF) || (flag == OFF && phases_in == ON))) {
    node0_printf("DUMMY: you fouled up the phases\n");
    terminate(1);
  }

  FORALLSITES(i, s) {
    for (dir = XUP; dir <= TUP; dir++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          s->link[dir].e[j][k].real *= s->phase[dir];
          s->link[dir].e[j][k].imag *= s->phase[dir];
        }
      }
    }
  }
  phases_in = flag;
}
// -----------------------------------------------------------------

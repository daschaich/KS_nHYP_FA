// -----------------------------------------------------------------
// Calculate the connected part of hadronic vacuum polarization
// (two-point function of the conserved staggered vector current)
#include "vacpol_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return CG iteration count
int vacuum_polarization() {
  int i, j, k, icol, cgTot = 0, cgs, mu, nu, tmp_parity = 0;
  int x = x_src, y = y_src, z = z_src, t = t_src;
  char src_parity;
  double invtime;
  complex tc, tempvec, tempaxi, contact[4];
  Real neg = -1.0, tol = sqrt(rsqmin);
  site* s;
  msg_tag *tag, *tag2;
  matrix tmat, tmat2, sourcelink[4];

  // Initialize link at source point
  FORALLUPDIR(nu) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        sourcelink[nu].e[j][k].real = 0;
        sourcelink[nu].e[j][k].imag = 0;
      }
    }
  }

  // Set link at source point, needed for point-split correlation function
  if (node_number(x, y, z, t) == mynode()) {
    i = node_index(x, y, z, t);
    sourcelink[XUP] = lattice[i].link[XUP];
    sourcelink[YUP] = lattice[i].link[YUP];
    sourcelink[ZUP] = lattice[i].link[ZUP];
    sourcelink[TUP] = lattice[i].link[TUP];
    if (lattice[i].parity == ODD)
      tmp_parity = 1;     // Default above is even
  }
  // Sum parity and sourcelink to "broadcast" them to all nodes
  g_intsum(&tmp_parity);
  if (tmp_parity == 1)
    src_parity = ODD;
  else
    src_parity = EVEN;
  FORALLUPDIR(nu)
    g_veccomplexsum((complex*)&sourcelink[nu], 9);

  // Go: invert on point source for each color
  for (icol = 0; icol < 3; icol++) {
    invtime = -dclock();
    // Set point source in chi
    clear_latvec(F_OFFSET(chi), EVENANDODD);
    if (node_number(x, y, z, t) == mynode()) {
      i = node_index(x, y, z, t);
      lattice[i].chi.c[icol].real = 1.0;
    }

    // Invert: g_rand = M chi, psi = (MM^dag)^{-1} g_rand = M^{-1} chi
    clear_latvec(F_OFFSET(g_rand), EVENANDODD);
    clear_latvec(F_OFFSET(psi), EVENANDODD);
    cgs = mat_invert_uml(F_OFFSET(chi), F_OFFSET(psi),
                         F_OFFSET(g_rand), mass);
    cgTot += cgs;

    // Copy solution into the propagator
    FORALLSITES(i, s) {
      for (j = 0; j < 3; j++)
        s->prop.e[j][icol] = s->psi.c[j];
    }
    invtime += dclock();
    node0_printf("color %d at (%d %d %d %d): iters = %d secs = %.4g\n",
                 icol, x, y, z, t, cgs, invtime);
  }

  // Invert on point sources for each color at each y+nu
  FORALLUPDIR(nu) {
    // Move point source by nu
    x = x_src;  y = y_src;  z = z_src;  t = t_src;
    switch(nu) {
      case XUP : x = (x_src + 1) % nx; break;
      case YUP : y = (y_src + 1) % ny; break;
      case ZUP : z = (z_src + 1) % nz; break;
      case TUP : t = (t_src + 1) % nt; break;
      default  : node0_printf("ERROR in vacpol.c, aborting\n"); terminate(1);
    }

    for (icol = 0; icol < 3; icol++) {
      invtime = -dclock();
      // Set point source in chi
      clear_latvec(F_OFFSET(chi), EVENANDODD);
      if (node_number(x, y, z, t) == mynode()) {
        i = node_index(x, y, z, t);
        lattice[i].chi.c[icol].real = 1.0;
      }

      // Invert: g_rand = M chi, psi = (MM^dag)^{-1} g_rand = M^{-1} chi
      clear_latvec(F_OFFSET(g_rand), EVENANDODD);
      clear_latvec(F_OFFSET(psi), EVENANDODD);
      cgs = mat_invert_uml(F_OFFSET(chi), F_OFFSET(psi),
                           F_OFFSET(g_rand), mass);
      cgTot += cgs;

      // Copy solution into the propnu array
      FORALLSITES(i, s) {
        for (j = 0; j < 3; j++)
          s->propnu[nu].e[j][icol] = s->psi.c[j];
      }
      invtime += dclock();
      node0_printf("color %d at (%d %d %d %d): iters = %d secs = %.4g\n",
                   icol, x, y, z, t, cgs, invtime);
    }
  }

  FORALLUPDIR(nu) {               // Source direction -- call src pt y
    FORALLUPDIR(mu) {             // Sink direction   -- call sink pt x
      // gen_pt[0] is D_{x + mu; y + nu}
      // gen_pt[1] is D_{x + mu; y}
      tag = start_gather_site(F_OFFSET(propnu[nu]), sizeof(matrix),
                               mu, EVENANDODD, gen_pt[0]);
      tag2 = start_gather_site(F_OFFSET(prop), sizeof(matrix),
                               mu, EVENANDODD, gen_pt[1]);
      wait_gather(tag);
      wait_gather(tag2);

      // Construct the vector and axial two-point functions
      // NB the staggered phases and BCs are already in the links!
      FORALLSITES(i, s) {
        // [U_{mu}(x) D_{x + mu; y + nu}]^{dag} D_{x; y} U_{nu}(y)
        mult_nn(&(s->link[mu]), (matrix*)gen_pt[0][i], &tmat);
        mult_an(&tmat, &(s->prop), &tmat2);
        mult_nn(&tmat2, &(sourcelink[nu]), &tmat);
        tc = trace(&tmat);
        // Extra (-1)^{x + y} from Ddag
        // Axial correlator has additional (-1)^{x + y} -- cancels negative
        // Fixing sign in D&D negates axial correlator again
        if (s->parity != src_parity)
          CMULREAL(tc, neg, tc);
        tempaxi = tc;
        tempvec = tc;

        // Ddag_{x; y + nu} U_{mu}(x) D_{x + mu; y} U_{nu}(y)
        mult_an(&(s->propnu[nu]), &(s->link[mu]), &tmat);
        mult_nn(&tmat, (matrix*)gen_pt[1][i], &tmat2);
        mult_nn(&tmat2, &(sourcelink[nu]), &tmat);
        tc = trace(&tmat);
        // Extra (-1)^{x + y + nu} from Ddag
        if (s->parity == src_parity)
          CMULREAL(tc, neg, tc);
        CSUM(tempvec, tc);
        // Additional (-1)^{x + y} for axial correlator -- always negative
        if (s->parity != src_parity)
          CMULREAL(tc, neg, tc);
        CSUM(tempaxi, tc);

        // [U_{mu}(x) D_{x + mu; y}]^{dag} D_{x, y + nu} Udag_{nu}
        mult_nn(&(s->link[mu]), (matrix*)gen_pt[1][i], &tmat);
        mult_an(&tmat, &(s->propnu[nu]), &tmat2);
        mult_na(&tmat2, &(sourcelink[nu]), &tmat);
        tc = trace(&tmat);
        // Extra (-1)^{x + mu + y} from Ddag
        if (s->parity == src_parity)
          CMULREAL(tc, neg, tc);
        CSUM(tempvec, tc);
        // Additional (-1)^{x + y} for axial correlator -- always negative
        if (s->parity != src_parity)
          CMULREAL(tc, neg, tc);
        CSUM(tempaxi, tc);

        // Ddag_{x; y} U_{mu}(x) D_{x + mu; y + nu} Udag_{nu}
        mult_an(&(s->prop), &(s->link[mu]), &tmat);
        mult_nn(&tmat, (matrix*)gen_pt[0][i], &tmat2);
        mult_na(&tmat2, &(sourcelink[nu]), &tmat);
        tc = trace(&tmat);
        // Extra (-1)^{x + mu + y + nu} from Ddag
        // Axial correlator has additional (-1)^{x + y} -- cancels negative
        // Fixing sign in D&D negates axial correlator again
        if (s->parity != src_parity)
          CMULREAL(tc, neg, tc);
        CSUM(tempaxi, tc);
        CSUM(tempvec, tc);

        // Factor of 2 for each propagator from MILC conventions
        // Factor of 1/2 for each current, so do nothing
        // Make sure correlator is really real
        if (abs(tempvec.imag) > rsqmin) {
          printf("WARNING: Im[VEC^{%d%d}(%d %d %d %d)] = %.4g > %.4g\n",
                 mu, nu, s->x, s->y, s->z, s->t, tc.imag, rsqmin);
        }
        if (abs(tempaxi.imag) > rsqmin) {
          printf("WARNING: Im[AXI^{%d%d}(%d %d %d %d)] = %.4g > %.4g\n",
                 mu, nu, s->x, s->y, s->z, s->t, tc.imag, rsqmin);
        }
        s->vacpol[mu][nu] = tempvec.real;
        s->axial[mu][nu] = tempaxi.real;
      }
      cleanup_gather(tag);
      cleanup_gather(tag2);
    } // Sink direction mu

    // Contact terms
    // According to Y. Aoki, these should be the same for V and A

    // First contact term: x = src + nu
    x = x_src;  y = y_src;  z = z_src;  t = t_src;
    switch(nu) {
      case XUP : x = (x_src + 1) % nx; break;
      case YUP : y = (y_src + 1) % ny; break;
      case ZUP : z = (z_src + 1) % nz; break;
      case TUP : t = (t_src + 1) % nt; break;
      default  : node0_printf("ERROR in vacpol.c, aborting\n"); terminate(1);
    }
    contact[nu] = cmplx(0.0, 0.0);
    if (node_number(x, y, z, t) == mynode()) {
      i = node_index(x, y, z, t);
      mult_nn(&(lattice[i].prop), &(sourcelink[nu]), &tmat);
      contact[nu] = trace(&tmat);
      // Factor of 2 for each propagator from MILC conventions
      // Factor of 1/2 for each current, so do nothing
    }
    // Sum contact[nu] to broadcast it to all nodes
    g_complexsum(&(contact[nu]));

    // Second contact term: x = src
    if (node_number(x_src, y_src, z_src, t_src) == mynode()) {
      i = node_index(x_src, y_src, z_src, t_src);
      mult_na(&(lattice[i].propnu[nu]), &(lattice[i].link[nu]), &tmat);
      tc = trace(&tmat);
      // Factor of 2 for each prop from MILC conventions
      // Factor of 1/2 for current, so do nothing
      CSUB(contact[nu], tc, contact[nu]);

      // Make sure contact term is really real
      if (fabs(contact[nu].imag) > tol) {
        printf("WARNING: Im[contact^{%d}(%d %d %d %d)] = %.4g > %.4g\n",
               nu, x_src, y_src, z_src, t_src, contact[nu].imag, tol);
      }
      lattice[i].vacpol[nu][nu] -= contact[nu].real;
      lattice[i].axial[nu][nu] -= contact[nu].real;
    }
  } // Source direction nu
  return cgTot;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute and print the vector current divergence f_mu(x) - f_mu(x-mu)
void divergence() {
  int i, mu, nu;
  site *s;
  Real div, tol = sqrt(rsqmin);
  msg_tag* tag[4];

  FORALLUPDIR(mu) {
    tag[mu] = start_gather_site(F_OFFSET(vacpol[mu]), 4 * sizeof(Real),
                                OPP_DIR(mu), EVENANDODD, gen_pt[mu]);
  }
  FORALLUPDIR(mu)
    wait_gather(tag[mu]);

  FORALLSITES(i, s) {
    FORALLUPDIR(nu) {
      div = 0.0;
      FORALLUPDIR(mu)
        div += s->vacpol[mu][nu] - ((Real *)gen_pt[mu][i])[nu];
      if (fabs(div) > tol)
        printf("DIV^%d(%d %d %d %d) %.4g\n",
               nu, s->x, s->y, s->z, s->t, div);
    }
  }

  FORALLUPDIR(mu)
    cleanup_gather(tag[mu]);
}
// -----------------------------------------------------------------

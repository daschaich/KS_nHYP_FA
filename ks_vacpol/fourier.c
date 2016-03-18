// -----------------------------------------------------------------
// Calculate the Fourier transform of the two-point function
#include "vacpol_includes.h"

// Threshold to print warning about non-zero Im(VEC) and Im(Axi)
#define IMAG_TOL 1e-10
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Helper function getQ calculates Qhat^2 from spatial point
Real getQ(int l, int k, int j, int i) {
  Real tr, QSq = 0;

  if (l <= nt / 2)
    tr = l * 2 * PI / nt;
  else
    tr = (nt - l) * 2 * PI / nt;
  QSq += tr * tr;

  if (k <= nx / 2)
    tr = k * 2 * PI / nx;
  else
    tr = (nx - k) * 2 * PI / nx;
  QSq += tr * tr;

  if (j <= ny / 2)
    tr = j * 2 * PI / ny;
  else
    tr = (ny - j) * 2 * PI / ny;
  QSq += tr * tr;

  if (i <= nx / 2)
    tr = i * 2 * PI / nx;
  else
    tr = (nx - i) * 2 * PI / nx;
  QSq += tr * tr;

  return QSq;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function map matches a given Qsq to the appropriate component
// in the given array, adding it if necessary
int map(double Qsq, double *moms, int *lenpt) {
  int i, len = *lenpt;
  int test, key = floor(1e4 * Qsq);     // Avoid roundoff shenanigans

  // If Qsq already in array, return appropriate index
  if (len > 0) {
    for (i = 0; i < len; i++) {
      test = floor(1e4 * moms[i]);
      if (key == test)
        return i;
    }
  }

  // Otherwise add new Qsq to the array
  moms[len] = Qsq;
  *lenpt = len + 1;
  return len;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fourier transform and decomposition
// Return number of distinct Qsq
int fourier(double *longiv, double *transv,
            double *longia, double *transa) {

  register int ii, i, j, k, l, m, n;
  int label, len = 0;
  Real p1, p2, p3, p4, px, py, pz, pt, q[4];
  Real max_x = nx, max_y = ny, max_z = nz, max_t = nt;   // Do not like
  Real pdotx, theta, theta2, Qsq;
  double dtime = -dclock();
  complex cc, phase;
  dcomplex sum, ax_sum;
  site *s;

  // Vector and axial transverse and longitudinal components
  for (i = 0; i < volume; i++) {    // Initialize to zero
    longiv[i] = 0.0;
    transv[i] = 0.0;
    longia[i] = 0.0;
    transa[i] = 0.0;
  }

  // Loop over momenta
  for (l = 0; l < max_t; l++) {
    p4 = l * 2 * PI / nt;
    for (k = 0; k < max_z; k++) {
      p3 = k * 2 * PI / nz;
      for (j = 0; j < max_y; j++) {
        p2 = j * 2 * PI / ny;
        for (i = 0; i < max_x; i++) {
          p1 = i * 2 * PI / nx;
          Qsq = getQ(l, k, j, i);
          if (Qsq > Qmax)
            continue;
          for (m = 0; m < 4; m++) {
            // Extra phases for point-split currents:
            //   theta  = q_mu / 2 = pi n_mu / L_mu
            //   theta2 = q_nu / 2 = pi n_nu / L_nu
            q[XUP] = (m == XUP) ? p1 : 0;
            q[YUP] = (m == YUP) ? p2 : 0;
            q[ZUP] = (m == ZUP) ? p3 : 0;
            q[TUP] = (m == TUP) ? p4 : 0;
            theta = 0.5 * (q[XUP] + q[YUP] + q[ZUP] + q[TUP]);
            for (n = 0; n < 4; n++) {
              q[XUP] = (n == XUP) ? p1 : 0;
              q[YUP] = (n == YUP) ? p2 : 0;
              q[ZUP] = (n == ZUP) ? p3 : 0;
              q[TUP] = (n == TUP) ? p4 : 0;
              theta2 = 0.5 * (q[XUP] + q[YUP] + q[ZUP] + q[TUP]);

              // Loop over sites on node, accounting for source location
              sum = cmplx(0, 0);
              ax_sum = cmplx(0, 0);
              FORALLSITES(ii, s) {
                px = (s->x - x_src) * p1;
                py = (s->y - y_src) * p2;
                pz = (s->z - z_src) * p3;
                pt = (s->t - t_src) * p4;
                pdotx = px + py + pz + pt;

                // Negative sign is convention of restrict_fourier
                phase = ce_itheta(-pdotx - theta + theta2);
                CMULREAL(phase, lattice[ii].vacpol[m][n], cc);
                CSUM(sum, cc);
                CMULREAL(phase, lattice[ii].axial[m][n], cc);
                CSUM(ax_sum, cc);
              }
              g_dcomplexsum(&sum);
              g_dcomplexsum(&ax_sum);

              // Check that sum and ax_sum are really real
              if (abs(sum.imag) > IMAG_TOL) {
                printf("WARNING: sum(%.4g) = %.4g > %.4g)\n",
                       Qsq, sum.imag, IMAG_TOL);
              }
              if (abs(ax_sum.imag) > IMAG_TOL) {
                printf("WARNING: ax_sum(%.4g) = %.4g > %.4g)\n",
                       Qsq, ax_sum.imag, IMAG_TOL);
              }

              // Decompose into transverse and longitudinal parts
              // labelled by Qsq
              if (Qsq != 0) {
                label = map(Qsq, moms, &len);
                count[label]++;               // For averaging
                longiv[label] -= q[m] * sum.real * q[n] / Qsq;
                longia[label] -= q[m] * ax_sum.real * q[n] / Qsq;
                if (m == n) {                 // Trace
                  transv[label] += sum.real / 3;
                  transa[label] += ax_sum.real / 3;
                }
              }
//              node0_printf("VACPOL0 %d %d mom %d %d %d %d %.4g %.4g\n",
//                           m, n, i, j, k, l,
//                           (double)sum.real, (double)sum.imag);
            }
          }
        }
      }
    }
  }

  // Normalize output
  for (i = 0; i < len; i++) {
    transv[i] += longiv[i] / 3;
    transa[i] += longia[i] / 3;
    longiv[i] /= count[i];
    longia[i] /= count[i];
    transv[i] /= count[i];
    transa[i] /= count[i];
  }

  dtime += dclock();
  node0_printf("Fourier secs = %.4g\n", dtime);
  return len;
}
// -----------------------------------------------------------------

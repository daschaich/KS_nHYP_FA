// -----------------------------------------------------------------
// Helper routines for nHYP blocking

// The general calculation in compute_RS, which follows hep-lat/0702028,
// crashes on the following conditions
// 1) Q has three equal eigenvalues (e.g. ordered start)
// 2) Q has zero eigenvalue

// Here we solve this by
// 1) add IR regulator,   Q = Omega^dag Omega + IR_STAB
// 2) use series approx for small S (and R)
// If the eigenvalues are degenerate to O(epsilon),
// then S=O(epsilon^2) and R=O(epsilon^3).

// IR_STAB and EPS_SQ must be set by the application!
// (Typically in defines.h)
#include "nhyp_includes.h"

static inline double sqr(double x) {return x * x;}
static inline double pow3(double x) {return x * x * x;}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void compute_RS(matrix *Q, double *R, double *S) {
  double dtmp;
  complex ctmp;
  double e22 = Q->e[1][1].real - Q->e[0][0].real;
  double e33 = Q->e[2][2].real - Q->e[0][0].real;
  double q012, q022, q122;

  q012 = cabs_sq(&Q->e[0][1]);
  q022 = cabs_sq(&Q->e[0][2]);
  q122 = cabs_sq(&Q->e[1][2]);

  dtmp = q012 + q022 + q122 +(sqr(e22 - e33) + e22 * e33) / 3;
  *S = dtmp / 3;

  CMUL(Q->e[0][1], Q->e[1][2], ctmp);
  dtmp = ctmp.real * Q->e[2][0].real - ctmp.imag * Q->e[2][0].imag;

  dtmp = 6 * dtmp
       + q012 * (e22-2.*e33) + q022 * (e33 - 2 * e22) + q122 * (e33 + e22)
       + (e22 * e22 * (2 * e22 - 3 * e33)
          + e33 * e33 * (2 * e33 - 3 * e22)) / 9;

  *R = dtmp / 6;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the Cayley Hamilton coefficients for the inverse square root
//   1 / sqrt(Q) = f[0] I + f[1] Q +f[2] Q^2
// The b[][] matrix contains the derivatives of f wrt to the traces of Q
//   b[i][j] = d f[i] / d c[j] with c[j] = 1/(j+1) trace Q^(j+1)
// The flag compute_b switches on/off the computation of these derivatives
#ifndef NHYP_DEBUG
void compute_fhb(matrix *Q, Real *f, Real b[3][3], int compute_b)
#else
void compute_fhb(matrix *Omega, matrix *Q,
                 Real *f, Real b[3][3], int compute_b)
#endif
{
  double S, R;
  double g_avr, g_avr_sq, g_avr_cube, two_sqrt_S, theta;
  double g[3];
  double sg0,sg1,sg2;
  double pi23 = 2.09439510239319549; // 2pi / 3
  double costheta, costheta_p_23pi, costheta_m_23pi;
  double u0, u1, p, den;
  double u02, u03, u04, u12, u13, u14, u16;
  double acos_input;

  compute_RS(Q, &R, &S);
  // Average eigenvalue, strictly positive for any positive IR_STAB
  g_avr = (Q->e[0][0].real + Q->e[1][1].real + Q->e[2][2].real) / 3;
  g_avr_sq = g_avr * g_avr;
  g_avr_cube = g_avr * g_avr_sq;

#ifdef NHYP_DEBUG
  int nhyp_debug_flag = 0;

  // Check that S is non-negative
  // This condition could be encountered on any node
  // Presumably it is rare enough that no two nodes will attempt
  // to print out simultaneously
  if (S < 0) {
    printf("NHYP_DEBUG(1): S<0\n  S = %.12e, R = %.12e, g_avr = %.12e\n",
           S, R, g_avr);
    nhyp_debug_flag = 1;
  }
  if (g_avr < 0) {
    printf("NHYP_DEBUG(2): g_avr<0\n  S = %.12e, R = %.12e, g_avr = %.12e\n",
           S, R, g_avr);
    nhyp_debug_flag = 1;
  }
#endif

  // If the three eigenvalues are exactly degenerate, R = S = 0
  // When S / g_avr_sq (and R / g_avr_cube) are small we apply an approximation.
  if (S / g_avr_sq > EPS_SQ) {
    // General case
    two_sqrt_S = 2 * sqrt(S);
    acos_input = R / (S * sqrt(S));

#ifdef NHYP_DEBUG
    if (S == 0) {
      printf("NHYP_DEBUG(3): S=0 in general case; ");
      printf("check EPS_SQ\n  S = %.12e, R = %.12e, acos_input = %.12e\n",
             S, R, acos_input);
      nhyp_debug_flag = 1;
    }
#endif

    // Absolute value of R/(S*sqrt(S) could be slightly above 1 due
    // to single/double precision issues.  For tiny deviation, set it
    // to +-1.  Issue warning for deviations large than 1e-5
    if (acos_input > 1) {
#ifdef NHYP_DEBUG
      if (acos_input > 1.00001) {
        printf("NHYP_DEBUG_WARNING(4): acos_input > 1.00001\n");
        printf("  S = %.12e, R = %.12e, acos_input = %.12e\n",
               S, R, acos_input);
        nhyp_debug_flag = 1;
      }
#endif
      acos_input = 1;
    }

    if (acos_input < -1) {
#ifdef NHYP_DEBUG
      if (acos_input < -1.00001) {
        printf("NHYP_DEBUG_WARNING(4): acos_input < -1.00001\n");
        printf("  S = %.12e, R = %.12e, acos_input = %.12e\n",
               S, R, acos_input);
        nhyp_debug_flag=1;
      }
#endif
      acos_input = -1;
    }

    // Now we can safely call acos
    theta = acos(acos_input) / 3.0;
    costheta = cos(theta);
    costheta_p_23pi = cos(theta + pi23);
    costheta_m_23pi = cos(theta - pi23);

    // Eigenvalues of Q
    g[0] = g_avr + two_sqrt_S * costheta;
    g[1] = g_avr + two_sqrt_S * costheta_m_23pi;
    g[2] = g_avr + two_sqrt_S * costheta_p_23pi;

#ifdef NHYP_DEBUG
    if (g[0] < 0 || g[1] < 0 || g[2] < 0) {
      printf("NHYP_DEBUG(5): Negative eigenvalue for Q\n");
      printf("  S = %.12e, R = %.12e, acos_input = %.12e, ",
             S, R, acos_input);
      printf("g[0] = %.12e, g[1] = %.12e, g[2] = %.12e\n",
             g[0], g[1], g[2]);
      nhyp_debug_flag = 1;
    }
#endif

    sg0 = sqrt(g[0]);
    sg1 = sqrt(g[1]);
    sg2 = sqrt(g[2]);

    // The coefficients and derivatives are symmetric under permutations
    // of g0, g1, g2. Since they are also polynomial in the sqrt's, they can be
    // reduced to elementary symmetric functions u0, u1, p
    u0 = sg0+sg1+sg2;
    u1 = sg0*sg1+sg0*sg2+sg1*sg2;
    p  = sg0*sg1*sg2;
  }

  // Otherwise, S / g_avr_sq < EPS_SQ
  // Use a series approximation for u0, u1, p, in terms of R and S
  else {
    S /= g_avr_sq;
    R /= g_avr_cube;
    u0 = sqrt(g_avr) * (3 - 3 * S / 4 + 3 * R / 8);
    u1 = g_avr*(3 - 9 * S / 4 + 7 * R / 8);
    p  = sqrt(g_avr_cube) * (1 - 3 * S / 2 - R / 2);
  }

  // Some powers
  u02 = u0 * u0;
  u12 = u1 * u1;

  den = p * (u0 * u1 - p);

#ifdef NHYP_DEBUG
  if (den == 0) {
    printf("NHYP_DEBUG(6): den=0\n");
    printf("  S = %.12e, R = %.12e, g_avr = %.12e, den = %.12e\n",
           S, R, g_avr, den);
    nhyp_debug_flag = 1;
  }
#endif

  // Den is never zero for any positive IR_STAB
  f[0] = (-p * (u02 + u1) + u0 * u12) / den;
  f[1] = (-p + u0* (2 * u1 - u02)) / den;
  f[2] = u0 / den;

  // Most of the time, the coefficients are all we need
  if (compute_b == 0) {
#ifdef NHYP_DEBUG
    if (nhyp_debug_flag == 1) {
      printf("NHYP_DEBUG ENCOUNTERED.  TERMINATING.\n");
      terminate(1);
    }
#endif
    return;
  }

  u03 = u02 * u0;
  u04 = u03 * u0;

  u13 = u12 * u1;
  u14 = u12 * u12;
  u16 = u12 * u14;

  den = 2 * pow3(p * (u0 * u1 - p));

#ifdef NHYP_DEBUG
  if (den == 0) {
    printf("NHYP_DEBUG(7): den=0\n");
    printf("  S = %.12e, R = %.12e, g_avr = %.12e, den = %.12e\n",
           S, R, g_avr, den);
    nhyp_debug_flag = 1;
  }
#endif

  b[0][0] = (-u03 * u16
             + p * (3 * u02 * (u02 + u1) * u14
                    + p * (-3 * u0 * (4 * u02 + u1) * u13
                           +p * (u02 * (-u04 + 3 * u02 * u1 + 16 * u12) + u13
                                 + p * (-4 * u0 * (u02 + 2 * u1) + p))))) / den;

  b[0][1] = (u03 * (u02 - 2 * u1) * u14
             + p * (-(u02 * (u04 + u02 * u1 - 6 * u12) * u12)
                    + p * (-(u03 * (u04 - 6 * u02 * u1 + 6 * u12)
                           + 6 * u13 * u0)
                           + p * (5 * u02 * (-u02 + 2 * u1)
                                  + 2 * u12 - 3 * u0 * p)))) / den;
  b[1][0] = b[0][1];

  b[0][2] = (-u03 * u14
             + p * (u02 * (u02 + 3 * u1) * u12
                    + p * (u0 * (u04 - 4 * u02 * u1 - 3 * u12)
                           + p * (4 * u02 + u1)))) / den;
  b[2][0] = b[0][2];

  b[1][1] = (-u03 * sqr(u02 - 2 * u1) * u12
             + p * (-((u02 - 3 * u1) * sqr(u03 - 2 * u0 * u1))
                    + p * (-(u0 * (5 * u02 - 6 * u1) * (u02 - 2 * u1))
                           + p * (-3 * u02 + 3 * u1)))) / den;

  b[1][2] = (u03 * (u02 - 2 * u1) * u12
             + p * (u02 * (u04 - 5 * u02 * u1 + 6 * u12)
                    + p * (4 * u03 - 6 * u0 * u1 + p))) / den;
  b[2][1] = b[1][2];

  b[2][2] = (u0 * (-u02 * u12 + p * (-u03 + 3 * u0 * u1 - 3 * p))) / den;

#ifdef NHYP_DEBUG
  if (nhyp_debug_flag == 1) {
    printf("NHYP_DEBUG ENCOUNTERED.  TERMINATING.\n");
    terminate(1);
  }
#endif
}
// -----------------------------------------------------------------

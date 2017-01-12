// -----------------------------------------------------------------
// Check unitarity of the link matrices, terminate if not unitary
#include "generic_includes.h"

#define TOLERANCE 0.0001
#define STRONG    // Check row orthogonality as well as norms
//#define UNIDEBUG
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_su3(su3_matrix *c) {
  register int i, j, k;
  register Real ar, ai, ari, max = 0.0;

  // Check normalization of each row
  for (i = 0; i < 3; i++) {
    ar = (*c).e[i][0].real * (*c).e[i][0].real
       + (*c).e[i][0].imag * (*c).e[i][0].imag
       + (*c).e[i][1].real * (*c).e[i][1].real
       + (*c).e[i][1].imag * (*c).e[i][1].imag
       + (*c).e[i][2].real * (*c).e[i][2].real
       + (*c).e[i][2].imag * (*c).e[i][2].imag;
    ar = fabs(sqrt((double)ar) - 1.0);
    if (max < ar)   // Maximum normalization deviation from 1
      max = ar;
  }

#ifdef STRONG
  // Test orthogonality of row i and row j
  for (i = 0; i < 3; i++) {
    for (j = i + 1; j < 3; j++) {
      ar = 0.0;   // Real part of i dot j
      ai = 0.0;   // Imag part of i dot j
      for (k = 0; k < 3; k++) {
        ar += (*c).e[i][k].real * (*c).e[j][k].real
            + (*c).e[i][k].imag * (*c).e[j][k].imag;
        ai += (*c).e[i][k].real * (*c).e[j][k].imag
            - (*c).e[i][k].imag * (*c).e[j][k].real;
      }
      ari = sqrt((double)(ar * ar + ai * ai));
      if (max < ari)
        max = ari;
    }
  }
#endif

  return max;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_unitarity() {
  register int i,dir;
  int ii, jj, status = 0;
  register site *s;
  register su3_matrix *mat;
  Real deviation, max_deviation = 0.0;
  double av_deviation = 0.0;
  union {
    Real fval;
    int ival;
  } ifval;

  FORALLSITES(i, s) {
#ifdef SCHROED_FUN
    for (dir=XUP; dir<=TUP; dir++) if (dir==TUP || s->t>0){
#else
    for (dir=XUP; dir<=TUP; dir++) {
#endif
      mat = (su3_matrix *)&(s->link[dir]);
      deviation = check_su3(mat);
      if (deviation > TOLERANCE) {
        printf("Unitarity problem on node %d, site %d, dir %d, deviation=%f\n",
               mynode(), i, dir, deviation);
        printf("SU3 matrix:\n");
        for (ii = 0; ii < 3; ii++) {
          for (jj = 0; jj < 3; jj++) {
            printf("%f ", (*mat).e[ii][jj].real);
            printf("%f ", (*mat).e[ii][jj].imag);
          }
          printf("\n");
        }
        printf("repeat in hex:\n");
        for (ii = 0; ii < 3; ii++) {
          for (jj = 0; jj < 3; jj++) {
            ifval.fval = (*mat).e[ii][jj].real;
            printf("%08x ", ifval.ival);
            ifval.fval = (*mat).e[ii][jj].imag;
            printf("%08x ", ifval.ival);
          }
          printf("\n");
        }
        printf("  \n\n");
        fflush(stdout);
        status++;
        break;
      }
      if (status)
        break;
      if (max_deviation < deviation)
        max_deviation = deviation;
      av_deviation += deviation * deviation;
    }
    if (status)
      break;
  }

  // Poll nodes for problems
  g_intsum(&status);
  if (status > 0) {
    node0_printf("Terminated due to unacceptable unitarity violation(s)\n");
    terminate(1);
  }

  av_deviation = sqrt(av_deviation / (4.0 * volume));
#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %.4g, ave %.4g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE)
    printf("Unitarity problem on node %d, maximum deviation = %.4g\n",
           mynode(), max_deviation);
  return max_deviation;
}
// -----------------------------------------------------------------

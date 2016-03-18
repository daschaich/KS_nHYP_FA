// -----------------------------------------------------------------
// Check unitarity of the link matrices, terminate if not unitary
#include "generic_includes.h"

#define TOLERANCE (0.0001)
#define STRONG    // Check row orthogonality as well as norms
/*#define UNIDEBUG */
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_su3(su3_matrix *c) {
  register int i;
  register Real ar, ai, ari, max;

  /* first normalize row */
  for (i = 0, max = 0.0; i < 3; ++i) {
    ar = (*c).e[i][0].real * (*c).e[i][0].real +    /* sum of squares of row */
      (*c).e[i][0].imag * (*c).e[i][0].imag +
      (*c).e[i][1].real * (*c).e[i][1].real +
      (*c).e[i][1].imag * (*c).e[i][1].imag +
      (*c).e[i][2].real * (*c).e[i][2].real +
      (*c).e[i][2].imag * (*c).e[i][2].imag;
    ar = fabs(sqrt((double)ar) - 1.0);
    if (max < ar)
      max = ar;
  }

#ifdef STRONG
  /* Test orthogonality of row 0 and row 1 */
  ar = (*c).e[0][0].real * (*c).e[1][0].real +     /* real part of 0 dot 1 */
    (*c).e[0][0].imag * (*c).e[1][0].imag +
    (*c).e[0][1].real * (*c).e[1][1].real +
    (*c).e[0][1].imag * (*c).e[1][1].imag +
    (*c).e[0][2].real * (*c).e[1][2].real +
    (*c).e[0][2].imag * (*c).e[1][2].imag;
  ai = (*c).e[0][0].real * (*c).e[1][0].imag -     /* imag part of 0 dot 1 */
    (*c).e[0][0].imag * (*c).e[1][0].real +
    (*c).e[0][1].real * (*c).e[1][1].imag -
    (*c).e[0][1].imag * (*c).e[1][1].real +
    (*c).e[0][2].real * (*c).e[1][2].imag -
    (*c).e[0][2].imag * (*c).e[1][2].real;

  ari = sqrt((double)(ar*ar + ai*ai));
  if (max<ari) max=ari;

  /* Test orthogonality of row 0 and row 2 */
  ar = (*c).e[0][0].real * (*c).e[2][0].real +     /* real part of 0 dot 1 */
    (*c).e[0][0].imag * (*c).e[2][0].imag +
    (*c).e[0][1].real * (*c).e[2][1].real +
    (*c).e[0][1].imag * (*c).e[2][1].imag +
    (*c).e[0][2].real * (*c).e[2][2].real +
    (*c).e[0][2].imag * (*c).e[2][2].imag;
  ai = (*c).e[0][0].real * (*c).e[2][0].imag -     /* imag part of 0 dot 1 */
    (*c).e[0][0].imag * (*c).e[2][0].real +
    (*c).e[0][1].real * (*c).e[2][1].imag -
    (*c).e[0][1].imag * (*c).e[2][1].real +
    (*c).e[0][2].real * (*c).e[2][2].imag -
    (*c).e[0][2].imag * (*c).e[2][2].real;

  ari = sqrt((double)(ar*ar + ai*ai));
  if (max<ari) max=ari;

  /* Test orthogonality of row 1 and row 2 */
  ar = (*c).e[1][0].real * (*c).e[2][0].real +     /* real part of 0 dot 1 */
    (*c).e[1][0].imag * (*c).e[2][0].imag +
    (*c).e[1][1].real * (*c).e[2][1].real +
    (*c).e[1][1].imag * (*c).e[2][1].imag +
    (*c).e[1][2].real * (*c).e[2][2].real +
    (*c).e[1][2].imag * (*c).e[2][2].imag;
  ai = (*c).e[1][0].real * (*c).e[2][0].imag -     /* imag part of 0 dot 1 */
    (*c).e[1][0].imag * (*c).e[2][0].real +
    (*c).e[1][1].real * (*c).e[2][1].imag -
    (*c).e[1][1].imag * (*c).e[2][1].real +
    (*c).e[1][2].real * (*c).e[2][2].imag -
    (*c).e[1][2].imag * (*c).e[2][2].real;

  ari = sqrt((double)(ar*ar + ai*ai));
  if (max < ari)
    max = ari;
#endif

  return max;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_unitarity() {
  register int i,dir;
  register site *s;
  register su3_matrix *mat;
  int ii, jj, status = 0;
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
        printf("  \n \n");
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

  av_deviation = sqrt(av_deviation / (4.0 * i));
#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %g, ave %g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE)
    printf("Unitarity problem on node %d, maximum deviation=%f\n",
           mynode(), max_deviation);
  return max_deviation;
}
// -----------------------------------------------------------------

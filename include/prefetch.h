#ifndef _PREFETCH_H
#define _PREFETCH_H

/***************************************************************************/
/*          Cache manipulation for a variety of architectures              */
/***************************************************************************/

/* Abbreviations for data types:
   M  SU(3) matrix
   V  SU(3) vector
*/

#if (defined P3 || defined P4) && defined __GNUC__

#include "../include/prefetch_asm.h"

#else

#ifndef PREFETCH

/***************************************************************************/
/*              Ignore Cache Manipulation Macros                           */
/***************************************************************************/

#define prefetch_M(a0)
#define prefetch_V(a0)
#define prefetch_VV(a0, a1)
#define prefetch_VVV(a0, a1, a2)
#define prefetch_VVVV(a0, a1, a2, a3)
#define prefetch_VVVVV(a0, a1, a2, a3, a4)
#define prefetch_4MVVVV(a0, a1, a2, a3, a4)
#define prefetch_4MV4V(a0, a1, a2)

#else

/***************************************************************************/
/*                  Cache Manipulation Macros                              */
/*                Prefetch via subroutine calls                            */
/***************************************************************************/

void _prefetch_M(matrix *m);
void _prefetch_V(vector *v);
void _prefetch_VV(vector *v, vector *v);
void _prefetch_VVV(vector *v, vector *v, vector *v);
void _prefetch_VVVV(vector *v, vector *v, vector *v, vector *v);
void _prefetch_VVVVV(vector *v, vector *v, vector *v, vector *v, vector *v);
void _prefetch_4MVVVV(matrix *m, vector *v, vector *v, vector *v, vector *v);
void _prefetch_4MV4V(matrix *m, vector *v, vector *v);

#define prefetch_M(a0)                    _prefetch_M(a0)
#define prefetch_V(a0)                    _prefetch_V(a0)
#define prefetch_VV(a0, a1)                _prefetch_VV(a0, a1)
#define prefetch_VVV(a0, a1, a2)            _prefetch_VVV(a0, a1, a2)
#define prefetch_VVVV(a0, a1, a2, a3)        _prefetch_VVVV(a0, a1, a2, a3)
#define prefetch_VVVVV(a0, a1, a2, a3, a4)    _prefetch_VVVVV(a0, a1, a2, a3, a4)
#define prefetch_4MVVVV(a0, a1, a2, a3, a4)   _prefetch_4MVVVV(a0, a1, a2, a3, a4)
#define prefetch_4MV4V(a0, a1, a2)          _prefetch_4MV4V(a0, a1, a2)

#endif /* NOPREFETCH */

#endif /* P4 or P3 */

#endif /* _PREFETCH_H */

#ifndef _PREFETCH_H
#define _PREFETCH_H

#include "../include/su3.h"  // Define vector and matrix

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

#if ! defined PREFETCH

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

void _prefetch_M(matrix *);
void _prefetch_V(vector *);
void _prefetch_VV(vector *, vector *);
void _prefetch_VVV(vector *, vector *, vector *);
void _prefetch_VVVV(vector *, vector *, vector *, vector *);
void _prefetch_VVVVV(vector *, vector *, vector *, vector *, vector *);
void _prefetch_4MVVVV(matrix *, vector *, vector *, vector *, vector *);
void _prefetch_4MV4V(matrix *, vector *, vector *);

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

// -----------------------------------------------------------------
// Cache prefetching with 8-byte alignment and 64-byte cache line
// This is the vanilla C version for use when there is no assembly code
// Some compilers will not do these null fetches
// Compiling with the -g option may help with that
// Otherwise you may have to get the compiler-generated assembly code
// and edit by hand

// Fetch addr, addr + cache, ..., addr + size - align
// to get datum of size "size" starting at address "addr"
// with alignment "align" and cache line "cache"

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/prefetch.h"

// 16 floats in cache line
//  2 floats per alignment boundary

// matrix: 18 Reals
#define _pftch_M(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16);

// vector: 6 Reals
#define _pftch_V(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 4 );

// 4 vectors: 24 Reals
#define _pftch_4V(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16); \
  dummy = *((Real *)(a) + 22);

// 4 matrices: 72 Reals
#define _pftch_4M(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16); \
  dummy = *((Real *)(a) + 32); \
  dummy = *((Real *)(a) + 48); \
  dummy = *((Real *)(a) + 64); \
  dummy = *((Real *)(a) + 70);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void _prefetch_M(matrix *a0) {
  register Real dummy;

  _pftch_M(a0);
}

void _prefetch_V(vector *a0) {
  register Real dummy;

  _pftch_V(a0);
}

void _prefetch_VV(vector *a0, vector *a1) {
  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
}

void _prefetch_VVV(vector *a0, vector *a1, vector *a2) {
  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
  _pftch_V(a2);
}

void _prefetch_VVVV(vector *a0, vector *a1, vector *a2,
                    vector *a3) {

  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
  _pftch_V(a2);
  _pftch_V(a3);
}

void _prefetch_VVVVV(vector *a0, vector *a1, vector *a2,
         vector *a3, vector *a4) {
  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
  _pftch_V(a2);
  _pftch_V(a3);
  _pftch_V(a4);
}

void _prefetch_4MVVVV(matrix *a0, vector *a1, vector *a2,
                      vector *a3, vector *a4) {

  register Real dummy;

  _pftch_4M(a0);
  _pftch_V(a1);
  _pftch_V(a2);
  _pftch_V(a3);
  _pftch_V(a4);
}

void _prefetch_4MV4V(matrix *a0, vector *a1, vector *a2) {
  register Real dummy;

  _pftch_4M(a0);
  _pftch_V(a1);
  _pftch_4V(a2);
}

#undef _pftch_M
#undef _pftch_V
#undef _pftch_4V
#undef _pftch_4M
// -----------------------------------------------------------------

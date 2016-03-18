/********************** show_opts.c *************************************/
/* List options selected in the compilation */

#include "generic_includes.h"

void show_generic_opts( void ){

#if (PRECISION==1)
  node0_printf("Generic single precision\n");
#else
  node0_printf("Generic double precision\n");
#endif

#ifdef C_GLOBAL_INLINE
  node0_printf("C_GLOBAL_INLINE\n");
#endif

#ifdef SSE_GLOBAL_INLINE
  node0_printf("SSE_GLOBAL_INLINE\n");
#endif
}

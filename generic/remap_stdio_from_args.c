// -----------------------------------------------------------------
// Reopen the first three command-line arguments
// as stdin, stdout and stderr, respectively
// Call only after consuming and removing all other command-line arguments
// This is in case input and output redirection is not an option
// Has been made fancier 'upstream', which may be worth merging...

#include "generic_includes.h"

int remap_stdio_from_args(int argc, char *argv[]) {
  FILE *fp;

  // stdin is remapped only on node 0 on any machine
  if (argc > 1 && mynode() == 0) {
    fp = freopen(argv[1], "r", stdin);
    if (fp == NULL) {
      node0_printf("Can't open stdin file %s for reading\n", argv[1]);
      return 1;
    }
  }
  if (argc > 2) {
    fp = freopen(argv[2], "w", stdout);
    if (fp == NULL) {
      node0_printf("Can't open stdout file %s for writing\n", argv[2]);
      return 1;
    }
  }
  if (argc > 3) {
    fp = freopen(argv[3], "w", stderr);
    if (fp == NULL) {
      node0_printf("Can't open stderr file %s for writing\n", argv[3]);
      return 1;
    }
  }
  return 0;
}
// -----------------------------------------------------------------

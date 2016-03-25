// -----------------------------------------------------------------
// Main procedure for SU(3) lattice loading check
#define CONTROL
#include "load_includes.h"

int main(int argc, char *argv[]) {
  int prompt;

  // Setup
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc, &argv);
  g_sync();
  prompt = setup();

  // Load input and run (loop removed)
  readin(prompt);
  return 0;
}
// -----------------------------------------------------------------

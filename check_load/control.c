// -----------------------------------------------------------------
// Main procedure for SU(3) lattice loading check
#define CONTROL
#include "load_includes.h"

int main(int argc, char *argv[]) {
  int prompt;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run (loop removed)
  readin(prompt);
  return 0;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Main procedure for lattice doubling
#define CONTROL
#include "double_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[])  {
  register int i, dir;
  register site *s;
  int j, prompt, *index;
  int ox, oy, oz, ot;   // Original nx, ny, nz, nt, to avoid lots of " / 2"
  double ssplaq, stplaq, dtime;
  matrix **bak;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();

  // Load input and run (loop removed)
  if (readin(prompt) == 0) {
    dtime = -dclock();

    // Check initial plaquette
    plaquette(&ssplaq, &stplaq);
    node0_printf("Initial plaquettes: %.8g %.8g\n", ssplaq, stplaq);

    // Allocate backup space and save links along with corresponding index
    bak = (matrix **)malloc(volume * sizeof(matrix *));
    index = (int *)malloc(volume * sizeof(int));
    for (j = 0; j <= volume; j++)
      bak[j] = (matrix *)malloc(4 * sizeof(matrix));

    FORALLSITES(i, s) {
      for (dir = XUP; dir <= TUP; dir++)
        mat_copy(&(s->link[dir]), &(bak[i][dir]));

      j = s->x * nt * nz * ny + s->y * nt * nz + s->z * nt + s->t;
      index[j] = i;
    }

    // Reset with all lattice dimensions doubled
    free_lattice();
    ox = nx;  oy = ny;  oz = nz;  ot = nt;
    nx *= 2;  ny *= 2;  nz *= 2;  nt *= 2;
    volume *= 16;
    setup_layout();
    make_lattice();
    reset_machine();    // Clear old gathers
    make_nn_gathers();

    // Create dummy lattice of appropriate size
    startflag = FRESH;    // Makes IO ignore startfile
    sprintf(startfile, "NULL"); // Check that this doesn't matter
    startlat_p = reload_lattice(startflag, startfile);

    // Copy saved links into new lattice
    FORALLSITES(i, s) {
      // Where we are in the original volume -- copy from this saved index
      j = (s->x % ox) * ot * oz * oy + (s->y % oy) * ot * oz
                  + (s->z % oz) * ot + (s->t % ot);
      for (dir = XUP; dir <= TUP; dir++)
        mat_copy(&(bak[index[j]][dir]), &(s->link[dir]));
    }

    // Check final plaquette
    plaquette(&ssplaq, &stplaq);
    node0_printf("Final plaquettes: %.8g %.8g\n", ssplaq, stplaq);

    // Save lattice
    if (saveflag != FORGET)
      save_lattice(saveflag, savefile, stringLFN);

    node0_printf("RUNNING COMPLETED\n");
    dtime += dclock();
    node0_printf("Time = %.4g seconds\n", dtime);
    fflush(stdout);
  } // readin(prompt) == 0
  return 0;
}
// -----------------------------------------------------------------

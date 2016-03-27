// -----------------------------------------------------------------
// Main procedure for staggered vacuum polarization calculation
#define CONTROL
#include "vacpol_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char **argv) {
  int i, prompt, pnt[4], idx, mu, nu;
  double ssplaq, stplaq, dtime;
  char vecName[80], axiName[80], md[80];
  site *s;

  // Setup
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc, &argv);
  g_sync();
  prompt = setup();

  // Load input and run (loop removed)
  // Does readin do a rephase(ON) on the read-in lattice?
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

#ifndef PBOUND
  node0_printf("\nANTIPERIODIC BC in time direction\n");
#else
  node0_printf("\nPERIODIC BC in time direction\n");
#endif

  // Smear, then calculate vector and axial vacuum polarization tensors
  block_and_fatten();
  rephase(OFF);
  plaquette(&ssplaq, &stplaq);
  rephase(ON);
  node0_printf("Plaquettes after smearing: %.8g %.8g\n", ssplaq, stplaq);
  total_iters = vacuum_polarization();
  divergence();       // Check divergence -- only prints if non-zero
  leanlinks();

  // Now we have both the vector vacuum polarization in s->vacpol[mu][nu]
  // and the axial vacuum polarization in s->axial[mu][nu]
  // Save them in QDP format for offline Fourier transform analysis
  QDP_initialize(&argc, &argv);
  QDP_profcontrol(0);           // No profiling
  QDP_set_latsize(ndim, latsize);
  QDP_create_layout();

  QLA_Real *tempaxi, *tempvec;
  QDP_Real *axi[4], *vec[4];
  QDP_Writer *vecFile, *axiFile;
  QDP_String *meta;

  // Open files for printing
  meta = QDP_string_create();
  QDP_string_set(meta, "Staggered axial correlators");
  sprintf(axiName, "KSaxi%s", outpat);
  axiFile = QDP_open_write(meta, axiName, QDP_SINGLEFILE);

  QDP_string_set(meta, "Staggered vector correlators");
  sprintf(vecName, "KSvec%s", outpat);
  vecFile = QDP_open_write(meta, vecName, QDP_SINGLEFILE);

  // Translate from MILC layout to QDP format
  FORALLUPDIR(mu) {
    axi[mu] = QDP_create_R();
    vec[mu] = QDP_create_R();
  }
  FORALLUPDIR(mu) {
    FORALLUPDIR(nu) {
      QDP_R_eq_zero(axi[nu], QDP_all);
      QDP_R_eq_zero(vec[nu], QDP_all);
      tempaxi = QDP_expose_R(axi[nu]);
      tempvec = QDP_expose_R(vec[nu]);
      FORALLSITES(i, s) {
        pnt[0] = s->x;    pnt[1] = s->y;    pnt[2] = s->z;    pnt[3] = s->t;
        this_node = mynode();
        if (QDP_node_number(pnt) != this_node) {
          printf("PROBLEM: node %d vs. %d...\n",
                 QDP_node_number(pnt), this_node);
        }
        idx = QDP_index(pnt);
        tempaxi[idx] = lattice[i].axial[mu][nu];
        tempvec[idx] = lattice[i].vacpol[mu][nu];
      }
      QDP_reset_R(axi[nu]);
      QDP_reset_R(vec[nu]);
    }

    sprintf(md, "Axial mu = %i", mu);
    QDP_string_set(meta, md);
    QDP_vwrite_R(axiFile, meta, axi, ndim);

    sprintf(md, "Vector mu = %i", mu);
    QDP_string_set(meta, md);
    QDP_vwrite_R(vecFile, meta, vec, ndim);
  }
  QDP_close_write(axiFile);
  QDP_close_write(vecFile);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Includes and definitions
#include <ctype.h>
#include <math.h>
#include <qcdlib.h>
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print usage information if necessary
void print_usage(char *argv0) {
  char *fmt = " %-9s %s\n";

  printf("%s [options]\n", argv0);
  printf("options:\n");
  printf(fmt, "o pat", "input/output file pattern (required)");
  printf(fmt, "x #", "source point (four components)");
  printf("\n");
  QDP_abort(1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up source point
void getPoint(int j, int ndim, int *pnt, int argc, char *argv[]) {
  int i;

  if (j == 0) {                 // No point given in input, use default
    for (i = 0; i < ndim; ++i)
      pnt[i] = 0;
  }
  else {
    if (isdigit(argv[j][0]))
      pnt[0] = atoi(argv[j]);
    else
      pnt[0] = atoi(argv[j] + 1);
    for (i = 1; i < ndim; ++i) {
      if ((++j < argc) && ((argv[j][0] == '-') || isdigit(argv[j][0])))
        pnt[i] = atoi(argv[j]);
      else
        pnt[i] = 0;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Load correlators and print out zero-momentum-projected functions
int main(int argc, char *argv[]) {
  int i, j, mu, ndim, *latsize, Nt;
  int pnt[4] = {0, 0, 0, 0}, src[4] = {0, 0, 0, 0};
  char prec;
  QLA_Real *data, *axi[4], *vec[4];
  QDP_Real *corr[4];
  char *outpat = NULL, fn[80];
  QDP_Reader *qr;
  QDP_String *md;
  FILE *out = NULL;

  if (argc < 2)
    print_usage(argv[0]);

  QDP_initialize(&argc, &argv);
  QDP_profcontrol(0);           // No profiling

  j = 0;    // Check for source point

  for(i = 1; i < argc; i++) {
    char c, *argp = &argv[i][1];
    c = argv[i][0];
    if(argv[i][1] == '\0') {
      if(i + 1 < argc) {
        argp = argv[i + 1];
        i++;
      }
      else argp = NULL;
    }
    switch (c) {
    case 'o': outpat = argp; break;
    case 'x' :
      j = i;
      while (i + 1 < argc
            && ((argv[i + 1][0] == '-') || isdigit(argv[i + 1][0])))
        ++i;
      break;
    default: print_usage(argv[0]);
    }
  }

  if (outpat == NULL) {
    printf("Error: output pattern not specified!\n");
    print_usage(argv[0]);
  }

  sprintf(fn, "KSaxi%s", outpat);
  QCD_latticePrecGetFromFile(&latsize, &ndim, &prec, fn);
  QCD_latticeSet(latsize, ndim);
  QDP_create_layout();
  Nt = latsize[ndim - 1];

  // Set up source point -- default zero
  getPoint(j, ndim, src, argc, argv);

  QCD_latticePrint();
  printf("source point =");
  for (i = 0; i < ndim; i++)
    printf(" %d", src[i]);
  printf("\noutpat = %s\n", outpat);

  // Arrays
  for (mu = 0; mu < ndim; mu++) {
    corr[mu] = QDP_create_R();
    axi[mu] = malloc(Nt * sizeof(QLA_Real));
    vec[mu] = malloc(Nt * sizeof(QLA_Real));
    for (j = 0; j < Nt; j++) {
      axi[mu][j] = 0.0;
      vec[mu][j] = 0.0;
    }
  }

  // Load diagonal axial correlators and project to zero momentum
  printf("Loading correlator files %s and KSvec%s\n", fn, outpat);
  md = QDP_string_create();
  qr = QDP_open_read(md, fn);   // Set to KSaxi above
  for (mu = 0; mu < 4; mu++) {
    QDP_vread_R(qr, md, corr, 4);
    data = QDP_expose_R(corr[mu]);
    for (i = 0; i < QDP_sites_on_node; i++) {
      QDP_get_coords(pnt, QDP_this_node, i);

      // Account for non-zero source point
      j = (pnt[3] + latsize[3] - src[3]) % latsize[3];
      axi[mu][j] += data[i];
    }
    QDP_reset_R(corr[mu]);
  }
  QDP_close_read(qr);

  // Load vector correlators and project to zero momentum
  // Only consider first three; checked that fourth is useless
  sprintf(fn, "KSvec%s", outpat);
  qr = QDP_open_read(md, fn);
  for (mu = 0; mu < 3; mu++) {
    QDP_vread_R(qr, md, corr, 4);
    data = QDP_expose_R(corr[mu]);
    for (i = 0; i < QDP_sites_on_node; i++) {
      QDP_get_coords(pnt, QDP_this_node, i);

      // Account for non-zero source point
      j = (pnt[3] + latsize[3] - src[3]) % latsize[3];
      vec[mu][j] += data[i];
    }
    QDP_reset_R(corr[mu]);
  }
  QDP_close_read(qr);

  // Print out zero-momentum-projected correlators
  sprintf(fn, "KSdecay%s", outpat);
  printf("Writing output file %s\n", fn);
  out = fopen(fn, "w");
  for (j = 0; j < Nt; j++) {
    fprintf(out, "%d\t", j);
    for (mu = 0; mu < 4; mu++)
      fprintf(out, "%.4g\t", axi[mu][j]);
    for (mu = 0; mu < 2; mu++)
      fprintf(out, "%.4g\t", vec[mu][j]);
    fprintf(out, "%.4g\n", vec[2][j]);
  }
  fclose(out);

  // Clean up
  for (mu = 0; mu < ndim; mu++) {
    QDP_destroy_R(corr[mu]);
    free(axi[mu]);
    free(vec[mu]);
  }
  QDP_finalize();
  return 0;
}
// -----------------------------------------------------------------

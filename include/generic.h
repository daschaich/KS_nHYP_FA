// -----------------------------------------------------------------
// Macros and declarations for miscellaneous generic routines
#ifndef _GENERIC_H
#define _GENERIC_H

// Other generic directory declarations are elsewhere:
//   See comdefs.h for communications
//   See io_lat.h for I/O
#include <stdio.h>
#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
/* ax_gauge.c */
void ax_gauge();

/* bsd_sum.c */
int32type bsd_sum (char *data,int32type total_bytes);

// check_unitarity.c
Real check_unitarity();

// reunitarize2.c
void reunitarize();
int reunit_su3(su3_matrix *c);

/* d_linktrsum */
void d_linktrsum(double_complex *linktrsum);

/* d_plaq?.c */
void d_plaquette(double *ss_plaq, double *st_plaq);

// field_strength.c
// link_src is offset for su3_matrix link[4] in site struct
// field_dest is offset for su3_matrix fieldstrength[6] in site struct
void make_field_strength(field_offset link_src, field_offset field_dest);

/* gaugefix.c and gaugefix2.c */
void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
        Real gauge_fix_tol, field_offset diffmat, field_offset sumvec,
        int nvector, field_offset vector_offset[], int vector_parity[],
        int nantiherm, field_offset antiherm_offset[], int antiherm_parity[]);

/* gauge_force_imp.c and gauge_force_symzk1_qop.c */
void imp_gauge_force(Real eps, field_offset mom_off );

/* gauge_stuff.c */
double imp_gauge_action();
void g_measure();
void make_loop_table();
void dsdu_qhb_subl(int dir, int subl);
int get_max_length();
int get_nloop();
int get_nreps();
int *get_loop_length();
int *get_loop_num();
int ***get_loop_table();
Real **get_loop_coeff();

/* glueball_op.c */
void make_glueball_ops();
void measure_glueball_ops();

/* hvy_pot.c */
void hvy_pot(field_offset links );

/* io_detect.c */
int get_file_type(char *filename);
int io_detect(char *filename, file_table ft[], int ntypes);
int io_detect_fm(char *filename);

/* io_helpers.c */
gauge_file *save_lattice(int flag, char *filename, char *stringLFN );
gauge_file *reload_lattice(int flag, char *filename);
int ask_corr_file(FILE *fp, int prompt, int *flag, char* filename);
int ask_starting_lattice(FILE *fp, int prompt, int *flag, char *filename );
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename );
int ask_gauge_fix(FILE *fp, int prompt, int *flag);
int ask_ildg_LFN(FILE *fp, int prompt, int flag, char *stringLFN);
void coldlat();
void funnylat();
int get_check_tag(FILE *fp, char *tag, char *myname);
int get_f(FILE *fp, int prompt, char *variable_name_string, Real *value );
int get_i(FILE *fp, int prompt, char *variable_name_string, int *value );
char *get_next_tag(FILE *fp, char *tag, char *myname);
int get_vi(FILE *fp, int prompt, char *variable_name_string,
      int *value, int nvalues );
int get_vf(FILE *fp, int prompt, char *variable_name_string,
      Real *value, int nvalues );
int get_s(FILE *fp, int prompt, char *variable_name_string, char *value );
int get_sn(FILE *fp, int prompt, char *variable_name_string, char *value );
int get_vs(FILE *fp, int prompt, char *tag, char *value[], int nvalues );
int get_prompt(FILE *fp, int *value );

/* io_source_cmplx_fm.c */
void r_source_cmplx_fm_to_site(char *filename, field_offset dest_site,
             int t0, int source_type);
void r_source_cmplx_fm_to_field(char *filename, complex *dest_field,
        int t0, int source_type);

/* layout_*.c */
int io_node(const int node);
void setup_layout();
int node_number(int x, int y, int z, int t);
int node_index(int x, int y, int z, int t);
size_t num_sites(int node);
const int *get_logical_dimensions();
const int *get_logical_coordinate();
void get_coords(int coords[], int node, int index);

/* make_lattice.c */
void make_lattice();
void free_lattice();

/* nersc_cksum.c */
u_int32type nersc_cksum();

/* make_global_fields.c */
void make_global_fields();

/* path_product.c */
void path_product(const int *dir, const int length, su3_matrix *tempmat1);
void path_product_fields(su3_matrix *src[4], const int *dir,
                         const int length, su3_matrix *tempmat1);
void path_prod_subl(const int *dir, const int length, const int subl,
                    su3_matrix *tempmat1);

/* plaquette4.c */
void plaquette(Real *ss_plaq, Real *st_plaq);

/* ploop.c */
complex ploop(int dir);

/* project_su3_hit.c */
void project_su3(
   su3_matrix *w,         /* input initial guess. output resulting
                             SU(3) matrix */
   su3_matrix *q,         /* starting 3 x 3 complex matrix */
   int Nhit,              /* number of SU(2) hits. 0 for no projection */
   Real tol              /* tolerance for SU(3) projection.
           If nonzero, treat Nhit as a maximum
           number of hits.  If zero, treat Nhit
           as a prescribed number of hits. */
   );

/* rand_gauge.c */
void rand_gauge(field_offset G);

/* ranmom.c */
void ranmom();

/* remap standard I/O */
int remap_stdio_from_args(int argc, char *argv[]);

/* ranstuff.c */
void initialize_prn(double_prn *prn_pt, int seed, int index);
Real myrand(double_prn *prn_pt);

/* restrict_fourier.c */
void setup_restrict_fourier(int *key, int *slice);
void restrict_fourier_site(
     field_offset src,   /* src is field to be transformed */
     int size,     /* Size of field in bytes.  The field must
          consist of size/sizeof(complex) consecutive
          complex numbers.  For example, an su3_vector
          is 3 complex numbers. */
     int isign);   /* 1 for x -> k, -1 for k -> x */

void restrict_fourier_field(
     complex *src,       /* src is field to be transformed */
     int size,     /* Size of field in bytes.  The field must
          consist of size/sizeof(complex) consecutive
          complex numbers.  For example, an su3_vector
          is 3 complex numbers. */
     int isign);   /* 1 for x -> k, -1 for k -> x */
void cleanup_restrict_fourier();

/* show_generic_opts.c */
void show_generic_opts();

/* Do Morninstar-Peardon stout smearing to construct unitary W from
   smeared link V and unsmeared link U */
void stout_smear(su3_matrix *W, su3_matrix *V, su3_matrix *U);

/* For quark source and sink routines - both Wilson and KS */
/* The Weyl representation types are included for w_source_h */
enum source_type {
  UNKNOWN = 0,
  COMPLEX_FIELD_FILE,
  COMPLEX_FIELD_FM_FILE,
  COMPLEX_FIELD_STORE,
  CORNER_WALL,
  CUTOFF_GAUSSIAN,
  CUTOFF_GAUSSIAN_WEYL,
  COVARIANT_GAUSSIAN,
  DERIV1,
  DERIV2_D,
  DERIV2_B,
  DERIV3_A,
  DIRAC_FIELD_FILE,
  DIRAC_FIELD_FM_FILE,
  DIRAC_FIELD_STORE,
  EVEN_WALL,
  EVENANDODD_WALL,
  FAT_COVARIANT_GAUSSIAN,
  FAT_COVARIANT_GAUSSIAN_DERIV1,
  FAT_COVARIANT_GAUSSIAN_DERIV2_D,
  FAT_COVARIANT_GAUSSIAN_DERIV2_B,
  GAUSSIAN,
  POINT,
  POINT_WEYL,
  RANDOM_VECTOR_WALL,
  ROTATE_3D,
  WAVEFUNCTION_FILE,
  VECTOR_FIELD_FILE,
  VECTOR_FIELD_STORE
};
#endif
// -----------------------------------------------------------------

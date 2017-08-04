#ifndef _GENERIC_KS_H
// -----------------------------------------------------------------
// Macros and declarations for generic KS routines
#define _GENERIC_KS_H

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/generic_quark_types.h"
#include "../include/comdefs.h"
#ifdef HAVE_QOP
#include <qop.h>
#endif
#ifdef HAVE_QIO
#include <qio.h>
#endif

#define ALL_T_SLICES -1
#define MAXDESCRP 128
#define MAXSRCLABEL 8
// -----------------------------------------------------------------



// -----------------------------------------------------------------
/* Structure defining staggered quark source */
/* There must be a color member */
/* Add other members to suit the generic_ks code
   that builds the source.  Ignore the members you don't need. */
typedef struct {
  int type;           /* source type for most source builders */
  char descrp[MAXDESCRP]; /* alpha description for most */
  char label[MAXSRCLABEL]; /* Abbreviation of description */
  int color;          /* source color */
  int x0,y0,z0,t0;    /* source coordinates for most */
  char source_file[MAXFILENAME]; /* file name for some sources */
  int flag;           /* mode of reading or writing for some sources */
  int file_initialized;  /* has source file been initialized? */
  int mom[3];         /* momentum insertion for some actions */
#ifdef HAVE_QIO
  QIO_Reader *infile;
  QIO_Writer *outfile;
#endif
  complex *c_src;     /* Pointer for complex source field storage */
  vector *cv_src; /* Pointer for SU3 vector source field storage */
  int ksource;        /* Counter for a list of sources */
} ks_quark_source;

/* Structure specifying each rotation and reflection of each kind of
  path.  */
#define MAX_PATH_LENGTH 16
typedef struct {
  int dir[MAX_PATH_LENGTH]; /* directions in path */
  int length;   /* length of path */
  Real coeff;         /* coefficient, including minus sign if backwards */
  Real forwback;  /* +1 if in forward Dslash, -1 if in backward */
} Q_path;

/* Structure defining the fermion action using paths or optimized
   coefficients */

#ifdef HISQ
#define MAX_NAIK 4 // max number of quarks that require Naik epsilon correction
                   // without c quark this is normally 1, however this
                   // constant should be set to at least 2
typedef struct {
  Real *act_path_coeff;    /* For optimized Asqtad action */
  int num_q_paths;         /* For all actions */
  Q_path *q_paths;         /* For all actions */
//  Real naik_mass;             /* The mass last used in the Naik term */
} ks_component_paths;

typedef struct {
  ks_component_paths p1, p2, p3;
  int umethod;
  int ugroup;
  int constructed;       /* Boolean */
} ks_action_paths;

typedef struct {
  // Flags: 1 if the corresponding link field is valid
  // (valid means it corresponds to the current links in the site structure)
  int valid_U_links, valid_V_links, valid_W_links, valid_Y_links,
    valid_X_links;
  int valid_all_links; // should be 1 if ALL links are valid
  Real valid_Xfat_mass, valid_Xlong_mass;
  // phases...in = 1 if KS and antiperiodic BC signs are absorbed in the links
  int phases_in_U, phases_in_V, phases_in_W, phases_in_Y,
    phases_in_Xfat, phases_in_Xlong;
  matrix *U_link[4]; // original gauge matrices, stored as four fields
  matrix *V_link[4]; // first iteration of fattening
  matrix *Y_unitlink[4]; // unitary projection of V_link, U(3)
  matrix *W_unitlink[4]; // special unitary projection of Y_link, SU(3)
  matrix *XX_fatlink[MAX_NAIK][4];
  matrix *XX_longlink[MAX_NAIK][4];
  matrix *X_fatlink[4]; // these arrays are only pointers,
  matrix *X_longlink[4];// they are not malloc'ed
  int last_used_X_set; // set of X links to which X_fat/long points
  int current_X_set;   // set of X links to which X_fat/long should switch
} hisq_links_t;

#else  /* Non-HISQ actions */
typedef struct {
  Real *act_path_coeff;    /* For optimized Asqtad action */
  int num_q_paths;         /* For all actions */
  Q_path *q_paths;         /* For all actions */
  int constructed;         /* Boolean */
} ks_action_paths;
#endif

/* Structure defining the precomputed links for the FN actions */

typedef struct {
  int valid;
  matrix *fat;
  matrix *lng;
  matrix *fatback;
  matrix *lngback;
  ks_action_paths *ap;  /* For EO actions */
#ifdef HISQ
  Real mass;    /* The mass last used in the coefficients */
  hisq_links_t hl;
#endif
#ifdef HAVE_QOP
  int valid_qop_F;
  int valid_qop_D;
  QOP_F3_FermionLinksAsqtad *qop_F_l;
  QOP_D3_FermionLinksAsqtad *qop_D_l;
#endif
} ferm_links_t;

int congrad(int niter, int nrestart, Real rsqmin, int parity, Real *rsq);
void copy_latvec(field_offset src, field_offset dest, int parity);
void dslash_site(field_offset src, field_offset dest, int parity);
void dslash_site_special(field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start);
void dslash_field(vector *src, vector *dest, int parity);
void dslash_field_special(vector *src, vector *dest,
    int parity, msg_tag **tag, int start);

void checkmul();
void phaseset();
void rephase(int flag);

void prefetch_vector(vector *);
void prefetch_matrix(matrix *);

// Replaced in KS_nHYP_FA
//int ks_congrad(field_offset src, field_offset dest, Real mass,
//    int niter, int nrestart, Real rsqmin, int prec,
//    int parity, Real *rsq, ferm_links_t *fn);

int ks_congrad_field(vector *src, vector *dest,
          quark_invert_control *qic, Real mass,
          ferm_links_t *fn);

// Replaced in KS_nHYP_FA
//int ks_congrad_site(field_offset src, field_offset dest,
//         quark_invert_control *qic, Real mass,
//         ferm_links_t *fn);

int ks_congrad_two_src(/* Return value is number of iterations taken */
    field_offset src1,    /* source vector (type vector) */
    field_offset src2,
    field_offset dest1, /* solution vectors */
    field_offset dest2,
    Real mass1,
    Real mass2,
    int niter,    /* maximal number of CG interations */
    int nrestart,       /* maximum number of restarts */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* internal precision for the inversion */
    int parity,   /* parity to be worked on */
    Real  *final_rsq_ptr,   /* final residue squared */
    ferm_links_t *fn       /* Storage for fermion links */
   );

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_temps();

void dslash_fn_site(field_offset src, field_offset dest, int parity,
         ferm_links_t *fn);
void dslash_fn_site_special(field_offset src, field_offset dest,
           int parity, msg_tag **tag, int start,
           ferm_links_t *fn);
void ddslash_fn_du0_site(field_offset src, field_offset dest, int parity,
        ferm_links_t *fn, ferm_links_t *fn_dmdu0);

void dslash_fn_field(vector *src, vector *dest, int parity,
          ferm_links_t *fn);
void dslash_fn_field_special(vector *src, vector *dest,
           int parity, msg_tag **tag, int start,
           ferm_links_t *fn);
void ddslash_fn_du0_field(vector *src, vector *dest, int parity,
         ferm_links_t *fn, ferm_links_t *fn_dmdu0);

void dslash_eo_site(field_offset src, field_offset dest, int parity,
         ferm_links_t *fn);

/* The following three do not exist yet (3/05 -CD) */
void dslash_eo_site_special(field_offset src, field_offset dest,
           int parity, msg_tag **tag, int start,
           ferm_links_t *fn);
void dslash_eo_field(vector *src, vector *dest, int parity,
          ferm_links_t *fn);
void dslash_eo_field_special(vector *src, vector *dest,
            int parity, msg_tag **tag, int start,
            ferm_links_t *fn);

int congrad_ks(           /* Return value is number of iterations taken */
     field_offset src,       /* type vector* (preloaded source) */
     field_offset dest,      /* type vector*  (answer and initial guess) */
     quark_invert_control *qic, /* inverter control */
     void *dmp               /* parameters defining the Dirac matrix */
    );

int ks_invert(/* Return value is number of iterations taken */
    field_offset src,   /* type vector or multi_vector
         (preloaded source) */
    field_offset dest,  /* type vector or multi_vector
         (answer and initial guess) */
    int (*invert_func)(field_offset src, field_offset dest,
           quark_invert_control *qic,
           Real mass, ferm_links_t *fn),
    quark_invert_control *qic, /* inverter control */
    Real mass,
    ferm_links_t *fn
   );

int ks_invert_ksqs(/* Return value is number of iterations taken */
    ks_quark_source *ksqs, /* source parameters */
    int (*source_func_field)(vector *src,
            ks_quark_source *ksqs),  /* source function */
    vector *dest,  /* answer and initial guess */
    int (*invert_func)(vector *src, vector *dest,
           quark_invert_control *qic, Real mass,
           ferm_links_t *fn),
    quark_invert_control *qic, /* inverter control */
    Real mass,
    ferm_links_t *fn
       );

/* in ks_multicg.c */
enum ks_multicg_opt_t {OFFSET, HYBRID, FAKE, REVERSE, REVHYB};
const char *ks_multicg_opt_chr();

int ks_multicg(        /* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *offsets,  /* the offsets */
    int num_offsets,  /* number of offsets */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr, /* final residue squared */
    ferm_links_t *fn       /* Storage for fat and Naik links */
   );

int ks_multicg_p(      /* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *offsets,  /* the offsets */
    int num_offsets,  /* number of offsets */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn       /* Storage for fat and Naik links */
   );

int ks_multicg_offset( /* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *offsets,  /* the offsets */
    int num_offsets,  /* number of offsets */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn      /* Storage for fat and Naik links */
   );

int ks_multicg_mass( /* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *masses, /* the masses */
    int num_masses, /* number of masses */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn       /* Storage for fat and Naik links */
   );

int ks_multicg_hybrid( /* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *offsets,  /* the offsets */
    int num_offsets,  /* number of offsets */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn       /* Storage for fat and Naik links */
   );

int ks_multicg_reverse(/* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *masses, /* the masses */
    int num_masses, /* number of masses */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn      /* Storage for fat and Naik links */
   );

int ks_multicg_fake( /* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *offsets,  /* the offsets */
    int num_offsets,  /* number of offsets */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired internal precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn      /* Storage for fat and Naik links */
   );

int ks_multicg_revhyb( /* Return value is number of iterations taken */
    field_offset src, /* source vector (type vector) */
    vector **psim,  /* solution vectors */
    Real *offsets,  /* the offsets */
    int num_offsets,  /* number of offsets */
    int niter,    /* maximal number of CG interations */
    Real rsqmin,  /* desired residue squared */
    int prec,           /* desired intermediate precision */
    int parity,   /* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    ferm_links_t *fn      /* Storage for fat and Naik links */
   );

/* d_congrad_opt.c */

void clear_latvec(field_offset v,int parity);

void scalar_mult_add_latvec(field_offset src1, field_offset src2,
          Real scalar, field_offset dest, int parity);
void scalar2_mult_add_vector(vector *a, Real s1, vector *b,
         Real s2, vector *c);
void scalar2_mult_add_latvec(field_offset src1,Real scalar1,
           field_offset src2,Real scalar2,
           field_offset dest,int parity);

// f_meas.c
// Removed ferm_links_t, now returns total iterations
int f_meas_imp(field_offset chi_off, field_offset psi_off, Real mass);

// fpi_2.c -- stripped ferm_links_t and multiple masses
int fpi_2(Real *masses);

/* fermion_force_asqtad*.c */
void eo_fermion_force_oneterm(Real eps, Real weight, field_offset x_off,
             int prec, ferm_links_t *fn,
             ks_action_paths *ap);
void eo_fermion_force_twoterms(Real eps, Real weight1, Real weight2,
        field_offset x1_off, field_offset x2_off,
        int prec, ferm_links_t *fn,
        ks_action_paths *ap);
void fermion_force_asqtad_block(Real eps, Real *residues,
         vector **xxx, int nterms, int veclength,
         int prec, ferm_links_t *fn,
         ks_action_paths *ap);
void fermion_force_asqtad_multi(Real eps, Real *residues,
         vector **xxx, int nterms, int prec,
         ferm_links_t *fn, ks_action_paths *ap);

/* fermion_force_fn_multi.c */

enum ks_multiff_opt_t {ASVEC, FNMAT, FNMATREV};

const char *ks_multiff_opt_chr();

int eo_fermion_force_set_opt(char opt_string[]);
void eo_fermion_force_multi(Real eps, Real *residues, vector **xxx,
           int nterms, int prec, ferm_links_t *fn,
           ks_action_paths *ap);
void fermion_force_asqtad_block(Real eps, Real *residues,
         vector **xxx, int nterms, int veclength,
         int prec, ferm_links_t *fn,
         ks_action_paths *ap);
void fermion_force_asqtad_multi(Real eps, Real *residues, vector **xxx,
         int nterms, int prec,
         ferm_links_t *fn, ks_action_paths *ap);
void fermion_force_fn_multi(Real eps, Real *residues, vector **multi_x,
           int nterms, int prec, ferm_links_t *fn,
             ks_action_paths *ap);
void fermion_force_fn_multi_reverse(Real eps, Real *residues,
             vector **multi_x, int nterms,
             ferm_links_t *fn, ks_action_paths *ap);
void fermion_force_fn_multi_june05(Real eps, Real *residues,
            vector **multi_x, int nterms,
            ferm_links_t *fn, ks_action_paths *ap);

/* fermion_links_fn.c and fermion_links_hisq.c */
void init_ferm_links(ferm_links_t *fn);
void load_ferm_links(ferm_links_t *fn, ks_action_paths *ap);
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap);
void invalidate_all_ferm_links(ferm_links_t *fn);
void invalidate_fn_links(ferm_links_t *fn);

/* fermion_links_helpers.c and fermion_links_hisq_helpers.c */
void load_longlinks(ferm_links_t *fn, ks_action_paths *ap);
void load_fatlinks(ferm_links_t *fn, ks_action_paths *ap);
void load_longbacklinks(ferm_links_t *fn);
void load_fatbacklinks(ferm_links_t *fn);
void free_fn_links(ferm_links_t *fn);
void free_fn_links_dmdu0(ferm_links_t *fn);
#ifdef HISQ
void load_fatlinks_hisq(matrix **Src, ks_component_paths *app,
       matrix **Dest);
void load_longlinks_hisq(matrix **Src, ks_component_paths *app,
        matrix **Dest);
#endif
void custom_rephase(matrix **internal_links, int flag, int *status_now);

/* ff_opt.c */
#ifndef VECLENGTH
#define VECLENGTH 4
#endif

typedef struct { vector v[VECLENGTH]; } veclist;
void mult_adj_su3_fieldlink_latveclist(matrix *link,
               veclist **src_pt, veclist *dest, int listlength);
void mult_su3_sitelink_latveclist(int dir,
   veclist **src_pt, veclist *dest, int listlength);
void scalar_mult_add_latveclist_proj(anti_hermitmat *mom,
           veclist *back, veclist *forw, Real *coeff, int listlength);
void scalar_mult_add_latveclist(veclist *dest,
            veclist *src, Real *s, int listlength);

/* flavor_ops.c */
void sym_shift(int dir, field_offset src,field_offset dest);
void zeta_shift(int n, int *d, field_offset src, field_offset dest);
void eta_shift(int n, int *d, field_offset src, field_offset dest);

void mult_flavor_vector(int mu, field_offset src, field_offset dest);
void mult_flavor_tensor(int mu, int nu, field_offset src, field_offset dest);
void mult_flavor_pseudovector(int mu, field_offset src, field_offset dest);
void mult_flavor_pseudoscalar(field_offset src, field_offset dest);

void mult_spin_vector(int mu, field_offset src, field_offset dest);
void mult_spin_tensor(int mu, int nu, field_offset src, field_offset dest);
void mult_spin_pseudovector(int mu, field_offset src, field_offset dest);
void mult_spin_pseudoscalar(field_offset src, field_offset dest);

/* grsource.c */
void grsource(int parity);

/* fermion_links_fn.c */
void init_ferm_links(ferm_links_t *fn);
void load_ferm_links(ferm_links_t *fn, ks_action_paths *ap);
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap);
void invalidate_all_ferm_links(ferm_links_t *fn);

/* ks_source.c */
void alloc_ksqs_c_src(ks_quark_source *ksqs);
void alloc_ksqs_cv_src(ks_quark_source *ksqs);
void clear_ksqs(ks_quark_source *ksqs);
int ask_ks_quark_source(FILE *fp, int prompt, int *type, char *descrp);
int ask_output_ks_quark_source_file(FILE *fp, int prompt,
             int *flag, int *source_type,
             int *t0, char *descrp, char *filename);
int choose_usqcd_ks_file_type(int source_type);
int get_ks_quark_source(FILE *fp, int prompt, ks_quark_source *ksqs);
int get_ks_quark_sink(FILE *fp, int prompt, ks_quark_source *ksqs);
void init_ksqs(ks_quark_source *ksqs);
int ks_source_site(field_offset src, ks_quark_source *ksqs);
int ks_source_field(vector *src, ks_quark_source *ksqs);
void ks_sink_site(field_offset snk, ks_quark_source *ksqs);
void ks_sink_field(complex *snk, ks_quark_source *ksqs);
void ks_sink_scalar(field_offset snk, ks_quark_source *ksqs);
int ks_source_write(vector *src, ks_quark_source *ksqs);
void r_close_ks_source(ks_quark_source *ksqs);
void r_open_ks_source(ks_quark_source *ksqs);
void w_close_ks_source(ks_quark_source *ksqs);
void w_open_ks_source(ks_quark_source *ksqs, char *fileinfo);

/* ks_utilities.c */
void clear_v_field(vector *v);
void clear_ksp_field(ks_prop_field ksp);
vector *create_v_field();
ks_prop_field create_ksp_field();
void copy_ksp_field(ks_prop_field kspcopy, ks_prop_field ksp);
vector *create_v_field();
ks_prop_field create_ksp_field_copy(ks_prop_field k);
void destroy_v_field(vector *v);
void destroy_ksp_field(ks_prop_field ksp);

// mat_invert.c
int mat_invert_uml(field_offset src, field_offset dest,
                   field_offset temp, Real mass);

/* mu.c and mu_fast.c */
void M_derivatives(field_offset phi_off, field_offset xxx_off,
       field_offset xxx1_off, Real mass,
       ferm_links_t *fn, ferm_links_t *fn_dmdu0);
void Deriv_O6(field_offset phi_off, field_offset xxx_off,
        field_offset xxx1_off, Real mass,
        ferm_links_t *fn, ferm_links_t *fn_dmdu0);

/* multimass_inverter.c */
#define MAX_MMINV_NMASSES 32
#define MAX_MMINV_SOURCES 16
typedef struct {
  Real masses[MAX_MMINV_NMASSES];
  int nmasses;
  int n_sources;
  int r0[MAX_MMINV_SOURCES][4];
  Real tol;
  Real rsqprop;
} params_mminv;

int multimass_inverter(params_mminv *mminv, ferm_links_t *fn);

// nl_spectrum.c -- stripped ferm_links_t
int nl_spectrum(Real vmass, field_offset tempvec1, field_offset tempvec2,
     field_offset tempmat1, field_offset tempmat2);

/* path_transport.c */
void link_gather_connection_hisq(matrix *src,
          matrix *dest, matrix *work, int dir);
void path_transport_site(field_offset src, field_offset dest, int parity,
        int *dir, int length);
void path_transport_field(vector * src, vector * dest, int parity,
         int *dir, int length);
void path_transport_hwv_site(field_offset src, field_offset dest, int parity,
            int *dir, int length);
void path_transport_connection(matrix * src, matrix * dest, int parity, int *dir, int length);
void link_transport_connection(matrix * src, matrix * dest, matrix * work, int dir);
void path_transport_connection_hisq(matrix * src, matrix **links, matrix * dest,
    int parity, int *dir, int length);
void link_transport_connection_hisq(matrix * src, matrix **links, matrix * dest,
    matrix * work, int dir);

/* quark_stuff.c and quark_stuff_hisq.c */
void init_path_table(ks_action_paths *ap);
int make_path_table(ks_action_paths *ap, ks_action_paths *ap_dmdu0);

/* show_generic_ks_opts.c */
void show_generic_ks_opts();

/* spectrum.c */
int spectrum(ferm_links_t *fn);

// spectrum2.c -- stripped ferm_links_t
int spectrum2(Real vmass, field_offset temp1, field_offset temp2);

/* spectrum_fzw.c */
int spectrum_fzw(Real vmass, field_offset temp1, field_offset temp2,
      ferm_links_t *fn);

/* spectrum_hybrids.c */
int spectrum_hybrids(Real mass, field_offset temp, Real tol,
          ferm_links_t *fn);

/* spectrum_mom.c */
int spectrum_mom(Real qmass, Real amass, field_offset temp, Real tol,
      ferm_links_t *fn);

/* spectrum_multimom.c */
int spectrum_multimom(Real dyn_mass, Real low_mass, Real mass_inc,
           int nmasses, Real tol, ferm_links_t *fn);

/* spectrum_nd.c */
int spectrum_nd(Real mass1, Real mass2, Real tol,
     ferm_links_t *fn);

// spectrum_nlpi2.c -- stripped ferm_links_t
int spectrum_nlpi2(Real fmass, Real amass);

/* spectrum_singlets */
int spectrum_singlets(Real mass, Real tol, field_offset temp_offset,
           ferm_links_t *fn);
#endif
// -----------------------------------------------------------------

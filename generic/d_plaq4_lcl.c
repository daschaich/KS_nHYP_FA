/************************** d_plaq4_lcl.c **************************/
/* This version mallocs the temporary matrix */

/* Double precision version of "plaquette4.c" including optional
   Schroedinger functional - UMH - 1/27/00 */

/* All matrices declared to be NCOLxNCOL.  -bqs 12/06	*/

/* Measure the average plaquette of the space-space and
   space-time plaquettes */

/* Added plaq measurement on hypersurfaces; turn on with:          */
/* #define LOCAL_PLAQ */

/* Added measurement of min plaq per config for tuning nhyp as in
Hasenfratz & Knechtli, hep-lat/0103029; turn on with:          */
#define MIN_PLAQ

/* Added local plaq measurement for plotting plaq distribution
 (frep only); turn on with:  */
/* #define ALL_PLAQ */
/*** CAUTION: Do not run ALL_PLAQ with MPI!! ***/

/* The following should be defined self-consistently.
   The value of MY_DIR is printed out on the first time
   the routine is called.                                          */

#ifdef LOCAL_PLAQ
#define MY_X x
#define MY_DIR XUP
#define MY_N nx
static int print_dir=0;
#endif

/* plaq_prll[xx] gives the average of the plaquettes lying in a
   fixed MY_X hypersurface, as a function of the MY_X coordinate.
   Similarly the other plaquettes are given by plaq_prll[xx].      */

#include "generic_includes.h"

void d_plaquette_lcl(double *ss_plaq,double *st_plaq) {
/* mat is scratch space of size matrix */
matrix_f *mat;
register int i,dir1,dir2;
register site *s;
register matrix_f *m1,*m4;
matrix_f mtmp;
double ss_sum,st_sum,cur_plaq;
#ifdef MIN_PLAQ
double min_plaq=NCOL;
#endif
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;

#ifdef LOCAL_PLAQ
int xx;
double *plaq_perp,*plaq_prll;

    plaq_perp = (double *)malloc(MY_N*sizeof(double));
    plaq_prll = (double *)malloc(MY_N*sizeof(double));

    for(xx=0;xx<MY_N;xx++){
       plaq_perp[xx] = 0.0;
       plaq_prll[xx] = 0.0;
    }
#endif

    mat = (matrix_f *)malloc(sizeof(matrix_f)*sites_on_node);
    if(mat == NULL)
      {
	printf("plaquette: can't malloc mat\n");
	fflush(stdout); terminate(1);
      }

    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(matrix_f),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(matrix_f),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(s->link[dir1]);
		m4 = &(s->link[dir2]);
		mult_su3_an_f(m4,m1,&mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);
/* SF: supersede the gather when necessary                            */
            FORALLSITES(i,s){
                gen_pt[0][i]=CHOOSE_NBR(i,s,dir1,link_bndr_up[dir2],0);
            }

/* SF: we need the space-time plaquettes associated with t=0 as well  */
            if (dir1==TUP ) {
               FORALLSITES(i,s){
	          mult_su3_nn_f( &(mat[i]), (matrix_f *)(gen_pt[0][i]),
                        &mtmp);
                  cur_plaq = (double)
		      realtrace_su3_f((matrix_f *)(gen_pt[1][i]),&mtmp);
#ifdef MIN_PLAQ
                  if(cur_plaq<min_plaq) min_plaq=cur_plaq;
#endif
                  st_sum += cur_plaq;
#ifdef LOCAL_PLAQ
                  if(dir1==MY_DIR || dir2==MY_DIR){
                      plaq_perp[s->MY_X] += cur_plaq;
                  }
                  else{
                      plaq_prll[s->MY_X] += cur_plaq;
                  }
#endif
               }
            }
/* SF: space-space plaquettes only for t>0  */
            else {
               FORALLSITESDOMAIN(i,s){
		  mult_su3_nn_f( &(mat[i]), (matrix_f *)(gen_pt[0][i]),
			&mtmp);
                  cur_plaq = (double)
		        realtrace_su3_f((matrix_f *)(gen_pt[1][i]),&mtmp);
#ifdef MIN_PLAQ
                  if(cur_plaq<min_plaq) min_plaq=cur_plaq;
#endif
                  ss_sum += cur_plaq;
#ifdef LOCAL_PLAQ
                  if(dir1==MY_DIR || dir2==MY_DIR){
                      plaq_perp[s->MY_X] += cur_plaq;
                  }
                  else{
                      plaq_prll[s->MY_X] += cur_plaq;
                  }
#endif
               }
            }

/*	    FORALLSITES(i,s){
#ifdef SCHROED_FUN
*		if(dir1==TUP ){
*		    if(s->t==(nt-1)){
*			mult_su3_nn_f( &mat[i],
*			    &(s->boundary[dir2]), &mtmp);
*		    }
*		    else{
*			mult_su3_nn_f( &mat[i],
*			    (matrix_f *)(gen_pt[0][i]), &mtmp);
*		    }
*		    st_sum +=
*			realtrace_su3_f((matrix_f *)(gen_pt[1][i]), &mtmp);
*		}
*		else if(s->t > 0){
*		    mult_su3_nn_f( &mat[i], (matrix_f *)(gen_pt[0][i]),
*			&mtmp);
*		    ss_sum +=
*			realtrace_su3_f((matrix_f *)(gen_pt[1][i]), &mtmp);
*		}
#else
*		mult_su3_nn_f( &mat[i], (matrix_f *)(gen_pt[0][i]),
*		    &mtmp);
*
*		if(dir1==TUP )st_sum += (double)
*		    realtrace_su3_f((matrix_f *)(gen_pt[1][i]),&mtmp);
*		else          ss_sum += (double)
*		    realtrace_su3_f((matrix_f *)(gen_pt[1][i]),&mtmp);
#endif
*	    }
*/
	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	} /* dir2 */
    } /* dir1 */
    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );
/* there are three space-space and three time-space plaquettes */
#ifdef SF
    *ss_plaq = ss_sum /((double)(3*nx*ny*nz*(nt-1)));
#else
    *ss_plaq = ss_sum /((double)(3*nx*ny*nz*nt));
#endif
    *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));

    free(mat);

#ifdef LOCAL_PLAQ
    for(xx=0;xx<MY_N;xx++){
       g_doublesum( &(plaq_perp[xx]) );
       g_doublesum( &(plaq_prll[xx]) );
    }
    /* normalization */
    for(xx=0;xx<MY_N;xx++){
#ifndef SF
       plaq_perp[xx] *= ((double)(MY_N))/(3*nx*ny*nz*nt);
       plaq_prll[xx] *= ((double)(MY_N))/(3*nx*ny*nz*nt);
#else /* SF */
       if(MY_DIR==TUP){
          plaq_perp[xx] /= (double)(3*nx*ny*nz);
          plaq_prll[xx] /= (double)(3*nx*ny*nz);
       }
       else{
          plaq_perp[xx] *= ((double)(MY_N))/(nx*ny*nz*(3*nt-2));
          plaq_prll[xx] *= ((double)(MY_N))/(nx*ny*nz*(3*nt-1));
       }
#endif /* SF */
    }
    /* printout */

   if(this_node==0) {
      if(print_dir==0){
         printf("LOCAL_PLAQ [0=XUP,..., 3=TUP] dir=%d\n", MY_DIR);
         print_dir=1;
      }
      printf("THIN_PLAQ_PERP");
      for(xx=0;xx<MY_N;xx++){
         printf(" %e", (double)plaq_perp[xx]);
      }
      printf("\n");
      printf("THIN_PLAQ_PRLL");
      for(xx=0;xx<MY_N;xx++){
         printf(" %e", (double)plaq_prll[xx]);
      }
      printf("\n");
   }
#endif /* LOCAL_PLAQ */

#ifdef MIN_PLAQ
   min_plaq = -min_plaq;
   g_doublemax(&min_plaq);
   min_plaq = -min_plaq;
   node0_printf("MIN_PLAQ_FUND %e\n",min_plaq);
#endif

} /* d_plaquette4 */


/*** fermion irrep's plaquette  ***/
void d_plaquette_frep_lcl(double *ss_plaq_frep, double *st_plaq_frep) {
/* mat is scratch space of size matrix */
matrix *mat;
register int i,dir1,dir2;
register site *s;
register matrix *m1,*m4;
matrix mtmp;
double ss_sum,st_sum,cur_plaq;
#ifdef MIN_PLAQ
double min_plaq=DIMF;
#endif
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;

#ifdef LOCAL_PLAQ
int xx;
double *plaq_perp,*plaq_prll;

    plaq_perp = (double *)malloc(MY_N*sizeof(double));
    plaq_prll = (double *)malloc(MY_N*sizeof(double));

    for(xx=0;xx<MY_N;xx++){
       plaq_perp[xx] = 0.0;
       plaq_prll[xx] = 0.0;
    }
#endif

    mat = (matrix *)malloc(sizeof(matrix)*sites_on_node);
    if(mat == NULL)
      {
	printf("plaquette: can't malloc mat\n");
	fflush(stdout); terminate(1);
      }

    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(s->link[dir1]);
		m4 = &(s->link[dir2]);
		mult_su3_an(m4,m1,&mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);
/* SF: supersede the gather when necessary                            */
            FORALLSITES(i,s){
                gen_pt[0][i]=CHOOSE_NBR(i,s,dir1,link_bndr_up[dir2],0);
            }

/* SF: we need the space-time plaquettes associated with t=0 as well  */
            if (dir1==TUP ) {
               FORALLSITES(i,s){
	          mult_su3_nn( &(mat[i]), (matrix *)(gen_pt[0][i]),
                        &mtmp);
                  cur_plaq = (double)
		      realtrace_su3((matrix *)(gen_pt[1][i]),&mtmp);
#ifdef MIN_PLAQ
                  if(cur_plaq<min_plaq) min_plaq=cur_plaq;
#endif
#ifdef ALL_PLAQ
                  printf("ALL_PLAQ %d %d %d %d %d %d %e\n",
                         s->x,s->y,s->z,s->t,dir1,dir2,cur_plaq);
#endif
                  st_sum += cur_plaq;
#ifdef LOCAL_PLAQ
                  if(dir1==MY_DIR || dir2==MY_DIR){
                      plaq_perp[s->MY_X] += cur_plaq;
                  }
                  else{
                      plaq_prll[s->MY_X] += cur_plaq;
                  }
#endif
               }
            }
/* SF: space-space plaquettes only for t>0  */
            else {
               FORALLSITESDOMAIN(i,s){
		  mult_su3_nn( &(mat[i]), (matrix *)(gen_pt[0][i]),
			&mtmp);
                  cur_plaq = (double)
		        realtrace_su3((matrix *)(gen_pt[1][i]),&mtmp);
#ifdef MIN_PLAQ
                  if(cur_plaq<min_plaq) min_plaq=cur_plaq;
#endif
#ifdef ALL_PLAQ
                  printf("ALL_PLAQ %d %d %d %d %d %d %e\n",
                         s->x,s->y,s->z,s->t,dir1,dir2,cur_plaq);
#endif
                  ss_sum += cur_plaq;
#ifdef LOCAL_PLAQ
                  if(dir1==MY_DIR || dir2==MY_DIR){
                      plaq_perp[s->MY_X] += cur_plaq;
                  }
                  else{
                      plaq_prll[s->MY_X] += cur_plaq;
                  }
#endif
               }
            }


	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	} /* dir2 */
    } /* dir1 */
    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );
/* there are three space-space and three time-space plaquettes */
#ifdef SF
    *ss_plaq_frep = ss_sum /((double)(3*nx*ny*nz*(nt-1)));
#else
    *ss_plaq_frep = ss_sum /((double)(3*nx*ny*nz*nt));
#endif
    *st_plaq_frep = st_sum /((double)(3*nx*ny*nz*nt));

    free(mat);

#ifdef LOCAL_PLAQ
    for(xx=0;xx<MY_N;xx++){
       g_doublesum( &(plaq_perp[xx]) );
       g_doublesum( &(plaq_prll[xx]) );
    }
    /* normalization */
    for(xx=0;xx<MY_N;xx++){
#ifndef SF
       plaq_perp[xx] *= ((double)(MY_N))/(3*nx*ny*nz*nt);
       plaq_prll[xx] *= ((double)(MY_N))/(3*nx*ny*nz*nt);
#else /* SF */
       if(MY_DIR==TUP){
          plaq_perp[xx] /= (double)(3*nx*ny*nz);
          plaq_prll[xx] /= (double)(3*nx*ny*nz);
       }
       else{
          plaq_perp[xx] *= ((double)(MY_N))/(nx*ny*nz*(3*nt-2));
          plaq_prll[xx] *= ((double)(MY_N))/(nx*ny*nz*(3*nt-1));
       }
#endif /* SF */
    }
    /* printout */

   if(this_node==0) {
/*    printf("LOCAL_PLAQ dir=%d\n", MY_DIR); */
      printf("FAT_PLAQ_PERP");
      for(xx=0;xx<MY_N;xx++){
         printf(" %e", (double)plaq_perp[xx]);
      }
      printf("\n");
      printf("FAT_PLAQ_PRLL");
      for(xx=0;xx<MY_N;xx++){
         printf(" %e", (double)plaq_prll[xx]);
      }
      printf("\n");
   }
#endif /* LOCAL_PLAQ */

#ifdef MIN_PLAQ
   min_plaq = -min_plaq;
   g_doublemax(&min_plaq);
   min_plaq = -min_plaq;
   node0_printf("MIN_PLAQ_FERM %e\n",min_plaq);
#endif

} /* d_plaquette_frep */


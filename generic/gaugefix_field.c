/************************** gaugefix_field.c *******************************/
/* Fix Coulomb or Lorentz gauge by doing successive SU(2) gauge hits */
/* Uses double precision global sums */
/* This version does automatic reunitarization at preset intervals */
/* C. DeTar 10-22-90 */
/* T. DeGrand 1993 */
/* U.M. Heller 8-31-95 */
/* C. DeTar 10-11-97 converted to generic */
/* C. DeTar 12-26-97 added automatic reunitarization */
/* C. DeTar 11-24-98 remove superfluous references to p2 (was for ks phases) */
/* B. Svetitsky 12-18-07 adjust for general NCOL, DIMF */
/* YS Jan 2010 adapt to SF */
/* YS Aug 2011 "field" version for fundamental rep links only */

/* Prototype...

   void gaugefix_field(int gauge_dir,Real relax_boost,int max_gauge_iter,
	      Real gauge_fix_tol )

   -------------------------------------------------------------------
   EXAMPLE:  Fixing the link matrices to Coulomb gauge:

   gaugefix_field(TUP,(Real)1.5,500,(Real)1.0e-7);

   -------------------------------------------------------------------

   gauge_dir     specifies the direction of the "time"-like hyperplane
                 for the purposes of defining Coulomb or Lorentz gauge
      TUP    for evaluating propagators in the time-like direction
      ZUP    for screening lengths.
      8      for Lorentz gauge
   relax_boost	   Overrelaxation parameter
   max_gauge_iter  Maximum number of iterations
   gauge_fix_tol   Stop if change is less than this
   diffmat         Scratch space for an su3 matrix
   sumvec          Scratch space for an su3 vector
   NOTE: if diffmat or sumvec are negative, gaugefix mallocs its own
   scratch space. */

#include "generic_includes.h"
#define REUNIT_INTERVAL 20

/*    CDIF(a,b)         a -= b						      */
								/*  a -= b    */
#define CDIF(a,b) { (a).real -= (b).real; (a).imag -= (b).imag; }

/* Scratch space */

static matrix_f *diffmatp;               /* malloced diffmat pointer */
static vector_f *sumvecp;                /* malloced sumvec pointer */
void accum_gauge_hit(int gauge_dir,int parity);
void do_hit(int gauge_dir, int parity, int p, int q, Real relax_boost);
double get_gauge_fix_action(int gauge_dir,int parity);
void gaugefixstep(int gauge_dir,double *av_gauge_fix_action,Real relax_boost);
void gaugefixscratch();

/***************************************************************************/
void accum_gauge_hit(int gauge_dir,int parity)
{

/* Accumulates sums and differences of link matrices for determining optimum */
/* hit for gauge fixing */
/* Differences are kept in diffmat and the diagonal elements of the sums */
/* in sumvec  */

  register int j;
  register matrix_f *m1,*m2;
  register int dir,i;
  register site *s;

  /* Clear sumvec and diffmat */

  FORSOMEPARITY(i,s,parity)
    {
      clear_mat_f(&diffmatp[i]);
      clearvec_f(&sumvecp[i]);
    }

  /* Subtract upward link contributions */

  FORSOMEPARITYDOMAIN(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
        {
          /* Upward link matrix */
          m1 = &(s->link[dir]);
          sub_matrix_f( &diffmatp[i], m1, &diffmatp[i]);
          for(j=0;j<NCOL;j++)CSUM( sumvecp[i].c[j],m1->e[j][j]);
        }
    }

  /* Add downward link contributions */

  FORSOMEPARITYDOMAIN(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
        {
          /* Downward link matrix */
          m2 = (matrix_f *)gen_pt[dir][i];
          add_matrix_f( &diffmatp[i], m2, &diffmatp[i]);
          /* Add diagonal elements to sumvec  */
          for(j=0;j<NCOL;j++)CSUM( sumvecp[i].c[j], m2->e[j][j]);
        }
    }
} /* accum_gauge_hit */


void do_hit(int gauge_dir, int parity, int p, int q, Real relax_boost)
{
  /* Do optimum SU(2) gauge hit for p, q subspace */

  Real a0,a1,a2,a3,asq,a0sq,x,r,xdr;
  register int dir,i;
  register site *s;
  su2_matrix u;

  /* Accumulate sums for determining optimum gauge hit */

  accum_gauge_hit(gauge_dir,parity);

  FORSOMEPARITYDOMAIN(i,s,parity)
    {
      /* The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)*/
      /* The locally optimum unnormalized components a0, aj are determined */
      /* from the current link in direction dir and the link downlink */
      /* in the same direction on the neighbor in the direction opposite dir */
      /* The expression is */
      /* a0 = Sum dir Tr Re 1       * (downlink dir + link dir) */
      /* aj = Sum dir Tr Im sigma j * (downlink dir - link dir)  j = 1,2, 3 */
      /*   where 1, sigma j are unit and Pauli matrices on the p,q subspace */
 /*
      a0 =  s->sumvec.c[p].real + s->sumvec.c[q].real;
      a1 =  s->diffmat.e[q][p].imag + s->diffmat.e[p][q].imag;
      a2 = -s->diffmat.e[q][p].real + s->diffmat.e[p][q].real;
      a3 =  s->diffmat.e[p][p].imag - s->diffmat.e[q][q].imag;
*/
      a0 =     sumvecp[i].c[p].real +  sumvecp[i].c[q].real;

      a1 =     diffmatp[i].e[q][p].imag + diffmatp[i].e[p][q].imag;
      a2 =    -diffmatp[i].e[q][p].real + diffmatp[i].e[p][q].real;
      a3 =     diffmatp[i].e[p][p].imag - diffmatp[i].e[q][q].imag;

      /* Over-relaxation boost */

      /* This algorithm is designed to give little change for large |a| */
      /* and to scale up the gauge transformation by a factor of relax_boost*/
      /* for small |a| */

      asq = a1*a1 + a2*a2 + a3*a3;
      a0sq = a0*a0;
      x = (relax_boost*a0sq + asq)/(a0sq + asq);
      r = sqrt((double)(a0sq + x*x*asq));
      xdr = x/r;
      /* Normalize and boost */
      a0 = a0/r; a1 = a1*xdr; a2 = a2*xdr; a3 = a3*xdr;

      /* Elements of SU(2) matrix */

      u.e[0][0] = cmplx( a0, a3);
      u.e[0][1] = cmplx( a2, a1);
      u.e[1][0] = cmplx(-a2, a1);
      u.e[1][1] = cmplx( a0,-a3);


      /* Do SU(2) hit on all upward links */

      FORALLUPDIR(dir)
	left_su2_hit_n_f(&u,p,q,&(s->link[dir]));

      /* Do SU(2) hit on all downward links */

      FORALLUPDIR(dir)
	right_su2_hit_a_f(&u,p,q,(matrix_f *)gen_pt[dir][i]);

    }

  /* Exit with modified downward links left in communications buffer */
} /* do_hit */


double get_gauge_fix_action(int gauge_dir,int parity)
{
  /* Adds up the gauge fixing action for sites of given parity */
  /* Returns average over these sites */
  /* The average is normalized to a maximum of 1 when all */
  /* links are unit matrices */

  register int dir,i,ndir;
  register site *s;
  register matrix_f *m1, *m2;
  double gauge_fix_action;
  complex trace;

  gauge_fix_action = 0.0;

  FORSOMEPARITYDOMAIN(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
	{
	  m1 = &(s->link[dir]);
	  m2 = (matrix_f *)gen_pt[dir][i];

	  trace = trace_su3_f(m1);
	  gauge_fix_action += (double)trace.real;

	  trace = trace_su3_f(m2);
 	  gauge_fix_action += (double)trace.real;
	}
    }

  /* Count number of terms to average */
  ndir = 0; FORALLUPDIRBUT(gauge_dir,dir)ndir++;

  /* Sum over all sites of this parity */
  g_doublesum( &gauge_fix_action);

  /* Average is normalized to max of 1/2 on sites of one parity */
  return(gauge_fix_action /((double)(2*NCOL*ndir*nx*ny*nz*nt)));
} /* get_gauge_fix_action */


void gaugefixstep(int gauge_dir,double *av_gauge_fix_action,Real relax_boost)
{
  /* Carry out one iteration in the gauge-fixing process */

  int parity;
  msg_tag *mtag[8];
  Real gauge_fix_action;
  register int dir,i,c1,c2;
  register site *s;

  /* Alternate parity to prevent interactions during gauge transformation */
  *av_gauge_fix_action = 0.;
  g_sync();
  fflush(stdout);

  for(parity = ODD; parity <= EVEN; parity++)
    {
      /* Start gathers of downward links */

      FORALLUPDIR(dir)
	{
	  mtag[dir] = start_gather_site( F_OFFSET(link[dir]), sizeof(matrix_f),
			   OPP_DIR(dir), parity, gen_pt[dir] );
	}

      /* Wait for gathers */

      FORALLUPDIR(dir)
         {
	  wait_gather(mtag[dir]);
	 }

      /* Total gauge fixing action for sites of this parity: Before */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);

      /* Do optimum gauge hit on all subspaces */

      for (c1=0;c1<NCOL-1;c1++)for(c2=c1+1;c2<NCOL;c2++){
        do_hit(gauge_dir,parity,c1,c2, relax_boost);
      }

      /* Total gauge fixing action for sites of this parity: After */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);

      *av_gauge_fix_action += gauge_fix_action;

      /* Scatter downward link matrices by gathering to sites of */
      /* opposite parity */

      FORALLUPDIR(dir)
	{
	  /* Synchronize before scattering to be sure the new modified link */
	  /* matrices are all ready to be scattered and diffmat is not */
	  /* overwritten before it is used */
	  g_sync();

/* SF:
dir=t
The "scattering back" must be done at t=0 as well: U_t(t=0) is modified
by gauge transformations on sites with t=1, but must be saved back to t=0.
So there should be no 'DOMAIN' for this operation.

At t=nt-1 the copied-back U_t are never modified
since we don't perform gauge transformations at t=0=nt.

dir=x,y,z
For the same reason, U_dir(t=0) is always copied back unmodified.
*/

	  /* First copy modified link for this dir */
	  /* from comm buffer or node to diffmat */

	  FORSOMEPARITY(i,s,parity)
	  {
	      mat_copy_f((matrix_f *)(gen_pt[dir][i]), &diffmatp[i]);
	  }

	  /* Now we are finished with gen_pt[dir] */
	  cleanup_gather(mtag[dir]);

	  /* Synchronize to make sure the previous copy happens before the */
	  /* subsequent gather below  */
	  g_sync();

	  /* Gather diffmat onto sites of opposite parity */
	  mtag[dir] = start_gather_field( diffmatp, sizeof(matrix_f),
				      dir, OPP_PAR(parity), gen_pt[dir] );

	  wait_gather(mtag[dir]);

         /* Copy modified matrices into proper location */

         FORSOMEPARITY(i,s,OPP_PAR(parity))
	      mat_copy_f((matrix_f *)(gen_pt[dir][i]),&(s->link[dir]));

	  cleanup_gather(mtag[dir]);
	}

    }
} /* gaugefixstep */

void gaugefixscratch()
{
      diffmatp = (matrix_f *)malloc(sizeof(matrix_f)*sites_on_node);
      if(diffmatp == NULL)
	{
	  node0_printf("gaugefix: Can't malloc diffmat\n");
	  fflush(stdout);terminate(1);
	}

      sumvecp = (vector_f *)malloc(sizeof(vector_f)*sites_on_node);
      if(sumvecp == NULL)
	{
	  node0_printf("gaugefix: Can't malloc sumvec\n");
	  fflush(stdout);terminate(1);
	}
} /* gaugefixscratch */

void gaugefix_field(int gauge_dir,Real relax_boost,int max_gauge_iter,
	      Real gauge_fix_tol)
{
  int gauge_iter;
  double current_av, old_av = 0, del_av = 0;

  /* We require at least 8 gen_pt values for gauge fixing */
  if(N_POINTERS < 8)
    {
      printf("gaugefix: N_POINTERS must be at least %d.  Fix the code.\n",
	     N_POINTERS);
      fflush(stdout); terminate(1);
    }


  /* Set up work space */
  gaugefixscratch();

  /* Do at most max_gauge_iter iterations, but stop after the second step if */
  /* the change in the avg gauge fixing action is smaller than gauge_fix_tol */

  for (gauge_iter=0; gauge_iter < max_gauge_iter; gauge_iter++)
    {
      gaugefixstep(gauge_dir,&current_av,relax_boost);

      if(gauge_iter != 0)
	{
	  del_av = current_av - old_av;
	  if (fabs(del_av) < gauge_fix_tol) break;
	}
      old_av = current_av;

      /* Reunitarize when iteration count is a multiple of REUNIT_INTERVAL */
      if((gauge_iter % REUNIT_INTERVAL) == (REUNIT_INTERVAL - 1))
	{
/**	  node0_printf("step %d av gf action %.8e, delta %.3e\n",
		       gauge_iter,current_av,del_av); **/
	  reunitarize();
	}
    }
  /* Reunitarize at the end, unless we just did it in the loop */
  if((gauge_iter % REUNIT_INTERVAL) != 0)
    reunitarize();

  /* Free workspace */
  free(diffmatp);
  free(sumvecp);

  if(this_node==0)
    printf("GFIX: Ended at step %d. Av gf action %.8e, delta %.3e\n",
	   gauge_iter,(double)current_av,(double)del_av);
}

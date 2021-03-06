/********** ape_smear.c*********************************************/

/* Does APE smearing */
/* Doug Toussaint 1/30/96 
   CD 11/25/01 Consolidated with other smearing routines 
   CD 11/15/01 malloc temp field

   Entry point ape_smear_dir for smearing links in one specified direction.
   Entry point ape_smear for smearing links in all four directions.

   Construct smeared links.   Smeared link is sum
   of link plus all staples:

                ---------               -------> dir1
                |       |
                |       |
                |       |               ^
                X--------               | dir2
                |       |               |
                |       |
                |       |
                ---------

   staple_weight is the weight of a staple relative to the original link

   That is we compute

   U_smeared = w_link U_link + w_staple * sum U_staple

   where w_staple/w_link = staple_weight

   w_link + w_staple is normalized so that

   w_link + nstaples * link_u0 * link_u0 * w_staple = 1

   where nstaples = 6 usually, but nstaples = 4 for spacelike
   links when called with space_only = 1.

   If an SU(3) projection is done, only the relative weight matters.

   This scheme takes care of four conventions in MILC use:

   (1) smooth_inst application

       with link_u0 = 1 and 
       staple_weight = ape_weight/(nstaples*(1 - ape_weight))
       we have
       w_link = 1 - ape_weight
       w_staple = ape_weight/nstaples

   (2) spectrum_hybrids4.c

       with link_u0 = u0
       we have
       simple_weight = 1/staple_weight
       norm_factor = 1/(nstaples * u0 * u0 + simple_weight)
       w_link = simple_weight * norm_factor
       w_staple = norm_factor

   (3) hvy_qpot and string_break application

       use 
       link_u0 = 1  (doesn't matter, because we project)
       staple_weight = 1/smear_fac
       space_only = 1

   (4) wilson_hybrids and clover_hybrids application
  
       use
       link_u0 
         = sqrt((1.0 - norm_factor*simple_weight)/(norm_factor*nstaples))
       staple_weight = 1/simple_weight
       space_only = 1
       nhits = 0
*/

#include "generic_includes.h"

/* Smear in a specified source direction. */
void ape_smear_dir(
  field_offset src,       /* field offset for matrix[4] type 
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  field_offset dest,      /* field offset for matrix type 
			     pointing to a specific direction 
			     output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  register int i,dir2;
  register site *s;
  matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  Real w_link, w_staple, norm_factor;
  int nstaples;
  matrix *temp;
  
  /* Allocate temporary space for staple calculation */
  temp = (matrix *)malloc(sites_on_node*sizeof(matrix));
  if(temp == NULL){
    printf("ape_smear: No room for temp\n");
    terminate(1);
  }
  
  nstaples = (space_only==1 && dir1 != TUP) ? 4 : 6;
  norm_factor = 1.0/(staple_weight*nstaples*link_u0*link_u0 + 1.0);
  w_staple = norm_factor*staple_weight;
  w_link = norm_factor;
  
  /* dest <- src w_link */ 
  FORALLSITES(i,s)
    {
      scalar_mult_matrix( 
		     &(((matrix *)F_PT(s,src))[dir1]), w_link, 
		     (matrix *)F_PT(s,dest) );
    }
  for(dir2=XUP;dir2<=(space_only==1?ZUP:TUP);dir2++)if(dir2!=dir1){
    
    /* Upper staple, and simple link */
    mtag0 = start_gather_site( src+dir2*sizeof(matrix),
			  sizeof(matrix), dir1, EVENANDODD, gen_pt[0] );
    mtag1 = start_gather_site( src+dir1*sizeof(matrix),
			  sizeof(matrix), dir2, EVENANDODD, gen_pt[1] );
    wait_gather(mtag0);
    wait_gather(mtag1);
    
    /* dest += w_staple * upper staple */ 
    FORALLSITES(i,s)
      {
	mult_na( (matrix *)gen_pt[1][i],
		     (matrix *)gen_pt[0][i], &tmat1 );
	mult_nn( &(((matrix *)F_PT(s,src))[dir2]), 
		     &tmat1, &tmat2 );
	scalar_mult_add_matrix( 
			   (matrix *)F_PT(s,dest),
			   &tmat2, w_staple,
			   (matrix *)F_PT(s,dest) );
      }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    /* lower staple */
    mtag0 = start_gather_site( src+dir2*sizeof(matrix),
			  sizeof(matrix), dir1,
			  EVENANDODD, gen_pt[0] );
    wait_gather(mtag0);
    FORALLSITES(i,s)
      {
	mult_nn( &(((matrix *)F_PT(s,src))[dir1]),
		     (matrix *)gen_pt[0][i], &tmat1 );
	mult_an( &(((matrix *)F_PT(s,src))[dir2]),
		     &tmat1, &temp[i] );
      }
    cleanup_gather(mtag0);
    mtag1 = start_gather_field( temp, sizeof(matrix),
				    OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    
    /* dest += w_staple * lower staple */ 
    FORALLSITES(i,s){
      scalar_mult_add_matrix( 
			 (matrix *)F_PT(s,dest),
			 (matrix *)gen_pt[1][i], w_staple, 
			 (matrix *)F_PT(s,dest) );
    }
    cleanup_gather(mtag1);
    
  } /* dir2 loop */
  
  /* project links onto SU(3) if nhits > 0 */
  if(nhits > 0){
    FORALLSITES(i,s){
      /* Use partially reunitarized link for guess */
      tmat1 = *((matrix *)F_PT(s,dest));
      reunit_su3(&tmat1);
      project_su3(&tmat1,
		  (matrix *)F_PT(s,dest),nhits,tol);
      /* Copy projected matrix to dest */
      *((matrix *)F_PT(s,dest)) = tmat1;
    }
  }
  
  free(temp);
  
} /* ape_smear_dir */


void ape_smear(
  field_offset src,       /* field offset for matrix type 
			     input unsmeared links */
  field_offset dest,      /* field offset for matrix type 
			     output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  register int dir1;
  
  for(dir1=XUP;dir1<=TUP;dir1++){
    ape_smear_dir(src,dir1,dest+dir1*sizeof(matrix),
		  staple_weight,link_u0,space_only,nhits,tol);
  }
} /* ape_smear */

/* Smear in a specified source direction. 
   Gauge matrices are stored four per site. They are taken at
   stride 4 (matrices) starting from the matrix src[dir] 
   The result is stored likewise in dest */

/* NOTE: This code can be unified with the ape_smear_dir code if we
   specify the start and stride and use declare_strided_gather
   throughout - CD */

void ape_smear_field_dir(
  matrix *src,        /* matrix[4] type 
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  matrix *dest,       /* matrix[4] type smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  register int i,dir2;
  register site *s;
  matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  Real w_link, w_staple, norm_factor;
  int nstaples;
  matrix *temp;
  
  /* Allocate temporary space for staple calculation */
  temp = (matrix *)malloc(sites_on_node*sizeof(matrix));
  if(temp == NULL){
    printf("ape_smear: No room for temp\n");
    terminate(1);
  }
  
  nstaples = (space_only==1 && dir1 != TUP) ? 4 : 6;
  norm_factor = 1.0/(staple_weight*nstaples*link_u0*link_u0 + 1.0);
  w_staple = norm_factor*staple_weight;
  w_link = norm_factor;
  
  /* dest <- src w_link */ 
  FORALLSITES(i,s)
    {
      scalar_mult_matrix( &src[4*i+dir1], w_link, &dest[4*i+dir1] );
    }
  for(dir2=XUP;dir2<=(space_only==1?ZUP:TUP);dir2++)if(dir2!=dir1){
    
    /* Upper staple, and simple link */
    mtag0 = declare_strided_gather( (char *)&src[dir2], 4*sizeof(matrix),
				    sizeof(matrix), dir1, 
				    EVENANDODD, gen_pt[0] );
    prepare_gather(mtag0);
    do_gather(mtag0);

    mtag1 = declare_strided_gather( (char *)&src[dir1], 4*sizeof(matrix),
				    sizeof(matrix), dir2, 
				    EVENANDODD, gen_pt[1] );
    prepare_gather(mtag1);
    do_gather(mtag1);

    wait_gather(mtag0);
    wait_gather(mtag1);
    
    /* dest += w_staple * upper staple */ 
    FORALLSITES(i,s)
      {
	mult_na( (matrix *)gen_pt[1][i],
		     (matrix *)gen_pt[0][i], &tmat1 );
	mult_nn( &src[4*i+dir2], &tmat1, &tmat2 );
	scalar_mult_add_matrix( &dest[4*i+dir1], &tmat2, w_staple,
				    &dest[4*i+dir1] );
      }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    /* lower staple */
    mtag0 = declare_strided_gather( (char *)&src[dir2], 4*sizeof(matrix),
				    sizeof(matrix), dir1,
				    EVENANDODD, gen_pt[0] );
    prepare_gather(mtag0);
    do_gather(mtag0);

    wait_gather(mtag0);
    FORALLSITES(i,s)
      {
	mult_nn( &src[4*i+dir1], (matrix *)gen_pt[0][i], &tmat1 );
	mult_an( &src[4*i+dir2], &tmat1, &temp[i] );
      }
    cleanup_gather(mtag0);
    mtag1 = start_gather_field( temp, sizeof(matrix),
				    OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    
    /* dest += w_staple * lower staple */ 
    FORALLSITES(i,s){
      scalar_mult_add_matrix( &dest[4*i+dir1], 
			 (matrix *)gen_pt[1][i], w_staple, 
			 &dest[4*i+dir1] );
    }
    cleanup_gather(mtag1);
    
  } /* dir2 loop */
  
  /* project links onto SU(3) if nhits > 0 */
  if(nhits > 0){
    FORALLSITES(i,s){
      /* Use partially reunitarized link for guess */
      tmat1 = dest[4*i+dir1];
      reunit_su3(&tmat1);
      project_su3(&tmat1, &dest[4*i+dir1], nhits, tol);
      /* Copy projected matrix to dest */
      dest[4*i+dir1] = tmat1;
    }
  }
  
  free(temp);
  
} /* ape_smear_dir */

/* Input field has four contigous SU(3) matrices per site */

void ape_smear_field(
  matrix *src,       /* Gauge field input unsmeared */
  matrix *dest,      /* Gauge field output smeared */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  register int dir1;
  
  for(dir1=XUP;dir1<=TUP;dir1++){
    ape_smear_field_dir(src,dir1,dest,staple_weight,
			link_u0,space_only,nhits,tol);
  }
} /* ape_smear */


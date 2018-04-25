/*********************** rand_gauge.c ***************************/
/* original code by UMH */
/* 2/19/98 Version 5 port CD */

/* Makes a random gauge transformation on the gauge fields
   and then reunitarizes them */
/* Warning KS fermion applications: Must be done with KS phases OUT! */
/* Requires site workspace G */

#include "generic_includes.h"

void randomize(field_offset G, Real radius);
void gauge_trans(field_offset G);

void rand_gauge(field_offset G)
	/* G holds the gauge transformation matrices */
{
    randomize(G, 1.0);
    gauge_trans(G);
    reunitarize();
}

void randomize(field_offset G, Real radius)
{
  register int a, b, i;
  site *s;

  FORALLSITES(i,s) {
    for(a=0; a<3; a++) for(b=0; b<3; b++)
      (*(matrix *)F_PT(s,G)).e[a][b] 
             = cmplx(radius*((Real)drand48()-0.5),
                     radius*((Real)drand48()-0.5));
    reunit_su3((matrix *)F_PT(s,G));
  }
}

void gauge_trans(field_offset G)
{
  register int i,mu;
  site *s;
  matrix tmp;
  msg_tag *tag[4];

  FORALLUPDIR(mu) 
    tag[mu] = start_gather_site(G,sizeof(matrix),mu,EVENANDODD,
		       gen_pt[mu]);

  FORALLUPDIR(mu) {
    wait_gather(tag[mu]);
    FORALLSITES(i,s) {

       mult_an((matrix *)F_PT(s,G), &(s->link[mu]), &tmp);
       mult_nn(&tmp, (matrix *)gen_pt[mu][i],
		       &(s->link[mu]));

    }
    cleanup_gather(tag[mu]);
  }
}

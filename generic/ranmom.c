/*************************** ranmom.c *******************************/
/* MIMD version 7 */

/* Produce Gaussian random momenta for the gauge fields. */

#include "generic_includes.h"
#include <defines.h>                 /* For SITERAND */

void ranmom(){
register int i,dir;
register site *s;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
#ifdef SCHROED_FUN
	    if(dir==TUP || s->t>0){
#endif
#ifdef SITERAND
		random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
		    &(s->site_prn) );
#else
		random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
		    &node_prn );
#endif
/*  It would not matter if we did not set the spatial momenta at t=0
    to zero, because we never attempt to update the t=0 spatial links
*/
#ifdef SCHROED_FUN
	    }
	    else{
		s->mom[dir].m00im = 0.0;
		s->mom[dir].m11im = 0.0;
		s->mom[dir].m01.real = 0.0;
		s->mom[dir].m01.imag = 0.0;
#if (NCOL>2)
		s->mom[dir].m22im = 0.0;
		s->mom[dir].m02.real = 0.0;
		s->mom[dir].m02.imag = 0.0;
		s->mom[dir].m12.real = 0.0;
		s->mom[dir].m12.imag = 0.0;
#if (NCOL>3)
		s->mom[dir].m33im = 0.0;
		s->mom[dir].m03.real = 0.0;
		s->mom[dir].m03.imag = 0.0;
		s->mom[dir].m13.real = 0.0;
		s->mom[dir].m13.imag = 0.0;
		s->mom[dir].m23.real = 0.0;
		s->mom[dir].m23.imag = 0.0;
#endif
#endif
	    }
#endif
	}
    }
}


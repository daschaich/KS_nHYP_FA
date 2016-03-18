/******** layout_fixedmap.c *********/
/* MIMD version 6 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

#define XBLOCKS 1
#define YBLOCKS 8
#define ZBLOCKS 8
#define TBLOCKS 16

/* This version divides the lattice any of the four directions.
   To fit on a machine with fixed topology, the divisors of each
   of the directions are fixed.
*/

/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
     This routine sets the global variables "sites_on_node",
     "even_sites_on_node" and "odd_sites_on_node".
   num_sites(node) returns the number of sites on a node
   node_number(x,y,z,t) returns the node number on which a site lives.
   node_index(x,y,z,t) returns the index of the site on the node - ie the
     site is lattice[node_index(x,y,z,t)].
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"

int squaresize[4];	/* dimensions of hypercubes */
int nsquares[4];	/* number of hypercubes in each direction */

void setup_layout(){
register int i,j,k,dir;
    node0_printf("LAYOUT = Fixed:  %d x %d x %d x %d\n",XBLOCKS,YBLOCKS,ZBLOCKS,TBLOCKS);

    /* Figure out dimensions of rectangle */
    if( XBLOCKS*YBLOCKS*ZBLOCKS*TBLOCKS != numnodes() ){
       node0_printf("DUMMY, wrong number of nodes for this layout\n");
       terminate(0);
    }
    if( nx%XBLOCKS != 0 || ny%YBLOCKS !=0 || nz%ZBLOCKS != 0 || nt%TBLOCKS != 0){
       node0_printf("DUMMY, can't lay out this lattice with this grid\n");
       terminate(0);
    }
    nsquares[XUP] = XBLOCKS;
    nsquares[YUP] = YBLOCKS;
    nsquares[ZUP] = ZBLOCKS;
    nsquares[TUP] = TBLOCKS;
    squaresize[XUP] = nx/XBLOCKS; squaresize[YUP] = ny/YBLOCKS;
    squaresize[ZUP] = nz/ZBLOCKS; squaresize[TUP] = nt/TBLOCKS;

    sites_on_node =
	    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];
    /* Need even number of sites per hypercube */
    if( mynode()==0)if( sites_on_node%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
	terminate(0);
    }
if( mynode()==0)
  printf("ON EACH NODE %d x %d x %d x %d\n",squaresize[XUP],squaresize[YUP],
                squaresize[ZUP],squaresize[TUP]);
if( mynode()==0 && sites_on_node%2 != 0)
	printf("WATCH OUT FOR EVEN/ODD SITES ON NODE BUG!!!\n");
    even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

int node_number(int x,int y,int z,int t) {
register int i;
    x /= squaresize[XUP]; y /= squaresize[YUP];
    z /= squaresize[ZUP]; t /= squaresize[TUP];
    i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
    return( i );
}

int node_index(int x,int y,int z,int t) {
register int i,xr,yr,zr,tr;
    xr = x%squaresize[XUP]; yr = y%squaresize[YUP];
    zr = z%squaresize[ZUP]; tr = t%squaresize[TUP];
    i = xr + squaresize[XUP]*( yr + squaresize[YUP]*( zr + squaresize[ZUP]*tr));
    if( (x+y+z+t)%2==0 ){	/* even site */
	return( i/2 );
    }
    else {
	return( (i + sites_on_node)/2 );
    }
}

unsigned int num_sites(int node) {
    return( sites_on_node );
}


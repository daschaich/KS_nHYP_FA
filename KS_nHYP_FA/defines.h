// -----------------------------------------------------------------
// Compiler macros common to all targets in this application
#ifndef _DEFINES_H
#define _DEFINES_H

#define SITERAND          // Use site-based random number generators
//#define TIMING
//#define CGTIME
//#define CG_DEBUG
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fix maximum number of staggered field pairs compile time
#define MAX_FIELDS 2
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// nHYP stuff
#define NHYP_DEBUG

// Hard-coded smearing parameters
// Uncomment the following line to hard-code smearing parameters
// Otherwise they are read from the infile
//#define HARD_CODE_SMEAR
//#ifdef CONTROL
//Real alpha_smear[3] = {0.5, 0.5, 0.4};
//#else
//extern Real alpha_smear[3];
//#endif

// IR regulator:  Q = Omega^dag Omega + IR_STAB
// This slightly changes the definition of the nHYP link. Fine.
// Since we're adding a constant, any derivative of Q is unchanged.
// Used in block_nhyp.c, force_nhyp.c
#define IR_STAB 1.0e-6

// Neighborhood of 0 where we use approximate R and S for u0, u1, p
// Bypass R / (-S)^{3 / 2} = 0 / 0 in calculation of Q^{-1/2}
// Used in nhyp.c
#define EPS_SQ 1.0e-5
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Integrator stuff
// Maximum number of Hasenbusch masses
#define MAX_MASSES 2

// Omelyan lambda, 2lambda and 1 - 2lambda
#define LAMBDA 0.193
#define TWO_LAMBDA 0.386
#define LAMBDA_MID 0.614
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#endif
// -----------------------------------------------------------------

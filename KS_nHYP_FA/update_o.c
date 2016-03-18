// -----------------------------------------------------------------
// Update lattice, with Omelyan integrator and Hasenbusch masses

// Note that for the final accept/reject,
// we already have a good solution to the CG:
// The last update was of the momenta

// The input file gave us the number of levels, trajectory length, and
// nsteps[MAX_MASSES + 1] (the gauge level is the last)
// num_masses <= MAX_MASSES is the total number of
// real dynamical fermion doublets plus the Hasenbusch preconditioners

// This routine begins at "integral" time
// with H and U evaluated at same time

// Includes and definitions
#include "ks_dyn_includes.h"

// Local function prototypes
void predict_next_psi(Real *oldtime, Real *newtime,
                      Real *nexttime, int level);
int update_step_two(Real *oldtime, Real *newtime, Real *nexttime,
                    double *fnorm, double *gnorm);
int update_step_three(Real *oldtime, Real *newtime, Real *nexttime,
                      double *fnorm, double *gnorm);
double update_gauge_step(Real eps);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Helper wrapper for staggered conjugate gradient
void doCG(int level, Real M, int *iters) {
#ifdef CG_DEBUG
  node0_printf("UPDATE: calling congrad for level %d\n", level);
#endif
  int j;

  block_N_fatten(Nsmear);
  for (j = 0; j < full_fields; j++)
    *iters += ks_congrad(F_OFFSET(chi[j][level]), F_OFFSET(psi[j][level]), M,
                         EVENANDODD);

  if (half_fields == 1) {
    j = full_fields;
    *iters += ks_congrad(F_OFFSET(chi[j][level]), F_OFFSET(psi[j][level]), M,
                         EVEN);
  }
  leanlinks();
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int iters = 0, j;
  Real cg_time[2], old_cg_time[2], next_cg_time[2];
#ifdef HMC_ALGORITHM
  double startaction, endaction, change;
  Real xrandom;
#endif
  double gnorm = 0.0, fnorm[2];
  fnorm[0] = 0.0;
  fnorm[1] = 0.0;
  max_gf = 0.0;
  max_ff[0] = 0.0;
  max_ff[1] = 0.0;

  // Refresh the momenta
  ranmom();

  // Generate a pseudofermion configuration only at start
  // Also clear psi fields -- zero is best guess with new source
  block_N_fatten(Nsmear);
  for (j = 0; j <= full_fields; j++)
    clear_latvec(F_OFFSET(psi[j][0]), EVENANDODD);

  if (num_masses == 1) {                    // The simple case
    for (j = 0; j < full_fields; j++)
      grsource_imp(F_OFFSET(chi[j][0]), mass, EVENANDODD);
    if (half_fields == 1) {
      j = full_fields;
      grsource_imp(F_OFFSET(chi[j][0]), mass, EVEN);
    }
  }
  else {  // num_masses = 2, the more complicated case
    for (j = 0; j <= full_fields; j++)
      clear_latvec(F_OFFSET(psi[j][1]), EVENANDODD);

    for (j = 0; j < full_fields; j++) {
      grsource_Hasen(F_OFFSET(chi[j][0]), &iters, EVENANDODD);
      grsource_imp(F_OFFSET(chi[j][1]), MH, EVENANDODD);
    }
    if (half_fields == 1) {
      j = full_fields;
      grsource_Hasen(F_OFFSET(chi[j][0]), &iters, EVEN);
      grsource_imp(F_OFFSET(chi[j][1]), MH, EVEN);
    }
  }

  leanlinks();
  old_cg_time[0] = -1;
  old_cg_time[1] = -1;
  cg_time[0] = -1;
  cg_time[1] = -1;

  // Do conjugate gradient to get psi = (M^dag M)^(-1) chi
  doCG(0, mass, &iters);
  if (num_masses == 2)
    doCG(1, MH, &iters);

  cg_time[0] = 0;
  cg_time[1] = 0;

#ifdef HMC_ALGORITHM    // Find action
  startaction = d_action();
  // node0_printf("startaction= %g\n", startaction);

  // Copy link field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));
#endif

  // Do microcanonical updating
  if (num_masses == 2)
    iters += update_step_three(old_cg_time, cg_time, next_cg_time,
                               fnorm, &gnorm);
  else
    iters += update_step_two(old_cg_time, cg_time, next_cg_time,
                             fnorm, &gnorm);

#ifdef HMC_ALGORITHM    // Find action, then accept or reject
  // Do conjugate gradient to get (M^dag M)^(-1) chi
  doCG(0, mass, &iters);
  if (num_masses == 2)
    doCG(1, MH, &iters);

  endaction = d_action();
  // printf("endaction= %g\n", endaction);

  change = endaction - startaction;

  // Reject configurations giving overflow
  if (fabs((double)change) > 1e20) {
    node0_printf("WARNING: Correcting apparent overflow: Delta S = %.4g\n",
                 change);
    change = 1e20;
  }

  // Decide whether to accept.  If not, copy old link field back
  // Careful: must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);

  broadcast_float(&xrandom);
  if (exp(-change) < (double)xrandom) {
    gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));
    node0_printf("REJECT: ");
  }
  else
    node0_printf("ACCEPT: ");

  node0_printf("delta S = %.4g start %.12g end %.12g\n",
               change, startaction, endaction);
#endif // HMC

  node0_printf("MONITOR_FORCE_FERMION0 %.4g %.4g\n",
               fnorm[0] / (double)(2 * nsteps[0]), max_ff[0]);
  node0_printf("MONITOR_FORCE_FERMION1 %.4g %.4g\n",
               fnorm[1] / (double)(4 * nsteps[0] * nsteps[1]), max_ff[1]);
  // gnorm divided by nsteps_gauge every time gauge_update_step called
  node0_printf("MONITOR_FORCE_GAUGE    %.4g %.4g\n",
               gnorm / (double)(4 * nsteps[0] * nsteps[1]), max_gf);

  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use linear extrapolation to predict next conjugate gradient solution
// Need even and odd sites depending on Nf
void predict_next_psi(Real *oldtime, Real *newtime, Real *nexttime,
                      int level) {

  register int i;
  register site *s;
  register Real x;
  int j;
  su3_vector tvec;

  if (newtime[level] != oldtime[level])
    x = (nexttime[level] - newtime[level])
      / (newtime[level] - oldtime[level]);
  else
    x = 0;

  for (j = 0; j < full_fields; j++) {
    if (oldtime[level] < 0) {
      FORALLSITES(i, s)
        s->old_psi[j][level] = s->psi[j][level];
    }
    else {
      FORALLSITES(i, s) {
        sub_su3_vector(&(s->psi[j][level]), &(s->old_psi[j][level]), &tvec);
        s->old_psi[j][level] = s->psi[j][level];
        scalar_mult_add_su3_vector(&(s->psi[j][level]), &tvec, x,
                                   &(s->psi[j][level]));
      }
    }
  }
  if (half_fields == 1) {
    j = full_fields;
    if (oldtime[level] < 0) {
      FOREVENSITES(i, s)
        s->old_psi[j][level] = s->psi[j][level];
    }
    else {
      FOREVENSITES(i, s) {
        sub_su3_vector(&(s->psi[j][level]), &(s->old_psi[j][level]), &tvec);
        s->old_psi[j][level] = s->psi[j][level];
        scalar_mult_add_su3_vector(&(s->psi[j][level]), &tvec, x,
                                   &(s->psi[j][level]));
      }
    }
  }
  oldtime[level] = newtime[level];
  newtime[level] = nexttime[level];
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Two-level Omelyan integrator for single-pseudofermion code
int update_step_two(Real *old_cg_time, Real *cg_time, Real *next_cg_time,
                    double *fnorm, double *gnorm) {

  if (num_masses == 2) {
    node0_printf("Warning: running two-level integrator with 1 PF\n");
    exit(1);
  }
  int iters = 0, n, tot_steps = 2 * nsteps[0] * nsteps[1];
  Real f_eps, g_eps, tr;

  f_eps = traj_length / ((Real) tot_steps);
  g_eps = f_eps / (2 * (Real)nsteps[MAX_MASSES]);

  // Initial lambda fermion update step
  // (Already did CG to get starting action)
  tr = fermion_force(0, f_eps * LAMBDA);
  fnorm[0] += tr;
  if (tr > max_ff[0])
    max_ff[0] = tr;

  for (n = 1; n <= tot_steps; n++) {
    tr = update_gauge_step(g_eps);
    *gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

    // (1 - 2lambda) fermion force update step
    next_cg_time[0] = cg_time[0] + f_eps;
    predict_next_psi(old_cg_time, cg_time, next_cg_time, 0);
    doCG(0, mass, &iters);
    tr = fermion_force(0, f_eps * LAMBDA_MID);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;

    tr = update_gauge_step(g_eps);
    *gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

    // 2lambda fermion force update step except on last iteration
    // On last iteration, just lambda fermion force update step
    next_cg_time[0] = cg_time[0] + f_eps;
    predict_next_psi(old_cg_time, cg_time, next_cg_time, 0);
    doCG(0, mass, &iters);
    if (n < tot_steps)
      tr = fermion_force(0, f_eps * TWO_LAMBDA);
    else
      tr = fermion_force(0, f_eps * LAMBDA);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;
  }
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Three-level Omelyan integrator for Hasenbusch-preconditioned code
// Set level to use either two or (inefficiently) one pseudofermion
int update_step_three(Real *old_cg_time, Real *cg_time, Real *next_cg_time,
                      double *fnorm, double *gnorm) {

  if (num_masses != 2) {
    node0_printf("Error: need 2 PF for three-level integrator\n");
    exit(1);
  }
  int iters = 0, outer, inner;
  Real f_eps0, f_eps1, g_eps, tr;

  f_eps0 = traj_length / (Real)nsteps[0];
  f_eps1 = f_eps0 / (2 * (Real)nsteps[1]);
  g_eps = f_eps1 / (2 * (Real)nsteps[MAX_MASSES]);

  // Initial lambda update steps for both inner and outer force loops
  // (Already did CG on both levels to get starting action)
  tr = fermion_force(0, f_eps0 * LAMBDA);
  fnorm[0] += tr;
  if (tr > max_ff[0])
    max_ff[0] = tr;
  tr = fermion_force(1, f_eps1 * LAMBDA);
  fnorm[1] += tr;
  if (tr > max_ff[1])
    max_ff[1] = tr;
  for (outer = 1; outer <= nsteps[0]; outer++) {
    for (inner = 1; inner <= nsteps[1]; inner++) {
      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // (1 - 2lambda) inner force update step
      next_cg_time[1] = cg_time[1] + f_eps1;
      predict_next_psi(old_cg_time, cg_time, next_cg_time, 1);
      doCG(1, MH, &iters);
      tr = fermion_force(1, f_eps1 * LAMBDA_MID);
      fnorm[1] += tr;
      if (tr > max_ff[1])
        max_ff[1] = tr;

      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // 2lambda inner force update step, since combined
      // with initial update step of next loop
      next_cg_time[1] = cg_time[1] + f_eps1;
      predict_next_psi(old_cg_time, cg_time, next_cg_time, 1);
      doCG(1, MH, &iters);
      tr = fermion_force(1, f_eps1 * TWO_LAMBDA);
      fnorm[1] += tr;
      if (tr > max_ff[1])
        max_ff[1] = tr;
    } // Inner update steps

    // (1 - 2lambda) outer force update step
    next_cg_time[0] = cg_time[0] + f_eps0;
    predict_next_psi(old_cg_time, cg_time, next_cg_time, 0);
    doCG(0, mass, &iters);
    tr = fermion_force(0, f_eps0 * LAMBDA_MID);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;

    // Initial lambda update step done above
    for (inner = 1; inner <= nsteps[1]; inner++) {
      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // (1 - 2lambda) inner force update step
      next_cg_time[1] = cg_time[1] + f_eps1;
      predict_next_psi(old_cg_time, cg_time, next_cg_time, 1);
      doCG(1, MH, &iters);
      tr = fermion_force(1, f_eps1 * LAMBDA_MID);
      fnorm[1] += tr;
      if (tr > max_ff[1])
        max_ff[1] = tr;

      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // 2lambda inner force update step, combined with
      // initial update step of next outer,
      // except on very last last iteration
      next_cg_time[1] = cg_time[1] + f_eps1;
      predict_next_psi(old_cg_time, cg_time, next_cg_time, 1);
      doCG(1, MH, &iters);
      if (inner < nsteps[1] || outer < nsteps[0])
        tr = fermion_force(1, f_eps1 * TWO_LAMBDA);
      else
        tr = fermion_force(1, f_eps1 * LAMBDA);
      fnorm[1] += tr;
      if (tr > max_ff[1])
        max_ff[1] = tr;
    } // Inner update steps

    // Except on last iteration, 2lambda outer force update step
    // On last iteration, just lambda outer force update step
    next_cg_time[0] = cg_time[0] + f_eps0;
    predict_next_psi(old_cg_time, cg_time, next_cg_time, 0);
    doCG(0, mass, &iters);
    if (outer < nsteps[0])
      tr = fermion_force(0, f_eps0 * TWO_LAMBDA);
    else
      tr = fermion_force(0, f_eps0 * LAMBDA);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;
  } // Outer loop
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Omelyan version; ``dirty'' speeded-up version
// Returns norm of force for averaging
double update_gauge_step(Real eps) {
  int n = nsteps[MAX_MASSES], i;
  double norm = 0;

#ifdef CG_DEBUG
  node0_printf("UPDATE: %d gauge steps of %.4g\n", n, eps);
#endif

  norm += gauge_force(eps * LAMBDA);
  for (i = 1; i <= n; i++) {
    update_u(0.5 * eps);
    norm += gauge_force(eps * LAMBDA_MID);
    update_u(0.5 * eps);

    // 2lambda step except on last iteration; then only lambda
    if (i < n)
      norm += gauge_force(eps * TWO_LAMBDA);
    else
      norm += gauge_force(eps * LAMBDA);
  }

  // Reunitarize the gauge field
  rephase(OFF);
  reunitarize();
  rephase(ON);
  return norm / n;
}
// -----------------------------------------------------------------

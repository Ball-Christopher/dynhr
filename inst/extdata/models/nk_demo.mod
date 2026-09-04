// @dynhr-model
// name: Three-equation New Keynesian demonstration model
// short: nk_demo
// tier: 2
// endo: 8
// exo: 3
// params: 11
// has_lead: true
// has_measurement: true
// has_estimation: true
// has_analytical_ss: true
// features: perturbation, kalman_filter, estimation, stoch_simul
// refs: Standard textbook three-equation New Keynesian model. See e.g.
//       Gali (2015) "Monetary Policy, Inflation, and the Business Cycle",
//       2nd ed., ch. 3, or Woodford (2003) "Interest and Prices", ch. 4,
//       for the derivation. Cited for the ECONOMICS only -- see
//       nk_demo_SOURCE.md: this file and its dataset were written for dynhr
//       and are covered by the package's MIT licence.
// notes: The bundled worked example for the README quick start and the
//        estimation/diagnostics vignettes. Deliberately small so a full
//        parse -> compile -> steady state -> mode-finding -> sampling pass
//        runs in a few seconds from a clean install.
//
//        Three shocks and three observables, so the observation block is
//        non-singular and the exact Gaussian likelihood is available
//        without a measurement-error floor.
// @dynhr-model-end

// ---------------------------------------------------------------------------
// Linearised around the zero-inflation steady state; y, pi, r, g, u are
// log-deviations (SS = 0). The three observables carry the steady-state
// levels, so the data look like conventional macro series.
//
//   IS curve      y_t  = E_t y_{t+1} - sigma (r_t - E_t pi_{t+1}) + g_t
//   Phillips      pi_t = beta E_t pi_{t+1} + kappa y_t + u_t
//   Taylor rule   r_t  = rho_r r_{t-1}
//                        + (1 - rho_r)(phi_pi pi_t + phi_y y_t) + e_m,t
//   Demand        g_t  = rho_g g_{t-1} + e_g,t
//   Cost-push     u_t  = rho_u u_{t-1} + e_u,t
//
//   ygr_t  = gam_bar + y_t - y_{t-1}     quarterly output growth, %
//   infl_t = pi_bar  + 4 pi_t            annualised inflation, %
//   intr_t = r_bar   + 4 r_t             annualised nominal rate, %
// ---------------------------------------------------------------------------

var y pi r g u ygr infl intr;

varexo e_g e_u e_m;

parameters beta kappa sigma rho_r rho_g rho_u phi_pi phi_y
           gam_bar pi_bar r_bar;

beta    = 0.99;    // discount factor
kappa   = 0.10;    // slope of the Phillips curve
sigma   = 1.00;    // intertemporal elasticity of substitution
rho_r   = 0.75;    // interest-rate smoothing
rho_g   = 0.85;    // demand-shock persistence
rho_u   = 0.50;    // cost-push persistence
phi_pi  = 1.50;    // Taylor coefficient on inflation
phi_y   = 0.25;    // Taylor coefficient on the output gap
gam_bar = 0.50;    // mean quarterly output growth, %
pi_bar  = 2.00;    // mean annualised inflation, %
r_bar   = 4.00;    // mean annualised nominal rate, %

model;
  y    = y(+1) - sigma * (r - pi(+1)) + g;
  pi   = beta * pi(+1) + kappa * y + u;
  r    = rho_r * r(-1) + (1 - rho_r) * (phi_pi * pi + phi_y * y) + e_m;
  g    = rho_g * g(-1) + e_g;
  u    = rho_u * u(-1) + e_u;
  ygr  = gam_bar + y - y(-1);
  infl = pi_bar + 4 * pi;
  intr = r_bar + 4 * r;
end;

initval;
  y    = 0;
  pi   = 0;
  r    = 0;
  g    = 0;
  u    = 0;
  ygr  = gam_bar;
  infl = pi_bar;
  intr = r_bar;
end;

shocks;
  var e_g; stderr 0.30;
  var e_u; stderr 0.15;
  var e_m; stderr 0.20;
end;

varobs ygr infl intr;

estimated_params;
// name,        distribution,  p1,    p2,    lower, upper
  kappa,        beta_pdf,      0.10,  0.05,  0.001, 0.95;
  phi_pi,       normal_pdf,    1.50,  0.25,  1.010, 4.00;
  phi_y,        normal_pdf,    0.25,  0.10,  0.000, 2.00;
  rho_r,        beta_pdf,      0.75,  0.10,  0.010, 0.99;
  rho_g,        beta_pdf,      0.85,  0.08,  0.010, 0.99;
  rho_u,        beta_pdf,      0.50,  0.15,  0.010, 0.99;
  stderr e_g,   inv_gamma_pdf, 0.30,  2.00;
  stderr e_u,   inv_gamma_pdf, 0.15,  2.00;
  stderr e_m,   inv_gamma_pdf, 0.20,  2.00;
end;

steady;
check;
stoch_simul(order=1, irf=0, periods=0);

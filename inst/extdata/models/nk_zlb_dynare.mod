// @dynhr-model
// name: NK model with ZLB (Dynare OccBin syntax)
// short: nk_zlb_dynare
// tier: 2
// endo: 3
// exo: 2
// params: 8
// has_lead: true
// has_measurement: false
// has_estimation: false
// has_analytical_ss: true
// features: obc, dynare_occbin_parity
// refs: Guerrieri-Iacoviello (2015); Giovannini-Pfeiffer-Ratto (2021)
// notes: Three-equation linearised NK model with ZLB on nominal rate.
//        Parameterisation identical to nk_2obc.mod (ZLB spec only).
//        Written in Dynare's occbin_constraints block syntax so Dynare 7.0
//        can solve it with its OccBin PKF filter.  Used as the parity
//        reference for kalman_filter_obc_pkf() in test-obc-pkf-dynare.R.
// @dynhr-model-end

// Three-equation linearised NK (log-deviations from SS).
// All variables are in deviation form; SS = 0 for all.
//
//   IS curve:      y_t = E[y_{t+1}] - sigma*(r_t - E[pi_{t+1}]) + eps_d
//   Phillips:      pi_t = beta*E[pi_{t+1}] + kappa*y_t
//   Taylor rule:   r_t = rho_r*r_{t-1} + (1-rho_r)*(phi_pi*pi_t + phi_y*y_t) + eps_m
//
// ZLB: r_t = max(r_lb, r_t^natural)   (r_lb < 0 in deviation form)

var y pi r;

varexo eps_d eps_m;

parameters sigma beta kappa rho_r phi_pi phi_y r_lb;

sigma  = 1.0;
beta   = 0.99;
kappa  = 0.1;
rho_r  = 0.7;
phi_pi = 1.5;
phi_y  = 0.5;
r_lb   = -0.01;   // ZLB: 1 ppt below SS (matching nk_2obc.mod)

model(linear);
// IS curve — no OBC
y = y(+1) - sigma*(r - pi(+1)) + eps_d;

// New Keynesian Phillips curve — no OBC
pi = beta*pi(+1) + kappa*y;

// Taylor rule — SLACK regime (standard)
[name='Taylor rule', relax='ZLB']
r = rho_r*r(-1) + (1 - rho_r)*(phi_pi*pi + phi_y*y) + eps_m;

// Taylor rule — BINDING regime (ZLB active, r pegged at r_lb)
[name='Taylor rule', bind='ZLB']
r = r_lb;

end;

occbin_constraints;
name 'ZLB'; bind r < r_lb; relax r > r_lb;
end;

initval;
y  = 0;
pi = 0;
r  = 0;
end;

shocks;
var eps_d;  stderr 0.01;
var eps_m;  stderr 0.005;
end;

steady;
check;
stoch_simul(order=1, irf=0, periods=0);

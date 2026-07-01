// @dynhr-model
// name: NK model with two OBCs (ZLB + output floor)
// short: nk_2obc
// tier: 2
// endo: 3
// exo: 2
// params: 7
// has_lead: true
// has_measurement: false
// has_estimation: true
// has_analytical_ss: true
// features: obc, multi_constraint
// refs: Guerrieri-Iacoviello (2015, OccBin); Phase O1 multi-constraint test fixture
// notes: Minimal linearised NK model with two simultaneously-bindable OBCs.
//        Constraint 1: ZLB on nominal rate r > r_lb (equation: Taylor rule).
//        Constraint 2: Output floor y > y_lb (equation: IS curve).
//        Both constraints can bind independently or simultaneously.
//        Designed to exercise the 2^k = 4 regime logic in Phase O1.
// @dynhr-model-end

// Three-equation linearised NK model (log-deviations from SS).
// All variables are in deviation form; SS = 0 for all.
//
//   IS curve:     y_t = E[y_{t+1}] - sigma*(r_t - E[pi_{t+1}]) + eps_d
//   Phillips:     pi_t = beta*E[pi_{t+1}] + kappa*y_t
//   Taylor rule:  r_t = rho_r*r_{t-1} + (1-rho_r)*(phi_pi*pi_t + phi_y*y_t) + eps_m
//
// OBC 1 (ZLB):     r > r_lb   (r_lb < 0 in deviation form, i.e. -r_ss)
// OBC 2 (y floor): y > y_lb   (y_lb < 0, e.g. -0.05 = 5 ppt below SS)

var y pi r;

varexo eps_d eps_m;

parameters sigma beta kappa rho_r phi_pi phi_y r_lb y_lb;

sigma  = 1.0;
beta   = 0.99;
kappa  = 0.1;
rho_r  = 0.7;
phi_pi = 1.5;
phi_y  = 0.5;
r_lb   = -0.01;   // ZLB: 1 ppt below SS rate (in deviation form)
y_lb   = -0.05;   // output floor: 5 ppt below SS

model(linear);
// IS curve — OBC on y (output floor)
[mcp = 'y > -0.05']  y = y(+1) - sigma * (r - pi(+1)) + eps_d;

// New Keynesian Phillips curve (no OBC)
pi = beta * pi(+1) + kappa * y;

// Taylor rule — OBC on r (zero lower bound)
[mcp = 'r > -0.01']  r = rho_r * r(-1) + (1 - rho_r) * (phi_pi * pi + phi_y * y) + eps_m;
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

estimated_params;
// name,      distribution,  p1 (mean/lo), p2 (std/hi), lower,  upper
rho_r,        beta_pdf,      0.70,          0.10,         0.10,   0.99;
phi_pi,       normal_pdf,    1.50,          0.20,         1.01,   4.00;
sigma,        normal_pdf,    1.00,          0.30,         0.10,   5.00;
end;

steady;
check;
stoch_simul(order=1, irf=0, periods=0);

// @dynhr-model
// name: Real Business Cycle (baseline)
// short: rbc
// tier: 1
// endo: 7
// exo: 1
// params: 7
// has_lead: true
// has_measurement: false
// has_estimation: false
// has_analytical_ss: false
// features: perturbation, kalman_filter, stoch_simul
// refs: King, Plosser, Rebelo (1988); standard textbook RBC
// notes: Cleanest smoke-test model. Single TFP shock, Cobb-Douglas production,
//        nonlinear Euler equation. Used for Stage C/D/E/F parity tests.
// @dynhr-model-end
// Reference model: King-Plosser-Rebelo RBC
// Endogenous: 7, Exogenous: 1, Parameters: 7
// Textbook real business cycle with Cobb-Douglas production and technology shock

var c y k l w r a;

varexo eps_a;

parameters beta alpha delta sigma eta rho_a sigma_a;

beta = 0.99;
alpha = 0.33;
delta = 0.025;
sigma = 1;
eta = 1;
rho_a = 0.95;
sigma_a = 0.01;

model;
% Euler equation (consumption-savings)
c^(-sigma) = beta * c(+1)^(-sigma) * (1 + r(+1) - delta);

% Production function
y = a * k(-1)^alpha * l^(1 - alpha);

% Capital law of motion
k = (1 - delta) * k(-1) + y - c;

% Labor supply (static)
w = eta * c^sigma * l;

% Marginal product of capital
r = alpha * a * k(-1)^(alpha - 1) * l^(1 - alpha);

% Marginal product of labor
w = (1 - alpha) * a * k(-1)^alpha * l^(-alpha);

% Technology process
log(a) = rho_a * log(a(-1)) + eps_a;
end;

// Initval values close to the analytical RBC steady state
// (computed offline for beta=0.99, alpha=0.33, delta=0.025, sigma=eta=1).
// Dynare 0.10's preprocessor does NOT evaluate expressions in initval —
// every value must be a literal numeric.
initval;
a = 1.0;
k = 24.7;
l = 0.928;
y = 2.78;
c = 2.16;
w = 2.00;
r = 0.0351;
end;

shocks;
var eps_a;
stderr 0.01;
end;

// @dynhr:expectations
// rho_a_hi:   type="param_range", variable="rho_a",   min=0.5,   max=0.999, description="Tech shock AR(1) should be persistent"
// sigma_a_lo: type="param_range", variable="sigma_a", min=0.001, max=0.1,   description="Tech shock SD should be small"
// @dynhr:end

steady;
check;
stoch_simul(order=1, irf=20, periods=0);

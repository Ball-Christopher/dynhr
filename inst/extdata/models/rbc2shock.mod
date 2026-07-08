// @dynhr-model
// name: Real Business Cycle (2-shock: TFP + discount-factor)
// short: rbc2shock
// tier: 1
// endo: 8
// exo: 2
// params: 9
// has_lead: true
// has_measurement: false
// has_estimation: false
// has_analytical_ss: false
// features: perturbation
// refs: King, Plosser, Rebelo (1988); extended with AR(1) preference shock
// notes: Two-shock RBC used for multi-shock order-3 perturbation parity tests.
//        TFP shock (eps_a) + AR(1) discount-factor shock (eps_b).
//        No auxiliary variables: all leads ≤ +1, all lags ≤ -1.
//        States: k, a, b (3 states).  Dynare golden: other_ignore/dynare_golden_rbc2shock/.
// @dynhr-model-end
// Reference model: King-Plosser-Rebelo RBC + AR(1) discount-factor shock
// Endogenous: 8, Exogenous: 2, Parameters: 9
// Real business cycle with Cobb-Douglas production, technology shock,
// and a stochastic preference (discount factor) shock.

var c y k l w r a b;

varexo eps_a eps_b;

parameters beta alpha delta sigma eta rho_a sigma_a rho_b sigma_b;

beta    = 0.99;
alpha   = 0.33;
delta   = 0.025;
sigma   = 1;
eta     = 1;
rho_a   = 0.95;
sigma_a = 0.01;
rho_b   = 0.90;
sigma_b = 0.01;

model;
% Euler equation (with discount-factor shock b multiplying beta)
c^(-sigma) = beta * exp(b) * c(+1)^(-sigma) * (1 + r(+1) - delta);

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

% Discount-factor shock (log-deviation from zero)
b = rho_b * b(-1) + eps_b;
end;

// Initval: same as rbc.mod for variables c,y,k,l,w,r,a; b=0 at steady state.
initval;
a = 1.0;
b = 0.0;
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
var eps_b;
stderr 0.01;
end;

steady;
check;
stoch_simul(order=3, irf=0, periods=0);

function [f_m, f_c, f_g, f_fv, f_fc, solver] = dyn_pid_init(Ts)
%[f_m, f_c, f_g, f_fv, f_fc, solver] = dyn_init()
% This function returns the dynamic model of irb120
% tau = f_m*ddq + f_c*dq + f_g + f_fv*dq + f_fc*sign(dq)
% f_m, f_c, f_g are casadi Functions. 
% f_fv and f_fc are diagonal matrices
% solver() is a casadi.integrator() for robot simulation
import casadi.*

% Sign = @(x) 2/(1+exp(-2000*x))-1;% smooth approximation of sign(X)
Sign = @(x) 0*x;

% Base parameters from identification experiments
fv_index = [2,11,21,31,41,51];
fc_index = fv_index + 1;
Theta = [0.8257, 14.5172, 12.9110, 0.4266, -0.0352, -0.0799, 0.2086, 0.6680, 2.2940, -0.0468, 17.7907, 9.9381, 0.2132, -0.0368, -0.0574, 0.2467, -0.0206, 0.4034, 1.0592, 0.2941, 6.7992, 5.2652, 0.1004, 0.0258, 0.0355, -0.0385, 0.0124, 0.0038, 0.0759, 0.0639, 0.7472, 1.0598, -0.1588, 0.0376, -0.0080, 0.1119, 0.0046, 0.0296, 0.0157, 0.0206, 1.2398, 1.5322, -0.0015, -0.0629, 0.0138, -0.0171, 0.0123, 0.0004, -0.0350, 0.0255, 0.5774, 1.2584].';
% Isn't the firction coefficient too large?

Theta(fv_index) = 0; % Negelet viscous friction for now 
% Theta(fc_index) = 0;
% Symbolic variables
q = SX.sym('q',6);
qd = SX.sym('qd',6);
qdd = SX.sym('qdd',6);
% States for integral action of PID controllers
nu = SX.sym('nu',6);
% Reference + feedforward
qref = SX.sym('qref',6);
qdref = SX.sym('qdref',6);
tau = SX.sym('tau',6);

% Evaluate regressor matrix (actual model)
Hb = regr(q,qd,qdd,Sign);
HbTh = Hb*Theta;

f_fv = diag(Theta(fv_index));
f_fc = diag(Theta(fc_index));

Fv = f_fv*qd;
Fc = f_fc*Sign(qd);
M = jacobian(HbTh, qdd);
C = jacobian(HbTh-Fv-Fc, qd)/2;
G = HbTh - Fv - Fc - M*qdd - C*qd;

% Functions for evaluation
f_m = casadi.Function('M', {q}, {M});
f_c = casadi.Function('C', {q, qd}, {C});
% G should be a function of q only. But it needs q, qd and qdd due to
% numerical error. The last two terms have no obvious contibution to the
% value of f_g hence they could be simply set to 0.
f_g = casadi.Function('G', {q, qd, qdd}, {G});

% PID parameters and model
% [121 121 101 50 51 50] : Transmission ratio between motors and joints
Kp = [10; 10; 10; 10; 10; 10];
Ki = [1*121^2; 1*121^2; 0.3*101^2; 0.1*50^2; 0.1*51^2; 0.1*50^2];
Kv = [0.1*121^2; 0.1*121^2; 0.03*101^2; 0.01*50^2; 0.01*51^2; 0.01*50^2];
pid_e1 = qref - q;
pid_e2 = Kp.*pid_e1 + qdref - qd;
pid_u = nu + Kv.*pid_e2 + tau;
% State vector: q, qdot, nu
x = vertcat(q,qd,nu);
% Algebraic states: qdotdot
z = qdd;
% Ordinary differential equations
rhs = vertcat(qd,qdd,Ki.*pid_e2);
% Algebraic equations
fz = pid_u - HbTh;

% Construct differential algebraic equations (DAEs)
dae.x = x;
dae.z = z;
dae.p = vertcat(tau,qref,qdref);
dae.ode = rhs;
dae.alg = fz;
opts.tf = Ts;
solver = integrator('integrator','idas',dae,opts);
% This DAE has the following form:
% ODE: d([q;qd;nu])/dt = [qd;qdd;Ki*pid_e2]
% AE: 0 = pid_u - phi(q, qd, qdd)*Theta
% It's written like this because casadi don't know dq/dt = qd, d(qd)/dt =
% qdd, etc. We have to write this relationship explictly.
end


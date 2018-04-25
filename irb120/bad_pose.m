% This is one of the bad pose
q_bad = [0.3933;0.3933;0.3933;0.3933;0.3933;0.3933];

% Solver
import casadi.*

% The coefficient of Sign can't be too large otherwise we get nan in C.
Sign = @(x) 2./(1+exp(-100*x))-1;% smooth approximation of sign(X)
% Sign = @(x) 0*x;

% Base parameters from identification experiments
fv_index = [2,11,21,31,41,51];
fc_index = fv_index + 1;
Theta = [0.8257, 14.5172, 12.9110, 0.4266, -0.0352, -0.0799, 0.2086, 0.6680, 2.2940, -0.0468, 17.7907, 9.9381, 0.2132, -0.0368, -0.0574, 0.2467, -0.0206, 0.4034, 1.0592, 0.2941, 6.7992, 5.2652, 0.1004, 0.0258, 0.0355, -0.0385, 0.0124, 0.0038, 0.0759, 0.0639, 0.7472, 1.0598, -0.1588, 0.0376, -0.0080, 0.1119, 0.0046, 0.0296, 0.0157, 0.0206, 1.2398, 1.5322, -0.0015, -0.0629, 0.0138, -0.0171, 0.0123, 0.0004, -0.0350, 0.0255, 0.5774, 1.2584].';
% Theta = [5.0796669673e-01 , 4.4568082460e+00 , 1.3787950000e+01 , 6.7745800133e-01 , 3.1212315804e-03 , 1.7898565256e-02 , 6.2402396070e-04 , 4.9529718097e-01 , 2.3508014800e+00 , 1.7878490000e-02 , 4.4568082460e+00 , 1.3787950000e+01 , 2.7761448193e-01 , 7.7343764000e-05 , 4.0935039840e-03 , 2.2414331697e-01 , 1.8299624120e-04 , 4.0197360000e-01 , 9.9362446000e-01 , 2.2513635000e-01 , 2.2581545660e+00 , 1.0389365000e+01 , 2.1507750870e-03 , -3.7504812860e-04 , 2.7741073260e-04 , 4.0212860522e-03 , 3.6945008786e-03 , -7.9412200000e-03 , 3.2607490000e-02 , 1.9425000000e-02 , 6.2243750000e-01 , 2.1890000000e+00 , 2.2019101148e-03 , -4.5413173000e-06 , -4.8526850500e-05 , 3.0755032090e-03 , 4.6311668500e-05 , 1.2837100000e-03 , 2.8373164000e-02 , 3.1836620000e-02 , 6.4758397500e-01 , 2.2327800000e+00 , 1.4461740000e-09 , 1.8080000000e-11 , 2.5788860000e-09 , -3.0631088000e-07 , 6.8001446626e-05 , 2.2600000000e-07 , -1.8080000000e-05 , 1.9425000000e-02 , 3.8679000000e-01 , 1.0698250000e+00]';
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
% M = SX.sym('M', 6,6);
% for i = 1:6
%     temp = HbTh(i);
%     M(i, :) = temp.gradient(qdd);
% end
% C = SX.sym('C', 6,6);
% for i = 1:6
%     temp = HbTh(i) - Fv(i) - Fc(i);
%     C(i, :) = temp.gradient(qd)/2;
% end
M = jacobian(HbTh, qdd);
C = jacobian(HbTh-Fv-Fc, qd)/2;
% M = HbTh.jacobian(qdd);
% tempC = HbTh - Fv - Fc;
% C = tempC.jacobian(qd)/2;
G = HbTh - Fv - Fc - M*qdd - C*qd;

% Functions for evaluation
f_m = casadi.Function('M', {q}, {M});
f_c = casadi.Function('C', {q, qd}, {C});
% G should be a function of q only. But it needs q, qd and qdd due to
% numerical error. The last two terms have no obvious contibution to the
% value of f_g hence they could be simply set to 0.
f_g = casadi.Function('G', {q, qd, qdd}, {G});
f_hb = casadi.Function('Hb', {q, qd, qdd}, {HbTh});

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
% Formulating the optimization problem using casadi.
% This script won't work without necessary variables

%% Optimal control problem
opti = casadi.Opti();
%decision variables
% a is the control input -> a is piecewise constant -> b is piecewise linear
% -> tau is piecewise nonlinear.
a = opti.variable(1,N0);
b = opti.variable(1,N0+1);
tau = opti.variable(6,N0);
%objective function
%   % reshape a and b for tau calculation
b_temp = (b(1:end-1)+b(2:end))/2;
aa = repmat(a, 6, 1);
bb = repmat(b_temp, 6, 1);
%   % m, c, g are calculated inside the while loop
m = opti.parameter(6, N0);
c = opti.parameter(6, N0);
g = opti.parameter(6, N0);
bk = ( (b(2:end)).^0.5 + (b(1:end-1)).^0.5 )/2 ;
bk_inv = 1./bk; % 2/(sqrt(bk)+sqrt(bk+1)), numerical approximation of 1/sqrt(b)
tau_sq = (tau./repmat(tau_max, 1, N0)).^2; % tau^2/tau_max^2
tau_sum = sum(tau_sq, 1); % summation over 6 joints
ds = opti.parameter(1,N0); % grid size
grid_control = opti.parameter(1,N0); % grid_control is either 0 or 1, for Filip's algorithm, not used
T = ds*bk_inv'; % Time spent
% E = ds*((bk_inv.*grid_control).*tau_sum)'; % Thermal energy
E = ds*(bk_inv.*tau_sum)'; % Thermal energy, integration over s
%   % set objective function
opti.minimize(E);
% joint torque & path info
% See 'Time-Optimal Control of Robotic Manipulators Along Specified Paths'
% for a clear explanation
opti.subject_to(tau == m.*aa + c.*bb + g); 
% dynamic constraints b'(s) = 2*a
for kk=1:N0
    opti.subject_to(b(kk+1)-b(kk)-2*a(kk)*ds(kk)==0);
end
% why should we impose torque constraint on all of them? We only apply the
% first one? 
% If we don't limit all the torques and if the path will be regenerated,
% we'll have some wild path that is ok in the first step. But if the path
% is regenerated again and again (which is very likely to happen), we will
% eventually get a path that is not feasible even at the beginning.
for i = 1:N0
    for j = 1:6
        opti.subject_to(-tau_max(j)<=tau(j,i)<=tau_max(j)); % torque is limited
    end
end
% % set initial value
b_init = opti.parameter(1,1);
T_remain = opti.parameter(1,1);
opti.subject_to(b(1) == b_init);
opti.subject_to(b>=0);
% Warm start
opti.set_initial(b, b_opt);
opti.set_initial(a, a_opt);
opti.set_initial(tau, tau_opt);
% Duration constraint
opti.subject_to(T <= T_remain);
% Choose the solver as ipopt
opti.solver('ipopt');

%% Optimal problem for si estimation
optis = casadi.Opti();
q_meas = optis.parameter(6,1); % measured position
dqdt_meas = optis.parameter(6,1); % measured velocity [rad/s]
b_optimized = optis.parameter(1, 2); % [b_opt(1), b_opt(2)] for interpolation
si_old = optis.parameter(1,1); % previous si
si_orig = optis.parameter(1,1); % s0, the first value of s sequence
grid_size = optis.parameter(1,1); % ds(1), grid size = (1-s0)/N
c0op = optis.parameter(6,1); % path coefficients
c1op = optis.parameter(6,1); % q(s)=c0s^4+...+c3s+c4
c2op = optis.parameter(6,1);
c3op = optis.parameter(6,1);
c4op = optis.parameter(6,1);
si_new = optis.variable(1,1); % decision variable, our new si0

b_sinew = (b_optimized(2)-b_optimized(1))/grid_size * (si_new-si_orig) + b_optimized(1);
dqds_meas = dqdt_meas/sqrt(b_sinew);
q_try = c0op*si_new^4+c1op*si_new^3+c2op*si_new^2+c3op*si_new+c4op; % q(si_new)
dq_try = 4*c0op*si_new*3+3*c1op*si_new^2+2*c2op*si_new+c3op; % dq(si_new)
ep_s = (q_try - q_meas)'*(q_try - q_meas); % position error, squared
ev_s = (dq_try - dqds_meas)'*(dq_try - dqds_meas); % velocity error, squared

optis.minimize(ep_s+0*ev_s); % If the robot won't move, tune the factor before ev_s
optis.subject_to(si_new <= si_old+0.1); % smaller searching range
optis.subject_to(si_new >= si_old-0.1);
optis.subject_to(0 <= si_new <= 1);
optis.solver('ipopt');
% From experience, if the robot wouldn't move, there was something wrong with
% finding si.

%% Optimization problem for new path generation
sz = 5; % 4th order polynomial
optip = casadi.Opti();
bp = optip.variable(sz,1);
sp = optip.parameter(sz, N0);
dsp = optip.parameter(sz-1, N0);
asp = optip.parameter(1, N0);
qsp_init = optip.parameter(1,1);
qsp_end = optip.parameter(1,1);
dqdsp_init = optip.parameter(1,1);

bsp = bp'*sp;
dqdsp = bp(2:end)'*dsp;
esp = asp - bsp;
ep = esp*esp';

optip.minimize(ep);
optip.subject_to(bsp(1) == qsp_init);
optip.subject_to(bsp(end) == qsp_end);
optip.subject_to(dqdsp(1) == dqdsp_init);
optip.subject_to(dqdsp(end) == 0);

optip.solver('ipopt');
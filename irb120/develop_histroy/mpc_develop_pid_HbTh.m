% This is worse than mcg approach. The invalid point is still there and
% other results are not good either. Even new singular points appears
% E.g. pi/20*ones(6,1)

import casadi.*

NMAX = 100;
N0 = 10;

fs = 250;
Ts = 1/fs;
% The interception is hard coded here.
q_itc = pi/10*ones(6,1);
T_itc = 2;
% Max gear torque
tg_max = [0.633; 0.552; 0.325; 0.146; 0.128; 0.201];
% transmission ratio
trans_ratio = [121; 121; 101; 50; 51; 50];
% Max joint torque
tau_max = tg_max.*trans_ratio;
%% System dynamics
% Sign = @(x) sign(x);
% % Sign = @(x) 2/(1+exp(-2000*x))-1;% smooth approximation of sign(X)
% fv_index = [2,11,21,31,41,51];
% fc_index = fv_index + 1;
% Theta = [0.8257, 14.5172, 12.9110, 0.4266, -0.0352, -0.0799, 0.2086, 0.6680, 2.2940, -0.0468, 17.7907, 9.9381, 0.2132, -0.0368, -0.0574, 0.2467, -0.0206, 0.4034, 1.0592, 0.2941, 6.7992, 5.2652, 0.1004, 0.0258, 0.0355, -0.0385, 0.0124, 0.0038, 0.0759, 0.0639, 0.7472, 1.0598, -0.1588, 0.0376, -0.0080, 0.1119, 0.0046, 0.0296, 0.0157, 0.0206, 1.2398, 1.5322, -0.0015, -0.0629, 0.0138, -0.0171, 0.0123, 0.0004, -0.0350, 0.0255, 0.5774, 1.2584].';
% Theta(fv_index) = 0;

[solver, Sign, f_hbth] = dyn_HbTh_init(Ts);

%% Initializing
x = zeros(12,NMAX+1); % actual states of the robot
u_act = zeros(6,NMAX); % control histroy
t = zeros(1,NMAX+1); % time
b_init_k = 1e-6;
dq_pre = 0.001*ones(6, 1); % dq/ds after previous execution
ddq_pre = 0.001*ones(6, 1); % dq^2/d^2s after previous execution
q_ref = zeros(6,NMAX+1);
dq_ref = zeros(6, NMAX+1);
a_opt = ones(1, N0);
b_opt = ones(1, N0+1);
tau_opt = zeros(6, N0);
err_init = abs(x(1:6,1) - q_itc);
err = zeros(1,NMAX);
% nr should be an odd number so that ceil(s) and floor(s) won't give the
% same value;
nr = 999;
k = 1;
counter = 1;
ER = 100;
flag_close = 0;
flag_precise = -1;
flag_init = 1;
flag_Nchange = 0;

%Test variables
imprecise_counter = -1;
ei = zeros(6,NMAX);
si_int = 0;
flag_integratorF = 0;

%% Fomulating optimization problems
% Define a path
N = N0;
      
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
aa = ones(6,1)*a;
bb = ones(6,1)*b_temp;

% read path info q(s), dq/ds, ddq/dss
q_op = opti.parameter(6, N0);
dq_op = opti.parameter(6, N0);
ddq_op = opti.parameter(6, N0);

bk = ( (b(2:end)).^0.5 + (b(1:end-1)).^0.5 )/2 ;
bk_inv = 1./bk; % Numerical approximation of 1/sqrt(b)

tau_sq = (tau./(tau_max*ones(1, N0))).^2;
tau_sum = sum(tau_sq, 1);


ds = opti.parameter(1,N0);
grid_control = opti.parameter(1,N0); % grid_control is either 0 or 1
T = ds*bk_inv'; % Time spent
% E = ds*((bk_inv.*grid_control).*tau_sum)'; % Thermal energy
E = ds*(bk_inv.*tau_sum)'; % Thermal energy
%   % set objective function
opti.minimize(E);

% joint torque & path info
% See 'Time-Optimal Control of Robotic Manipulators Along Specified Paths'
% for a clear explanation
for i = 1:N0
    HbTh = f_hbth(q_op(:, i),dq_op(:, i)*sqrt(b(i)),ddq_op(:, i)*b(i)+dq_op(:, i)*a(i));
    opti.subject_to(tau(:, i) == HbTh);
end
% dynamic constraints -
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

%% Optimization problem to find si
q_now = zeros(6,1);
% c1= dq_pre + 2*q_now - 2*q_itc;
% c2 = 3*q_itc - 3*q_now - 2*dq_pre;
% c3 = dq_pre;
% c4 = q_now;
% qs = @(s) c1*s.^3 + c2*s.^2 + c3*s + c4;
% dqs = @(s) 3*c1*s.^2 + 2*c2*s + c3;
% ddqs = @(s) 6*c1*s + 2*c2;
c1 = q_itc - q_now;
c2 = 100;
sx = casadi.SX.sym('sx', 1, 1);
geo_path = 2/(1+exp(-c2*(sx-1)));
qs = casadi.Function('qs', {sx}, {geo_path});
dqs = casadi.gradient(qs, sx);
ddqs = casadi.gradient(dqs, sx);
qs = c1*qs + q_now;
dqs = c1*dqs;
ddqs = c1*ddqs;

optis = casadi.Opti();
q_meas = optis.parameter(6,1);
si_old = optis.parameter(1,1);
si_new = optis.variable(1,1);
e_s = (qs(si_new) - q_meas)'*(qs(si_new) - q_meas);
optis.minimize(e_s);
optis.subject_to(si_new <= si_old+0.1);
optis.subject_to(si_new >= si_old-0.1);
optis.solver('ipopt');



%% Start running
qh = zeros(6, N0);
dqh = zeros(6, N0);
ddqh = zeros(6, N0);
temp_m = zeros(6, 6);
temp_c = zeros(6, 6);
temp_g = zeros(6, 1);
while k < 3 || norm(x(:, k) - x(:,k-1))>1e-6
    
    if counter > NMAX
        break
    end
    
    err(k) = norm(x(1:6,k) - q_itc);
    if (err(k)<1e-2) && (err(k)>err(k-1))
        break
    end
    
    %% Solve the optimization problem
    if flag_precise ~= 1
        % first step or imprecise execution, regenerate path and s
        s = linspace(0,1,N0+1);
        si = 0;
        q_now = x(1:6,k);
        
        c1= dq_pre + 2*q_now - 2*q_itc;
        c2 = 3*q_itc - 3*q_now - 2*dq_pre;
        c3 = dq_pre;
        c4 = q_now; 
        qs = @(s) c1*s.^3 + c2*s.^2 + c3*s + c4;
        dqs = @(s) 3*c1*s.^2 + 2*c2*s + c3;
        ddqs = @(s) 6*c1*s + 2*c2;

        flag_init = 0;
        flag_precise = 1;
        imprecise_counter = imprecise_counter + 1;
    else
        % keep tracking the same path
        s = linspace(si, 1, N0+1);
    end
    
    ds_k = s(2:end) - s(1:end-1);
    sh = s(1:end-1) + 0.5*ds_k;
    qh = qs(sh);
    dqh = dqs(sh);
    ddqh = ddqs(sh);
    %*******************m, c, g********************    
%     m_k = zeros(size(qh));
%     c_k = zeros(size(qh));
%     g_k = zeros(size(qh));
%     % Let's ignore friction for now
%     for i = 1:N0
%         temp_m = full( f_m(qh(:, i)) );
%         temp_c = full( f_c(qh(:, i), dqh(:, i)) );
%         temp_g = full( f_g(qh(:, i), [0;0;0;0;0;0], [0;0;0;0;0;0]) );
%         
%         m_k(:, i) = temp_m*dqh(:, i);
%         c_k(:, i) = temp_m*ddqh(:, i) + temp_c*dqh(:,i);
%         g_k(:, i) = temp_g + f_fc*sign(dqh(:, i));
%     end
    
    gcp = ones(1,N0);
%     if flag_Nchange == 1
%         gcp(1, N:end) = 0;
%     end
    
    opti.set_value(grid_control, gcp);
    opti.set_value(q_op, qh);
    opti.set_value(dq_op, dqh);
    opti.set_value(ddq_op, ddqh);
    opti.set_value(b_init, b_init_k);
    opti.set_value(T_remain, T_itc);
    opti.set_value(ds, ds_k);
    
    % solve NLP
    sol = opti.solve();   % actual solve
    
    counter = counter + 1;
    
    %% Post-processing
    % retrieving variables
    bk_inv_opt = sol.value(bk_inv);
    b_opt = sol.value(b);
    a_opt = sol.value(a);
    tau_opt = sol.value(tau);
    dt = ds_k.*bk_inv_opt;
    q_ref_k= qs(si+ds_k(1));
    dq_ref_k = dqs(si+ds_k(1))*sqrt(b_opt(2));
    
    %% update
    number_of_execution = round(dt(1)/Ts);
    
    if dt(1) < Ts
        number_of_execution = 1;
    end
    x_inter = zeros(12,number_of_execution+1);
    x_inter(:,1) = x(:,k);
    
    % tau is piecewise nonlinear. It should be calculated one by one
    % The calculation of tau is a feedback procedure
    
    % Integral state
    nu = zeros(6,1);
    slope_b = (b_opt(2) - b_opt(1))/ds_k(1);
    for i = 1:number_of_execution
        % q, dq, ddq
        q_inter = x_inter(1:6,i);
        dq_inter = dqs(si);
        ddq_inter = ddqs(si);
        b_inter = slope_b*( si - s(1) ) + b_opt(1);% b is piecewise linear
        % reference
        ds_inter = Ts*Ts*slope_b/4 + Ts*sqrt(b_inter);
        b_inter_next = slope_b*ds_inter + b_inter;
        q_ref_inter = qs(si+ds_inter);
        dq_ref_inter = dqs(si+ds_inter)*sqrt(b_inter_next);
%         % m c g
%         temp_m = full( f_m(q_inter) );
%         temp_c = full( f_c(q_inter, dq_inter) );
%         temp_g = full( f_g(q_inter, [0;0;0;0;0;0], [0;0;0;0;0;0]) );
%         m_inter = temp_m*dq_inter;
%         c_inter = temp_m*ddq_inter + temp_c*dq_inter;
%         g_inter = temp_g + f_fc*sign(dq_inter);
        % get tau
        u = f_hbth(q_inter,dq_inter*sqrt(b_inter),ddq_inter*b_inter+dq_inter*a_opt(1));
%         u = m_inter*a_opt(1) + c_inter*b_inter + g_inter;
        % Integrate
%         
%         dq_ref_inter = zeros(6,1);
%         
        out = solver('x0',[x_inter(:, i);nu],'p',[u;q_ref_inter;dq_ref_inter]);          
        temp_state = full(out.xf);
        x_inter(:, i+1) = temp_state(1:12);
        nu = temp_state(13:end);
        % update si
        optis.set_value(si_old, si);
        optis.set_value(q_meas, x_inter(1:6, i+1));
        sol_s = optis.solve();
        si = sol_s.value(si_new);
        
    end
    
    if flag_integratorF == 1
        for ii = 1:i-1
            x(:,k+ii) = x_inter(:, ii+1);
            t(k+ii) = t(k+ii-1) + Ts;
        end
        break
    end
    x(:,k+1) = x_inter(:,end);
    
    b_init_k = slope_b*(si - s(1)) + b_opt(1);
    dq_pre = x(7:12,k+1)/sqrt(b_init_k);
    
    t(k+1) = t(k) + number_of_execution*Ts;
    T_itc = T_itc - number_of_execution*Ts;
    if T_itc < 0
        warning('time out');
        break
    end
    q_ref(:,k+1) = q_ref_k;
    dq_ref(:, k+1) = dq_ref_k;
    
    ei(:, k) = q_ref(:, k) - x(1:6, k);
    % Let's disable the precision check for now
    q_now = x(1:6,k+1);
    %     if q_ref_k(1) > 1/k || q_ref_k(2) > 1/k
    %         flag_precise = 0;
    %     end

    if si >= 1-1/nr
        break
    end
    
    k = k + 1;
end

arm_plot;

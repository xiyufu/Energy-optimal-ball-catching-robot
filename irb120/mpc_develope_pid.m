
import casadi.*

NMAX = 100;
N0 = 15;

% Solving the optimization problem of a,b,tau is fast
% But solving the optimization problem of si for every sample period is too
% much of work. Might have to reduce fs (below the robot's fs)
fs = 250;
Ts = 1/fs;

% Joint Range
joint_range = [-165, 165;...
               -110, 110;...
               -110, 70;...
               -160, 160;...
               -120, 120;...
               -400, 400];
% q_itc = joint_range(:, 1) + (joint_range(:,2) - joint_range(:, 1)).*rand(6,1);
% q_itc = q_itc*pi/180;
% while max(abs(q_itc)) > pi % make sure q is within the range of [-pi, pi]
%     qindex = (abs(q_itc) > pi);
%     q_itc(qindex) = q_itc(qindex) - sign(q_itc(qindex))*pi;
% end
[test_traj, tc, idx] = ball_traj_gen([2;0;0], [-1;0;9], [0.512;0;0.63], 0.01, 3);
q_itc = end_configuration(test_traj(:,idx), irb120);
q_itc = q_itc';
T_itc = tc;
% Max gear torque
tg_max = [0.633; 0.552; 0.325; 0.146; 0.128; 0.201];
% transmission ratio
trans_ratio = [121; 121; 101; 50; 51; 50];
% Max joint torque
tau_max = tg_max.*trans_ratio;

%% System dynamics
[f_m, f_c, f_g, f_fv, f_fc, solver, Sign, f_hb] = dyn_pid_init(Ts);

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
tau_opt = ones(6, N0);
% err_init = abs(x(1:6,1) - q_itc);
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
s_sequence = zeros(1, NMAX);

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

%   % m, c, g are calculated inside the while loop
m = opti.parameter(6, N0);
c = opti.parameter(6, N0);
g = opti.parameter(6, N0);

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
opti.subject_to(tau == m.*aa + c.*bb + g); 
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
c1= dq_pre + 2*q_now - 2*q_itc;
c2 = 3*q_itc - 3*q_now - 2*dq_pre;
c3 = dq_pre;
c4 = q_now;
qs = @(s) c1*s.^3 + c2*s.^2 + c3*s + c4;
dqs = @(s) 3*c1*s.^2 + 2*c2*s + c3;
ddqs = @(s) 6*c1*s + 2*c2;

optis = casadi.Opti();
q_meas = optis.parameter(6,1);
dq_meas = optis.parameter(6,1);
si_old = optis.parameter(1,1);
si_new = optis.variable(1,1);
ep_s = (qs(si_new) - q_meas)'*(qs(si_new) - q_meas);
ev_s = (dqs(si_new) - dq_meas)'*(dqs(si_new) - dq_meas);
optis.minimize(ep_s+1e-6*ev_s);
optis.subject_to(si_new <= si_old+0.1);
optis.subject_to(si_new >= si_old-0.1);
optis.solver('ipopt');
% From experience, if the robot would move, there was something wrong with
% finding si.


%% Start running
qh = zeros(6, N0);
dqh = zeros(6, N0);
ddqh = zeros(6, N0);
temp_m = zeros(6, 6);
temp_c = zeros(6, 6);
temp_g = zeros(6, 1);
% nu = zeros(6,1);
% while k < 3 || norm(x(:, k) - x(:,k-1))>1e-6
while counter < NMAX
    
%     if counter > NMAX
%         break
%     end
    
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
    
    % %test
    s_sequence(k) = si;
    
    % %end test
    
    ds_k = s(2:end) - s(1:end-1);
    sh = s(1:end-1) + 0.5*ds_k;
    qh = qs(sh);
    dqh = dqs(sh);
    ddqh = ddqs(sh);
    %*******************m, c, g********************    
    m_k = zeros(size(qh));
    c_k = zeros(size(qh));
    g_k = zeros(size(qh));
    % Let's ignore friction for now
    for i = 1:N0
        temp_m = full( f_m(qh(:, i)) );
        temp_c = full( f_c(qh(:, i), dqh(:, i)) );
        temp_g = full( f_g(qh(:, i), [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        
        m_k(:, i) = temp_m*dqh(:, i);
        c_k(:, i) = temp_m*ddqh(:, i) + temp_c*dqh(:,i);
        g_k(:, i) = temp_g + f_fc*Sign(dqh(:, i));
    end
    
    gcp = ones(1,N0);
    if flag_Nchange == 1
        % We have to control the grid that matters
        gcp(1, N:end) = 0;
    end
    
    opti.set_value(grid_control, gcp);
    opti.set_value(m, m_k);
    opti.set_value(c, c_k);
    opti.set_value(g, g_k);
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
    for i = 1:number_of_execution
        slope_b = (b_opt(2) - b_opt(1))/ds_k(1);
        b_inter = slope_b*( si - s(1) ) + b_opt(1);% b is piecewise linear
        % q, dq, ddq
        q_inter = x_inter(1:6,i);
%         dq_inter = dqs(si);
        dq_inter = x_inter(7:12,i)/sqrt(b_inter);
        ddq_inter = ddqs(si);        
        % reference
        ds_inter = Ts*Ts*slope_b/4 + Ts*sqrt(b_inter);
        b_inter_next = slope_b*ds_inter + b_inter;
        q_ref_inter = qs(si+ds_inter);
        dq_ref_inter = dqs(si+ds_inter)*sqrt(b_inter_next);
        % m c g
        temp_m = full( f_m(q_inter) );
        temp_c = full( f_c(q_inter, dq_inter) );
        temp_g = full( f_g(q_inter, [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        m_inter = temp_m*dq_inter;
        c_inter = temp_m*ddq_inter + temp_c*dq_inter;
        g_inter = temp_g + f_fc*Sign(dq_inter);
        % get tau
        u = m_inter*a_opt(1) + c_inter*b_inter + g_inter;
        % Integrate
%         scale_big = blkdiag(scale,scale,eye(6));
%         out = solver('x0',scale_big\[x_inter(:, i);nu],'p',[u;q_ref_inter;dq_ref_inter]);          
%         temp_state = scale_big*full(out.xf);
        try
            out = solver('x0',[x_inter(:, i);nu],'p',[u;q_ref_inter;dq_ref_inter]);          
        catch
            q_bad = x_inter(1:6,i);
            warning('IDAS failed');
            flag_integratorF = 1;
            break;
        end
        temp_state = full(out.xf);
        x_inter(:, i+1) = temp_state(1:12);
        nu = temp_state(13:end);
        % update si
        optis.set_value(si_old, si);
        optis.set_value(q_meas, x_inter(1:6, i+1));
        optis.set_value(dq_meas, x_inter(7:12, i+1)/sqrt(b_inter_next));
        sol_s = optis.solve();
        si = sol_s.value(si_new);
        
    end
    
    if flag_integratorF == 1
%         for ii = 1:i-1
%             x(:,k+ii) = x_inter(:, ii+1);
%             t(k+ii) = t(k+ii-1) + Ts;
%         end
        break
    end
    x(:,k+1) = x_inter(:,end);
    
    b_init_k = (b_opt(2) - b_opt(1))/ds_k(1)*(si - s(1)) + b_opt(1);
    dq_pre_k = x(7:12,k+1)/sqrt(b_init_k);
    dq_pre = dq_pre_k;
    
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
    dq_now = x(7:12, k+1);
    %     if q_ref_k(1) > 1/k || q_ref_k(2) > 1/k
    %         flag_precise = 0;
    %     end

    if si >= 1-1/1000
        break
    end
    
    k = k + 1;
end

% Termination
dt_sum = 0;
u_sum = zeros(6,1);
n_sum = 0;
dt_total = sum(dt(2:end));
nu = zeros(6,1);
for i = 2:N
    dt_sum = dt_sum + dt(i);
    u_sum = u_sum + tau_opt(:,i);
    n_sum = n_sum + 1;
    if dt_sum >= Ts
        u_sum = u_sum/n_sum;
        q_ref_inter = (q_itc - q_now)*(i-1)/(N-1) + q_now; %note this q_now doesn't change
        dq_ref_inter = (0 - dq_now)*(i-1)/(N-1) + dq_now;
        out = solver('x0',[x(:, k);nu],'p',[u_sum;q_ref_inter;dq_ref_inter]);
        temp_state = full( out.xf );
        x(:,k+1) = temp_state(1:12);
        t(k+1) = t(k) + Ts;
        k = k+1;
        q_ref(:,k) = q_ref_inter;
        dq_ref(:, k) = dq_ref_inter;
        u_sum = 0;
        n_sum = 0;
        dt_sum = 0;
    end
end
%% Show the results
arm_plot;

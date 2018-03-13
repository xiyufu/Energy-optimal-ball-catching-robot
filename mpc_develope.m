% TODO:
% The maximum torque in manual is not enough for gravity compensation...

% Strange, the calcuation doesn't give the 4th joint enough torque (for
% acceleration) and it will fall back. The estimation of si is influenced
% by the 4th joint and goes back (si_old > si_new, because the 4th joint
% is going back) Maybe I should check the Jacobian at this point? A
% singular pose?
import casadi.*

NMAX = 500;
N0 = 20;

fs = 250;
Ts = 1/fs;
% The interception is hard coded here.
q_itc = (pi/3)*ones(6,1);
T_itc = 2;
tau_max = 50;
% [0.633104; 0.552033; 0.324992; 0.146561; 0.128488; 0.200659]; % Nm, from gear side

%% System dynamics
[f_m, f_c, f_g, f_fv, f_fc, solver] = dyn_init();

% % Forward kinematics
% fkine = @ (q) [0.3*cos(q(1))+0.3*cos(q(1)+q(2)), 0.3*sin(q(1))+0.3*sin(q(1)+q(2))];
% 
% %****************** M(q), C(q,dq), G(q)********************
% M = @(x) [4.1575+0.09*cos(x(2)), 1.0225+0.045*cos(x(2));
%     1.0225+0.045*cos(x(2)), 1.0225];
% 
% C = @(x, dx) [0.1, -0.045*(dx(2) + 2*dx(1))*sin(x(2));
%     0.045*dx(1)*sin(x(2)), 0.1];
% 
% G = @(x) [5.886*cos(x(1))+1.4715*cos(x(1)+x(2));
%     1.4715*cos(x(1)+x(2))];
% 
% Mi = @(x) pinv(M(x));
% 
% % x_dot = f(x,u)
% zero22 = zeros(2,2);
% zero21 = zeros(2,1);
% eye22 = eye(2,2);
% 
% f = @(x, u) [zero22, eye22;
%     zero22, -Mi(x(1:2))*C(x(1:2),x(3:4))]*x + [zero22; Mi(x(1:2))]*u + [zero21; -Mi(x(1:2))*G(x(1:2))];

%% Initializing
x = zeros(12,NMAX+1); % actual states of the robot
u_act = zeros(6,NMAX); % control histroy
t = zeros(1,NMAX+1); % time
b_init_k = 0.001;
dq_pre = 0.001*ones(6, 1); % dq/ds after previous execution
ddq_pre = 0.001*ones(6, 1); % dq^2/d^2s after previous execution
q_ref = zeros(6,NMAX+1);
a_opt = ones(1, N0);
b_opt = ones(1, N0+1);
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
ei = zeros(2,NMAX);
si_int = 0;

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
bb = [b_temp;b_temp;b_temp;b_temp;b_temp;b_temp];
aa = [a;a;a;a;a;a];

%   % m, c, g are calculated inside the while loop
m = opti.parameter(6, N0);
c = opti.parameter(6, N0);
g = opti.parameter(6, N0);

bk = ( (b(2:end)).^0.5 + (b(1:end-1)).^0.5 )/2 ;
bk_inv = 1./bk; % Numerical approximation of 1/sqrt(b)

tau_sq = tau.^2;
tau_sum = sum(tau_sq, 1)/(6*(tau_max'*tau_max));

ds = opti.parameter(1,N0);
grid_control = opti.parameter(1,N0); % grid_control is either 0 or 1
T = ds*bk_inv'; % Time spent
% E = ds*((bk_inv.*grid_control).*tau_sum)'; % Thermal energy
E = ds*(bk_inv.*tau_sum)'; % Thermal energy
%   % set objective function
opti.minimize(E);

opti.subject_to(tau == m.*aa + c.*bb + g); % joint torque & path info
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
    opti.subject_to(-tau_max<=tau(:,i)<=tau_max); % torque is limited
end


b_init = opti.parameter(1,1);
T_remain = opti.parameter(1,1);

opti.subject_to(b(1) == b_init);
opti.subject_to(b>=0);

% Warm start
opti.set_initial(b, b_opt);
opti.set_initial(a, a_opt);

% Duration constraint
opti.subject_to(T <= T_remain);

% Choose the solver as ipopt
opti.solver('ipopt');

%% Optimization problem to find si
optis = casadi.Opti();
q_meas = optis.parameter(6,1);
spre = optis.parameter(1,1);
snext = optis.parameter(1,1);
q_spre = optis.parameter(6,1);
q_snext = optis.parameter(6,1);
si_mid = optis.parameter(1,1);
q_smid = optis.parameter(6,1);

si_new = optis.variable(1,1);
qsi_new = (q_snext-q_spre)*(si_new-si_mid)/(snext-spre) + q_smid;

e_s = (qsi_new-q_meas)'*(qsi_new-q_meas);
optis.minimize(e_s);
optis.subject_to(si_new <= snext);
optis.subject_to(si_new >= spre);
optis.solver('ipopt');



%% Start running
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
        
        qs = zeros(6, nr+1);
        dqs = zeros(6, nr+1);
        ddqs = zeros(6, nr+1);
        for traj_iter = 1:6
            [qs(traj_iter, :), dqs(traj_iter, :), ddqs(traj_iter, :), ~, ~] = polytraj(q_itc(traj_iter), nr, q_now(traj_iter), dq_pre(traj_iter), ddq_pre(traj_iter), 0);
        end
%         [q1s, dq1s, ddq1s, ~, ~] = polytraj(q_itc(1), nr, q_now(1), dq_pre(1), ddq_pre(1), 0);
%         [q2s, dq2s, ddq2s, ~, ~] = polytraj(q_itc(2), nr, q_now(2), dq_pre(2), ddq_pre(2), 0);
%         [q3s, dq3s, ddq3s, ~, ~] = polytraj(q_itc(3), nr, q_now(3), dq_pre(3), ddq_pre(3), 0);
%         [q4s, dq4s, ddq4s, ~, ~] = polytraj(q_itc(4), nr, q_now(4), dq_pre(4), ddq_pre(4), 0);
%         [q5s, dq5s, ddq5s, ~, ~] = polytraj(q_itc(5), nr, q_now(5), dq_pre(5), ddq_pre(5), 0);
%         [q6s, dq6s, ddq6s, ~, ~] = polytraj(q_itc(6), nr, q_now(6), dq_pre(6), ddq_pre(6), 0);
        
        flag_init = 0;
        flag_precise = 1;
        imprecise_counter = imprecise_counter + 1;
    else
        % keep tracking the same path
        s = linspace(si, 1, N0+1);
    end
    
    ds_k = s(2:end) - s(1:end-1);
    sh = s(1:end-1) + 0.5*ds_k;
    sh_int = sh*nr;
    % Find q(sh) by interpolation
    sh_low = floor(sh_int)+1;
    sh_high = ceil(sh_int)+1;
    
    % sh_high - sh_low == 1, so we can remove the denominator
    qh = zeros(6, N0);
    dqh = zeros(6, N0);
    ddqh = zeros(6, N0);
    for mid_iter = 1:6
        qh(mid_iter, :) = (qs(mid_iter, sh_high)-qs(mid_iter, sh_low)).*(sh_int-sh_low)+qs(mid_iter, sh_low);
        dqh(mid_iter, :) = (dqs(mid_iter, sh_high)-dqs(mid_iter, sh_low)).*(sh_int-sh_low)+dqs(mid_iter, sh_low);
        ddqh(mid_iter, :) = (ddqs(mid_iter, sh_high)-ddqs(mid_iter, sh_low)).*(sh_int-sh_low)+ddqs(mid_iter, sh_low);
    end
%     qh1 = (q1s(sh_high)-q1s(sh_low)).*(sh_int-sh_low)+q1s(sh_low);
%     qh2 = (q2s(sh_high)-q2s(sh_low)).*(sh_int-sh_low)+q2s(sh_low);
%     dqh1 = (dq1s(sh_high)-dq1s(sh_low)).*(sh_int-sh_low)+dq1s(sh_low);
%     dqh2 = (dq2s(sh_high)-dq2s(sh_low)).*(sh_int-sh_low)+dq2s(sh_low);
%     ddqh1 = (ddq1s(sh_high)-ddq1s(sh_low)).*(sh_int-sh_low)+ddq1s(sh_low);
%     ddqh2 = (ddq2s(sh_high)-ddq2s(sh_low)).*(sh_int-sh_low)+ddq2s(sh_low);
% 
%     qh = [qh1;qh2];
%     dq = [dqh1;dqh2];
%     ddq = [ddqh1;ddqh2];
    
    %*******************m, c, g*********************
    % memory allocation
    temp_m = zeros(6, 6);
    temp_c = zeros(6, 6);
    temp_g = zeros(6, 1);
    
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
        g_k(:, i) = temp_g + f_fc*sign(dqh(:, i));
    end
    
    gcp = ones(1,N0);
    if flag_Nchange == 1
        % We have to control the grid that matters
        gcp(1, N:end) = 0;
    end
    
%     % % tset
%     
%     if k < 6
%         mcg_tracker(:, k) = [m_k(:,1);c_k(:,1);g_k(:,1)];
%         qh_tracker(:,k) = qh(:, 1);
%         dq_tracker(:, k) = dq(:, 1);
%     else
%         k;
%     end
%     % %
    opti.set_value(grid_control, gcp);
    opti.set_value(m, m_k);
    opti.set_value(c, c_k);
    opti.set_value(g, g_k);
    opti.set_value(b_init, b_init_k);
    opti.set_value(T_remain, T_itc);
    opti.set_value(ds, ds_k);
%     opti.set_value(a_warm, a_opt);
%     opti.set_value(b_warm, b_opt);
    
    % solve NLP
    sol = opti.solve();   % actual solve
    
    counter = counter + 1;
    
    %% Post-processing
    % retrieving variables
    bk_inv_opt = sol.value(bk_inv);
    b_opt = sol.value(b);
    a_opt = sol.value(a);
    
    dt = ds_k.*bk_inv_opt;
    
    
    %% update
    number_of_execution = floor(dt(1)/Ts);
    
    if dt(1) < Ts
        number_of_execution = 1;
%         break
%         N = N - 1;
%         if N == 0
%             break
%         end
%         flag_Nchange = 1;
%         continue;
    end
    x_inter = zeros(12,number_of_execution+1);
    x_inter(:,1) = x(:,k);
    
    % tau is piecewise nonlinear. It should be calculated one by one
    % The calculation of tau is a feedback procedure
    for i = 1:number_of_execution
        % q, dq, ddq
        q_inter = x_inter(1:6,i);
        
        b_inter = ( (b_opt(2) - b_opt(1))/ds_k(1) )*( si - s(1) ) + b_opt(1);% b is piecewise linear
        % linear interpolation
        si_scaled = si*nr;
        si_low = floor(si_scaled) + 1;
        si_high = ceil(si_scaled) + 1;
        si_scaled = si_scaled + 1;
        if si_low == si_high
            dq_inter = dqs(:, si_low);
            ddq_inter = ddqs(:, si_low);
        else
            dq_inter = (dqs(:, si_high)-dqs(:, si_low))*(si_scaled-si_low)/(si_high-si_low)+dqs(:, si_low);
%            dq_inter2 = (dq2s(si_high)-dq2s(si_low))*(si_scaled-si_low)/(si_high-si_low)+dq2s(si_low);
            ddq_inter = (ddqs(:, si_high)-ddqs(:, si_low))*(si_scaled-si_low)/(si_high-si_low)+ddqs(:, si_low);
%             ddq_inter2 = (ddq2s(si_high)-ddq2s(si_low))*(si_scaled-si_low)/(si_high-si_low)+ddq2s(si_low);
%             dq_inter = [dq_inter1;dq_inter2];
%             ddq_inter = [ddq_inter1;ddq_inter2];
        end
        % m c g
        temp_m = full( f_m(q_inter) );
        temp_c = full( f_c(q_inter, dq_inter) );
        temp_g = full( f_g(q_inter, [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        m_inter = temp_m*dq_inter;
        c_inter = temp_m*ddq_inter + temp_c*dq_inter;
        g_inter = temp_g;
        % get tau
        u = m_inter*a_opt(1) + c_inter*b_inter + g_inter;
        
        out = solver('x0',x_inter(:, i),'p',u);
        x_inter(:, i+1) = full( out.xf );
%         k1 = f(x_inter(:, i), u);
%         k2 = f(x_inter(:,i)+Ts*k1/2, u);
%         k3 = f(x_inter(:,i)+Ts*k2/2, u);
%         k4 = f(x_inter(:,i)+Ts*k3,   u);
%         x_inter(:,i+1) = x_inter(:,i) + Ts*(k1+2*k2+2*k3+k4)/6;
        
        % update si
        % First step, search si in [si_low-50, si_low+50], get the index (a
        % integer) that is closest to si. It's called si_int
        search_low = max(1, si_low-100);
        search_high = min(nr+1, si_low+100);
        e_inter = x_inter(1:6, i+1) - qs(:, search_low:search_high);
        e_inter_sq = zeros(1,length(e_inter));
        for e_inter_i = 1:length(e_inter)
            e_inter_sq(e_inter_i) = e_inter(:,e_inter_i)'*e_inter(:,e_inter_i);
        end
        si_int_test = si_int;
        [~, si_int] = min(e_inter_sq);
        si_int = si_int + search_low - 1;
        
        if si_int_test > si_int
            warning('si reversed');
        end
        % Second step, linearize q(s) around si_int ([si_int-1, si_int+1])
        % q(si) = (q(si_int+1)-q(si_int-1))/2 * (si - si_int) + q(si_int)
        % Solve min{|q(si)-q_now|^2}, s.t. si_int-1<si<si_int+1
        si_int_low = max(1,si_int-1);
        si_int_high = min(nr+1, si_int+1);
        optis.set_value(spre, si_int_low);
        optis.set_value(snext, si_int_high);
        optis.set_value(q_meas, x_inter(1:6, i+1));
        optis.set_value(q_spre, qs(:, si_int_low));
        optis.set_value(q_snext, qs(:, si_int_high));
        optis.set_value(si_mid, si_int);
        optis.set_value(q_smid, qs(:, si_int));
        
        sols = optis.solve();
        
        si = sols.value(si_new);
        si = (si-1)/nr;

    end
    
    x(:,k+1) = x_inter(:,end);
    
    % Let's disable the precision check for now
    q_ref_k= (qs(:, si_int_high)-qs(:, si_int_low))*(si*nr+1-si_int)/(si_int_high-si_int_low) + qs(:, si_int);
    q_now = x(1:6,k+1);
%     if q_ref_k(1) > 1/k || q_ref_k(2) > 1/k
%         flag_precise = 0;
%     end
    
    b_init_k = (b_opt(2) - b_opt(1))/ds_k(1)*(si - s(1)) + b_opt(1);
    dq_pre_k = x(7:12,k+1)/sqrt(b_init_k);
    dq_pre = dq_pre_k;
%     ddq_pre1 = (ddq1s(si_int_high)-ddq1s(si_int_low))*(si*nr+1-si_int)/(si_int_high-si_int_low)+ddq1s(si_int);
%     ddq_pre2 = (ddq2s(si_int_high)-ddq2s(si_int_low))*(si*nr+1-si_int)/(si_int_high-si_int_low)+ddq2s(si_int);
    ddq_pre = (ddqs(:, si_int_high)-ddqs(:, si_int_low))*(si*nr+1-si_int)/(si_int_high-si_int_low)+ddqs(:, si_int);
    
    t(k+1) = t(k) + number_of_execution*Ts;
    T_itc = T_itc - number_of_execution*Ts;
    if T_itc < 0
        warning('time out');
        break
    end
    
    q_ref(:,k+1) = q_ref_k;
%     ei(:, k) = q_ref(:, k) - x(1:6, k);
    k=k+1;
end

figure; hold
sc1 = scatter(t(1:k), x(1,1:k));
sc1.Marker = 'x';
sc2 = scatter(t(1:k), x(2,1:k));
sc2.Marker = 'x';
sc3 = scatter(t(1:k), x(3,1:k));
sc3.Marker = 'x';
sc4 = scatter(t(1:k), x(4,1:k));
sc4.Marker = 'x';
sc5 = scatter(t(1:k), x(5,1:k));
sc5.Marker = 'x';
sc6 = scatter(t(1:k), x(6,1:k));
sc6.Marker = 'x';

plot(t(1:k), q_ref(:,1:k));
plot(t(1:k), q_itc(1)*ones(1,k),t(1:k), q_itc(2)*ones(1,k));


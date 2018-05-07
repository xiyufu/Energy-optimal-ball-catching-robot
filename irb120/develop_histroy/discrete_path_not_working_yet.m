% TODO:
% The maximum torque in manual is not enough for gravity compensation...

import casadi.*

NMAX = 100;
N0 = 10;

fs = 250;
Ts = 1/fs;
% The interception is hard coded here.
q_itc = pi/3*ones(6,1);
T_itc = 5;
tau_max = 50;
% [0.633104; 0.552033; 0.324992; 0.146561; 0.128488; 0.200659]; % Nm, from gear side

%% System dynamics
[f_m, f_c, f_g, f_fv, f_fc, solver] = dyn_init(Ts);
Sign = @(x) 2/(1+exp(-2000*x))-1;

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
ei = zeros(6,NMAX);
si_int = 0;
s_sequence = zeros(1, NMAX);
mcg_tracker = zeros(18, NMAX);
u_tracker = zeros(6, NMAX);
a_tracker = zeros(1, NMAX);
b_tracker = zeros(1, NMAX);
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
        
        qs = zeros(6, nr+1);
        dqs = zeros(6, nr+1);
        ddqs = zeros(6, nr+1);
        for traj_iter = 1:6
            [qs(traj_iter, :), dqs(traj_iter, :), ddqs(traj_iter, :), ~, ~] = polytraj(q_itc(traj_iter), nr, q_now(traj_iter), dq_pre(traj_iter), ddq_pre(traj_iter), 0);
        end      
        flag_init = 0;
        flag_precise = 1;
        imprecise_counter = imprecise_counter + 1;
    else
        % keep tracking the same path
        s = linspace(si, 1, N0+1);
    end
    
    ds_k = s(2:end) - s(1:end-1);
    sh = s(1:end-1) + 0.5*ds_k;
    sh_scaled = sh*nr;
    % Find q(sh) by interpolation
    sh_low = floor(sh_scaled);
    sh_high = ceil(sh_scaled);
    % sh_high - sh_low == 1 (because nr = 999, si_scaled~=int), so we can remove the denominator
    qh = (qs(:, sh_high)-qs(:, sh_low)).*(sh_scaled-sh_low)+qs(:, sh_low);
    dqh = (dqs(:, sh_high)-dqs(:, sh_low)).*(sh_scaled-sh_low)+dqs(:, sh_low);
    ddqh = (ddqs(:, sh_high)-ddqs(:, sh_low)).*(sh_scaled-sh_low)+ddqs(:, sh_low);
    
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
        g_k(:, i) = temp_g + f_fc*sign(dqh(:, i));
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
    
    % retrieving variables
    bk_inv_opt = sol.value(bk_inv);
    b_opt = sol.value(b);
    a_opt = sol.value(a);
    dt = ds_k.*bk_inv_opt;
    
    
    %% update
    number_of_execution = round(dt(1)/Ts);
    
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
        % linear interpolation of q(s), dq(s) and ddq(s)
        si_scaled = si*nr;
        si_low = floor(si_scaled) + 1;
        si_high = ceil(si_scaled) + 1;
        si_scaled = si_scaled + 1;
        if si_low == si_high % si_scaled is an integer
%             q_inter = dqs(:, si_low);
            dq_inter = dqs(:, si_low);
            ddq_inter = ddqs(:, si_low);
        else
%             q_inter = (qs(:, si_high)-qs(:, si_low))*(si_scaled-si_low)/(si_high-si_low)+qs(:, si_low);
            dq_inter = (dqs(:, si_high)-dqs(:, si_low))*(si_scaled-si_low)/(si_high-si_low)+dqs(:, si_low);
            ddq_inter = (ddqs(:, si_high)-ddqs(:, si_low))*(si_scaled-si_low)/(si_high-si_low)+ddqs(:, si_low);
        end
        % m c g
        temp_m = full( f_m(q_inter) );
        temp_c = full( f_c(q_inter, dq_inter) );
        temp_g = full( f_g(q_inter, [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        m_inter = temp_m*dq_inter;
        c_inter = temp_m*ddq_inter + temp_c*dq_inter;
        g_inter = temp_g + f_fc*sign(dq_inter);
        % get tau
        u = m_inter*a_opt(1) + c_inter*b_inter + g_inter;
        % Integrate
        out = solver('x0',x_inter(:, i),'p',u);
        x_inter(:, i+1) = full( out.xf );
        
        % update si
        % First step, search si in [si_low-50, si_low+50], get the index (a
        % integer) that is closest to si. It's called si_int
        search_low = max(1, si_low-20);
        search_high = min(nr+1, si_low+20);
        e_inter = x_inter(1:6, i+1) - qs(:, search_low:search_high);
        e_inter_sq = zeros(1,length(e_inter));
        for e_inter_i = 1:length(e_inter)
            e_inter_sq(e_inter_i) = e_inter(:,e_inter_i)'*e_inter(:,e_inter_i);
        end
        si_int_test = si_int;
        [~, si_int] = min(e_inter_sq);
        si_int = si_int + search_low - 1;
        
        % test
        mcg_tracker(:, k) = [m_inter;c_inter;g_inter];
        a_tracker(k) = a_opt(1);
        b_tracker(k) = b_opt(1);
        if i == 1
            u_tracker(:, k) = u;
        end
        if si_int_test > si_int
            warning('si reversed');
        end
        % end test %
        
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
    
    s_sequence(k) = si_int_test;
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
    ddq_pre = (ddqs(:, si_int_high)-ddqs(:, si_int_low))*(si*nr+1-si_int)/(si_int_high-si_int_low)+ddqs(:, si_int);
    
    t(k+1) = t(k) + number_of_execution*Ts;
    T_itc = T_itc - number_of_execution*Ts;
    if T_itc < 0
        warning('time out');
        break
    end
    
    q_ref(:,k+1) = q_ref_k;
    ei(:, k) = q_ref(:, k) - x(1:6, k);
    k=k+1;
    
    if si >= 1-1/nr
        break
    end
end

arm_plot;


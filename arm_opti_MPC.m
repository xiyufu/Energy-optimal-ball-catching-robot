% Simulation of a 2 dof robot arm using MPC

NN = 500;

fs = 250;
Ts = 1/fs;
% The interception is hard coded here.
q_itc = [0.1260; 2.1567];
T_itc = 1.1497;
tau_max = 100; % Nm

%% System dynamics
% Forward kinematics
fkine = @ (q) [0.3*cos(q(1))+0.3*cos(q(1)+q(2)), 0.3*sin(q(1))+0.3*sin(q(1)+q(2))];

%****************** M(q), C(q,dq), G(q)********************
M = @(x) [4.1575+0.09*cos(x(2)), 1.0225+0.045*cos(x(2));
    1.0225+0.045*cos(x(2)), 1.0225];

C = @(x, dx) [0.1, -0.045*(dx(2) + 2*dx(1))*sin(x(2));
    0.045*dx(1)*sin(x(2)), 0.1];

G = @(x) [5.886*cos(x(1))+1.4715*cos(x(1)+x(2));
    1.4715*cos(x(1)+x(2))];

Mi = @(x) pinv(M(x));

% x_dot = f(x,u)
zero22 = zeros(2,2);
zero21 = zeros(2,1);
eye22 = eye(2,2);

f = @(x, u) [zero22, eye22;
    zero22, -Mi(x(1:2))*C(x(1:2),x(3:4))]*x + [zero22; Mi(x(1:2))]*u + [zero21; -Mi(x(1:2))*G(x(1:2))];

%% Initializing
x = zeros(4,NN+1); % actual states of the robot
u_act = zeros(2,NN); % control histroy
t = zeros(1,NN+1); % time
b_init = 0.001; % initial value of b, calculated in previous step£¬ set to 0.001 at first
dq_pre = [0.1;0.1]; % dq/ds after previous execution
ddq_pre = 0.1; % dq^2/d^2s after previous execution
q_ref = zeros(2,NN+1);
err_init = abs(x(1:2,1) - q_itc);

N0 = 20;
k = 1;
counter = 1;
ER = 100;
flag_close = 0;
flag_precise = 1;
flag_init = 1;
si = 0;

%Test variables
imprecise_counter = 0;
q_near = zeros(NN,N0+1); %path generated at each step
dt_near = zeros(NN,N0); % time spent on each step
different_si = zeros(2,NN+1);

sd_sim = zeros(NN,1);
ds_sim = zeros(NN,1);
%% Start running
% while norm(x(1:2,k)-q_itc)>norm(err_init)/ER
while k < 3 || norm(x(:, k) - x(:,k-1))>1e-6

    if counter > NN
        break
    end
 
    if norm(x(1:2,k)-q_itc)<norm(err_init)/20
        flag_close = flag_close + 1;
    end
    N = N0; % horizon

    %% Solve the optimization problem
    if flag_init == 1 || flag_precise == 0 
        % first step or imprecise execution, regenerate path and s
        s = linspace(0,1,N+1);
        si = 0;
        q_now = x(1:2,k);
        
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
        s = linspace(si, 1, N+1);
    end
    
    ds = s(2:end) - s(1:end-1);
    sh = s(1:end-1) + 0.5*ds;
    
    qh = qs(sh);
    dq = dqs(sh);
    ddq = ddqs(sh);
    
    %*******************m, c, g*********************
    m = zeros(size(qh));
    for i = 1:N
        m(:, i) = M(qh(:,i))*dq(:, i);
    end
    
    c = zeros(size(qh));
    for i = 1:N
        c(:, i) = M(qh(:,i))*ddq(:, i) + C(qh(:,i), dq(:,i))*dq(:,i);
    end
    
    g = zeros(size(qh));
    for i = 1:N
        g(:, i) = G(qh(:,i));
    end
    
    %***************formulate the problem****************
    opti = casadi.Opti();
    
    %decision variables
    % a is the control input -> a is piecewise constant -> b is piecewise linear
    % -> tau is piecewise nonlinear.
    a = opti.variable(1,N);
    b = opti.variable(1,N+1);
    
    %objective function
    bb = [(b(1:end-1)+b(2:end))/2;(b(2:end)+b(1:end-1))/2];
    aa = [a;a];
    tau = m.*aa + c.*bb + g; % joint torque
    
    bk = ( (b(2:end)).^0.5 + (b(1:end-1)).^0.5 )/2 ;
    bk_inv = 1./bk;
    
    tau_sq = tau.^2;
    tau_sum = (tau_sq(1,:) + tau_sq(2,:))/(2*tau_max^2);
    
    T = ds*bk_inv';
    E = ds*(bk_inv.*tau_sum)';
    
    opti.minimize(E);
    
    % dynamic constraints -
    for kk=1:N
        opti.subject_to(b(kk+1)-b(kk)-2*a(kk)*ds(kk)==0);
    end
    
    for i = 1:N
        opti.subject_to(-tau_max<=tau(:,i)<=tau_max); % torque is limited
    end
       
    
    opti.subject_to(b(1) == b_init);
    
    opti.subject_to(b>=0);
    opti.set_initial(b, b_init);
    
    % Duration constraint
    opti.subject_to(T <= T_itc);
    
    % solve NLP
    opti.solver('ipopt'); % set numerical backend
    sol = opti.solve();   % actual solve
    
    counter = counter + 1;

    
    %% Post-processing
    % retrieving variables
    bk_inv_opt = sol.value(bk_inv);
    b_opt = sol.value(b);
    a_opt = sol.value(a);
    
    dt = ds.*bk_inv_opt;
    
    sd_sim(k) =mean(sqrt(b_opt(1:2)));
    ds_sim(k) = ds(1);
    
    
    %% update
    number_of_execution = floor(dt(1)/Ts); 
    
    x_inter = zeros(4,number_of_execution+1);
    x_inter(:,1) = x(:,k);

        % It's meaningless to have a time interval smaller than the sample
    % period. We have to decrease N0
    if min(dt) < Ts
        N0 = N0 - 1;
        if N0 == 0
            break
        end
        
        continue
    end
    
    % tau is piecewise nonlinear. It should be calculated one by one
    % The calculation of tau is a feedback procedure
    
    for i = 1:number_of_execution
        % q, dq, ddq
        q_inter = x_inter(1:2,i);

        b_inter = ( (b_opt(2) - b_opt(1))/ds(1) )*( si - s(1) ) + b_opt(1);% b is piecewise linear
        % dq_inter = x_inter(3:4,i)/sqrt(b_inter); 
        dq_inter = dqs(si);
        ddq_inter = ddqs(si); 
        % m c g
        m_inter = M(q_inter)*dq_inter;
        c_inter = M(q_inter)*ddq_inter + C(q_inter, dq_inter)*dq_inter;
        g_inter = G(q_inter);
        % get tau
        u = m_inter*a_opt(1) + c_inter*b_inter + g_inter;
        
        % Matlab solver should be used here instead of doing it manually
        k1 = f(x_inter(:, i), u);
        k2 = f(x_inter(:,i)+Ts*k1/2, u);
        k3 = f(x_inter(:,i)+Ts*k2/2, u);
        k4 = f(x_inter(:,i)+Ts*k3,   u);
        x_inter(:,i+1) = x_inter(:,i) + Ts*(k1+2*k2+2*k3+k4)/6;
        
        % update si
        q_meas = x_inter(1:2,i+1);
        optis = casadi.Opti();
        si_new = optis.variable(1,1);
        e_s = (qs(si_new)-q_meas)'*(qs(si_new)-q_meas);
        optis.minimize(e_s);
        optis.subject_to(0 <= si_new <= 1);
        optis.solver('ipopt');
        sol_s = optis.solve();
        si = sol_s.value(si_new);
        
    end
    

    different_si(1,k) = si;
    
    x(:,k+1) = x_inter(:,end);
    
    q_now = x(1:2,k+1);
    if norm(qs(si)-q_now) > 5e-3
%          flag_precise = 0;
    end
    % store trajectory at every step
    temp = qs(s);
    q_near(k,1:length(s)) = temp(2,:);
    dt_near(k,1:length(dt)) = dt;
    
    b_init = (b_opt(2) - b_opt(1))/ds(1)*(si - s(1)) + b_opt(1);
    dq_pre = x(3:4,k+1)/sqrt(b_init);

    t(k+1) = t(k) + number_of_execution*Ts;
    T_itc = T_itc - number_of_execution*Ts;
    if T_itc < 0
        warning('time out');
        break
    end
    
    q_ref(:,k+1) = qs(s(2));
    k=k+1;
end

figure; hold
sc1 = scatter(t(1:k), x(1,1:k));
sc1.Marker = 'x';
sc2 = scatter(t(1:k), x(2,1:k));
sc2.Marker = 'x';

plot(t(1:k), q_ref(:,1:k));
plot(t(1:k), q_itc(1)*ones(1,k),t(1:k), q_itc(2)*ones(1,k));
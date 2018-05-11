% Energy optimal ball catching robot
% main function
eopt_init; % Initialize
opti_formulation; % Formulate optimization problems in casadi
[f_m, f_c, f_g, f_fv, f_fc, solver, Sign, f_hb] = dyn_pid_init(Ts); % dynamic model
q_itc_track = zeros(6, NMAX); % check how q_itc changes
k = 1;

% % %%%%%%%%%%%%%%%%%%%%
% % Joint Range
% joint_range = [-165, 165;...
%                -110, 110;...
%                -110, 70;...
%                -160, 160;...
%                -120, 120;...
%                -400, 400];
% q_itc = joint_range(:, 1) + (joint_range(:,2) - joint_range(:, 1)).*rand(6,1);
% q_itc = q_itc*pi/180;
% while max(abs(q_itc)) > pi % make sure q is within the range of [-pi, pi]
%     qindex = (abs(q_itc) > pi);
%     q_itc(qindex) = q_itc(qindex) - sign(q_itc(qindex))*pi;
% end
% T_itc = 0.5 + 0.1*rand(1,1);
% % %%%%%%%%%%%%%%%%%%%%%
q_itc = [-1.5610;1.5871;-1.4412;1.8197;0.1606;0.6442];
T_itc = 0.54;
while k < NMAX
    
%     % %%%%%%%%%%%%%%%%%%%%
%     if k == 2
%         error('feasible');
%     end
%     %   %%%%%%%%%%%%%%%%%%%
    q_now = x(1:6, k);
    dq_now = x(7:12, k);
    % ball trajectory estimation
%     sigma = 0.05/k;
%     % ball_init_pos = [2;0;0];% + sigma*rand(3,1);
%     % ball_init_vel = [-1;0;9];% + sigma*rand(3,1);
%     % ball_init_vel = [-3.15;0;2.86];
%     ee_init_pos = [0.512;0;0.63];
%     [test_traj, tc, idx] = ball_traj_gen(ball_init_pos, ball_init_vel, ee_init_pos, 0.005, 3);
%     
%     % catch-point determination & inverse kinematics
%     q_itc = end_configuration(test_traj(:,idx), irb120);
%     if isempty(q_itc)
%         error('The ball will not go inside the workspace');
%     end
%     q_itc = q_itc';
%     q_itc_track(:,k) = q_itc;
%     T_itc = tc;

    if k == 1
        q_itc_old = q_itc;
    end
    if norm(q_itc_old-q_itc) >= 0.05
        flag_precise = 0;
    else
        q_itc = q_itc_old;
    end
    
    % geometric path generation, 4th order polynomial
    if flag_init == 1 % 3rd order polynomial at the beginning
        c0 = zeros(6,1);
        c1= dq_pre + 2*q_now - 2*q_itc;
        c2 = 3*q_itc - 3*q_now - 2*dq_pre;
        c3 = dq_pre;
        c4 = q_now;
        si = 0;
        qs = @(s) c1*s.^3 + c2*s.^2 + c3*s + c4;
        dqs = @(s) 3*c1*s.^2 + 2*c2*s + c3;
        ddqs = @(s) 6*c1*s + 2*c2;
        
        flag_init = 0;
    elseif flag_precise == 0
        tempc = zeros(5, 6);
        for i = 1:6
            pathic = [c4(i);c3(i);c2(i);c1(i);c0(i)];
            % pathic, dqds_initpath, qs_initpath, qs_endpath, N0, si, optip
            %tempc(:, i) = new_path(pathic, x(6+i,k)/sqrt(b_init_k),q_now(i),q_itc(i),N0, si, optip);
            dqds_initpath = x(6+i,k)/sqrt(b_init_k);
            qs_initpath = q_now(i);
            qs_endpath = q_itc(i);
            new_path;
            tempc(:,i) = b_coefficients;
        end
        c4 = tempc(1, :)';
        c3 = tempc(2, :)';
        c2 = tempc(3, :)';
        c1 = tempc(4, :)';
        c0 = tempc(5, :)';
        qs = @(s) c0*s.^4 + c1*s.^3 + c2*s.^2 + c3*s + c4;
        dqs = @(s) 4*c0*s.^3 + 3*c1*s.^2 + 2*c2*s + c3;
        ddqs = @(s) 12*c0*s.^2 + 6*c1*s + 2*c2; 
        
        flag_precise = 1;
        imprecise_counter = imprecise_counter + 1;
    end
    
    % get the collocation points
    s = linspace(si, 1, N0+1);    
    ds_k = s(2:end) - s(1:end-1);
    sh = s(1:end-1) + 0.5*ds_k;
    qh = qs(sh);
    dqh = dqs(sh);
    ddqh = ddqs(sh);
    
    % Evaluate m,c,g at current s 
    m_k = zeros(size(qh));
    c_k = zeros(size(qh));
    g_k = zeros(size(qh));
    for i = 1:N0
        temp_m = full( f_m(qh(:, i)) );
        temp_c = full( f_c(qh(:, i), dqh(:, i)) );
        temp_g = full( f_g(qh(:, i), [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        
        m_k(:, i) = temp_m*dqh(:, i);
        c_k(:, i) = temp_m*ddqh(:, i) + temp_c*dqh(:,i);
        g_k(:, i) = temp_g + f_fc*Sign(dqh(:, i));
    end
    % feed data to the solver
    gcp = ones(1,N0);
    if flag_Nchange == 1 % for Filip's algorithm, not used
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
    sol = opti.solve();   % actual solve
    
    % get the optimal solution
    bk_inv_opt = sol.value(bk_inv);
    b_opt = sol.value(b);
    a_opt = sol.value(a);
    tau_opt = sol.value(tau);
    dt = ds_k.*bk_inv_opt;
    q_ref_k= qs(si+ds_k(1));
    dq_ref_k = dqs(si+ds_k(1))*sqrt(b_opt(2));
    
    % %%%%%%%test
    if k == 1
        dt_ref = zeros(1,NMAX);
        m_init = m_k(:, 1);
        c_init = c_k(:, 1);
        g_init = g_k(:, 1);
        a_opt_init = a_opt;
        b_opt_init = b_opt;
        tau_opt_init = tau_opt;
        E_init = sol.value(E);
        t_begining = sol.value(T);
    end
    dt_ref(k) = dt(1);
    % %%%%%%%%%end test
    
    % Send data to robot (simulation)
    number_of_execution = floor(dt(1)/Ts);
    if dt(1) < Ts
        number_of_execution = 1;
    end
    x_inter = zeros(12,number_of_execution+1); % states during internal loop
    x_inter(:,1) = x(:,k); % set initial values
    optis.set_value(b_optimized, [b_opt(1),b_opt(2)]); % Unchanged inside the loop ...
    optis.set_value(c0op, c0); % so we pass these values to optis now
    optis.set_value(c1op, c1);
    optis.set_value(c2op, c2);
    optis.set_value(c3op, c3);
    optis.set_value(c4op, c4);
    optis.set_value(grid_size, ds_k(1));
    optis.set_value(si_orig, s(1));
    nu = zeros(6,1); % state of integrator in PID control
    slope_b = (b_opt(2) - b_opt(1))/ds_k(1);
    for i = 1:number_of_execution
        b_inter = slope_b*( si - s(1) ) + b_opt(1);% b is piecewise linear
        % q, dqdt, ddqdtt
        q_inter = x_inter(1:6,i);
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
        try % simulation
            out = solver('x0',[x_inter(:, i);nu],'p',[u;q_ref_inter;dq_ref_inter]);          
        catch
            q_bad = x_inter(1:6,i);
            warning('IDAS failed');
            flag_integratorF = 1;
            break;
        end
        temp_state = full(out.xf); % update states
        x_inter(:, i+1) = temp_state(1:12);
        nu = temp_state(13:end);
        % update si
        optis.set_value(si_old, si);
        optis.set_value(q_meas, x_inter(1:6, i+1));
%         optis.set_value(dq_meas, x_inter(7:12, i+1)/sqrt(b_inter_next));
        optis.set_value(dqdt_meas, x_inter(7:12, i+1));
%         optis.set_value(b_optimized, [b_opt(1),b_opt(2)]);
%         optis.set_value(c0op, c0);
%         optis.set_value(c1op, c1);
%         optis.set_value(c2op, c2);
%         optis.set_value(c3op, c3);
%         optis.set_value(c4op, c4);
        optis.set_initial(si_new, si);
        sol_s = optis.solve();
        si = sol_s.value(si_new);
    
    end % end of execution
    
    if flag_integratorF == 1
        break % if IDAS failed, break
    end
    
    % update
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
%     q_ref(:, k+1) = q_ref_inter;
%     dq_ref(:, k+1) = dq_ref_inter;
    k = k + 1;
    
    % Finished?
    if si >= 1-1/50 || dt(1) < Ts
        break
    end    
    
    q_itc_old = q_itc;
end

%Termination
kt = k;
dt_sum = 0;
zero_vec = zeros(6,1);
u_sum = zeros(6,1);
n_sum = 0;
dt_total = sum(dt(2:end));
nu = zeros(6,1);
for i = 2:N0
    dt_sum = dt_sum + dt(i);
    u_sum = u_sum + tau_opt(:,i);
    n_sum = n_sum + 1;
    if dt_sum >= Ts || i == N0
        u_sum = u_sum/n_sum;
        q_ref_inter = (q_itc_old - q_now)*(i-1)/(N0-1) + q_now; %note this q_now doesn't change
        dq_ref_inter = (zero_vec - dq_now)*((i-1)/(N0-1))^4 + dq_now;
        out = solver('x0',[x(:, k);nu],'p',[u_sum;q_ref_inter;dq_ref_inter]);
        temp_state = full( out.xf );
        x(:,k+1) = temp_state(1:12);
        t(k+1) = t(k) + Ts;
        k = k+1;
        q_ref(:,k) = q_ref_inter;
        dq_ref(:, k) =  dq_ref_inter;
        u_sum = 0;
        n_sum = 0;
        dt_sum = 0;
    end
    
end

% plot results
arm_plot;
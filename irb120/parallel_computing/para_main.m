% % Parallel evaluating NLP and Linear predictor
% % The syntax of the following code is correct. 
% % However, CasADi doesn't support parfor for now.
% % Running this code will cause parallel pool shut down
% % due to the not supported SWIG serialization.

eopt_init;
opti_formulation;

[f_m, f_c, f_g, ~, f_fc, solver, Sign, ~] = dyn_pid_init(Ts);

% % %%%%%%%%%%%%%%%%%%%%
% %This commented section randomly generates interception point
% %in the workspace of IRB 120.
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
% T_itc = 2 + 0.5*rand(1,1);
% % %%%%%%%%%%%%%%%%%%%%%

% % %%%%%%%%%%%%%%%%%%%%%
% One almost feasible interception point that is infeasible
q_itc = [-1.5610;1.5871;-1.4412;1.8197;0.1606;0.6442];
% T_itc = 0.54;
T_itc = 2; % This is feasible
% % %%%%%%%%%%%%%%%%%%%%%

q_now = zeros(6,1);
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

% Group variables for functions running parallel 
% and try once in sequence
coefficients = [c0,c1,c2,c3,c4];
path = {qs, dqs, ddqs};
fs = {f_m, f_c, f_g, f_fc, Sign, solver};
opti_vari = {opti, a, b, tau, bk_inv };
opti_param = {grid_control, m, c, g, tau_holding, b_init, T_remain, ds};
[a_opt, b_opt, dt,s] = NLP_solver(si, T_itc, q_itc, b_init_k, N0, path, fs, opti_vari, opti_param);

optis_param = {q_meas, dqdt_meas, b_optimized, si_old, si_orig, grid_size, si_new};
path_param = {c0op, c1op, c2op, c3op, c4op};
[lp_result, ~] = linear_predictor(a_opt, b_opt, dt, si, s, coefficients, fs, optis, optis_param, path_param, path);

% Prepare for parfor
arg1 = {si, T_itc, q_itc, b_init_k, N0, path, fs, opti_vari, opti_param};
arg2 = {a_opt, b_opt, dt, si, s, coefficients, fs, optis, optis_param, path_param, path};
arguements = {arg1; arg2};
solutions = cell(1,2);
func_para = {@NLP_solver; @linear_predictor};

% Parallel computing.
% Better to start the parallel pool manually first since initializing it
% takes some time.
% Sadly this doesn't work, CasADi doesn't support parfor
parfor iter = 1:2
    solutions{iter} = func_para{iter}(arguements{iter, 1}{:});
end


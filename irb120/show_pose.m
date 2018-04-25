%% Kinematic model
import casadi.*
startup_rvc;

L(1) = Link('revolute', 'alpha', -pi/2, 'a', 0, 'd', 0.29, 'offset', 0);
L(2) = Link('revolute', 'alpha', 0, 'a', 0.27, 'd', 0, 'offset', -pi/2);
L(3) = Link('revolute', 'alpha', -pi/2, 'a', 0.07, 'd', 0, 'offset', 0);
L(4) = Link('revolute', 'alpha', pi/2, 'a', 0, 'd', 0.302, 'offset', 0);
L(5) = Link('revolute', 'alpha', -pi/2, 'a', 0, 'd', 0, 'offset', 0);
L(6) = Link('revolute', 'alpha', 0, 'a', 0, 'd', 0.21, 'offset', pi);

irb120 = SerialLink(L, 'name', 'irb120');

[f_m, f_c, f_g, f_fv, f_fc, solver, ~] = dyn_pid_init(0.04);

% Joint Range
joint_range = [-165, 165;...
               -110, 110;...
               -110, 70;...
               -160, 160;...
               -120, 120;...
               -400, 400];
           
% Max speed
vel_max = [250; 250; 250; 320; 320; 420];

joint_range = joint_range*pi/180;
vel_max = vel_max*pi/180;

% Max gear torque
tg_max = [0.633; 0.552; 0.325; 0.146; 0.128; 0.201];

% transmission ratio
trans_ratio = [121; 121; 101; 50; 51; 50];

% Max joint torque
tau_max = tg_max.*trans_ratio;

%% randomly check the condition number
N = 10000;
q = zeros(6, N);
m_cond = zeros(1, N);
for i = 1:N
    q(:, i) = joint_range(:, 1) + (joint_range(:,2) - joint_range(:, 1)).*rand(6,1);
    temp_m = full( f_m(q(:, i)) );
    m_cond(i) = cond(temp_m);
end
index = (m_cond > 2500);
q_ill = q(:, index);
for i = 1:length(q_ill)
    irb120.plot(q_ill(:,i)');
end

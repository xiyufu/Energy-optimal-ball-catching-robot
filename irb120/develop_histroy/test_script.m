import casadi.*
% %% check the integrator
% clear
% 
% Ts = 0.004;
% 
% q_now = zeros(6,1);
% dq_now = zeros(6,1);
% q_itc = [pi/6;pi/5;pi/4;pi/3;pi/2;pi];
% 
% c1 = -3*q_itc + 3*q_now + 2*dq_now;
% c2 = 4*q_itc - 4*q_now - 3*dq_now;
% c3 = zeros(6,1);
% c4 = dq_now;
% c5 = q_now;
% qs = @(s) c1*s.^4 + c2*s.^3 + c3*s.^2 + c4*s + c5;
% dqs = @(s) 4*c1*s.^3 + 3*c2*s.^2 + 2*c3*s + c4;
% ddqs = @(s) 12*c1*s.^2 + 6*c2*s + 2*c3;
% % c1= - 2*pi/4*ones(6,1);
% % c2 = 3*pi/4*ones(6,1);
% % qs = @(s) c1*s.^3 + c2*s.^2 ;
% % dqs = @(s) 3*c1*s.^2 + 2*c2*s;
% % ddqs = @(s) 6*c1*s + 2*c2;
% 
% t = linspace(0,1, 1/Ts);
% 
% [f_m, f_c, f_g, f_fv, f_fc, solver] = dyn_pid_init(Ts);
% 
% pos = qs(t);
% vel = dqs(t);
% 
% x = zeros(12, length(pos));
% nu = zeros(6, length(pos));
% u = zeros(6,1);
% for i = 1:length(pos)-1
%    out = solver('x0', [x(:, i); nu(:, i)], 'p', [u;pos(:, i); vel(:, i)]);
%     temp_state = full( out.xf );
%     x(:, i+1) = temp_state(1:12);
%     nu(:, i+1) = temp_state(13:end);
% end
% 
% % If the Sign function is set to 0,
% % We will encounter a 'IDA_TOO_MUCH_WORK' at the middle of our path. According
% % to the document, it means we have a 'singular point'. What is a singular
% % ponit? The one at which the Jacobian of the robot is singular?
% 
% % If the Sign function is set to 2/(1+exp(-2000*x))-1,
% % We get 'IDA_CONV_FAIL'. Casadi says that 'Calculating Jacobian failed'.
% 
% %This is true for mpc_develope_pid.m too. So the mpc itself might be ok.
% %This problem lies in the path?
% 
% % 
% % It's not likely that this Jacobian refers to the jacobian of DAE because
% % it is always singular...
% 
% figure; hold
% for i = 1:6
%     plot(t, pos(i, :));
%     plot(t, x(i, :),'x');
% end
% figure; hold
% for i = 7:12
%     plot(t, x(i, :), 'x');
%     plot(t, vel(i-6, :));
% end

%% Check f_m, f_c, f_g, f_fv
% %  % Problem: f_m*ddq+f_c*dq_f_g+f_fv*Sign(dq) != f_hbth
% %  % Solved: Sign should use ./ instead of /
% q_rand = rand(6,1);
% dq_rand = rand(6,1);
% ddq_rand = rand(6,1);
% 
% tmcg = full( f_m(q_rand)*ddq_rand + f_c(q_rand, dq_rand)*dq_rand + ...
%              f_g(q_rand,dq_rand,ddq_rand) + f_fc*Sign(dq_rand));
% thbth = full( f_hb(q_rand, dq_rand, ddq_rand) );

%% Check the details of f_c, f_m
% %  % Problem: Get NaN from f_c for some (q, dq). 
% %  % Solved: The Sign function was too steep, we have inf by using exp(.)
% %               in automatic differentiation. Change for smaller factor (100).
% [f_m, f_c, f_g, f_fv, f_fc, solver, Sign, f_hb] = dyn_pid_init(0.004);
% opts = struct('main', true);
% f_c.generate('fc.c', opts);
% % fc.c has a main function and could be compiled directly. But we need to
% % link the math library. So add -lm at the end of gcc command (order of -lm is
% % important, it's a bug of gcc that haven't be changed for many years)

%% Check the successful rate of the program
Nr = 100;
nf = 0;
for ir = 1:Nr
    try
        mpc_develope_pid;
    catch
        nf = nf + 1;
        continue;
    end
end
srate = (Nr-nf)/Nr;
%  % successful rate: 99%
%  % bad point: q_bad = [-0.0087;-0.0512;-0.0889;-0.0905;-0.2057;-0.4161]
%  % condition number of m at bad point: 124.4359
%  % eigen value: [2.4255;1.4355;0.4272;0.0346;0.0235;0.0195]
%  % If we take a look at the robot, we find that at this point, the robot
%  almost has two colinear revolute joints (joint 4 and 6), a really
%  singular point...
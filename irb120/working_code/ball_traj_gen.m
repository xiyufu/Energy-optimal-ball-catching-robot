function [traj_ball, tc, idx] = ball_traj_gen(init_pos, init_vel, ee_pose, Ts, total_time)
%[traj_ball, tc] = ball_traj_gen(init_pos, init_vel, Ts, total_time)
% This function takes the initial position and initial velocity of a ball
% (a tennis, specifically) and simulates its trajectory. Air drag is taken into account
% dv/dt = g - alpha*|v|*v
% dx/dt = v
import casadi.*

% Robot workspace range
range = [0.17^2; 0.58^2];

g = [0; 0; -9.81];
Cw = 0.45; % drag coefficient for a sphere
A = 0.002; %[m^2] section area of the ball (tennis)
rho = 1.293; % [kg/m^3] air density
m = 0.05; % [kg] mass of the ball (tennis)
alpha = Cw*A*rho/(2*m);

r = SX.sym('r', 3, 1);
dr = SX.sym('v', 3, 1);
ddr = SX.sym('a', 3, 1);

% state variable 
x = [r; dr];
% algebraic variable
z = ddr;
% ODE
rhs = [dr; ddr];
% AE
fz =  ddr - g + alpha*sqrt(dr'*dr)*dr;

dae.x = x;
dae.z = z;
dae.ode = rhs;
dae.alg = fz;
opts.tf = Ts;
solver = integrator('integrator','idas',dae,opts);

% generate trajectory
itc_index = [0; 0];
tc = 0;
x0 = [init_pos; init_vel];
N = total_time/Ts + 1;
traj_ball = zeros(6,N);
traj_ball(:,1) = x0;
for i = 1:N-1
    out = solver('x0',traj_ball(:,i));
    temp = full(out.xf);
    traj_ball(:,i+1) = temp;
    r = temp(1)^2 + temp(2)^2 + (temp(3)-0.29)^2;
    % r = 580 mm. The center locates at the 2nd joint.
    if itc_index(1) == 0
        if r < range(2)
            itc_index(1) = i;
            tc = temp(4)*(ee_pose(1)-temp(1))+temp(5)*(ee_pose(2)-temp(2))+temp(6)*(ee_pose(3)-temp(3));
            tc = tc/(temp(4)^2 + temp(5)^2 + temp(6)^2);
            break
        end
%     elseif intersection(2) == 0
%         if r < range(1)
%             intersection(2) = i;
%             break;
%         end
    end    
end

tc = tc + i*Ts;
idx = floor(tc/Ts);
for j = i:idx-1
    out = solver('x0',traj_ball(:,j));
    temp = full(out.xf);
    traj_ball(:,j+1) = temp;
end

end

%% test
% [test_traj, tc, idx] = ball_traj_gen([2;0;0], [-1;0;9], [0.512;0;0.63], 0.01, 3);
% plot3(test_traj(1,1:175),test_traj(2,1:175),test_traj(3,1:175))

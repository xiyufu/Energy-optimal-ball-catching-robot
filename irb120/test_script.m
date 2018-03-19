%% check the integrator
clear

Ts = 0.004;

q_now = zeros(6,1);
dq_now = zeros(6,1);
q_itc = [pi/6;pi/5;pi/4;pi/3;pi/2;pi];

c1 = -3*q_itc + 3*q_now + 2*dq_now;
c2 = 4*q_itc - 4*q_now - 3*dq_now;
c3 = zeros(6,1);
c4 = dq_now;
c5 = q_now;
qs = @(s) c1*s.^4 + c2*s.^3 + c3*s.^2 + c4*s + c5;
dqs = @(s) 4*c1*s.^3 + 3*c2*s.^2 + 2*c3*s + c4;
ddqs = @(s) 12*c1*s.^2 + 6*c2*s + 2*c3;
% c1= - 2*pi/4*ones(6,1);
% c2 = 3*pi/4*ones(6,1);
% qs = @(s) c1*s.^3 + c2*s.^2 ;
% dqs = @(s) 3*c1*s.^2 + 2*c2*s;
% ddqs = @(s) 6*c1*s + 2*c2;

t = linspace(0,1, 1/Ts);

[f_m, f_c, f_g, f_fv, f_fc, solver] = dyn_pid_init(Ts);

pos = qs(t);
vel = dqs(t);

x = zeros(12, length(pos));
nu = zeros(6, length(pos));
u = zeros(6,1);
for i = 1:length(pos)-1
   out = solver('x0', [x(:, i); nu(:, i)], 'p', [u;pos(:, i); vel(:, i)]);
    temp_state = full( out.xf );
    x(:, i+1) = temp_state(1:12);
    nu(:, i+1) = temp_state(13:end);
end

% We encountered a 'IDA_TOO_MUCH_WORK' at the middle of our path. According
% to the document, it means we have a 'singular point'. What is a singular
% ponit? The one at which the Jacobian of the robot is singular?
% 
% It's not likely that this Jacobian refers to the jacobian of DAE because
% it is always singular...

figure; hold
for i = 1:6
    plot(t, pos(i, :));
    plot(t, x(i, :),'x');
end
figure; hold
for i = 7:12
    plot(t, x(i, :), 'x');
    plot(t, vel(i-6, :));
end

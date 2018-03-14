%% check the integrator
clear

Ts = 0.004;

c1= - 2*pi/4*ones(6,1);
c2 = 3*pi/4*ones(6,1);
qs = @(s) c1*s.^3 + c2*s.^2 ;
dqs = @(s) 3*c1*s.^2 + 2*c2*s;
ddqs = @(s) 6*c1*s + 2*c2;

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

figure; hold
for i = 1:6
    plot(t, x(i, :),'--');
end
plot(t, pos(1,:));
figure; hold
for i = 7:12
    plot(t, x(i, :), '--');
end
plot(t, vel(1, :));
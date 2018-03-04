function [pos, vel, acc, px] = polytraj_posi(pos_end, vel_end, nr, pos_init, vel_init)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tm = 1;
a0 = pos_init;
a1 = vel_init;

optip = casadi.Opti();
a2 = optip.variable(1,1);
a3 = optip.variable(1,1);
a4 = optip.variable(1,1);
t = optip.variable(1,1);

optip.minimize(a4);
optip.subject_to(a4 + a3 + a2 + a1 + a0 == pos_end);
optip.subject_to(4*a4 + 3*a3 + 2*a2 + a1 == vel_end);
% If we have pos_end < pos_init && vel_init > 0, this problem has no
% solution. e.g. [pos,~,~,~] = polytraj_posi(10, 0, 1000, 15, 1);
% TOO BAD!
% TODO: If undershoot is inevitable, try to minimize the undershoot
if pos_end > pos_init
    optip.subject_to(4*a4*t^3 + 3*a3*t^2 + 2*a2*t + a1 >= 0);
else
    optip.subject_to(4*a4*t^3 + 3*a3*t^2 + 2*a2*t + a1 <= 0);
end
optip.subject_to(0<=t<=1);

optip.solver('ipopt');
sol = optip.solve();

px = [sol.value(a4) sol.value(a3) sol.value(a2) a1 a0];

t1 = linspace(0, tm, nr+1);

pos = polyval(px,t1);

vx = [4 3 2 1 0].*px;
vel = polyval(vx(1:end-1),t1);

ax = [4*3 3*2 2*1 0 0].*px;
acc = polyval(ax(1:end-2),t1);
end


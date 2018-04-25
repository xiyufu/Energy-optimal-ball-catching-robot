function [b_coefficients] = new_path(a, dqds_init, qs_init, qs_end, N, s0)
%[b] = new_path(a, dqds_init, qs_init, qs_end, N, s0)
%   takes current position and velocity, end position, current path
%   returns the closest new path 
%   a = [an; an-1; an-2; ... a0]
%   q_a(s) = a0*s^(n-1) + a1*s^(n-2) + ... + an
import casadi.*

sz = length(a);
opti = casadi.Opti();
b = opti.variable(sz,1);
if s0 > 1
    error('s0 should be smaller than 1');
end
s_temp = linspace(s0, 1, N);
s = zeros(sz, N);
for i = 1:length(a)
    s(i, :) = s_temp.^(i-1);
end
ds = s(1:end-1,:);
ds = repmat([1:sz-1]', 1, N).*ds;

bs = b'*s;
dqds = b(2:end)'*ds;
as = a'*s;
es = as - bs;
e = es*es';

opti.minimize(e);
opti.subject_to(bs(1) == qs_init);
opti.subject_to(bs(end) == qs_end);
opti.subject_to(dqds(1) == dqds_init);
opti.subject_to(dqds(end) == 0);

opti.solver('ipopt');
sol = opti.solve();
b_coefficients = sol.value(b);

end
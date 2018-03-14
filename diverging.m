% Variables used in this scripts come from mpc_develop.m
% Let's check the lth joint
l = 4;
pos = x(1:6, 1:k);
vel = x(7:12, 1:k);

q4 = pos(l, :);
dq4 = vel(l, :);

tau4_applied = u_tracker(l, 1:k);

s4 = zeros(1, k);
for i = 1:k
    [~, index4] = min(abs(q4(i) - qs(l,:)));
    s4(i) = index4;
end
es4 = s4 - [1, s_sequence(1:k-1)];
plot(es4);
% es4 decreased at the end because q4 is larger than the maximum of qs. The
% estimation of s can't be larger than 1000

tau4_needed = zeros(1, k);
m_tracker = mcg_tracker(1:6, :);
c_tracker = mcg_tracker(7:12, :);
g_tracker = mcg_tracker(13:end, :);

i = 1;
while i<=k 
        temp_m4 = full( f_m(pos(:, i)) );
        temp_c4 = full( f_c(pos(:, i), dqs(:, s4(i))) );
        temp_g4 = full( f_g(pos(:, i), [0;0;0;0;0;0], [0;0;0;0;0;0]) );
        m4 = temp_m4*dqs(:, s4(i));
        c4 = temp_m4*ddqs(:, s4(i)) + temp_c4*dqs(:, s4(i));
        g4 = temp_g4 + f_fc*sign(dqs(:, s4(i)));
        temp_u = m4*a_tracker(i) + c4*b_tracker(i) + g4;
        tau4_needed(i) = temp_u(l);
        i = i+1;
end

eu4 = (tau4_needed(1:i-1) - tau4_applied(1:i-1));
plot(eu4);
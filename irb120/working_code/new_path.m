if si > 1
    error('s0 should be smaller than 1');
end
s_temppath = linspace(si, 1, N0);
spath = zeros(sz, N0);
for ipath = 1:length(pathic)
    spath(ipath, :) = s_temppath.^(ipath-1);
end
dspath = spath(1:end-1,:);
dspath = repmat([1:sz-1]', 1, N0).*dspath;

% bs = b'*s;
% dqds = b(2:end)'*ds;
aspath = pathic'*spath;
optip.set_value(asp, aspath);
optip.set_value(dsp, dspath);
optip.set_value(sp, spath);
optip.set_value(qsp_init, qs_initpath);
optip.set_value(dqdsp_init, dqds_initpath);
optip.set_value(qsp_end, qs_endpath);
% es = as - bs;
% e = es*es';
% 
% opti.minimize(e);
% opti.subject_to(bs(1) == qs_init);
% opti.subject_to(bs(end) == qs_end);
% opti.subject_to(dqds(1) == dqds_init);
% opti.subject_to(dqds(end) == 0);

% opti.solver('ipopt');
sol = optip.solve();
b_coefficients = sol.value(bp);

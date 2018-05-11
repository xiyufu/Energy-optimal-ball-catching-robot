T_fail = opti.debug.value(T);
tau_fail = opti.debug.value(tau);
e_tau = repmat(tau_max, 1, N0) - tau_fail
e_T = T_itc - T_fail
q_itc
T_itc
% if T_fail >= T_itc
%     fprintf('Time constraint is active, T = %5.4f, T_itc = %5,4f\n', T_fail, T_itc);
% end
% index_active = (e_tau < 0);
% if ~all(e_tau(:))
%     fprintf('Torque constraint is active');
% end
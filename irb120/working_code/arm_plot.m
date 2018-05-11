%%
figure; hold
grid on
title('Joints position');
xlabel('time [s]', 'FontSize', 18);
ylabel('position [rad]', 'FontSize', 18);
if flag_integratorF == 1
    kF = k+ii;
else
    kF = k;
end
for i = 1:6
    plot(t(1:kF), x(i, 1:kF), 'x', 'LineWidth', 1.5, 'MarkerSize', 8);
end
plot(t(1:k), q_ref(:,1:k), 'LineWidth', 3);
legend({'pos 1','pos 2','pos 3','pos 4','pos 5','pos 6','ref 1','ref 2','ref 3','ref 4','ref 5','ref 6'},'FontSize',18);
%%
figure; hold
grid on
title('Joints velocity');
xlabel('time [s]', 'FontSize', 18);
ylabel('velocity [rad/s]', 'FontSize', 18);
for i = 7:12
    plot(t(1:kF), x(i,1:kF), 'x', 'LineWidth', 1.5, 'MarkerSize', 8);
end
plot(t(1:k), dq_ref(:,1:k), 'LineWidth', 3);
legend({'vel 1','vel 2','vel 3','vel 4','vel 5','vel 6','ref 1','ref 2','ref 3','ref 4','ref 5','ref 6'},'FontSize',18);
%% 
figure; hold
grid on
title('Velocity error');
ylabel('velocity error [rad/s]', 'FontSize', 18);
xlabel('time [s]', 'FontSize', 18);
e_vel = abs(dq_ref(:,1:k) - x(7:12, 1:k));
for i = 1:6
    plot(t(1:kF), e_vel(i, 1:kF), 'LineWidth', 3);
    plot(t(1:kF), e_vel(i, 1:kF), 'o', 'MarkerSize', 8);
end
legend({'joint 1', ' ','joint 2', ' ', 'joint 3', ' ', 'joint 4', ' ', 'joint 5', ' ', 'joint 6', ' '},'FontSize',18);
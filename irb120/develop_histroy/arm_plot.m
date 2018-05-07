figure; hold
if flag_integratorF == 1
    kF = k+ii;
else
    kF = k;
end
for i = 1:6
    plot(t(1:kF), x(i, 1:kF), 'x');
end
plot(t(1:k), q_ref(:,1:k));
% plot(t(1:kF), q_itc(1)*ones(1,kF),t(1:kF), q_itc(2)*ones(1,kF));

figure; hold
for i = 7:12
    plot(t(1:kF), x(i,1:kF), '--');
end
plot(t(1:k), dq_ref(:,1:k));

figure; hold
if flag_integratorF == 1
    kF = k+ii;
else
    kF = k;
end
sc1 = scatter(t(1:kF), x(1,1:kF));
sc1.Marker = 'x';
sc2 = scatter(t(1:kF), x(2,1:kF));
sc2.Marker = 'x';
sc3 = scatter(t(1:kF), x(3,1:kF));
sc3.Marker = 'x';
sc4 = scatter(t(1:kF), x(4,1:kF));
sc4.Marker = 'x';
sc5 = scatter(t(1:kF), x(5,1:kF));
sc5.Marker = 'x';
sc6 = scatter(t(1:kF), x(6,1:kF));
sc6.Marker = 'x';

plot(t(1:k), q_ref(:,1:k));
plot(t(1:kF), q_itc(1)*ones(1,kF),t(1:kF), q_itc(2)*ones(1,kF));
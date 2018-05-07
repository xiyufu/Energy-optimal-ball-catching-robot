function [pos,vel,acc,jerk,px]= polytraj(pos_end, nr, pos_init, vel_init, acc_init, jerk_init)

% [pos,vel,acc,jerk]= polytraj(d,Ts, nr, nr_tot);
% generation of trajectory with 9th order polynomial
% d = distance
% Ts = sample period
% nr : number of samples of period of motion
% nr_tot = total number of samples of trajectory
% pos = position
% vel = velocity
% acc = acceleration
% jerk = jerk
% Jan Swevers March 2006
% Modified by Xiyu Fu, Feb 2018

tm = 1;

% determine the coefficients of the polynomial trajectory
% p = a9 t^9 + a8 t^8 ... a1 t^1 + a0

a0 = pos_init; % p(0) 
a1 = vel_init; % v(0) 
a2 = acc_init; % a(0) 
a3 = jerk_init; % j(0) 
a4 = 0; % diff(j(0)) 

%p(tm) = d;
A = [tm^9 tm^8 tm^7 tm^6 tm^5]; 
b = pos_end - [tm^4 tm^3 tm^2 tm 1]*[a4;a3;a2;a1;a0]; 

%v(tm) = 0;
A = [A; 9*tm^8 8*tm^7 7*tm^6 6*tm^5 5*tm^4];
b = [b;0 - [4*tm^3 3*tm^2 2*tm 1 0]*[a4;a3;a2;a1;a0]];

%a(tm) = 0;
A = [A; 8*9*tm^7 7*8*tm^6 6*7*tm^5 5*6*tm^4 4*5*tm^3];
b = [b;0 - [3*4*tm^2 2*3*tm 1*2 0 0]*[a4;a3;a2;a1;a0]];

%j(tm) = 0;
A = [A; 7*8*9*tm^6 6*7*8*tm^5 5*6*7*tm^4 4*5*6*tm^3 3*4*5*tm^2];
b = [b;0 - [2*3*4*tm 1*2*3 0 0 0]*[a4;a3;a2;a1;a0]];

%diff(j(tm)) = 0;
A = [A; 6*7*8*9*tm^5 5*6*7*8*tm^4 4*5*6*7*tm^3 3*4*5*6*tm^2 2*3*4*5*tm^1];
b = [b;0 - [1*2*3*4 0 0 0 0]*[a4;a3;a2;a1;a0]];


x = A\b;
x = x';

% This is the vector of coefficents
px = [x a4 a3 a2 a1 a0];

t1 = linspace(0, tm, nr+1);

pos = polyval(px,t1);

vx = [9 8 7 6 5 4 3 2 1 0].*px;
vel = polyval(vx(1:end-1),t1);

ax = [9*8 8*7 7*6 6*5 5*4 4*3 3*2 2*1 0 0].*px;
acc = polyval(ax(1:end-2),t1);

jx = [9*8*7 8*7*6 7*6*5 6*5*4 5*4*3 4*3*2 3*2*1 0 0 0].*px;
jerk = polyval(jx(1:end-3),t1);


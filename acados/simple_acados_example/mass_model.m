close all
clear;
import casadi.*

% constants
M = 1; % [kg] mass of the point
Fmax = 10; % [N] upper limit of force
Fmin = -10; % [N] lower limit of force
distance = 10; % [m] distance to travel 
pmax = 20; % [W] the accuator is limited to 20W

% states
x1 = SX.sym('pos', 1,1); % position 1D
v1 = SX.sym('vel', 1,1); % velocity 1D
x = vertcat(x1, v1); % states

% control inputs
u = SX.sym('force', 1,1); % input force 1D

nx = length(x); % number of states
nu = length(u); % number of control inputs

% ODE, pull a mass point to 10m away, minimize time.
% xd = v
% xdd = u/M
f = [x(1);...
    u/M];

% casadi function of dynamics
dyn_mass = casadi.Function('dyn_mass', {x, u}, {f});

% constraints
% power limitation
p = u*x(2) - pmax;
plimt = casadi.Function('plimt', {x, u}, {p});

% generate c code
opts = struct('mex', false, 'with_header', true, 'with_export', false);
dyn_mass.generate('dyn_mass.c', opts);
plimt.generate('power_limit.c', opts);

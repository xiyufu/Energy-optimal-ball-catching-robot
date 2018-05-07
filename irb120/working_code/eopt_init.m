% Initialize 
% Define golbal variables

%% Define the robot's kinematics
startup_rvc;

L(1) = Link('revolute', 'alpha', -pi/2, 'a', 0, 'd', 0.29, 'offset', 0);
L(2) = Link('revolute', 'alpha', 0, 'a', 0.27, 'd', 0, 'offset', -pi/2);
L(3) = Link('revolute', 'alpha', -pi/2, 'a', 0.07, 'd', 0, 'offset', 0);
L(4) = Link('revolute', 'alpha', pi/2, 'a', 0, 'd', 0.302, 'offset', 0);
L(5) = Link('revolute', 'alpha', -pi/2, 'a', 0, 'd', 0, 'offset', 0);
L(6) = Link('revolute', 'alpha', 0, 'a', 0, 'd', 0.21, 'offset', pi);

irb120 = SerialLink(L, 'name', 'irb120');

%% Joint range and torque limits
% Joint Range
joint_range = [-165, 165;...
               -110, 110;...
               -110, 70;...
               -160, 160;...
               -120, 120;...
               -400, 400];           
% Max speed
vel_max = [250; 250; 250; 320; 320; 420];
joint_range = joint_range*pi/180;
vel_max = vel_max*pi/180;
% Max gear torque
tg_max = [0.633; 0.552; 0.325; 0.146; 0.128; 0.201];
% transmission ratio
trans_ratio = [121; 121; 101; 50; 51; 50];
% Max joint torque
tau_max = tg_max.*trans_ratio;

%% Variables in the algorithm
% Sample period and frequency
fs = 250; % Hz
Ts = 1/fs; % 4 ms
% Number of collocation points
N0 = 15;
% Maximum iteration times
NMAX = 100;
% Matrix for states storage
x = zeros(12,NMAX+1); % states
u_act = zeros(6,NMAX); % control histroy
t = zeros(1,NMAX+1); % time
b_init_k = 1e-6; % initial value of b, can't be 0.
dq_pre = 0.001*ones(6, 1); % dq/ds after previous execution, initial 0.001
ddq_pre = 0.001*ones(6, 1); % dq^2/d^2s after previous execution, initial 0.001
q_ref = zeros(6,NMAX+1); % reference position of joints
dq_ref = zeros(6, NMAX+1); % reference velocity of joints
a_opt = ones(1, N0); % optimal a (numeral), will be updated after every optimization
b_opt = ones(1, N0+1); % optimal b (numeral)
tau_opt = ones(6, N0); % optimal tau (numeral)
nr = 999; % nr is defined for discrete path, we are not using it 
% nr should be an odd number so that ceil(s) and floor(s) won't give the
% same value;
flag_precise = 1; % if our execution is finished precisely = 1 otherwise = 0
flag_init = 1; % if it's the first run = 1 otherwise = 0
flag_Nchange = 0; % This is for Filip's algorithm, it's not really working. Not used
imprecise_counter = 0; % How many times is our execution imprecise
ei = zeros(6,NMAX); % error between q_ref(k) and x(1:6,k)
si_int = 0; % for discrete path, not used here
flag_integratorF = 0; % test only
s_sequence = zeros(1, NMAX); % test only
qh = zeros(6, N0); % q( (sk + sk+1)/2 )
dqh = zeros(6, N0); % dq( (sk + sk+1)/2 )
ddqh = zeros(6, N0); % ddq( (sk + sk+1)/2 )
temp_m = zeros(6, 6); % temporal storage of m, c, g
temp_c = zeros(6, 6);
temp_g = zeros(6, 1);
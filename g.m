function [g_out] = g(parms, q)

g_out = zeros(6, 1);

x0 = cos(q(2));
x1 = 9.81*x0;
x2 = sin(q(3));
x3 = -x2;
x4 = sin(q(2));
x5 = -9.81*x4;
x6 = cos(q(3));
x7 = x1*x3 + x5*x6;
x8 = cos(q(4));
x9 = sin(q(4));
x10 = x1*x6 + x2*x5;
x11 = x10*x9;
x12 = cos(q(5));
x13 = x10*x8;
x14 = sin(q(5));
x15 = -x14;
x16 = x12*x7 + x13*x15;
x17 = sin(q(6));
x18 = x12*x13 + x14*x7;
x19 = cos(q(6));
x20 = -x19;
x21 = -x11;
x22 = -x21;
x23 = x17*x22 + x18*x20;
x24 = -x16;
x25 = parms(72)*x24 + parms(74)*x23;
x26 = x17*x18 + x20*x21;
x27 = parms(73)*x16 - parms(74)*x26;
x28 = parms(75)*x23;
x29 = 0.21*x17;
x30 = parms(75)*x26;
x31 = parms(60)*x16 + parms(61)*x22 + x17*x25 + 0.21*x19*x30 + x20*x27 + x28*x29;
x32 = -x7;
x33 = parms(72)*x26 - parms(73)*x23;
x34 = parms(59)*x21 - parms(60)*x18 + x33;
x35 = parms(47)*x11 + parms(48)*x32 + x12*x31 + x15*x34;
x36 = -x17;
x37 = -parms(59)*x24 - parms(61)*x18 + 0.21*x19*x28 - x20*x25 - x27*x36 - x29*x30;
x38 = parms(46)*x7 - parms(47)*x13 + x37;
x39 = parms(49)*x11 - parms(62)*x21 - x20*x30 - x28*x36;
x40 = 0.302*x39;
x41 = parms(62)*x18 + x17*x30 + x20*x28;
x42 = parms(62)*x16 + parms(75)*x16;
x43 = parms(49)*x13 + x12*x41 + x15*x42;
x44 = x43*x9;
x45 = parms(34)*x7 + x35*x8 + x38*x9 + x40*x8 - 0.302*x44;
x46 = parms(46)*x21 + parms(48)*x13 + x12*x34 + x14*x31;
x47 = -parms(34)*x10 + x46;
x48 = -x8;
x49 = x39*x48 + x44;
x50 = 0.07*x2^2 + 0.07*x6^2;
x51 = x49*x50;
x52 = x43*x8;
x53 = -parms(33)*x32 - parms(35)*x10 - x35*x9 - x38*x48 - x40*x9 - 0.302*x52;
x54 = parms(36)*x7 + parms(49)*x7 + x12*x42 + x14*x41;

g_out(1) = x0*(-parms(22)*x5 + x3*x47 + x3*x51 + x45*x6) + x4*x49*(-0.27*x0^2 - 0.27*x4^2) - x4*(parms(22)*x1 + x2*x45 + x47*x6 + x51*x6);
g_out(2) = parms(20)*x5 - parms(21)*x1 + 0.27*parms(23)*x5 + 0.27*x2*(parms(36)*x10 + x39*x9 + x52) + x50*x54 + x53 + 0.27*x54*x6;
g_out(3) = x53 + 0.07*x54;
g_out(4) = x46;
g_out(5) = x37;
g_out(6) = x33;

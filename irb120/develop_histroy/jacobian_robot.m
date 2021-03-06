function [jac_out] = jacobian(q)
import casadi.*

if isa(q,'casadi.SX')
    jac_out = SX.zeros(36,1);
elseif isa(q,'casadi.MX')
    jac_out = MX.zeros(36,1);
elseif isa(q,'casadi.DM')
    jac_out = DM.zeros(36,1);
else
    jac_out = zeros(36, 1);
end

x0 = sin(q(3));
x1 = sin(q(1));
x2 = cos(q(2));
x3 = x1*x2;
x4 = sin(q(2));
x5 = cos(q(3));
x6 = x4*x5;
x7 = x1*x6;
x8 = x0*x3 + x7;
x9 = cos(q(4));
x10 = sin(q(4));
x11 = cos(q(1));
x12 = -x11;
x13 = x10*x12 + x8*x9;
x14 = sin(q(5));
x15 = x13*x14;
x16 = cos(q(5));
x17 = x3*x5;
x18 = -x1;
x19 = x0*x4;
x20 = x17 + x18*x19;
x21 = x16*x20;
x22 = -0.07*x0;
x23 = x1*x4;
x24 = x0*x23;
x25 = 0.21*x16;
x26 = x0*x2;
x27 = -x26 - x6;
x28 = -0.21*x14;
x29 = x2*x5;
x30 = -x19 + x29;
x31 = x30*x9;
x32 = x25*x27 + x28*x31;
x33 = x22*x4 + 0.07*x29;
x34 = -0.302*x26 - 0.302*x6;
x35 = 0.27*x2 + x32 + x33 + x34;
x36 = x32 + x33 + x34;
x37 = x32 + x34;
x38 = x13*x28 + x20*x25;
x39 = 0.302*x17 - 0.302*x24;
x40 = x38 + x39;
x41 = x10*x8 + x11*x9;
x42 = x10*x30;
x43 = -x38;
x44 = -x15 + x21;
x45 = -x14;
x46 = x16*x27 + x31*x45;
x47 = x11*x29 + x12*x19;
x48 = x11*x26 + x11*x6;
x49 = x1*x10 + x48*x9;
x50 = x25*x47 + x28*x49;
x51 = 0.07*x11;
x52 = x26*x51 + x51*x6;
x53 = -0.302*x11*x19 + 0.302*x11*x29;
x54 = 0.27*x11*x4 + x50 + x52 + x53;
x55 = x50 + x53;
x56 = x10*x48 + x18*x9;
x57 = -x32;
x58 = x16*x47 + x45*x49;
x59 = 0.07*x0*x3 + 0.07*x7;
x60 = -x50;

jac_out(1) = -0.27*x1*x4 + 0.21*x15 - 0.302*x17 - 0.21*x21 + x22*x3 + 0.302*x24 - 0.07*x7;
jac_out(2) = x11*x35;
jac_out(3) = x11*x36;
jac_out(4) = x20*x37 - x27*x40;
jac_out(5) = x32*x41 + x42*x43;
jac_out(6) = x32*x44 + x43*x46;
jac_out(7) = x54;
jac_out(8) = x1*x35;
jac_out(9) = x1*x36;
jac_out(10) = x27*x55 - x37*x47;
jac_out(11) = x42*x50 + x56*x57;
jac_out(12) = x46*x50 + x57*x58;
jac_out(13) = 0;
jac_out(14) = x12*x54 + x18*(0.27*x23 + x38 + x39 + x59);
jac_out(15) = x12*(x50 + x52 + x53) + x18*(x38 + x39 + x59);
jac_out(16) = -x20*x55 + x40*x47;
jac_out(17) = x38*x56 + x41*x60;
jac_out(18) = x38*x58 + x44*x60;
jac_out(19) = 0;
jac_out(20) = x18;
jac_out(21) = x18;
jac_out(22) = x47;
jac_out(23) = x56;
jac_out(24) = x58;
jac_out(25) = 0;
jac_out(26) = x11;
jac_out(27) = x11;
jac_out(28) = x20;
jac_out(29) = x41;
jac_out(30) = x44;
jac_out(31) = 1;
jac_out(32) = 0;
jac_out(33) = 0;
jac_out(34) = x27;
jac_out(35) = x42;
jac_out(36) = x46;

jac_out = reshape(jac_out, 6,6);

end


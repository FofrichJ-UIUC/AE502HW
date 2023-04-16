clear all
close all
syms omega [3 1]
syms q [3 1]
syms p [3 1]

omega1 =0;omega2=0;p
Omega = skew(omega);

test = q.'*Omega*q


function M = skew(x)
M=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ].';

end
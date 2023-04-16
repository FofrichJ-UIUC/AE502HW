function [rVec,vVec] = oe2cart(a,e,i,RAAN,w,M,time,f,mu)
%UNTITLED3 Summary of this function goes here
%   outputs the vectors in 1x3 format
% inputs in degrees of angle
%time needed for solving kepler's equation

M = M*pi/180;

if abs(M)>0

    f = keplerOE(e,a,M,time,mu,false,1e-15);
end


i = i*pi/180;
RAAN = RAAN*pi/180;
w = w*pi/180;
f = f;

cf = cos(f);
sf = sin(f);

p = a*(1-e^2);
r = p/(1+e*cos(f));

r_focal = r.*[cf; sf; 0];

h = sqrt(mu*p);

v_focal = (mu/h).* [-sf; (e+cf); 0];
% v = sqrt(mu*(2/r -1/a));

cOm = cos(RAAN);
sOm = sin(RAAN);

ci = cos(i);
si = sin(i);



cw = cos(w);
sw = sin(w);

%using curtis Algorithm 4.5 because whatever I was doing before was wrong
Qxx = [-sOm*ci*sw+cOm*cw, -sOm*ci*cw-cOm*sw, sOm*si;
    cOm*ci*sw+sOm*cw, cOm*ci*cw-sOm*sw, -cOm*si;
    si*sw, si*cw, ci];

rVec = Qxx*r_focal;
vVec = Qxx*v_focal;

end
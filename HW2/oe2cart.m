function [rVec,vVec] = oe2cart(a,e,i,RAAN,w,f,mu)
%UNTITLED3 Summary of this function goes here
%   outputs the vectors in 1x3 format

i = i*pi/180;
RAAN = RAAN*pi/180;
w = w*pi/180;
f = f*pi/180;
r = a*(1-e^2)/(1+e*cos(f));

h = sqrt(mu*a*(1-e^2));

cOm = cos(RAAN);
sOm = sin(RAAN);

ci = cos(i);
si = sin(i);

f = f+w;

cf = cos(f);
sf = sin(f);

cw = cos(w);
sw = sin(w);

rVec = r.*[(cOm*cf-sOm*sf*ci) (sOm*cf+cOm*sf*ci) sf*si];
vVec = (-mu/h).*[(cOm*(sf+e*sw) + sOm*(cf+e*cw)*ci) (sOm*(sf+e*sw)-cOm*(cf+e*cw)*ci) -(cf+e*cw)*si];

end
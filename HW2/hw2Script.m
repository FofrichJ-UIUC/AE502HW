clear all
close all

pb = 3;


if pb ==1
% Problem 1
J2_Earth = 0.00108;
Re       = 6370;
mu = 398600;

rp = 600 + Re;
T0 = 24*3600*3;
a = (mu*(T0/(2*pi))^2)^(1/3);

Reasq = (Re/a)^2;

n = sqrt(mu/a^3);
e = 1-(rp/a);
i = acos(sqrt((3*n*J2_Earth*Reasq))/(20*((1-e^2)^2)));
om_dot = -3/2 * n * J2_Earth * Reasq* cos(i)/((1-e^2)^2);



% Problem 2
elseif pb==2

J2_Mars = 0.00196;
Re       = 3390;
mu = 42820;

rp = 400 + Re;
T0 = 24*3600 + 39*60 + 35;
a = (mu*(T0/(2*pi))^2)^(1/3);

Reasq = (Re/a)^2;

n = sqrt(mu/a^3);
e = 1-(rp/a);
i = acos(sqrt((3*n*J2_Mars*Reasq))/(20*((1-e^2)^2)));
om_dot = -3/2 * n * J2_Mars * Reasq* cos(i)/((1-e^2)^2);

% Problem 3
elseif pb==3

J2_Earth = 0.00108;
Re       = 6370;
mu = 398600;
tspan = linspace(0,100*86400,10000);
opts = odeset('RelTol',1e-11,'AbsTol',1e-13);

%km          deg                nondim    deg       deg     deg
a = 26000; i = 1.10654*180/pi; e = 0.74; w = 5; RAAN = 90; M = 10;

f = keplerOE(e,a,M,0,mu,false,1e-15);

[r1,v1] = oe2cart(a,e,i,RAAN,w,f,mu);


IC = [r1.';v1.'];

[t,y] = ode113(@(t,y) TwoBodywJ2(t,y,mu,J2_Earth,Re), tspan, IC,opts);

anew = zeros(length(t),1); 
inew = zeros(length(t),1); 
enew = zeros(length(t),1); 
wnew = zeros(length(t),1); 
RAANnew = zeros(length(t),1);

for j = 1:length(t)
    %[a,e,i,RAAN,w,f]
    output = rv2oe(y(j,1:3).',y(j,4:6).',mu);
    anew(j) = output(1);
    enew(j) = output(2);
    inew(j) = output(3);
    RAANnew(j) =output(4);
    wnew(j) = output(5);

end
figure()
plot(t/86400,inew*180/pi)
grid on 
grid minor
xlabel('Time (days)')
ylabel('Inclination (deg)')
figure()
plot(t/86400,RAANnew*180/pi)
grid on 
grid minor
xlabel('Time (days)')
ylabel('RAAN (deg)')
figure()
plot(t/86400,wnew*180/pi)
grid on 
grid minor
xlabel('Time (days)')
ylabel('Argument of Periapsis (deg)')
figure()
plot(t/86400,enew)
grid on 
grid minor
xlabel('Time (days)')
ylabel('Eccentricity')

figure()
plot(t/86400,anew)
grid on 
grid minor
xlabel('Time (days)')
ylabel('Semi-Major Axis (km)')


end






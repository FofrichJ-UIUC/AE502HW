clear all
close all



% Mean anomaly, arg peripsis, RAAN

omegaf = 0.01;
j = 45*pi/180;
a = 1;
e = 0.5;
L = 1;
G = (1-e^2)^0.5;
H = G*cos(j);

number = 100001;


%assume initial M ,omega, and Omega are zero

l = 0; %M
ll(1) =l;

g = 0;  %omega
gg(1) =g; 
h = 0;  %Omega
% h = g*cos(j);
hh(1) = h;
theta(1) = 0;
% ii(1) = j;



% a = a*ones(1,1001);
% e = e*ones(1,1001);
% j = j*ones(1,1001);



n = 1; %due to n % n = sqrt(mu/a^3)

time = linspace(0,100,number);

for i =2:length(time)
    delt = time(i) - time(i-1);
    theta(i) =theta(i-1) +omegaf*delt;
    
%     l = -delt* ((1/(L^3)) +omegaf*(1-e^2)^0.5 *cos(j));
     l = -delt*(1/(L^3));
    ll(i) = ll(i-1) +l;

%     g = -delt*(omegaf*cos(j));
%     g= -delt*(((1-e^2)/(G^3)) +omegaf *cos(j));
g = 0;
    gg(i) = gg(i-1) +g;

    h = -delt*omegaf;
%     h = -delt*(((1-e^2)*(cos(j)^2)/(H^3)) +omegaf );
    
    hh(i) = hh(i-1) +h;
%     testvar = acos(hh(i)/gg(i))
%     ii = sqrt(1-(gg(i)/ll(i))^2)

end

a = a*ones(1,number);
e = e*ones(1,number);
i = j*ones(1,number)*180/pi;
omega = gg*180/pi;
Omega =hh*180/pi;
M = ll*180/pi;

% i = j*ones(1,number);
% omega = gg;
% Omega =hh;
% M = ll;
for jj = 1:length(i)
%     f(jj) = keplerOE(e(jj),a(jj),M(jj),time(jj),1,false,1e-15);

    [r,v1] = oe2cart(a(jj),e(jj),i(jj),Omega(jj),omega(jj),M(jj),time(jj),0,1);


%     r = a(jj)*(1-e(jj)^2)/(1+e(jj)*cos(i(jj)*pi/180));

%     r1 = r*[cos(gg(jj)+f(jj))*cos(hh(jj))-sin(gg(jj)+f(jj))*sin(hh(jj))*cos(i(jj)*pi/180);
%         cos(gg(jj)+f(jj))*sin(hh(jj))-sin(gg(jj)+f(jj))*cos(hh(jj))*cos(i(jj)*pi/180);
%         sin(gg(jj)+f(jj))*sin(i(jj)*pi/180)];
    Q = [cos(theta(jj)) -sin(theta(jj)) 0;
        sin(theta(jj)) cos(theta(jj)) 0;
        0 0 1]; % https://en.wikipedia.org/wiki/Rotation_matrix
    R(:,jj) = r;

end

plot3(R(1,:),R(2,:),R(3,:))
grid on 
grid minor
xlabel('Inertial X-axis (ND)')
ylabel('Inertial Y-axis (ND)')
zlabel('Inertial Z-axis (ND)')
axis equal
title('Orbit Propagated for 100 time units')

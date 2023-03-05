function YY = TwoBodywJ2(t,XX,mu,J2,R)

rnorm       = norm(XX(1:3));


p = (3/2 *J2*R^2/rnorm^5) .*[XX(1)*(5*XX(3)^2/rnorm^2 -1);XX(2)*(5*XX(3)^2/rnorm^2 -1);XX(3)*(5*XX(3)^2/rnorm^2 -3)];



YY    = [XX(4); XX(5); XX(6); -mu*XX(1)/rnorm^3+p(1); -mu*XX(2)/rnorm^3+p(2); -mu*XX(3)/rnorm^3+p(3)];
end
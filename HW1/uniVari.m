function result1 = uniVari(r0,v0,dt,mu)
    %#universal variable calculation using curtis 3.3 and 3.4 algorithm
    r0norm = sqrt(r0(1)^2 + r0(2)^2 + r0(3)^2);
    vnorm = sqrt(v0(1)^2 + v0(2)^2 + v0(3)^2);
    vr    = dot(v0,r0)/r0norm;
    alpha = 2.0/r0norm - (vnorm^2)/mu;
    X = uniAnom(r0norm,vr,dt,mu,alpha);
    CS = stumpff(alpha*X^2);
    f = 1 - ((X^2)/r0norm)*CS(1);
    g = dt - (1/sqrt(mu))*(X^3)*CS(2);

    r = f*r0+g*v0;

    rnorm = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
    fdot = (sqrt(mu)/(rnorm*r0norm))*(CS(2)*alpha*(X^3)-X);
    gdot = (1-CS(1)*(X^2)/rnorm);

    v = fdot*r0 + gdot*v0;
    result1 = [r;v];
    
end 
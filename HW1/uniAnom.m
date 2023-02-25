function xstart = uniAnom(r0,vr,dt,mu,alpha)
    %#universal variable calculation using curtis 3.3 and 3.4 algorithm
    tol = 1e-8;
    xstart = sqrt(mu)*dt*abs(alpha);
    z = alpha*xstart^2;
    function ff =ratio(x,r0,vr,mu,dt,z)
        a = stumpff(z);
        smu = sqrt(mu);
        frac = r0*vr*(x)/smu;
        f = frac*x*a(1) + (1-alpha*r0)*(x^3)*a(2)+r0*x - smu*dt;
        df = frac*(1-alpha*(x^2)*a(2)) + (1-alpha*r0)*(x^2)*a(1)+r0;
        ff = f/df;
        
    end
    ff = ratio(xstart,r0,vr,mu,dt,z);
    while abs(ff)>tol
        xstart = xstart - ff;
        z = alpha*xstart^2;
        ff = ratio(xstart,r0,vr,mu,dt,z);
    end

    
end
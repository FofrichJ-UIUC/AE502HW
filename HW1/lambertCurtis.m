function output =lambertCurtis(r1,r2,dt,mu,prograde)
    %#Curtis Algorithm 5.2 that uses universal variable to calculate
    %#prograde = 1 then its prograde trajectory else its retrograde
     %   #r1 is 3x1 vector and so is r2
        tol = 1e-8;
    r1norm = sqrt(r1(1)^2 + r1(2)^2 + r1(3)^2);
    r2norm = sqrt(r2(1)^2 + r2(2)^2 + r2(3)^2);
    crossProd = cross(r1,r2);
    flag = crossProd(3);
    
    inVar = acos(dot(r1,r2)/(r1norm*r2norm));
    if prograde==1
        if flag>=0
            dtheta = inVar;
        else
            dtheta = 2*pi - inVar;
        end
    
    else
        if flag>=0
            
            dtheta = 2*pi - inVar;
        else
            dtheta = inVar;
        end
    end

    A = sin(dtheta)*sqrt((r1norm*r2norm)/(1-cos(dtheta)));

    z = -100;
    
    nmax = 500;
    n=0;

    
    ff = -1;
    
    while real(ff)<0 
        CS = stumpff(z);
        y = r1norm + r2norm +A*(z*CS(2)-1)/sqrt(CS(1));
        
        ff = (y/CS(1))^(1.5)*CS(2)+A*(y)^0.5-sqrt(mu)*dt;
       
        z = z +0.1;

        
        
        
    end

    ff =1;

    while abs(real(ff))>tol && n<nmax
        CS = stumpff(z);
        y = r1norm + r2norm +A*(z*CS(2)-1)/sqrt(CS(1));
        
        F = (y/CS(1))^(1.5)*CS(2)+A*(y)^0.5-sqrt(mu)*dt;
        if z ==0
            dF = (sqrt(2)/40)*(y)^(1.5) + (A/8)*((y)^0.5 + A*(1/(2*y) )^0.5);
        else
            dF = ((y/CS(1) )^1.5)*(    (1/(2*z))* (CS(1) - (3*CS(2)/(2*CS(1)))) + (3*CS(2)^2)/(4*CS(1))) + (A/8)*((y)^0.5*3*CS(2)/CS(1) + A*(CS(1)/y )^0.5);
        end
        ff = F/dF;
        z = z - ff;

        n = n+1;
        
        
    end
    z = real(z);
    CS = stumpff(z);
    y = r1norm + r2norm +A*(z*CS(2)-1)/sqrt(CS(1));
    

    f = 1 - y/r1norm;
    g = A*sqrt(y/mu);

    fdot = (sqrt(mu)/(r1norm*r2norm))*sqrt(y/CS(1))*(z*CS(2)-1);
    gdot = 1 - y/r2norm;
    
    v1 = (r2-f*r1)/g;
    v2 = (r2*gdot - r1)/g;
    output = [v1 v2];

    
end
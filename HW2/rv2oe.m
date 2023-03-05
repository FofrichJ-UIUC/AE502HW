function output = rv2oe(rVec,vVec,mu)
   % #rVec 3x1 and same for vVec
    
    r = norm(rVec);
    v = norm(vVec);
    
    I = [1; 0; 0];
    J = [0; 1; 0];
    K = [0; 0; 1];
    
    
    a = 1/(2/r - (v^2)/mu);

    eVec = ((v^2)/mu - 1/r)*rVec - 1/mu * (dot(rVec,vVec))*vVec;

    e = norm(eVec);
    
    hVec = cross(rVec,vVec);
    h = norm(hVec);
    
    i = acos(dot((hVec/h),K));
    
    nVec = cross(K, hVec);
    n = norm(nVec);
    
    if dot(nVec,J) < 0 
        RAAN = 2*pi - acos(dot(nVec/n,I));
    else
        RAAN = acos(dot(nVec/n,I));
    
    end
    
    if dot(eVec,K) < 0
        w = 2*pi - acos(dot(nVec/n,eVec/e));
    else
        w = acos(dot(nVec/n,eVec/e));
    end
    
    rdotv = dot(rVec,vVec);
    
    if rdotv < 0
    
        f = 2*pi - acos(dot(eVec/e,rVec/r));
    else
        f = acos(dot(eVec/e,rVec/r));
    end

    output = [a,e,i,RAAN,w,f];

    
end
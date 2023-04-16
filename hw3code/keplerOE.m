function f = keplerOE(e,a,M,t,mu,laguerre,tol)

%finds the True anamoly after a certain period of time
%This assumes that the delta time give is just after periapsis

% if t>0
% %t is time after Periapsis in seconds
% M = sqrt(mu/a^3)*t;     %Rads
% end


% %using the method from Orbital Mechanics by Prussing & Conway 2nd ed pg 31
% %eqn 2.16
u = M + e;              %Rads
eccAnomGuess = (M*(1-sin(u))+u*sin(M))/(1+sin(M)-sin(u));       %Rads


%order of the laguerre's method algorithm
num = 5;
err = 1;
iter = 0;

    if laguerre == true 
        while err>tol
            F0 = eccAnomGuess-e*sin(eccAnomGuess) -M;
            F1 = 1-e*cos(eccAnomGuess);

    
            pos = sign(F1);
            if pos >0

                eccAnom = eccAnomGuess - (num*F0)/(F1+sqrt(((num-1)^2)*F1^2 - num*(num-1)*F0*e*sin(eccAnomGuess)));
            else
 
                eccAnom = eccAnomGuess - (num*F0)/(F1-sqrt(((num-1)^2)*F1^2 - num*(num-1)*F0*e*sin(eccAnomGuess))); 
            end
    
            err = abs(abs(eccAnom)-abs(eccAnomGuess));
            eccAnomGuess = eccAnom;

            iter = iter+1;
    
        end
    else
        while err>1e-12
            F0 = eccAnomGuess-e*sin(eccAnomGuess) -M;
            F1 = 1-e*cos(eccAnomGuess);
            
    
            eccAnom = eccAnomGuess - (F0)/(F1);
    
            err = abs(abs(eccAnom)-abs(eccAnomGuess));
            eccAnomGuess = eccAnom;
            iter = iter+1;
             
    
        end

    end
% disp(iter)
eccAnom = eccAnomGuess;

%true anamoly Calc
f = 2*atan(sqrt((1+e)/(1-e))*tan(eccAnom/2));       %Rads



end
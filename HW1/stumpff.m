function result =stumpff(z)
    %#stumpff function to calculate the C and S functions
    
    
    

    if z>0
        S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
        C = (1 - cos(sqrt(z)))/z;
    elseif z==0
        S = 1.0/6.0;
        C = 1.0/2.0;
    else
        S = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z))^3;
        C = (cosh(sqrt(-z))-1)/(-z);
    end
    result = [C;S];
    
end
function [a1 b1]=sine_surrogate(c1,c2,p_t,u)
distance=norm(p_t-u);

a=abs(c1)/2;
b=-2*a*c2^2*distance+c1*c2*cos(c2*distance);

a=a*c2^2;

% b=b*c2;

if b>=0
    a1=a+(b/(distance*2));
    b1=-2*a1*u;
else
    a1=a;
    b1=b*(p_t-u)/distance-2*a1*u;
end




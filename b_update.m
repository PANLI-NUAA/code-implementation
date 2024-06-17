function [bd]=b_update(c_1,c_2,c_3,c_4)
%%It solves the minimization problem of c_1(b_1^2+b_2^2+1)^2+c_2(b_1+b_2)+c_3b_1+c_4b_2

a1=c_1/(c_3^2+c_4^2)^2;
b1=(2*c_1+c_2)/((c_3^2+c_4^2));
c=1;

a=4*a1;
b=2*b1;
y=(-c/(2*a)+sqrt(81*a^2*c^2+12*a*b^3)/(18*a^2))^(1/3);

solution=y-b/(3*a*y);
% a*solution^3+b*solution+1

A=[c_3 c_4;c_4 -c_3];
d=[solution;0];

bd=inv(A)*d;

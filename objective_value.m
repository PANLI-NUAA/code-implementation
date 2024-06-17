function [value]=objective_value(y,b_t,L,N,delta_f,p_t,u,c)      
z1=0;
      
for l=1:L
    for n=1:N
        q1=real(y(n,l))*real(b_t(l)*exp(-1i*2*pi*(n-1)*delta_f*norm(p_t-u(:,l))/c));
        q2=imag(y(n,l))*imag(b_t(l)*exp(-1i*2*pi*(n-1)*delta_f*norm(p_t-u(:,l))/c));
z1=-log(normcdf(q1))-log(normcdf(q2))+z1;
    end
end

value=z1;
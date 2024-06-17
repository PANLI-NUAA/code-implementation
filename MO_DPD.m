function [time,p_t,fz]=MO_DPD(B,N,b,noise,u,p,SNR,xx,yy)
delta_f=B/N;
L=length(u(1,:));
received_signal_power=abs(b).^2;


SIGMA=mean(received_signal_power)/10^(SNR/10);%%noise power
c=3*10^8;
%%channel response measurement
for l=1:L
    for n=1:N
        tau=norm(u(:,l)-p);
r_u(l,n)=b(l).*exp(-1i*2*pi*(n-1)*delta_f*tau/c);
    end
end
r=r_u+noise*(sqrt(SIGMA));
% r=r_u;
% mean(mean(abs((sqrt(SNR))*noise).^2))
DN=0:1:N-1;
DN=DN.';
%%quantization
y=sign(real(r))+1i*sign(imag(r));
y=y.';

y_R=real(y);
y_I=imag(y);
%%开始计时
tic;
%%initialization
   for i=1:length(xx)
        for j=1:length(yy)
            p_e=[xx(i);yy(j)];
 z=0;
            for l=1:L
                F2=exp(-1i*2*pi*DN*delta_f*norm(p_e-u(:,l))/c);
  z=z+abs(F2'*y(:,l))^2;
            end
            f(i,j)=z;
        end      
    end
    [I,J]=find(f==max(max(f)));
        p_0=[xx(I(1));yy(J(1))];

p_t=p_0;
p_tt=p_0+1;
for l=1:L
b_t(l)=exp(-1i*2*pi*DN*delta_f*norm(p_0-u(:,l))/c)'*y(:,l)/N;
end
iii=1;
while norm(p_t-p_tt)>=1e-4%1e3 is the stopping criterion 
%          fz(iii)=objective_value(y,b_t,L,N,delta_f,p_t,u,c);    
    for l=1:L
                F=b_t(l)*exp(-1i*2*pi*DN*delta_f*norm(p_t-u(:,l))/c);
                gc(:,l)=y_R(:,l).*g(y_R(:,l).*real(F))+1i*y_I(:,l).*g(y_I(:,l).*imag(F));
    end
    for i=1:length(xx)
        for j=1:length(yy)
            p_e=[xx(i);yy(j)];
 z=0;
            for l=1:L
                F2=exp(-1i*2*pi*DN*delta_f*norm(p_e-u(:,l))/c);
  z=z+abs(F2'*gc(:,l))^2;
            end
            f(i,j)=z;
        end      
    end
    [I,J]=find(f==max(max(f)));


    p_tt=p_t;
        p_t=[xx(I(1));yy(J(1))];
    for l=1:L
b_t(l)=exp(-1i*2*pi*DN*delta_f*norm(p_t-u(:,l))/c)'*gc(:,l)/N;
    end
        iii=iii+1;
      z1=0;
      
% for l=1:L
%     for n=1:N
%         q1=real(y(n,l))*real(b_t(l)*exp(-1i*2*pi*(n-1)*delta_f*norm(p_t-u(:,l))/c));
%         q2=imag(y(n,l))*imag(b_t(l)*exp(-1i*2*pi*(n-1)*delta_f*norm(p_t-u(:,l))/c));
% z1=-log( normcdf(q1))-log( normcdf(q2))+z1;
%     end
% end
%     fz(iii)=z1;    
end


time=toc;
%%结束计时
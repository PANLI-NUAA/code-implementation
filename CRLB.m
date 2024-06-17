function SPEB=CRLB(B,N,b,u,p,SNR,L)

Q=[(p(1)-u(1,:))./vecnorm(p-u);(p(2)-u(2,:))./vecnorm(p-u)];
delta_f=B/N;
received_signal_power=abs(b).^2;


sigma=sqrt(mean(received_signal_power)/10^(SNR/10));%%noise power

c=3*10^8;

FIM=zeros(2*L+2);
for l=1:L
    for n=1:N
        tau=norm(u(:,l)-p);
z=exp(-1i*2*pi*(n-1)*delta_f*tau/c);

g=real(b(l)*z);
h=imag(b(l)*z);

alpha=-g/(sigma/sqrt(2));
kappa=-h/(sigma/sqrt(2));

z11=-real(b(l))*sin(2*pi*(n-1)*delta_f*tau/c)+imag(b(l))*cos(2*pi*(n-1)*delta_f*tau/c);
z21=-real(b(l))*cos(2*pi*(n-1)*delta_f*tau/c)-imag(b(l))*sin(2*pi*(n-1)*delta_f*tau/c);

z11=z11*(2*pi*(n-1)*delta_f/c);
z21=z21*(2*pi*(n-1)*delta_f/c);
d1=[zeros(1,(l-1)*2) real(z) -imag(z) zeros(1,2*L-(l-1)*2-2) z11*Q(:,l).'];
d2=[zeros(1,(l-1)*2) imag(z)  real(z) zeros(1,2*L-(l-1)*2-2) z21*Q(:,l).'];


f1=4*((1/sqrt(2*pi))*exp(-alpha^2/2))^2/(erfc(alpha/sqrt(2))*erfc(-alpha/sqrt(2)));
f2=4*((1/sqrt(2*pi))*exp(-kappa^2/2))^2/(erfc(kappa/sqrt(2))*erfc(-kappa/sqrt(2)));


FIM=f1*d1.'*d1+f2*d2.'*d2+FIM;

    end
end

FIM=FIM/sigma^2*2;

CRLB=pinv(FIM);
SPEB=sqrt(trace(CRLB(end-1:end,end-1:end)));
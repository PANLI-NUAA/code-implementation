function [time,p_t,fz]=MO_DPD_new_SQUAREM(B,N,b,noise,u,p,SNR,xx,yy)
delta_f=B/N;
L=length(u(1,:));
received_signal_power=abs(b).^2;


SIGMA=mean(received_signal_power)/10^(SNR/10);%%noise power
% SNR=1/10^(SNR/20);


c=3*10^8;

A_r0=0.5*[zeros(2,2) eye(2);eye(2) zeros(2,2)];
Aio=[zeros(2,2) [0 -1;1 0];zeros(2,2) zeros(2,2)];

A_i0=(Aio+Aio.')/2;

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
%%begin timing
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
% p_0=[(max(xx)-min(xx))*rand(1)+min(xx)  (max(yy)-min(yy))*rand(1)+min(yy)].';

p_t=p_0;
p_tt=p_0+1;
for l=1:L
b_t(l)=exp(-1i*2*pi*DN*delta_f*norm(p_0-u(:,l))/c)'*y(:,l)/N;
end
iii=1;
while norm(p_t-p_tt)>=1e-4%1e3 is the stopping criterion 
%      fz(iii)=objective_value(y,b_t,L,N,delta_f,p_t,u,c);    
    p_tt=p_t;
    u_t0=[real(b_t) imag(b_t) p_t.'].';
    for iiz=1:2
    
    for l=1:L
                F=b_t(l)*exp(-1i*2*pi*DN*delta_f*norm(p_t-u(:,l))/c);
                gc_r(:,l)=y_R(:,l).*g(y_R(:,l).*real(F));
                gc_i(:,l)=y_I(:,l).*g(y_I(:,l).*imag(F));
    end
    
      a0=0;
      b0=zeros(2,1);
  for l=1:L
      c1=0;
      c2=0;
      c3=0;
      c4=0;
      for n=1:N
      theta=2*pi*(n-1)*delta_f*norm(p_t-u(:,l))/c;
            u2=[cos(theta) sin(theta)].';
            x_t=[real(b_t(l)) imag(b_t(l)) u2.'].';            
            A_r=A_r0-gc_r(n,l)*[zeros(2,2) zeros(2,2);zeros(2,2) eye(2)];
            A_i=A_i0-gc_i(n,l)*[zeros(2,2) zeros(2,2);zeros(2,2) eye(2)];
          [v1r, v2r, v3r]=Quartic_Surrogate(A_r,x_t);
          [v1i, v2i, v3i]=Quartic_Surrogate(A_i,x_t);
          
          c1=v1r+v1i+c1;
          c2=v2r+v2i+c2;
          c3=v3r(1)+v3i(1)+c3;
          c4=v3r(2)+v3i(2)+c4;
          
          [a1, b1]=cosine_surrogate(v3r(3)+v3i(3),2*pi*(n-1)*delta_f/c,p_t,u(:,l));
          [a2, b2]=sine_surrogate(v3r(4)+v3i(4),2*pi*(n-1)*delta_f/c,p_t,u(:,l));
          a0=a1+a2+a0;
          b0=b1+b2+b0;
          
      end
      [bd]=b_update(c1,c2,c3,c4);
      b_t(l)=bd.'*[1 1i].';
  end
        p_t=-b0/(2*a0);
        u_t1(:,iiz)=[real(b_t) imag(b_t) p_t.'].';
    end
        
        iii=iii+1;
 %%SQUAREM       
 rr=u_t1(:,1)-u_t0;
 vv=u_t1(:,2)-u_t1(:,1)-rr;
 if vv==0
     break;
 end
 aat=-norm(rr)/norm(vv);
 zzz=u_t0-2*aat*rr+aat^2*vv;
  while objective_value(y,zzz(1:L).'+1i*zzz(L+1:2*L).',L,N,delta_f,zzz(2*L+1:end),u,c)>objective_value(y,u_t0(1:L)+1i*u_t0(L+1:2*L),L,N,delta_f,u_t0(2*L+1:end),u,c)&&aat~=-1
      aat= (aat-1)/2;
     zzz=u_t0-2*aat*rr+aat^2*vv;
  end       
 
  p_t=zzz(2*L+1:end);
  b_t=zzz(1:L).'+1i*zzz(L+1:2*L).';

end

% find(fz(1:end-1)-fz(2:end)<=0)
% plot(fz)

time=toc;
%%end timing
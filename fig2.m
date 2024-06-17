clear all
close all

Mc=1000;%%number of MC 
L=3; 
B=5000000;
N=128;%%subcarriers

%%sensor position
u=[0 3000 0 ;0 0 3000];
%%source
p=[1498.7;3207.4];

%%noise
randn('state',1)
rn1=randn(L,N,2*Mc);
rn=(rn1(1:L,1:N,1:Mc)+1i*rn1(1:L,1:N,Mc+1:2*Mc))/sqrt(2);

%%path attenuation
b=[   0.3226 - 0.6887i;
   0.4409 - 0.4035i;
   0.4138 - 0.5091i];

z1=1500;
z2=3200;
xx=z1-100:10:z1+100;
yy=z2-100:10:z2+100;

SNRd=-20:5:25;
for i=1:10
SNR=SNRd(i);
parfor mc=1:1000
   [time, p_e]=MO_DPD(B,N,b,rn(:,1:N,mc),u,p,SNR,xx,yy);
   f1(mc)=norm(p_e-p)^2 ;
   Time1(mc)=time;
     [time, p_e]=MO_DPD_new_SQUAREM(B,N,b,rn(:,1:N,mc),u,p,SNR,xx,yy);
   f3(mc)=norm(p_e-p)^2 ; 
   Time3(mc)=time;
end
z1(i)= sqrt(mean(f1))

z3(i)= sqrt(mean(f3))
end

SNRd=-20:5:25;
for i=1:10


SNR=SNRd(i);

SPEB(i)=CRLB(B,N,b,u,p,SNR,L);
end

figure(1)

semilogy(SNRd,z1,'o','linewidth',1.5,'Color', [60/255, 128/255, 0],'markersize',8)
hold on
semilogy(SNRd,z3,'s','linewidth',1.5,'Color', [180/255, 80/255, 190/255],'markersize',8)
hold on
semilogy(SNRd,SPEB,'--k','linewidth',1.5)
legend('MO-DPD','Proposed ML-MM','One-bit CRLB')
xlabel('SNR [dB]')
ylabel('RMSE [m]')
grid on



clear all
close all

Mc=1000;%%MC
L=3; %%
B=5000000;%%bandwidth
N=128;%%subcarriers

%%sensor position
u=[0 3000 0 ;0 0 3000];


%%source
p=[1498.7;3207.4];

% noise
randn('state',1);
rn1=randn(L,N,2*Mc);
rn=(rn1(1:L,1:N,1:Mc)+1i*rn1(1:L,1:N,Mc+1:2*Mc))/sqrt(2);

%%path attenuation


b=[0.3226 - 0.6887i;
   0.4409 - 0.4035i;
   0.4138 - 0.5091i];


SNR=-5;


z1=1500;
z2=3200;


grid=0.2:0.4:2.6;


for i=1:7
N=128;
SNR=10;
xx=z1-200:10:z1+200;
yy=z2-200:10:z2+200;

xx2=z1-200:grid(i):z1+200;
yy2=z2-200:grid(i):z2+200;
for mc=1:1000
    
   [time, p_e]=MO_DPD(B,N,b,rn(:,1:N,mc),u,p,SNR,xx2,yy2);
   f1(mc)=norm(p_e-p)^2 ;
   Time1(mc)=time;
  
     [time, p_e]=MO_DPD_new_SQUAREM(B,N,b,rn(:,1:N,mc),u,p,SNR,xx,yy);
   f3(mc)=norm(p_e-p)^2 ; 
   Time3(mc)=time;

end
fd1(i)= sqrt(mean(f1))
% z2(i)= sqrt(mean(f2))
fd3(i)= sqrt(mean(f3))
Timed1(i)= mean(Time1);
Timed3(i)= mean(Time3);
end


for i=1:7
N=128;

SNR=10;

SPEB(i)=CRLB(B,N,b,u,p,SNR,L);
end

figure(1)

semilogy(grid,fd1,'-o','linewidth',1.5,'Color', [60/255, 128/255, 0])
hold on
semilogy(grid,fd3,'-s','linewidth',1.5,'Color', [180/255, 80/255, 190/255])
hold on
semilogy(grid,SPEB,'--k','linewidth',1.5)
legend('MO-DPD','Proposed ML-MM','One-bit CRLB')
xlabel('Gird size [m]', 'FontName', 'Arial', 'Interpreter', 'tex');
ylabel('RMSE [m]')
grid on

figure(2)

semilogy(grid,Timed1,'-o','linewidth',1.5,'Color', [60/255, 128/255, 0])
hold on
semilogy(grid,mean(Timed3)*ones(1,7),'-s','linewidth',1.5,'Color', [180/255, 80/255, 190/255])
% hold on
% semilogy(1:7,SPEB,'--k','linewidth',1.5)
legend('MO-DPD','Proposed ML-MM')
xlabel('Gird size [m]', 'FontName', 'Arial', 'Interpreter', 'tex');
ylabel('Average time [s]');


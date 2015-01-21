%例2.10.1 多维RLS算法仿真程序
clear;clc;Bushu=2000;
Q=[1 0;0 0.81];A=[-0.8 0;-0.3 0.3];B=[-0.3 0;0.3 0.5];
randn('seed',1);e=sqrt(Q)*randn(2,Bushu);
for t=1:Bushu
    u(1,t)=sin(2*pi*t/100);u(2,t)=cos(2*pi*t/100);
end
y(:,1)=e(:,1);
for i=2:Bushu
    y(:,i)=A*y(:,i-1)+B*u(:,i)+e(:,i);%生成AR()模型
end
%-----------------模型参数辨识------------------%
Na=1;%模型的待估参数Na
N=2;%预估参数的维数N
Nb=1;%手控u的参数个数Nb
NN=(Na+Nb)*N;
fai(:,1)=[0,0,0,0]';
sita(:,1:4)=[0 0 0 0;0 0 0 0];
p(:,1:4)=10^5*eye((Na+Nb)*N);
for i=2:Bushu
    fai(:,i)=[y(1,i-1),y(2,i-1),u(1,i-1),u(2,i-1)]';
    sita(:,(i-1)*NN+1:i*NN)=sita(:,(i-2)*NN+1:(i-1)*NN)+(y(:,i)-sita(:,(i-2)*NN+1:(i-1)*NN)*fai(:,i))*fai(:,i)'*p(:,NN*(i-2)+1:NN*(i-1))/(1+fai(:,i)'*p(:,NN*(i-2)+1:NN*(i-1))*fai(:,i));
    p(:,NN*(i-1)+1:NN*i)=p(:,NN*(i-2)+1:NN*(i-1))-p(:,NN*(i-2)+1:NN*(i-1))*fai(:,i)*fai(:,i)'*p(:,NN*(i-2)+1:NN*(i-1))/(1+fai(:,i)'*p(:,NN*(i-2)+1:NN*(i-1))*fai(:,i));
end 
%----------------对噪声方差的估计----------------%
for i=1:Bushu
    emixiu(:,i)=y(:,i)-sita(:,(i-1)*NN+1:i*NN)*fai(:,i);
end
taoe(:,1:N)=zeros(N,N);
for i=2:Bushu
    taoe(:,N*(i-1)+1:N*i)=taoe(:,N*(i-2)+1:N*(i-1))+[emixiu(:,i)*emixiu(:,i)'-taoe(:,N*(i-2)+1:N*(i-1))]/i;
end  
%--------------------作图部分--------------------%
t=1:Bushu;
subplot(2,2,1);
plot(t,sita(1,4*(t-1)+1),'r',t,sita(1,4*(t-1)+2),'r',t,sita(2,4*(t-1)+1),'r',t,sita(2,4*(t-1)+2),'r');
line([0,Bushu],[-0.8,-0.8]);line([0,Bushu],[0,0]);
line([0,Bushu],[-0.3,-0.3]);line([0,Bushu],[0.3,0.3]);axis([0 Bushu -1 1]);

subplot(2,2,2);
plot(t,sita(1,4*(t-1)+3),'r',t,sita(1,4*(t-1)+4),'r',t,sita(2,4*(t-1)+3),'r',t,sita(2,4*(t-1)+4),'r');
line([0,Bushu],[-0.3,-0.3]);line([0,Bushu],[0,0]);
line([0,Bushu],[0.3,0.3]);line([0,Bushu],[0.5,0.5]);axis([0 Bushu -1 1]);

subplot(2,2,3);
plot(t,taoe(1,2*(t-1)+1),'r',t,taoe(1,2*(t-1)+2),'r',t,taoe(2,2*(t-1)+2),'r');
line([0,Bushu],[1,1]);line([0,Bushu],[0,0]);
line([0,Bushu],[0.81,0.81]);axis([0,Bushu -0.5 1.5]);

%习题2.16 多维RLS算法仿真程序
clear;clc;Bushu=2000;
Q=[0.36 0;0 0.16];A=[0.9 0;0.4 0.6];
randn('seed',1);e=sqrt(Q)*randn(2,Bushu);
y(:,1)=e(:,1);
for i=2:Bushu
    y(:,i)=A*y(:,i-1)+e(:,i);%生成AR()模型
end
%-----------------模型参数辨识------------------%
Na=1;%模型的待估参数Na
N=2;%预估参数的维数N
NN=Na*N;
fai(:,1)=[0,0]';
sita(:,1:2)=[0 0;0 0];
p(:,1:2)=10^5*eye(NN);
for i=2:Bushu
    fai(:,i)=[y(1,i-1),y(2,i-1)]';
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
figure
subplot(2,2,1);
plot(t,sita(1,2*(t-1)+1),'r',t,sita(1,2*(t-1)+2),'r',t,sita(2,2*(t-1)+1),'r',t,sita(2,2*(t-1)+2),'r');
line([0,Bushu],[0.9,0.9]);line([0,Bushu],[0,0]);
line([0,Bushu],[0.4,0.4]);line([0,Bushu],[0.6,0.6]);axis([0 Bushu -0.5 1.2]);

subplot(2,2,2);
plot(t,taoe(1,2*(t-1)+1),'r',t,taoe(1,2*(t-1)+2),'r',t,taoe(2,2*(t-1)+2),'r');
line([0,Bushu],[0.36,0.36]);line([0,Bushu],[0,0]);
line([0,Bushu],[0.16,0.16]);axis([0,Bushu -0.2 0.5]);

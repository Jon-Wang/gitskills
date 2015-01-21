%习题2.5
clear;clc;Bushu=1000;
% * * * * * * 生成模型 * * * * * * %
randn('seed',1);e=sqrt(0.01)*randn(1,Bushu);
y(1)=e(1);y(2)=2.3*y(1)+e(2);y(3)=2.3*y(2)-1.85*y(1)+e(3);y(4)=2.3*y(3)-1.85*y(2)+0.601*y(1)+e(4);
for i=5:Bushu
    y(i)=2.3*y(i-1)-1.85*y(i-2)+0.601*y(i-3)-0.063*y(i-4)+e(i);
end
% * * * * * * 系统参数辨识 * * * * * * %
sita(:,Bushu)=zeros(4,1);
p(:,1:4)=10^5*eye(4);p(:,5:8)=10^5*eye(4);p(:,9:12)=10^5*eye(4);p(:,13:16)=10^5*eye(4);
for i=5:Bushu
    fai(:,i)=[y(i-1),y(i-2),y(i-3),y(i-4)]';
    sita(:,i)=sita(:,i-1)+p(:,4*(i-1)-3:4*(i-1))*fai(:,i)/(1+fai(:,i)'*p(:,4*(i-1)-3:4*(i-1))...
        *fai(:,i))*(y(i)-fai(:,i)'*sita(:,i-1));
p(:,4*i-3:4*i)=p(:,4*(i-1)-3:4*(i-1))-p(:,4*(i-1)-3:4*(i-1))*fai(:,i)*fai(:,i)'*p(:,4*(i-1)-3:4*(i-1)	)/(1+fai(:,i)'...
    *p(:,4*(i-1)-3:4*(i-1))*fai(:,i));
end
% * * * * * * 噪声方差的估计 * * * * * * %
for i=1:Bushu
    emixiu(i)=y(i)-fai(1:4,i)'*sita(1:4,i);
end
taoe(1)=emixiu(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*[emixiu(i)^2-taoe(i-1)];
end
% * * * * * * 作图部分 * * * * * * %
t=1:Bushu;
figure
subplot(2,2,1);plot(t,sita(1,t),'r');line([0,Bushu],[2.3,2.3]);axis([0 Bushu 0 3]);sita(1,Bushu)
subplot(2,2,2);plot(t,sita(2,t),'r');line([0,Bushu],[-1.85,-1.85]);axis([0 Bushu -2 0]);sita(2,Bushu)
subplot(2,2,3);plot(t,sita(3,t),'r');line([0,Bushu],[0.601,0.601]);axis([0 Bushu -2 4]);sita(3,Bushu)
subplot(2,2,4);plot(t,sita(4,t),'r');line([0,Bushu],[0.063,0.063]);axis([0 Bushu -4 2]);sita(4,Bushu)
figure
subplot(1,1,1);plot(t,taoe,'r');line([0,Bushu],[0.01,0.01]);axis([0 Bushu 0 0.03]);taoe(Bushu)





%例2.2.1 RLS算法仿真程序
clear;clc;Bushu=600;
randn('seed',1);e=sqrt(0.25^2)*randn(1,Bushu);
y(1)=e(1);y(2)=1.4*y(1)+e(2);
for i=3:Bushu
    y(i)=1.4*y(i-1)-0.45*y(i-2)+e(i);%生成AR(2)模型
end
%模型参数辨识
sita(:,Bushu)=zeros(2,1);p(:,1:2)=10^5*eye(2);p(:,3:4)=10^5*eye(2);
for i=3:Bushu
    fai(:,i)=[y(i-1),y(i-2)]';
    sita(:,i)=sita(:,i-1)+p(:,2*(i-1)-1:2*(i-1))*fai(:,i)/(1+fai(:,i)'*p(:,2*(i-1)-1:2*(i-1))*fai(:,i))*(y(i)-fai(:,i)'*sita(:,i-1));
    p(:,2*i-1:2*i)=p(:,2*(i-1)-1:2*(i-1))-p(:,2*(i-1)-1:2*(i-1))*fai(:,i)*fai(:,i)'*p(:,2*(i-1)-1:2*(i-1))/(1+fai(:,i)'*p(:,2*(i-1)-1:2*(i-1))*fai(:,i));
end
%对噪声方差的估计
for i=1:Bushu
    emixiu(i)=y(i)-fai(1:2,i)'*sita(1:2,i);
end
taoe(1)=emixiu(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*[emixiu(i)^2-taoe(i-1)];
end
%作图部分
t=1:Bushu;
subplot(2,2,1);plot(t,sita(1,t),'r');line([0,Bushu],[1.4,1.4])
subplot(2,2,2);plot(t,sita(2,t),'r');line([0,Bushu],[-0.45,-0.45])
subplot(2,2,3);plot(t,taoe,'r');line([0,Bushu],[0.0625,0.0625])
sita(1,600)
sita(2,600)
taoe(600)

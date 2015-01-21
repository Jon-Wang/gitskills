%习题2.6 RLS算法仿真程序
clear;clc;Bushu=800;
randn('seed',1);e=sqrt(0.5^2)*randn(1,Bushu);
y(1)=e(1);y(2)=2.8*y(1)+e(2);y(3)=2.8*y(1)-2.6*y(2)+e(3);
for i=4:Bushu
    y(i)=2.8*y(i-1)-2.6*y(i-2)+0.8*y(i-3)+e(i);%生成AR(3)模型
end
%模型参数辨识
sita(:,Bushu)=zeros(3,1);p(:,1:3)=10^5*eye(3);p(:,4:6)=10^5*eye(3);
for i=4:Bushu
    fai(:,i)=[y(i-1),y(i-2),y(i-3)]';
    sita(:,i)=sita(:,i-1)+p(:,3*(i-1)-5:3*(i-1)-3)*fai(:,i)/(1+fai(:,i)'*p(:,3*(i-1)-5:3*(i-1)-3)*fai(:,i))*(y(i)-fai(:,i)'*sita(:,i-1));
    p(:,3*i-5:3*i-3)=p(:,3*(i-1)-5:3*(i-1)-3)-p(:,3*(i-1)-5:3*(i-1)-3)*fai(:,i)*fai(:,i)'*p(:,3*(i-1)-5:3*(i-1)-3)/(1+fai(:,i)'*p(:,3*(i-1)-5:3*(i-1)-3)*fai(:,i));
end
p
%对噪声方差的估计
for i=1:Bushu
    emixiu(i)=y(i)-fai(1:3,i)'*sita(1:3,i);
end
taoe(1)=emixiu(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*[emixiu(i)^2-taoe(i-1)];
end
%作图部分
t=1:Bushu;
subplot(2,2,1);plot(t,sita(1,t),'r');line([0,Bushu],[2.8,2.8]);axis([0,Bushu,2,3]);
subplot(2,2,2);plot(t,sita(2,t),'r');line([0,Bushu],[-2.6,-2.6]);axis([0,Bushu,-3,-1.5]);
subplot(2,2,3);plot(t,sita(3,t),'r');line([0,Bushu],[0.8,0.8]);axis([0,Bushu,0,1]);
subplot(2,2,4);plot(t,taoe,'r');line([0,Bushu],[0.25,0.25]);axis([0,Bushu,0,1]);
sita(1,800)
sita(2,800)
sita(3,800)
taoe(800)


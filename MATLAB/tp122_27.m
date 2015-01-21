%习题2.7 带遗忘因子的RLS算法估计
clear;clc;Bushu=6000;lamda=0.99;pai=3.14159;
randn('seed',4);e=sqrt(1)*randn(1,Bushu);
for i=1:Bushu
u(i)=(-1)^(fix(i/100));
end
figure
t=1:Bushu
plot(t,u(t));axis([0,Bushu,-2,2]);
for i=1:Bushu
    a(i)=0.5-0.15*sin(2*pai*i/2000);
    b(i)=0.81-0.2*cos(2*pai*i/2000);
end
y(1)=e(1);
for i=2:Bushu
    y(i)=a(i)*y(i-1)+b(i)*u(i-1)+e(i);
end
fai(:,1)=[0 0]';sita(:,1)=zeros(2,1);
p(:,1:2)=10^5*eye(2);p(:,3:4)=10^5*eye(2);
for i=2:Bushu
    fai(:,i)=[y(i-1),u(i-1)];
    sita(:,i)=sita(:,i-1)+p(:,2*(i-1)-1:2*(i-1))*fai(:,i)/(lamda+fai(:,i)'*p(:,2*(i-1)-1:2*(i-1))...
        *fai(:,i))*(y(i)-fai(:,i)'*sita(:,i-1));
    p(:,2*i-1:2*i)=1/lamda*p(:,2*(i-1)-1:2*(i-1))-p(:,2*(i-1)-1:2*(i-1))*fai(:,i)*fai(:,i)'*p(:,2*(i-1)-1:2*(i-1))/(lamda+fai(:,i)'...
    *p(:,2*(i-1)-1:2*(i-1))*fai(:,i));
end
%对噪声方差的估计
for i=1:Bushu
    emixiu(i)=y(i)-fai(1:2,i)'*sita(1:2,i);
end
taoe(1)=emixiu(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*[emixiu(i)^2-taoe(i-1)];
end

t=1:Bushu;
figure
plot(t,a(t),t,sita(1,t),'r');axis([0,Bushu,-2,4]);
figure
plot(t,b(t),t,sita(2,t),'r');axis([0,Bushu,-2,4]);
figure
plot(t,taoe,'r');line([0,Bushu],[1,1])

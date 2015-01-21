%例3.8.1 状态分量解耦滤波器仿真源程序
clc;clear all;Bushu=300;
%---------------生成系统模型----------------%
Q=1;R=1;T=1;fai=[1 T;0 1];gama=[T^2/2 T]';H=[1 0];
randn('seed',1);w=sqrt(Q)*randn(1,Bushu);
randn('seed',2);v=sqrt(R)*randn(1,Bushu);
x(:,1)=zeros(2,1);
for t=2:Bushu
    x(:,t)=(fai*x(:,t-1)+gama*w(t));
end
t=1:Bushu;
x(:,t)=x(:,t)/2;
for t=1:Bushu
    y(t)=H*x(:,t)+v(t);
end
%-------------Riccati方程---------------%
n=2;
pp(:,:)=eye(n);
 for i=2:Bushu+1;
    temp=pp(:,n*(i-2)+1:n*(i-1));
    pp(:,n*(i-1)+1:n*i)=fai*[temp-temp*H'*inv(H*temp*H'+R)*H*temp]*fai'+gama*Q*gama';
    kf(:,i-1)=temp*H'*inv(H*temp*H'+R);
 end
 kf=kf(:,end);
 pusaif=(eye(n)-kf*H)*fai;
 %-------------解耦部分---------------%
 syms q;
 f=eye(2)-pusaif*q;
 phisym=det(f);phi=sym2poly(phisym);
 adj=simple(inv(f)*det(f));
 k1sym=[1 0]*adj*kf; k1=sym2poly(k1sym);
 k2sym=[0 1]*adj*kf; k2=sym2poly(k2sym);
 
 xjian(1,1)=-phi(2)*y(1);
 xjian(1,2)=-phi(2)*xjian(1,1)+k1(2)*y(2)+k1(1)*y(1);
 for i=3:Bushu
     xjian(1,i)=-phi(2)*xjian(1,i-1)-phi(1)*xjian(1,i-2)+k1(2)*y(i)+k1(1)*y(i-1);
 end
 
 xjian(2,1)=-phi(2)*y(1);
 xjian(2,2)=-phi(2)*xjian(1,1)+k2(2)*y(1)+k2(1)*y(2);
 for i=3:Bushu
     xjian(2,i)=-phi(2)*xjian(2,i-1)-phi(1)*xjian(2,i-2)+k2(2)*y(i)+k2(1)*y(i-1);
 end
 %-------------作图部分------------%
 t=1:Bushu;
 subplot(2,2,1);plot(t,x(1,t),'b',t,xjian(1,t),'r:');
 subplot(2,2,2);plot(t,x(2,t),'b',t,xjian(2,t),'r:');
 
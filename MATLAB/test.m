clc;clear;bushu=100;
q=1;r=1;
fai=[0.9 1;0 0.5];gama=[1 1]';h=[1 0];
randn('seed',1);w=sqrt(q)*randn(2,bushu+10);
randn('seed',1);v=sqrt(r)*randn(1,bushu+10);
x(:,1)=zeros(2,1);
y(1)=h*x(:,1)+v(1);
for i=2:bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=h*x(:,i)+v(i);
end
%生成riccati方程%
n=2;
pp(:,:)=eye(n);
for i=2:bushu+1
    temp=pp(:,n*(i-2)+1:n*(i-1));
    pp(:,n*(i-1)+1:n*(i))=fai*[temp-temp*h'*inv(h*temp*h'+r)*h*temp]*fai'+gama*q*gama';
    kf(:,i-1)=temp*h'*inv(h*temp*h'+r);
    kp(:,i-1)=fai*kf(:,i-1);
end
%滤波算法%
xjian(:,bushu)=zeros(n,1);pusaif=(eye(n)-kf(:,bushu)*h)*fai;
for i=2:bushu
    xjian(:,i)=pusaif*xjian(:,i-1)+kf(:,bushu)*y(:,i);
end
 xjian(:,i)
%预报算法%
xxjian(:,bushu)=zeros(n,1);pusaip=fai*(eye(n)-kf(:,bushu)*h);
for i=2:bushu
     xxjian(:,i)=pusaip*xxjian(:,i-1)+kp(:,bushu)*y(:,i);
end
 xxjian(:,bushu)
%作图部分%
t=1:bushu;
subplot(2,2,1);plot(t,x(1,t),'k',t,xjian(1,t),'r.');axis([0 bushu,-15 10]);
subplot(2,2,2);plot(t,x(2,t),'k',t,xjian(2,t),'r.');axis([0 bushu,-3 3]);
subplot(2,2,3);plot(t,x(1,t),'k',t,xxjian(1,t),'r.');axis([0 bushu,-15 10]);
subplot(2,2,4);plot(t,x(2,t),'k',t,xxjian(2,t),'r.');axis([0 bushu,-3 3]);
%习题3.3 Kalman滤波算法仿真程序
clear;clc;Bushu=300;
Q=0.64; %w(t)的方差
R=0.36;    %V(t)的方差
randn('seed',1);w=sqrt(Q)*randn(1,Bushu+2);
randn('seed',2);v=sqrt(R)*randn(1,Bushu+2);
fai=[0.6 -0.4;0.2 1.1];gama=[1 2]';H=[1 1];
x(:,1)=[0,0]';y(1)=H*x(:,1)+v(1);
for i=2:Bushu+2
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
%------------------Kalman滤波器----------------%
n=2; %维数
xjian(:,1)=zeros(2,1);p(:,:)=zeros(2);
for i=1:Bushu+1
   xxjian(:,i+1)=fai*xjian(:,i);%xxjian为x的预报
   xxxjian(:,i+1)=fai^2*xjian(:,i);%两步预报
   e(:,i+1)=y(i+1)-H*xxjian(:,i+1);%新息
   pp(:,n*(i-1)+1:n*i)=fai*p(:,n*(i-1)+1:n*i)*fai'+gama*Q*gama';
   kf(:,i+1)=pp(:,n*(i-1)+1:n*i)*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'+R);
   xjian(:,i+1)=xxjian(:,i+1)+kf(:,i+1)*e(:,i+1);
   p(:,n*i+1:n*(i+1))=[eye(2)-kf(:,i+1)*H]*pp(:,n*(i-1)+1:n*i);
   %mxjian(:,i+1)=xjian(:,i)+pp(:,n*(i-2)+1:n*(i-1))*[eye(2)-kf(:,i+1)*H]'*fai'*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'*R)*e(:,i+1);
end 
for i=2:Bushu
    mxjian(:,i+1)=xjian(:,i)+pp(:,n*(i-2)+1:n*(i-1))*[eye(2)-kf(:,i)*H]'*fai'*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'*R)*e(:,i+1);
    mmxjian(:,i+1)=mxjian(:,i)+pp(:,n*(i-2)+1:n*(i-1))*[eye(2)-kf(:,i)*H]'*fai'*[eye(2)-kf(:,i+2)*H]'*fai'*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'*R)*e(:,i+1);
end

%-------------------作图部分--------------------%
t=1:Bushu;
subplot(2,2,1);plot(t,x(1,t),'b',t,xjian(1,t),'r:');legend('滤波器','估值');%axis([10 100 -150 150]);legend('两步预报器');
subplot(2,2,2);plot(t,x(2,t),'b',t,xjian(2,t),'r:');legend('滤波器','估值');%axis([10 100 -150 150]);
subplot(2,2,3);plot(t,x(1,t),'b',t,xxxjian(1,t),'r:');legend('两步预报器','估值');
subplot(2,2,4);plot(t,x(2,t),'b',t,xxxjian(2,t),'r:');legend('两步预报器','估值');
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,mmxjian(1,t),'r:');legend('两步平滑器','估值');
subplot(2,2,2);plot(t,x(2,t),'b',t,mmxjian(2,t),'r:');legend('两步平滑器','估值');

